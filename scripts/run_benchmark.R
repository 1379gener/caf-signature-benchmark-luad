#!/usr/bin/env Rscript
# ==============================================================================
# run_benchmark.R —— 零分布基准（阶段 1 驱动脚本）
#
# 对每个「队列 × signature」组合：
#   1. 覆盖率检查（SAP §4.2，< 70% 判为不可评估）
#   2. 计算真实 signature 分数与单变量 Cox 的 Wald z
#   3. 生成零分布 A / B（A 类）或加权复合零分布（B 类），各 1000 次
#   4. 计算经验 p 值（SAP §4.3，双侧）
#
# 特性：
#   - 断点续跑：每个组合算完立刻落盘，重跑自动跳过已完成者
#   - 拒绝在 signature 清单锁定之前运行（SAP §3.6）
#   - 全部 seed 由 make_seed() 确定性生成，不得事后更换
#
# 用法：
#   Rscript run_benchmark.R --dry-run            # 只打印工作计划，不计算
#   Rscript run_benchmark.R --test               # 小规模冒烟测试（n_null=50，单队列）
#   Rscript run_benchmark.R                      # 全量
#   Rscript run_benchmark.R --cohort=TCGA_LUAD   # 限定队列
#   Rscript run_benchmark.R --sig=SIG-001,SIG-003
#   Rscript run_benchmark.R --workers=8          # GSVA 并行核数
#   Rscript run_benchmark.R --force              # 忽略已有结果，全部重算
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(k, d = NA) { h <- grep(paste0("^--", k, "="), args, value = TRUE)
                              if (length(h)) sub(paste0("^--", k, "="), "", h[1]) else d }
flag <- function(k) paste0("--", k) %in% args

DRY     <- flag("dry-run")
TEST    <- flag("test")
FORCE   <- flag("force")
WORKERS <- as.integer(opt("workers", "1"))
N_NULL  <- as.integer(opt("n-null", if (TEST) "50" else "1000"))

ROOT   <- path.expand("~/thesis_LUAD_CAF")
REGDIR <- file.path(ROOT, "03_signature_registry")
CHDIR  <- file.path(ROOT, "02_chapters/ch04_benchmark")
OUTDIR <- file.path(CHDIR, "data_derived/null_benchmark")
LOGDIR <- file.path(CHDIR, "logs")
for (d in c(OUTDIR, LOGDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

suppressMessages({
  source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
  source(file.path(CHDIR, "scripts/null_engine.R"))
})

LOG <- file.path(LOGDIR, format(Sys.time(), "benchmark_%Y%m%d_%H%M.log"))
sink(LOG, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
sh <- function(t) cat("\n  --- ", t, "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-28s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("零分布基准 · 阶段 1\n时间: ", format(Sys.time()), "\n日志: ", LOG, "\n", sep = "")
kv("模式", paste(c(if (DRY) "dry-run", if (TEST) "test", if (FORCE) "force",
                   if (!DRY && !TEST) "全量")[1], sep = ""))
kv("n_null", N_NULL)
kv("GSVA 并行核数", WORKERS)

# ============================================================ 1. 锁定闸
hr("1. signature 清单")
LOCK <- require_locked_registry(REGDIR)     # 未锁定则在此停止
kv("锁定文件", basename(LOCK$file))
kv("md5", LOCK$md5)
reg <- LOCK$table
need <- c("sig_id", "class")
if (!all(need %in% names(reg))) stop("锁定文件缺列: ", paste(setdiff(need, names(reg)), collapse = ", "))
kv("signature 数", sprintf("%d（A 类 %d，B 类 %d）", nrow(reg),
    sum(reg$class == "A"), sum(reg$class == "B")))

read_sig <- function(sid, cls) {
  fg <- file.path(REGDIR, "genes", paste0(sid, "_genes.txt"))
  if (!file.exists(fg)) stop("缺基因文件: ", fg)
  g <- trimws(readLines(fg)); g <- g[nzchar(g)]
  w <- NULL
  if (cls == "B") {
    fw <- file.path(REGDIR, "genes", paste0(sid, "_weights.csv"))
    if (!file.exists(fw)) stop("B 类缺权重文件: ", fw)
    d <- data.table::fread(fw)
    w <- setNames(as.numeric(d$weight), d$gene)[g]
    if (anyNA(w)) stop(sid, " 的权重与基因列表不匹配")
  }
  list(genes = g, weights = w)
}

# ============================================================ 2. 队列
hr("2. 队列")
cohort_filter <- opt("cohort"); sig_filter <- opt("sig")
CO <- list()
for (a in COHORTS) {
  x <- tryCatch(suppressMessages(load_cohort(a, verbose = FALSE)), error = function(e) NULL)
  if (is.null(x)) next
  if (!x$has_survival) { cat(sprintf("  %-12s 跳过（无生存数据）\n", a)); next }
  CO[[a]] <- x
  cat(sprintf("  %-12s %4d 样本  %3d 事件  %6d 基因  %s\n", a, x$n, x$n_events, x$n_genes,
      ifelse(x$n_events < 25, "!! 事件 < 25，降为敏感性分析队列（SAP §2.4）", "")))
}
if (!is.na(cohort_filter)) CO <- CO[intersect(names(CO), strsplit(cohort_filter, ",")[[1]])]
if (TEST && length(CO) > 1) CO <- CO[which.min(vapply(CO, function(x) x$n, 0L))]
kv("纳入队列", names(CO))

sig_ids <- reg$sig_id
if (!is.na(sig_filter)) sig_ids <- intersect(sig_ids, strsplit(sig_filter, ",")[[1]])
if (TEST) sig_ids <- head(sig_ids, 2L)
kv("纳入 signature", sig_ids)

# ============================================================ 3. 工作计划
hr("3. 工作计划")
plan <- data.table::CJ(cohort = names(CO), sig_id = sig_ids, sorted = FALSE)
plan <- merge(plan, reg[, c("sig_id", "class")], by = "sig_id", sort = FALSE)
plan[, outfile := file.path(OUTDIR, sprintf("%s__%s.rds", cohort, sig_id))]
plan[, done := file.exists(outfile) & !FORCE]
cat(sprintf("  组合总数 %d  已完成 %d  待算 %d\n", nrow(plan), sum(plan$done), sum(!plan$done)))

if (DRY) {
  print(plan[, .(cohort, sig_id, class, done)], row.names = FALSE)
  cat("\n  dry-run 结束，未做任何计算。\n"); quit(status = 0)
}

BP <- NULL
if (WORKERS > 1L && requireNamespace("BiocParallel", quietly = TRUE))
  BP <- BiocParallel::MulticoreParam(workers = WORKERS)

# ============================================================ 4. 主循环
hr("4. 计算")
t_all <- Sys.time()
POOLS <- list(); RZ <- list()

for (r in seq_len(nrow(plan))) {
  cohort <- plan$cohort[r]; sid <- plan$sig_id[r]; cls <- plan$class[r]
  tag <- sprintf("%s x %s (%s)", cohort, sid, cls)
  if (plan$done[r]) { cat(sprintf("  [%2d/%2d] %-34s 已存在，跳过\n", r, nrow(plan), tag)); next }

  cat(sprintf("\n  [%2d/%2d] %s\n", r, nrow(plan), tag))
  t0 <- Sys.time()
  co <- CO[[cohort]]; sg <- read_sig(sid, cls)

  # ---- 覆盖率（SAP §4.2）
  cv <- gene_coverage(sg$genes, co$expr)
  cat(sprintf("        覆盖率 %d/%d = %.1f%%  %s\n", cv$n_present, cv$n_total,
              100 * cv$frac, ifelse(cv$evaluable, "", "-> 不可评估，跳过计算")))
  base_row <- list(cohort = cohort, sig_id = sid, class = cls,
    n_genes = cv$n_total, n_present = cv$n_present, coverage = cv$frac,
    evaluable = cv$evaluable, n = co$n, n_events = co$n_events,
    low_power = co$n_events < 25, missing_genes = paste(cv$missing, collapse = ";"))
  if (!cv$evaluable) { saveRDS(base_row, plan$outfile[r]); next }

  # ---- 候选池与秩矩阵（每队列只建一次）
  if (is.null(POOLS[[cohort]])) {
    POOLS[[cohort]] <- build_pool(co$expr)
    RZ[[cohort]]    <- build_rank_z(co$expr, POOLS[[cohort]])
    cat(sprintf("        候选池 %d 基因，%d 个联合箱\n",
                length(POOLS[[cohort]]$genes), length(unique(POOLS[[cohort]]$bin))))
  }
  pool <- POOLS[[cohort]]; rz <- RZ[[cohort]]
  gp <- intersect(sg$genes, pool$genes)
  if (length(gp) < 2L) { base_row$evaluable <- FALSE
    base_row$note <- "候选池中基因不足 2 个"; saveRDS(base_row, plan$outfile[r]); next }

  tt <- co$clin$OS_time; ee <- co$clin$OS_event

  # ---- 真实分数与 Cox
  if (cls == "A") {
    real_score <- score_sets_ssgsea(co$expr, list(gp), prefix = "real")[1, ]
  } else {
    real_score <- score_sets_weighted(co$expr, list(gp), sg$weights[gp])[1, ]
  }
  rc <- cox_one(real_score, tt, ee)
  cat(sprintf("        真实: z = %+.3f  HR %.3f (%.3f-%.3f)  p = %.4g\n",
              rc$z, rc$hr, rc$lo, rc$hi, rc$p))

  res <- c(base_row, list(z_real = rc$z, hr_real = rc$hr, hr_lo = rc$lo,
                          hr_hi = rc$hi, p_cox = rc$p, n_used = rc$n, ev_used = rc$events))

  # ---- 零分布
  run_null <- function(sets, null_type, weights = NULL) {
    sc <- if (is.null(weights)) {
      score_sets_ssgsea(co$expr, sets, prefix = null_type, bp = BP)
    } else {
      score_sets_weighted(co$expr, sets, weights)
    }
    zn <- cox_z_batch(sc, tt, ee)
    p  <- empirical_p(rc$z, zn)
    cat(sprintf("        零分布 %s: 有效 %d/%d  |z| 中位 %.2f  95%%分位 %.2f  p_emp = %.4f\n",
        null_type, sum(is.finite(zn)), length(zn), median(abs(zn), na.rm = TRUE),
        quantile(abs(zn), .95, na.rm = TRUE), p))
    list(z_null = zn, p_emp = p)
  }

  sA <- sample_null_A(pool, gp, n_null = N_NULL, seed = make_seed(cohort, sid, "A"))
  rA <- run_null(sA$sets, "A", weights = if (cls == "B") sg$weights[gp] else NULL)
  res$p_emp_A <- rA$p_emp; res$n_fallback_A <- sA$n_fallback
  res$seed_A <- make_seed(cohort, sid, "A")

  sB <- sample_null_B(pool, rz, gp, n_null = N_NULL, seed = make_seed(cohort, sid, "B"))
  rB <- run_null(sB$sets, "B", weights = if (cls == "B") sg$weights[gp] else NULL)
  res$p_emp_B <- rB$p_emp; res$seed_B <- make_seed(cohort, sid, "B")
  res$coexpr_real <- sB$real_coexpr; res$coexpr_null_med <- sB$null_coexpr_median
  res$coexpr_ratio <- sB$coexpr_ratio; res$coexpr_qc_pass <- sB$qc_pass
  res$n_relax_B <- sB$n_relax
  cat(sprintf("        零分布 B 质量验证: 真实 %.4f / 零分布 %.4f = %.3f  -> %s\n",
      sB$real_coexpr, sB$null_coexpr_median, sB$coexpr_ratio,
      ifelse(sB$qc_pass, "通过", "!! 未通过，本组合的 B 类结论作废（SAP §5.2）")))

  res$z_null_A <- rA$z_null; res$z_null_B <- rB$z_null
  res$elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  res$n_null <- N_NULL; res$engine_version <- ENGINE$version
  res$registry_md5 <- LOCK$md5; res$run_time <- format(Sys.time())
  saveRDS(res, plan$outfile[r])
  cat(sprintf("        完成，耗时 %.1f 分钟\n", res$elapsed_min))
}

# ============================================================ 5. 汇总
hr("5. 汇总")
fs <- list.files(OUTDIR, pattern = "\\.rds$", full.names = TRUE)
if (!length(fs)) { cat("  无结果文件\n"); quit(status = 0) }
rows <- lapply(fs, function(f) {
  x <- readRDS(f); x$z_null_A <- NULL; x$z_null_B <- NULL
  as.data.frame(x, stringsAsFactors = FALSE) })
S <- data.table::rbindlist(rows, fill = TRUE)
data.table::setorder(S, sig_id, cohort)

f_sum <- file.path(CHDIR, "data_derived",
                   sprintf("null_benchmark_summary_%s.csv", format(Sys.Date(), "%Y%m%d")))
data.table::fwrite(S, f_sum)
kv("汇总表", f_sum)
kv("行数", nrow(S))

show <- intersect(c("cohort","sig_id","class","coverage","n_events","z_real",
                    "hr_real","p_emp_A","p_emp_B","coexpr_ratio","coexpr_qc_pass",
                    "low_power","evaluable"), names(S))
print(S[, ..show], row.names = FALSE)

if ("p_emp_A" %in% names(S)) {
  sh("通过队列数（p_emp < 0.05，排除不可评估与低功效）")
  ok <- S[evaluable == TRUE & low_power == FALSE]
  if (nrow(ok)) {
    agg <- ok[, .(n_cohort = .N,
                  pass_A = sum(p_emp_A < 0.05, na.rm = TRUE),
                  pass_B = sum(p_emp_B < 0.05 & coexpr_qc_pass, na.rm = TRUE)),
              by = .(sig_id, class)]
    print(agg, row.names = FALSE)
    cat("\n  注：pass_B 仅计入共表达质量验证通过的组合。\n")
  }
}

hr("6. 环境")
cat(sprintf("  总耗时 %.1f 分钟\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
cat("  引擎 v", ENGINE$version, "  清单 md5 ", LOCK$md5, "\n", sep = "")
print(sessionInfo())
cat("\n日志: ", LOG, "\n", sep = "")
