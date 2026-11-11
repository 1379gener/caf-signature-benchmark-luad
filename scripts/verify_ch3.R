#!/usr/bin/env Rscript
# ==============================================================================
# verify_ch3.R —— 第 3 章 / SCI-2 阶段 1 头条数字的独立复算
#
#   放置：~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/verify_ch3.R
#   用法：Rscript verify_ch3.R
#
#   独立性说明
#   - 真实分数与 Cox z：从 Base A 原始表达矩阵重新打分（A 类直接调 GSVA ssGSEA，
#     B 类自写加权 Z 和），不读 RDS 里的 z_real；候选池沿用 null_engine::build_pool()
#     （池的定义属于方法本身）
#   - 经验 p：用 RDS 中存的零分布 z 向量按 SAP §4.3 公式重算，不读 p_emp_*
#   - 汇总数字（分类、通过率、角色、过拟合梯度）全部由复算值现算
#   - 「期望值」= 20260905 中间表所声称的值，只作对照，不作为输入
#
#   输出（data_derived/）：verify_ch3_combo_<日期>.csv、verify_ch3_summary_<日期>.csv、
#                          scoring_coverage_audit_<日期>.csv；日志写 logs/
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA) })
ROOT  <- path.expand("~/thesis_LUAD_CAF")
CH    <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD  <- file.path(ROOT, "03_signature_registry")
DD    <- file.path(CH, "data_derived")
STAMP <- format(Sys.Date(), "%Y%m%d")
dir.create(file.path(CH, "logs"), showWarnings = FALSE)
LOG <- file.path(CH, "logs", format(Sys.time(), "verify_ch3_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
CHK <- list()
check <- function(id, what, expected, observed, tol = 0) {
  ok <- if (is.numeric(expected) && is.numeric(observed)) isTRUE(all(abs(expected - observed) <= tol)) else
          identical(as.character(expected), as.character(observed))
  CHK[[length(CHK) + 1]] <<- data.table(id = id, what = what, expected = paste(expected, collapse = "/"),
                                       observed = paste(observed, collapse = "/"), status = ifelse(ok, "PASS", "FAIL"))
  cat(sprintf("  [%s] %-6s %-46s 期望 %-14s 复算 %s\n", ifelse(ok, "PASS", "FAIL"), id, what,
              paste(expected, collapse = "/"), paste(observed, collapse = "/")))
}
cat("verify_ch3  ", format(Sys.time()), "\n日志: ", LOG, "\n", sep = "")
cat("GSVA ", as.character(packageVersion("GSVA")), "\n", sep = "")

suppressMessages({
  source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
  source(file.path(CH, "scripts/null_engine.R"))
})

# ============================================================ 0. 锁定物与环境
hr("0. 锁定物与环境")
LOCK <- require_locked_registry(REGD); reg <- LOCK$table
check("V0.1", "签名注册表 md5", "4ad3218931d95dfb6ee67a7ba07b9a75", LOCK$md5)
check("V0.2", "GSVA 版本", "1.48.3", as.character(packageVersion("GSVA")))
check("V0.3", "注册签名数", 12, nrow(reg))

# ============================================================ 1. Base A
hr("1. Base A：逐队列读入、断言、建候选池")
CO <- list(); POOL <- list(); MAN <- list()
for (a in COHORTS) {
  co <- load_cohort(a, verbose = FALSE)
  m  <- co$man
  stopifnot(co$n == as.integer(m$n_final[1]))
  if (co$has_survival) stopifnot(co$n_events == as.integer(m$n_events[1]))
  MAN[[a]] <- data.table(cohort = a, n = co$n, events = co$n_events, n_genes = co$n_genes,
                         has_survival = co$has_survival, platform = m$platform[1])
  if (co$has_survival) { CO[[a]] <- co; POOL[[a]] <- build_pool(co$expr) }
}
MAN <- rbindlist(MAN); print(MAN, row.names = FALSE)
S <- MAN[has_survival == TRUE]
check("V1.1", "队列数（含 CPTAC）", 7, nrow(MAN))
check("V1.2", "生存队列数", 6, nrow(S))
check("V1.3", "生存队列样本合计", 1697, sum(S$n))
check("V1.4", "生存队列事件合计", 655, sum(S$events))
check("V1.5", "GSE31210 样本数（已排除 exclude）", 204, S[cohort == "GSE31210", n])

# ============================================================ 2. 逐组合复算
hr("2. 逐组合复算（真实分数 + Cox；经验 p 由零分布向量重算）")
read_sig <- function(sid, cls) {
  g <- trimws(readLines(file.path(REGD, "genes", paste0(sid, "_genes.txt")))); g <- g[nzchar(g)]
  w <- NULL
  if (cls == "B") { d <- fread(file.path(REGD, "genes", paste0(sid, "_weights.csv")))
                    w <- setNames(as.numeric(d$weight), d$gene)[g]; stopifnot(!anyNA(w)) }
  list(genes = g, weights = w)
}
wsum <- function(expr, genes, w) {                       # 自写：各基因队列内 Z 后按权重相加
  X <- expr[genes, , drop = FALSE]
  Z <- t(apply(X, 1, function(v) { s <- sd(v); if (!is.finite(s) || s == 0) v * 0 else (v - mean(v)) / s }))
  as.numeric(crossprod(w[genes], Z))
}
coxz <- function(sc, t, e) { ok <- !is.na(t) & t > 0 & !is.na(e) & is.finite(sc)
  s <- summary(coxph(Surv(t[ok], e[ok]) ~ as.numeric(scale(sc[ok]))))
  c(z = s$coefficients[1, "z"], hr = s$conf.int[1, 1], n = sum(ok), ev = sum(e[ok])) }
emp_p <- function(zr, zn) { zn <- zn[is.finite(zn)]; (1 + sum(abs(zn) >= abs(zr))) / (1 + length(zn)) }

rds <- list.files(file.path(DD, "null_benchmark"), pattern = "\\.rds$", full.names = TRUE)
check("V2.1", "组合 RDS 数", 72, length(rds))
R <- rbindlist(lapply(rds, function(f) {
  x <- readRDS(f); co <- CO[[x$cohort]]
  sg <- read_sig(x$sig_id, x$class)
  cv <- gene_coverage(sg$genes, co$expr)
  row <- data.table(cohort = x$cohort, sig_id = x$sig_id, class = x$class,
                    n_sig = cv$n_total, n_present = cv$n_present, coverage = cv$frac,
                    coverage_rds = x$coverage, evaluable = isTRUE(x$evaluable),
                    universe = nrow(co$expr), pool_size = length(POOL[[x$cohort]]$genes),
                    missing = paste(cv$missing, collapse = ";"), n_null = if (is.null(x$n_null)) NA_integer_ else x$n_null)
  if (!row$evaluable) return(row)
  gp <- intersect(sg$genes, POOL[[x$cohort]]$genes)
  sc <- if (x$class == "A") {
    as.numeric(GSVA::gsva(co$expr, list(real = gp), method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)[1, ])
  } else wsum(co$expr, gp, sg$weights)
  cz <- coxz(sc, co$clin$OS_time, co$clin$OS_event)
  zA <- x$z_null_A[is.finite(x$z_null_A)]; zB <- x$z_null_B[is.finite(x$z_null_B)]
  row[, `:=`(n_in_pool = length(gp), z_real = cz[["z"]], z_real_rds = x$z_real, hr = cz[["hr"]],
             n_used = cz[["n"]], ev_used = cz[["ev"]],
             p_emp_A = emp_p(cz[["z"]], zA), p_emp_A_rds = x$p_emp_A,
             p_emp_B = emp_p(cz[["z"]], zB), p_emp_B_rds = x$p_emp_B,
             zA_p95 = unname(quantile(abs(zA), .95)), coexpr_qc_pass = isTRUE(x$coexpr_qc_pass),
             low_power = isTRUE(x$low_power))]
  row
}), fill = TRUE)
E <- R[evaluable == TRUE]
cat(sprintf("  可评估 %d / 不可评估 %d\n", nrow(E), nrow(R) - nrow(E)))
cat(sprintf("  |z_real 复算 − RDS| 最大 %.2e；|p_emp_A 复算 − RDS| 最大 %.2e；|p_emp_B 复算 − RDS| 最大 %.2e\n",
            max(abs(E$z_real - E$z_real_rds)), max(abs(E$p_emp_A - E$p_emp_A_rds)), max(abs(E$p_emp_B - E$p_emp_B_rds))))
check("V2.2", "可评估组合数", 71, nrow(E))
check("V2.3", "z_real 复算与 RDS 一致（全部组合，容差 1e-6）", 0, sum(abs(E$z_real - E$z_real_rds) > 1e-6))
check("V2.4", "p_emp_A 复算与 RDS 一致（全部组合）", 0, sum(abs(E$p_emp_A - E$p_emp_A_rds) > 1e-9))
check("V2.5", "p_emp_B 复算与 RDS 一致（全部组合）", 0, sum(abs(E$p_emp_B - E$p_emp_B_rds) > 1e-9))
check("V2.6", "覆盖率复算与 RDS 一致", 0, sum(abs(R$coverage - R$coverage_rds) > 1e-9))
check("V2.7", "零分布次数 != 1000 的组合数（SAP 要求 1000）", 0, sum(E$n_null != 1000, na.rm = TRUE))
if (any(E$n_null != 1000, na.rm = TRUE)) print(E[n_null != 1000, .(cohort, sig_id, n_null)], row.names = FALSE)

# ============================================================ 3. 头条数字
hr("3. 阶段 1 头条数字（由复算值现算）")
X <- fread(sort(list.files(DD, "^results_extended_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1])
X <- X[evaluable == TRUE]
cls_of <- function(z, p) fifelse(abs(z) < 1.96, "无信号", fifelse(p >= 0.05, "不特异", "强于随机"))
E[, sig_class := cls_of(z_real, p_emp_A)]; X[, sig_class := cls_of(z_real, p_emp_A)]
tabE <- E[, .N, by = sig_class]; tabX <- X[, .N, by = sig_class]
for (k in c("无信号", "不特异", "强于随机"))
  check(paste0("V3.", match(k, c("无信号", "不特异", "强于随机"))), paste("信号分类：", k),
        tabX[sig_class == k, N], tabE[sig_class == k, N])
check("V3.4", "弱于自身零分布 p95 的组合数", sum(abs(X$z_real) < X$zA_p95), sum(abs(E$z_real) < E$zA_p95))

P <- fread(sort(list.files(DD, "^results_pass_counts_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1])
agg <- E[, .(pass_A = sum(p_emp_A < 0.05), pass_B_main = sum(p_emp_B < 0.05 & coexpr_qc_pass),
             pass_B_sens = sum(p_emp_B < 0.05)), by = sig_id][order(sig_id)]
P <- P[order(sig_id)]
check("V3.5", "逐签名 pass_A（12 个）", P$pass_A, agg$pass_A)
check("V3.6", "逐签名 pass_B_main", P$pass_B_main, agg$pass_B_main)
check("V3.7", "逐签名 pass_B_sens", P$pass_B_sens, agg$pass_B_sens)

# ============================================================ 4. 队列角色与过拟合梯度
hr("4. 队列角色通过率与过拟合梯度")
fa <- file.path(DD, "train_cohort_annotation.csv")
if (!file.exists(fa)) stop("缺 ", fa, "：先运行 make_train_annotation.R（人工标注，来源为各原文 Methods）")
TR <- fread(fa, na.strings = c("", "NA"))
stopifnot(nrow(TR) == 12, !anyDuplicated(TR$sig_id))
norm_co <- function(x) sub("^TCGALUAD$", "TCGA", toupper(gsub("[^A-Za-z0-9]", "", x)))
used <- setNames(lapply(reg$overlap_cohorts, function(s)
  if (is.na(s) || s %in% c("-", "待核实")) character(0) else norm_co(strsplit(s, ",")[[1]])), reg$sig_id)
E <- merge(E, TR[, .(sig_id, train_cohort)], by = "sig_id", all.x = TRUE)
E[, co_n := norm_co(cohort)]
E[, used_in_original := mapply(function(s, c0) c0 %in% used[[s]], sig_id, co_n)]
E[, is_train := !is.na(train_cohort) & train_cohort != "待核实" & norm_co(train_cohort) == co_n]
E[, role := fifelse(is_train, "开发队列", fifelse(used_in_original, "原文验证队列", "完全独立队列"))]
E[!is.na(train_cohort) & train_cohort == "待核实" & used_in_original, role := "原文用过（训练/验证未分）"]
RL <- E[, .(n = .N, pass_A = sum(p_emp_A < 0.05)), by = role]; print(RL, row.names = FALSE)
Xr <- fread(sort(list.files(DD, "^results_by_cohort_role_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1])
for (rr in c("开发队列", "原文验证队列", "完全独立队列"))
  check(paste0("V4.", match(rr, c("开发队列", "原文验证队列", "完全独立队列"))), paste("通过/组合：", rr),
        paste0(Xr[role == rr, sum(p_emp_A < 0.05)], "/", Xr[role == rr, .N]),
        paste0(RL[role == rr, pass_A], "/", RL[role == rr, n]))
G <- E[!is.na(train_cohort) & train_cohort != "待核实",
       .(z_train = z_real[role == "开发队列"][1], z_ext_median = median(z_real[role != "开发队列"])), by = sig_id]
G[, ratio := abs(z_train) / abs(z_ext_median)]
Go <- fread(sort(list.files(DD, "^overfit_gradient_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1])
GG <- merge(Go[, .(sig_id, ratio_old = ratio)], G[, .(sig_id, ratio)], by = "sig_id")
print(GG[, .(sig_id, ratio_old, ratio = round(ratio, 2))], row.names = FALSE)
check("V4.4", "过拟合梯度（有训练队列者，保留两位）", round(GG$ratio_old[!is.na(GG$ratio_old)], 2),
      round(GG$ratio[!is.na(GG$ratio_old)], 2), tol = 0.005)

# ============================================================ 5. 打分覆盖审计
hr("5. 打分覆盖审计（每个 签名 × 队列）")
AU <- R[, .(cohort, sig_id, class, platform = MAN$platform[match(cohort, MAN$cohort)], n_sig, n_present,
            coverage = round(coverage, 3), n_in_pool, universe, pool_size, evaluable, missing)]
print(AU[coverage < 1][order(coverage)], row.names = FALSE)
check("V5.1", "打分宇宙 < 10,000 基因的组合数", 0, sum(AU$universe < 10000))
fwrite(AU, file.path(DD, sprintf("scoring_coverage_audit_%s.csv", STAMP)))

# ============================================================ 6. 落盘与汇总
hr("6. 汇总")
fwrite(E[, !c("co_n"), with = FALSE], file.path(DD, sprintf("verify_ch3_combo_%s.csv", STAMP)))
CK <- rbindlist(CHK); fwrite(CK, file.path(DD, sprintf("verify_ch3_summary_%s.csv", STAMP)))
cat(sprintf("\n  PASS %d / FAIL %d\n", sum(CK$status == "PASS"), sum(CK$status == "FAIL")))
if (any(CK$status == "FAIL")) print(CK[status == "FAIL"], row.names = FALSE)
cat("\n  输出：verify_ch3_combo / verify_ch3_summary / scoring_coverage_audit _", STAMP, ".csv\n", sep = "")
print(sessionInfo())
