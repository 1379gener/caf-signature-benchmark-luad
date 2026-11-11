#!/usr/bin/env Rscript
# ==============================================================================
# sci2_sens_sig012_origweights.R —— SIG-012 原文权重敏感性分析（运行前登记：DEVIATIONS.md 2026-10-25 D）
#
#   主分析：SIG-012 按 ≥10 基因规则归 A 类，等权 ssGSEA。本脚本只做事先登记的敏感性分析，不改动主结果。
#   做法（全部事先规定）：
#     - 打分基因 = SIG-012_genes.txt ∩ 候选池（精确符号，与主分析相同）；权重取 genes/SIG-012_weights_original.csv
#     - 分数 = 各基因队列内 Z 标准化后按原文权重加权求和（与 B 类规则相同；缺失基因不重新分配权重）
#     - 单变量 Cox（per SD），Wald z
#     - 零分布：SAP §5.3 第 1 条加权复合零分布 = 逐基因联合箱匹配（sample_null_A）+ 同一组权重，1000 次；
#       seed = make_seed(队列, "SIG-012", "W")（新类型码，不复用主分析 seed）
#     - 不做零分布 B（与 G 的处理一致：零分布 B 只用于等权 A 类）
#     - 报告：z_w、HR_w、p_emp_W、零分布 |z| 95% 分位，与主分析 z、p_emp_A 并列；两种分数的 Spearman 相关
#   输出：data_derived/sens_SIG012_origweights_<日期>.csv（已存在即停止）；日志 logs/
#   放置：ch04_benchmark/scripts/   用法：Rscript sci2_sens_sig012_origweights.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD <- file.path(ROOT, "03_signature_registry"); DD <- file.path(CH, "data_derived")
OUT <- file.path(DD, sprintf("sens_SIG012_origweights_%s.csv", format(Sys.Date(), "%Y%m%d")))
if (file.exists(OUT)) stop("输出已存在，不覆盖：", OUT)
LOG <- file.path(CH, "logs", format(Sys.time(), "sens_SIG012_origweights_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
suppressMessages({
  source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
  source(file.path(CH, "scripts/null_engine.R"))
})
LOCK <- require_locked_registry(REGD); reg <- LOCK$table
stopifnot(LOCK$md5 == "4ad3218931d95dfb6ee67a7ba07b9a75", reg[sig_id == "SIG-012", class] == "A")

fw <- file.path(REGD, "genes/SIG-012_weights_original.csv")
cat("权重文件 ", fw, "  md5 ", unname(tools::md5sum(fw)), "\n", sep = "")
stopifnot(unname(tools::md5sum(fw)) == "c91c10cce0c0cd73f7d5808292f79a9b")
W <- fread(fw); stopifnot(!anyDuplicated(W$gene), !anyNA(W$weight))
g <- trimws(readLines(file.path(REGD, "genes/SIG-012_genes.txt"))); g <- g[nzchar(g)]
stopifnot(setequal(g, W$gene), length(g) == 11L)
w_all <- setNames(as.numeric(W$weight), W$gene)[g]

coxfit <- function(sc, t, e) { ok <- !is.na(t) & t > 0 & !is.na(e) & is.finite(sc)
  s <- summary(coxph(Surv(t[ok], e[ok]) ~ as.numeric(scale(sc[ok]))))
  c(z = s$coefficients[1, "z"], hr = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4]) }
emp_p <- function(zr, zn) { zn <- zn[is.finite(zn)]; (1 + sum(abs(zn) >= abs(zr))) / (1 + length(zn)) }
cls_of <- function(z, p) ifelse(abs(z) < 1.96, "无信号", ifelse(p >= 0.05, "不特异", "强于随机"))

rds <- list.files(file.path(DD, "null_benchmark"), "\\.rds$", full.names = TRUE)
RM <- rbindlist(lapply(rds, function(f) { x <- readRDS(f)
  if (x$sig_id != "SIG-012" || !isTRUE(x$evaluable)) return(NULL)
  data.table(cohort = x$cohort, z_main = x$z_real, p_emp_A_main = x$p_emp_A) }))
stopifnot(!anyDuplicated(RM$cohort))

OUTL <- list()
for (a in COHORTS) {
  co <- load_cohort(a, verbose = FALSE)
  stopifnot(!anyDuplicated(colnames(co$expr)), ncol(co$expr) == co$n, co$n == as.integer(co$man$n_final[1]))
  if (!co$has_survival || !a %in% RM$cohort) next
  stopifnot(nrow(co$clin) == co$n, co$n_events == as.integer(co$man$n_events[1]))
  pool <- build_pool(co$expr); gp <- intersect(g, pool$genes); w <- w_all[gp]
  real <- as.numeric(score_sets_weighted(co$expr, list(gp), w)[1, ])
  main <- as.numeric(GSVA::gsva(co$expr, list(m = gp), method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)[1, ])
  cr <- coxfit(real, co$clin$OS_time, co$clin$OS_event)
  cm <- coxfit(main, co$clin$OS_time, co$clin$OS_event)
  stopifnot(abs(cm[["z"]] - RM[cohort == a, z_main]) < 1e-6)        # 主分析闸门
  sd <- make_seed(a, "SIG-012", "W")
  nl <- sample_null_A(pool, gp, n_null = 1000L, seed = sd)
  S  <- score_sets_weighted(co$expr, nl$sets, w)
  zn <- apply(S, 1, function(s) coxfit(s, co$clin$OS_time, co$clin$OS_event)[["z"]])
  p  <- emp_p(cr[["z"]], zn)
  OUTL[[a]] <- data.table(cohort = a, n = co$n, events = co$n_events, n_genes_scored = length(gp),
    genes_missing = paste(setdiff(g, gp), collapse = ";"), seed_W = sd, n_fallback_W = nl$n_fallback,
    z_main = cm[["z"]], p_emp_A_main = RM[cohort == a, p_emp_A_main],
    z_w = cr[["z"]], hr_w = cr[["hr"]], hr_w_lo = cr[["lo"]], hr_w_hi = cr[["hi"]],
    p_emp_W = p, zW_null_p95 = unname(quantile(abs(zn[is.finite(zn)]), 0.95)), n_null_finite = sum(is.finite(zn)),
    rho_scores = cor(real, main, method = "spearman"))
  cat(sprintf("  %-10s 基因 %2d  z_main %6.2f  z_w %6.2f  p_emp_W %.3f  ρ(两分数) %.2f\n",
              a, length(gp), cm[["z"]], cr[["z"]], p, OUTL[[a]]$rho_scores))
}
X <- rbindlist(OUTL)
X[, `:=`(class_main = cls_of(z_main, p_emp_A_main), class_w = cls_of(z_w, p_emp_W))]
print(X[, .(cohort, n_genes_scored, z_main = round(z_main, 2), p_emp_A_main, class_main,
            z_w = round(z_w, 2), p_emp_W, zW_null_p95 = round(zW_null_p95, 2), class_w, rho_scores = round(rho_scores, 2))], row.names = FALSE)
cat(sprintf("\n分类一致的队列：%d / %d\n", sum(X$class_main == X$class_w), nrow(X)))
fwrite(X, OUT); cat("写出 ", OUT, "\n", sep = "")
cat("\n"); print(sessionInfo())
