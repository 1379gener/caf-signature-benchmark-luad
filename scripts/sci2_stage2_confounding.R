#!/usr/bin/env Rscript
# ==============================================================================
# sci2_stage2_confounding.R —— 阶段 2：混杂分解（SAP v1.2 §7，M1–M4 嵌套 Cox）
#
#   M1 = score；M2 = M1 + 纯度；M3 = M1 + 增殖；M4 = M1 + 纯度 + 增殖 + stage + age + sex
#   - 签名打分与阶段 1 完全同一路径：候选池 build_pool()，gp = intersect(签名基因, 池)，
#     A 类 score_sets_ssgsea()，B 类 score_sets_weighted()（null_engine.R）
#   - 闸门：全样本 M1 的 z 必须与阶段 1 RDS 的 z_real 逐个一致，否则停止
#   - 纯度只有 TCGA（Aran 2015 ABSOLUTE）；GEO 队列 M2/M4 按 SAP §7.2 标记为不可计算
#   - stage：二分类 stage_ge_II（2026-10-17 登记）；TCGA 另跑分期补充版 M4（敏感性）
#   - 衰减幅度在同一完整病例样本上计算：分母为在该样本上重新拟合的 M1（2026-10-17 登记）
#   - 分数与连续协变量在各模型自身的分析样本内 Z 标准化（per SD）
#
#   输入：data_derived/null_benchmark/*.rds（取可评估组合与 z_real 闸门）、
#         m2_purity_TCGA_*.csv、m3_prolif_metaPCNA_*.csv、m4_stage_binary_*.csv、m4_stage_binary_TCGA_sens_*.csv
#   输出：data_derived/stage2_models_<日期>.csv（逐模型）、stage2_attenuation_<日期>.csv（逐组合）；日志写 logs/
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_stage2_confounding.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD <- file.path(ROOT, "03_signature_registry"); DD <- file.path(CH, "data_derived")
suppressMessages({ source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
                   source(file.path(CH, "scripts/null_engine.R")) })
STAMP <- format(Sys.Date(), "%Y%m%d")
LOG <- file.path(CH, "logs", format(Sys.time(), "stage2_confounding_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
latest <- function(pat) { f <- sort(list.files(DD, pat, full.names = TRUE), decreasing = TRUE); stopifnot(length(f) >= 1); f[1] }
rd <- function(pat) { f <- latest(pat); cat("  读入 ", basename(f), "  md5 ", unname(tools::md5sum(f)), "\n", sep = "")
                      fread(f, na.strings = c("", "NA")) }

hr("0. 锁定物与输入")
LOCK <- require_locked_registry(REGD); stopifnot(LOCK$md5 == "4ad3218931d95dfb6ee67a7ba07b9a75")
PUR  <- rd("^m2_purity_TCGA_[0-9]{8}\\.csv$")
PRO  <- rd("^m3_prolif_metaPCNA_[0-9]{8}\\.csv$")
STG  <- rd("^m4_stage_binary_[0-9]{8}\\.csv$")
STS  <- rd("^m4_stage_binary_TCGA_sens_[0-9]{8}\\.csv$")
for (x in list(PUR, PRO, STG, STS)) stopifnot(!anyDuplicated(x[, .(cohort, sample_id)]))

read_sig <- function(sid, cls) {                         # 与 run_benchmark.R L73–86 相同
  g <- trimws(readLines(file.path(REGD, "genes", paste0(sid, "_genes.txt")))); g <- g[nzchar(g)]
  w <- NULL
  if (cls == "B") { d <- fread(file.path(REGD, "genes", paste0(sid, "_weights.csv")))
                    w <- setNames(as.numeric(d$weight), d$gene)[g]; stopifnot(!anyNA(w)) }
  list(genes = g, weights = w) }
zs <- function(x) as.numeric(scale(x))

#' 在给定样本上拟合一个 Cox 模型，返回 score 项的结果
fit_one <- function(d, covars, model, sample_name) {
  d <- copy(d); d[, s := zs(score)]
  for (v in intersect(covars, c("purity", "prolif", "age"))) set(d, j = v, value = zs(d[[v]]))
  rhs <- paste(c("s", covars), collapse = " + ")
  w <- character(0)
  m <- withCallingHandlers(coxph(as.formula(paste("Surv(t, e) ~", rhs)), data = d),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  sm <- summary(m); cf <- sm$coefficients["s", ]; ci <- sm$conf.int["s", ]
  data.table(model = model, sample = sample_name, n = m$n, events = m$nevent,
             logHR = unname(cf["coef"]), se = unname(cf["se(coef)"]), z = unname(cf["z"]), p = unname(cf["Pr(>|z|)"]),
             HR = unname(ci["exp(coef)"]), lo = unname(ci["lower .95"]), hi = unname(ci["upper .95"]),
             covariates = paste(covars, collapse = "+"), warnings = paste(unique(w), collapse = " | "))
}

hr("1. 逐组合拟合")
rds <- list.files(file.path(DD, "null_benchmark"), "\\.rds$", full.names = TRUE)
RES <- list(); CO <- list(); POOL <- list(); GATE <- list()
for (f in rds) {
  x <- readRDS(f); if (!isTRUE(x$evaluable)) next
  a <- x$cohort; sid <- x$sig_id; cls <- x$class
  if (is.null(CO[[a]])) {
    co <- load_cohort(a, verbose = FALSE); stopifnot(!anyDuplicated(colnames(co$expr)), ncol(co$expr) == co$n)
    CO[[a]] <- co; POOL[[a]] <- build_pool(co$expr) }
  co <- CO[[a]]; sg <- read_sig(sid, cls)
  gp <- intersect(sg$genes, POOL[[a]]$genes)
  sc <- if (cls == "A") score_sets_ssgsea(co$expr, list(gp), prefix = "real")[1, ] else
                        score_sets_weighted(co$expr, list(gp), sg$weights[gp])[1, ]
  id <- colnames(co$expr)
  D <- data.table(sample_id = id, score = as.numeric(sc), t = co$clin$OS_time, e = co$clin$OS_event,
                  age = suppressWarnings(as.numeric(co$clin$age)), sex = factor(co$clin$sex))
  D <- merge(D, PRO[cohort == a, .(sample_id, prolif = prolif_metaPCNA)], by = "sample_id", all.x = TRUE)
  D <- merge(D, STG[cohort == a, .(sample_id, stage = stage_ge_II)], by = "sample_id", all.x = TRUE)
  if (a == "TCGA_LUAD") {
    D <- merge(D, PUR[, .(sample_id, purity = purity_absolute)], by = "sample_id", all.x = TRUE)
    D <- merge(D, STS[, .(sample_id, stage_s = stage_ge_II_sens)], by = "sample_id", all.x = TRUE)
  } else D[, `:=`(purity = NA_real_, stage_s = NA_integer_)]
  stopifnot(nrow(D) == co$n)
  D <- D[!is.na(t) & t > 0 & !is.na(e) & is.finite(score)]          # 与 cox_one() 同一有效样本定义

  out <- list(fit_one(D, character(0), "M1", "full"))
  GATE[[length(GATE) + 1]] <- data.table(cohort = a, sig_id = sid, z_stage1 = x$z_real, z_now = out[[1]]$z)
  S3 <- D[!is.na(prolif)]
  out <- c(out, list(fit_one(S3, character(0), "M1", "S3"), fit_one(S3, "prolif", "M3", "S3")))
  if (a == "TCGA_LUAD") {
    S2 <- D[!is.na(purity)]
    S4 <- D[!is.na(purity) & !is.na(prolif) & !is.na(stage) & !is.na(age) & !is.na(sex)]
    S4s <- D[!is.na(purity) & !is.na(prolif) & !is.na(stage_s) & !is.na(age) & !is.na(sex)]
    S4s[, stage := stage_s]
    c4 <- c("purity", "prolif", "stage", "age", "sex")
    out <- c(out, list(fit_one(S2, character(0), "M1", "S2"), fit_one(S2, "purity", "M2", "S2"),
                       fit_one(S4, character(0), "M1", "S4"), fit_one(S4, c4, "M4", "S4"),
                       fit_one(S4s, character(0), "M1", "S4_sens"), fit_one(S4s, c4, "M4", "S4_sens")))
  }
  R <- rbindlist(out); R[, `:=`(cohort = a, sig_id = sid, class = cls, n_genes_scored = length(gp))]
  RES[[length(RES) + 1]] <- R
  cat(sprintf("  %-9s %-7s %s  M1 z %+.2f%s\n", a, sid, cls, out[[1]]$z, ifelse(any(nzchar(R$warnings)), "  [有警告]", "")))
}
M <- rbindlist(RES); setcolorder(M, c("cohort", "sig_id", "class", "n_genes_scored"))

hr("2. 闸门：全样本 M1 与阶段 1 的 z_real")
GT <- rbindlist(GATE); GT[, diff := abs(z_now - z_stage1)]
cat(sprintf("  组合 %d；最大 |差| %.2e\n", nrow(GT), max(GT$diff)))
stopifnot(nrow(GT) == 71, max(GT$diff) < 1e-6)

hr("3. 衰减幅度（同样本）")
att <- function(m1, mk) (m1 - mk) / m1
A <- M[, {
  g <- function(mo, sa) .SD[model == mo & sample == sa]
  r <- list()
  for (p in list(c("M2", "S2"), c("M3", "S3"), c("M4", "S4"), c("M4", "S4_sens"))) {
    b <- g("M1", p[2]); k <- g(p[1], p[2])
    if (nrow(b) && nrow(k)) r[[length(r) + 1]] <- data.table(comparison = paste0(p[1], "@", p[2]), n = k$n, events = k$events,
      z_M1 = b$z, HR_M1 = b$HR, HR_Mk = k$HR, p_Mk = k$p, attenuation = att(b$logHR, k$logHR),
      M1_signal = abs(b$z) >= 1.96)
  }
  rbindlist(r) }, by = .(cohort, sig_id, class)]
A[M1_signal == FALSE, attenuation_reported := NA_real_][M1_signal == TRUE, attenuation_reported := attenuation]
cat("\n  逐组合（attenuation_reported 仅在同样本 M1 |z| ≥ 1.96 时给出）：\n")
print(A[, .(cohort, sig_id, comparison, n, events, z_M1 = round(z_M1, 2), HR_M1 = round(HR_M1, 2),
            HR_Mk = round(HR_Mk, 2), att = round(attenuation_reported, 2))], row.names = FALSE)
cat("\n  M3（增殖）衰减：按队列汇总（仅 M1 有信号者）\n")
print(A[comparison == "M3@S3" & M1_signal == TRUE,
        .(n_combo = .N, att_median = round(median(attenuation), 2), att_min = round(min(attenuation), 2),
          att_max = round(max(attenuation), 2)), by = cohort], row.names = FALSE)
cat("\n  有拟合警告的模型：", M[nzchar(warnings), .N], "\n")
if (M[nzchar(warnings), .N]) print(M[nzchar(warnings), .(cohort, sig_id, model, sample, warnings)], row.names = FALSE)

hr("4. 写出")
f1 <- file.path(DD, sprintf("stage2_models_%s.csv", STAMP)); f2 <- file.path(DD, sprintf("stage2_attenuation_%s.csv", STAMP))
if (file.exists(f1) || file.exists(f2)) stop("已存在，不覆盖：", f1)
fwrite(M, f1); fwrite(A, f2)
cat("  ", f1, "\n  ", f2, "\n", sep = "")
print(sessionInfo())
sink(); close(zz)
