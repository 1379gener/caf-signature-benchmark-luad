#!/usr/bin/env Rscript
# ==============================================================================
# sci2_stage3_explore_fibro.R —— 事后探索性分析（2026-10-25 登记，见 DEVIATIONS.md）
#   问题：两个成纤维细胞读数预后方向相反；xCell_Fibroblasts 与增殖指数负相关。
#         校正增殖指数后，两者的 HR 如何变化？
#   模型（各队列，单变量 vs + 增殖）：Surv ~ scale(readout)；Surv ~ scale(readout) + scale(prolif)
#   读数只限预先指定的 2 个：MCP_Fibroblasts、xCell_Fibroblasts
#   数据：stage3_readouts_long_*.csv（读数）+ load_cohort()（生存）；有效样本定义同阶段 1–3
#   输出：data_derived/stage3_explore_fibro_prolif_<日期>.csv；日志写 logs/
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived")
suppressMessages(source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")))
STAMP <- format(Sys.Date(), "%Y%m%d")
LOG <- file.path(CH, "logs", format(Sys.time(), "stage3_explore_fibro_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)
FL <- sort(list.files(DD, "^stage3_readouts_long_[0-9]{8}\\.csv$", full.names = TRUE), decreasing = TRUE)[1]
cat("读数文件 ", basename(FL), "  md5 ", unname(tools::md5sum(FL)), "\n", sep = "")
L <- fread(FL); stopifnot(!anyDuplicated(L[, .(cohort, sample_id, readout)]))
RD <- c("MCP_Fibroblasts", "xCell_Fibroblasts")

fit <- function(d, f, model) {
  m <- coxph(f, data = d); s <- summary(m)
  data.table(model = model, n = m$n, events = m$nevent, HR = s$conf.int["r", 1], lo = s$conf.int["r", 3],
             hi = s$conf.int["r", 4], z = s$coefficients["r", "z"], p = s$coefficients["r", "Pr(>|z|)"]) }
OUT <- list()
for (a in sort(unique(L$cohort))) {
  co <- load_cohort(a, verbose = FALSE); cl <- co$clin
  stopifnot(!anyDuplicated(cl[[co$id_col]]), nrow(cl) == co$n)
  S <- data.table(sample_id = as.character(cl[[co$id_col]]), t = cl$OS_time, e = cl$OS_event)
  W <- dcast(L[cohort == a & readout %in% c(RD, "Prolif_metaPCNA")], sample_id ~ readout, value.var = "value")
  stopifnot(nrow(W) == co$n)
  D <- merge(S, W, by = "sample_id"); stopifnot(nrow(D) == co$n)
  D <- D[!is.na(t) & t > 0 & !is.na(e)]
  for (rv in RD) {
    d <- D[is.finite(D[[rv]]) & is.finite(Prolif_metaPCNA)]
    d[, r := as.numeric(scale(d[[rv]]))][, pr := as.numeric(scale(Prolif_metaPCNA))]
    rp <- cor(d$r, d$pr, method = "spearman")
    OUT[[length(OUT) + 1]] <- rbind(fit(d, Surv(t, e) ~ r, "univariable"), fit(d, Surv(t, e) ~ r + pr, "+prolif"))[,
      `:=`(cohort = a, readout = rv, rho_prolif = rp)]
  }
}
R <- rbindlist(OUT); setcolorder(R, c("cohort", "readout", "model"))
fo <- file.path(DD, sprintf("stage3_explore_fibro_prolif_%s.csv", STAMP)); if (file.exists(fo)) stop("已存在，不覆盖：", fo)
fwrite(R, fo)
options(width = 200, datatable.print.nrows = 100)
print(R[, .(cohort, readout, model, n, events, HR = round(HR, 2), CI = sprintf("%.2f–%.2f", lo, hi), z = round(z, 2),
            p = signif(p, 2), rho_prolif = round(rho_prolif, 2))][order(readout, cohort, model != "univariable")])
cat("\n汇总（逐读数）：\n")
print(dcast(R, readout + cohort ~ model, value.var = "z")[, .(n_cohort = .N,
  uni_dir = paste0(sum(univariable > 0), "+/", sum(univariable < 0), "−"), adj_dir = paste0(sum(`+prolif` > 0), "+/", sum(`+prolif` < 0), "−"),
  uni_sig = sum(abs(univariable) >= 1.96), adj_sig = sum(abs(`+prolif`) >= 1.96),
  z_uni_median = round(median(univariable), 2), z_adj_median = round(median(`+prolif`), 2)), by = readout])
cat("\n写出 ", fo, "\n", sep = "")
sink(); close(zz)
