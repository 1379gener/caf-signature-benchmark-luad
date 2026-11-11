# sci2_review_summaries_20261109.R —— 审稿准备：事后描述性汇总（只读）
# 登记：DEVIATIONS.md 2026-11-09「审稿准备：事后描述性汇总」（运行之前）
suppressPackageStartupMessages(library(data.table))
setwd(path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"))
DD <- "data_derived"; OD <- file.path(DD, "review_summaries_20261109")
if (dir.exists(OD)) stop("output dir exists, not overwriting: ", OD)
dir.create(OD); dir.create("logs", showWarnings = FALSE)
sink(file.path("logs", "review_summaries_20261109.log"), split = TRUE)
cat("run:", format(Sys.time()), "\n")

bh <- fread(file.path(DD, "sap6_bh_combo_20261025.csv"))
sc <- fread(file.path(DD, "classA_scani_bh_20261025.csv"))
B  <- bh[class == "B"]
stopifnot(nrow(bh) == 71, nrow(B) == 41, nrow(sc) == 36, sum(sc$sig_id == "SCANI_myCAF") == 6,
          sum(B$pass_A_bh) == 8, sum(sc$pass_A_bh) == 1)

cat("\n===== 1. 类 B 的 8 个通过组合在零分布 B（S3）下 =====\n")
s3 <- fread(file.path(DD, "supp_tables_20261025/S3_classB_nullB.csv"))
stopifnot(nrow(s3) == 41)
m1 <- merge(B[pass_A_bh == TRUE, .(cohort, sig_id, role, z_real, zA_p95, p_emp_A, q_A)],
            s3[, .(cohort, sig_id, zB_p95, p_emp_B, coexpr_ratio, coexpr_qc_pass)],
            by = c("cohort", "sig_id"), all.x = TRUE)
print(m1); fwrite(m1, file.path(OD, "01_classB_passes_nullB.csv"))
cat("8 个中共表达 QC 通过:", sum(m1$coexpr_qc_pass), "；其中 p_emp_B < 0.05:",
    m1[coexpr_qc_pass == TRUE, sum(p_emp_B < 0.05)], "\n")
cat("全部 41 个类 B 组合中 QC 通过:", sum(s3$coexpr_qc_pass), "；coexpr_ratio < 0:", sum(s3$coexpr_ratio < 0), "\n")

cat("\n===== 2/3. 零分布 A：p95、随机成员 |z|>=1.96 比例、最小可检测 HR =====\n")
files <- c(list.files(file.path(DD, "null_benchmark"), "\\.rds$", full.names = TRUE),
           list.files(file.path(DD, "null_benchmark_scani_supp"), "\\.rds$", full.names = TRUE))
x0 <- readRDS(files[1]); cat("rds fields:\n"); print(sapply(x0, function(v) paste(class(v)[1], length(v))))
pick <- function(x) {
  nm <- names(x); ok <- vapply(x, function(v) is.numeric(v) && length(v) == 1000, logical(1))
  c1 <- nm[ok & grepl("z", nm, ignore.case = TRUE) & grepl("A", nm)]
  if (length(c1) != 1) stop("cannot identify null-A z vector; length-1000 numeric fields: ", paste(nm[ok], collapse = ", "))
  c1 }
nul <- rbindlist(lapply(files, function(f) {
  x <- readRDS(f); if (!isTRUE(x$evaluable)) return(NULL)
  fld <- pick(x); zn <- x[[fld]]
  data.table(cohort = x$cohort, sig_id = x$sig_id, class = x$class, n_events = x$n_events, field = fld,
             p95 = quantile(abs(zn), 0.95, names = FALSE), frac_sig = mean(abs(zn) >= 1.96))
}))
cat("null-A field used:", unique(nul$field), "；evaluable combos:", nrow(nul), "\n"); stopifnot(nrow(nul) == 77)
chk <- merge(nul, rbind(bh[, .(cohort, sig_id, zA_p95)], sc[sig_id == "SCANI_myCAF", .(cohort, sig_id, zA_p95)]),
             by = c("cohort", "sig_id"))
cat("p95 与结果文件 zA_p95 的最大差:", max(abs(chk$p95 - chk$zA_p95)), "\n")
cx <- fread(file.path(DD, "stage3_cox_20261024.csv"))
se_c <- cx[is.finite(z) & abs(z) > 0.5 & HR > 0, .(se_med = median(abs(log(HR) / z)), n_used = .N), by = cohort]
print(se_c)
mde <- merge(nul, se_c[, .(cohort, se_med)], by = "cohort", all.x = TRUE)
mde[, hr_thr := exp(p95 * se_med)]
summ <- mde[, .(events = n_events[1], se_med = se_med[1], p95_med = median(p95), p95_min = min(p95), p95_max = max(p95),
                hr_thr_med = median(hr_thr), hr_thr_approx = exp(median(p95) / sqrt(n_events[1])),
                frac_sig_med = median(frac_sig), frac_sig_min = min(frac_sig), frac_sig_max = max(frac_sig)),
            by = .(cohort, class)][order(class, -events)]
print(summ); fwrite(mde, file.path(OD, "02_nullA_p95_mde_by_combo.csv")); fwrite(summ, file.path(OD, "02_nullA_summary_by_cohort.csv"))
cat("随机成员 |z|>=1.96 比例，按类别（跨组合中位数 [范围]）：\n")
print(nul[, .(median = median(frac_sig), min = min(frac_sig), max = max(frac_sig)), by = class])

cat("\n===== 4. 全局 BH 敏感性 =====\n")
allp <- rbind(sc[, .(cohort, sig_id, class = "A", z_real, p = p_emp_A, pass_family = pass_A_bh)],
              B[, .(cohort, sig_id, class, z_real, p = p_emp_A, pass_family = pass_A_bh)])
stopifnot(nrow(allp) == 77)
allp[, `:=`(q_global = p.adjust(p, "BH"))][, q_class := p.adjust(p, "BH"), by = class][, q_sig := p.adjust(p, "BH"), by = sig_id]
sig <- abs(allp$z_real) >= 1.96
res4 <- allp[, .(family_pass = sum(pass_family), global_pass = sum(q_global < 0.05 & abs(z_real) >= 1.96),
                 class_pass = sum(q_class < 0.05 & abs(z_real) >= 1.96), bysig_pass = sum(q_sig < 0.05 & abs(z_real) >= 1.96)), by = class]
print(res4); fwrite(allp, file.path(OD, "04_global_bh.csv"))
cat("全局 BH 下通过的组合：\n"); print(allp[q_global < 0.05 & abs(z_real) >= 1.96, .(cohort, sig_id, class, p, q_global)])

cat("\n===== 5. 签名与 meta-PCNA 重叠 =====\n")
print(fread(file.path(DD, "m3_prolif_metaPCNA_20261017_sig_overlap.csv")))

cat("\n===== 6. TCGA：签名与增殖、成纤维细胞读数的 Spearman rho =====\n")
cr <- fread(file.path(DD, "stage3_cor_TCGA_LUAD_20261024.csv"))
want <- grep("Prolif|MCP_Fibroblasts|EPIC_CAF|xCell_Fibroblasts", names(cr), value = TRUE); cat("columns:", want, "\n")
r6 <- cr[grepl("^SIG-", readout), c("readout", want), with = FALSE]
print(r6, digits = 2); fwrite(r6, file.path(OD, "06_TCGA_rho_prolif_fibro.csv"))

cat("\n===== 7. 成纤维细胞读数单变量显著数；增殖指数 z =====\n")
ex <- fread(file.path(DD, "stage3_explore_fibro_prolif_20261025.csv"))
print(ex[, .(n_cohorts = .N, n_sig = sum(abs(z) >= 1.96), z_min = min(z), z_max = max(z)), by = .(readout, model)])
pz <- cx[grepl("Prolif", readout)]; print(pz[, .(cohort, readout, events, HR, z)])
cat("补充表中含 Prolif 的文件：\n")
for (f in list.files(file.path(DD, "supp_tables_20261025"), "\\.csv$", full.names = TRUE)) {
  n <- length(grep("Prolif", readLines(f, warn = FALSE))); if (n > 0) cat("  ", basename(f), ":", n, "lines\n") }
cat("\ndone:", format(Sys.time()), "\n"); sink()
