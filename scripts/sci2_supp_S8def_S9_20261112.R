# sci2_supp_S8def_S9_20261112.R —— 由 11-11 三项事后分析的输出生成补充表 S8d/S8e/S8f、S9（参数），S1 增加 Nos. 35–37
# 登记：DEVIATIONS.md 2026-11-12「补充表 S8d–f、S9」（运行之前）。输入 supp_tables_20261111（不覆盖），输出 supp_tables_20261112
source("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_common_20261111.R")
suppressMessages(library(openxlsx))
inF  <- "data_derived/supp_tables_20261111/SCI-2_Supplementary_Tables_S1-S8_EN.xlsx"
outD <- "data_derived/supp_tables_20261112"; outF <- file.path(outD, "SCI-2_Supplementary_Tables_S1-S9_EN.xlsx")
stopifnot(file.exists(inF), !dir.exists(outD)); dir.create(outD)
sink("logs/supp_S8def_S9_20261112.log", split = TRUE); cat("run:", format(Sys.time()), "\n")
P  <- fread("data_derived/posctrl_20261111/posctrl_metaPCNA_20261111.csv")
MB <- fread("data_derived/mc_redraw_20261111/mc_redraw_borderline_20261111.csv")
MD <- fread("data_derived/mc_redraw_20261111/mc_redraw_decisions_20261111.csv")
MS <- fread("data_derived/mc_redraw_20261111/mc_se_all_20261111.csv")
G  <- fread("data_derived/gse68465_site_20261111/gse68465_site_strat_20261111.csv")
stopifnot(nrow(P) == 6, nrow(MB) == 16, nrow(MD) == 77, nrow(MS) == 77, nrow(G) == 12, all(G$null_reproduced))

## ---- 主分析判定（与 mc 脚本相同的 BH 族），及 GSE68465 分层后的判定
M <- read_main(); stopifnot(nrow(M) == 77)
bh <- function(D) { D <- copy(D); D[, q_A := p.adjust(p_emp_A, "BH"), by = .(cohort, class)]; D }
Q0 <- bh(M); cat("主分析通过（q_A < 0.05）：\n"); print(Q0[q_A < 0.05, .(cohort, sig_id, class, p_emp_A, q_A = round(q_A, 3))])
M2 <- copy(M); M2[G, on = .(cohort, sig_id), p_emp_A := i.p_emp_A_site_strat]
Q2 <- bh(M2)[cohort == "GSE68465"]; Q0g <- Q0[cohort == "GSE68465"]
site_cmp <- merge(Q0g[, .(sig_id, q_A_unstrat = q_A)], Q2[, .(sig_id, q_A_site_strat = q_A)], by = "sig_id")
n_pass_u <- sum(site_cmp$q_A_unstrat < 0.05); n_pass_s <- sum(site_cmp$q_A_site_strat < 0.05)
cat(sprintf("GSE68465 通过数（FDR）：未分层 %d，按中心分层 %d\n", n_pass_u, n_pass_s))

## ---- 摘要数字（供 S1 与正文）
dif <- MB$p_10000 - MB$p_1000
flip_raw <- MB[(p_1000 < 0.05) != (p_10000 < 0.05)]
n_dec_change <- MD[pass_A_1000 != pass_A_mix | pass_B_1000 != pass_B_mix, .N]
cat(sprintf("阳性对照：z %.2f–%.2f；零分布 A 通过 %d/6（p_emp < 0.05；q < 0.05 为 %d/6）；零分布 B 通过 %d/6（QC 通过 %d/6），p_emp_B %.3f–%.3f\n",
  min(P$z_real), max(P$z_real), sum(P$p_emp_A < 0.05), sum(P$q_A < 0.05), sum(P$qc_pass & P$p_emp_B < 0.05), sum(P$qc_pass), min(P$p_emp_B), max(P$p_emp_B)))
cat(sprintf("蒙特卡洛：77 组合 SE_A 中位 %.4f、最大 %.4f；边界 %d 组合，p 变化 %+.3f 至 %+.3f；未校正 0.05 判定改变 %d（%s）；FDR 判定改变 %d\n",
  median(MS$se_A), max(MS$se_A), nrow(MB), min(dif), max(dif), nrow(flip_raw), paste(flip_raw[, sprintf("%s %s %.3f→%.3f", cohort, sig_id, p_1000, p_10000)], collapse = "; "), n_dec_change))
cat(sprintf("GSE68465 分层：|z| ≥ 1.96 未分层 %d、分层 %d；p_emp_A < 0.05 未分层 %d（%s）、分层 %d（%s）；PH p < 0.05 未分层 %d、分层 %d；FDR 通过 %d → %d\n",
  sum(abs(G$z_unstrat) >= 1.96), sum(abs(G$z_site_strat) >= 1.96),
  sum(G$p_emp_A_unstrat < 0.05), paste(G[p_emp_A_unstrat < 0.05, sig_id], collapse = ", "),
  sum(G$p_emp_A_site_strat < 0.05), paste(G[p_emp_A_site_strat < 0.05, sig_id], collapse = ", "),
  sum(G$ph_p_unstrat < 0.05), sum(G$ph_p_site_strat < 0.05), n_pass_u, n_pass_s))

## ---- 工作簿
wb <- loadWorkbook(inF); hs <- createStyle(textDecoration = "bold")
add <- function(name, d) { addWorksheet(wb, name); writeData(wb, name, d, headerStyle = hs); setColWidths(wb, name, seq_along(d), "auto") }
r4 <- function(x) round(x, 4)
add("S8d_positive_control", P[, .(cohort, n_genes, n_genes_scored = n_used, z = r4(z_real), hr_per_sd = r4(hr_per_sd),
  abs_z_null_A_p95 = r4(zA_p95), p_emp_A = r4(p_emp_A), q_A_within_control = r4(q_A), abs_z_null_B_p95 = r4(zB_p95), p_emp_B = r4(p_emp_B),
  coexpr_ratio_B = r4(coexpr_ratio), qc_pass_B = qc_pass, n_relaxed_B = n_relax_B, seed_A, seed_B)])
add("S8e_monte_carlo", MB[, .(cohort, sig_id, class, null_distribution = null, p_emp_1000 = r4(p_1000), se_1000 = r4(sqrt(p_1000 * (1 - p_1000) / 1000)),
  p_emp_10000 = r4(p_10000), se_10000 = r4(se_10000), seed_10000, minutes)])
add("S8f_GSE68465_site", merge(G, site_cmp, by = "sig_id")[, .(sig_id, class, z_unstratified = r4(z_unstrat), p_emp_A_unstratified = r4(p_emp_A_unstrat),
  q_A_unstratified = r4(q_A_unstrat), z_site_stratified = r4(z_site_strat), p_emp_A_site_stratified = r4(p_emp_A_site_strat), q_A_site_stratified = r4(q_A_site_strat),
  abs_z_null_A_p95_stratified = r4(zA_p95_strat), ph_p_unstratified = signif(ph_p_unstrat, 3), ph_p_site_stratified = signif(ph_p_site_strat, 3), null_reproduced)])

## ---- S9 参数（取自代码常量、包版本与数据元信息）
se_rel <- tryCatch(S4Vectors::metadata(readRDS(path.expand("~/caf_spatial_project/data/TCGA_LUAD/TCGA_LUAD_SE.rds")))$data_release, error = function(e) NA_character_)
pv <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
S9 <- rbindlist(list(
 list("Software", "R", R.version$version.string, "sessionInfo"),
 list("Software", "Bioconductor", tryCatch(as.character(BiocManager::version()), error = function(e) pv("BiocVersion")), "BiocManager::version()"),
 list("Software", "Packages", paste(sprintf("%s %s", c("GSVA","survival","data.table","org.Hs.eg.db","MCPcounter","xCell","EPIC","estimate","dendextend","ggplot2","hgu133a.db","hgu133plus2.db","SummarizedExperiment","GEOquery"),
        vapply(c("GSVA","survival","data.table","org.Hs.eg.db","MCPcounter","xCell","EPIC","estimate","dendextend","ggplot2","hgu133a.db","hgu133plus2.db","SummarizedExperiment","GEOquery"), pv, "")), collapse = "; "), "packageVersion"),
 list("Random numbers", "RNG kind", paste(RNGkind(), collapse = ", "), "RNGkind()"),
 list("Random numbers", "Seed", sprintf("%d + (sum over i of i x UTF-8 code of character i of \"cohort|signature|null type\") mod 100000", ENGINE$seed_base), "null_engine.R make_seed"),
 list("Data", "TCGA-LUAD release", se_rel, "metadata of TCGA_LUAD_SE.rds"),
 list("Data", "TCGA-LUAD assay; OS", "tpm_unstrand, log2(TPM + 1); vital_status, days_to_death, days_to_last_follow_up", "03_baseA_tcga.R"),
 list("Data", "Linear-scale rule", "maximum >= 25 -> log2(x + 1)", "05_baseA_geo_rev3.R"),
 list("Null A", "Draws", as.character(ENGINE$n_null), "ENGINE$n_null"),
 list("Null A", "Candidate pool", sprintf("expressed (> 0) in >= %.0f%% of samples, variance > 0", 100 * ENGINE$detect_frac), "build_pool"),
 list("Null A", "Bins", sprintf("%d equal-frequency bins of mean log2 expression x %d of variance, within cohort", ENGINE$n_bins, ENGINE$n_bins), "build_pool"),
 list("Null B", "Draws", as.character(ENGINE$n_null), "ENGINE$n_null"),
 list("Null B", "Module", "seed gene from the joint bin of a random signature gene + k - 1 genes with highest Spearman correlation; signature genes excluded", "sample_null_B"),
 list("Null B", "Matching", sprintf("module mean expression and mean variance within +/-%.0f%% of the signature's; %d draws, then the next draw checked at +/-%.0f%% and kept regardless", 100 * ENGINE$match_tol, ENGINE$max_resample, 100 * ENGINE$match_tol_relax), "sample_null_B"),
 list("Null B", "Co-expression QC", sprintf("median over the first 200 modules of the mean pairwise Spearman correlation >= %.2f x the signature's", ENGINE$coexpr_min_ratio), "sample_null_B"),
 list("Class A score", "ssGSEA", sprintf("GSVA::gsva(method = \"ssgsea\", ssgsea.norm = FALSE); other arguments at the GSVA %s defaults (tau = 0.25, min.sz = 1, max.sz = Inf)", pv("GSVA")), "score_sets_ssgsea"),
 list("Class B score", "Weighted sum", "sum of published weight x within-cohort z score of gene; missing genes not reweighted", "score_sets_weighted"),
 list("Survival model", "Cox", "coxph(Surv(time, event) ~ scale(score)), Efron ties; Wald z", "cox_one, cox_z_batch"),
 list("Empirical p", "Formula", "(1 + #{|z_null| >= |z_real|}) / (1 + number of finite z_null), two-sided", "empirical_p"),
 list("FDR", "Benjamini-Hochberg", "within cohort x class (null A); within cohort for class A (null B)", "sci2_sap6_bh.R"),
 list("Readouts", "MCP-counter", "MCPcounter.estimate(featuresType = \"HUGO_symbols\")", "sci2_stage3_redundancy.R"),
 list("Readouts", "xCell", "xCellAnalysis(rnaseq = TRUE for TCGA-LUAD, FALSE otherwise; parallel.sz = 1)", "sci2_stage3_redundancy.R"),
 list("Readouts", "EPIC", "EPIC(bulk = TPM), default reference (TRef); TCGA-LUAD only", "sci2_stage3_redundancy.R"),
 list("Readouts", "ESTIMATE immune", "ImmuneSignature genes of estimate::SI_geneset scored by ssGSEA as a class A gene set", "sci2_stage3_redundancy.R"),
 list("Readouts", "Clustering", "hclust(as.dist(1 - |rho|), method = \"ward.D2\"); Spearman rho", "sci2_stage3_redundancy.R"),
 list("Confounding", "VIF", "1 / (1 - r^2), r = Pearson correlation of score and proliferation index", "sci2_stage2_confounding.R, sci2_stage2_diagnostics.R"),
 list("Proportional hazards", "Test", "cox.zph(transform = \"km\"); split at 3 years with survSplit", "sci2_ph_check_20261110.R"),
 list("Weights", "SIG-001", "recovered by least squares from the risk scores in the authors' public data (R^2 = 1; the published text gives CLEC3B with opposite signs)", "registry_notes.csv"),
 list("Weights", "SIG-004, -006, -007, -008, -011", "coefficients from the published formulas", "registry_notes.csv"),
 list("Weights", "SIG-002", "equal weights (MCP-counter fibroblast genes)", "registry_notes.csv")
))
setnames(S9, c("area", "item", "value", "source_in_code"))
stopifnot(!anyNA(S9$value)); add("S9_parameters", S9)

## ---- S1 Nos. 35–37
s1 <- read.xlsx(inF, sheet = "S1_deviation_log", check.names = FALSE); stopifnot(nrow(s1) == 34)
mk <- function(no, item, desc, eff) { r <- s1[1, ]; r[1, ] <- NA; r[1, 1] <- no; r[1, 2] <- "2026-11-11"; r[1, 3] <- item; r[1, 4] <- desc
  r[1, 5] <- "After"; r[1, 6] <- "Post hoc analysis registered before it was run; reported in Table S8; main analysis unchanged."; r[1, 7] <- eff; r }
new <- rbind(
 mk(35L, "Positive control (not in the plan)", "The 131 meta-PCNA proliferation genes were evaluated as a class A gene set against null distributions A and B (1,000 draws each).",
    sprintf("Exceeded null A in %d of 6 cohorts and null B in %d of 6 (Table S8d)." , sum(P$p_emp_A < 0.05), sum(P$qc_pass & P$p_emp_B < 0.05))),
 mk(36L, "Monte Carlo error (not in the plan)", sprintf("%d combinations with 0.02 < p_emp < 0.10 were re-evaluated with 10,000 new random gene sets.", nrow(MB)),
    sprintf("p_emp changed by %+.3f to %+.3f; no decision after false-discovery-rate control changed (Table S8e).", min(dif), max(dif))),
 mk(37L, "GSE68465 study site (not in the plan)", "Cox models in GSE68465 were stratified by contributing site (four sites), for the signature and for its 1,000 null-A gene sets.",
    sprintf("Passes after false-discovery-rate control in GSE68465: %d unstratified, %d stratified (Table S8f).", n_pass_u, n_pass_s)))
writeData(wb, "S1_deviation_log", new, startRow = nrow(s1) + 2, colNames = FALSE)
saveWorkbook(wb, outF)

## ---- 核对：旧 sheet 未变（S1 除外），新 sheet 无中文
for (n in getSheetNames(inF)) { a <- read.xlsx(inF, n, check.names = FALSE); b <- read.xlsx(outF, n, check.names = FALSE)
  cat(sprintf("%-28s %s\n", n, if (identical(a, b)) "identical" else sprintf("CHANGED (rows %d -> %d)", nrow(a), nrow(b)))) }
for (n in setdiff(getSheetNames(outF), getSheetNames(inF))) { b <- read.xlsx(outF, n); cat(sprintf("%-28s NEW %d rows; CJK: %s\n", n, nrow(b), any(grepl("[一-鿿]", unlist(b))))) }
cat("md5", tools::md5sum(outF), "\n"); cat("done:", format(Sys.time()), "\n"); sink()
