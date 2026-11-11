# sci2_supp_S8_20261109.R —— 在定稿 S1–S7 EN 工作簿基础上：S1 追加 2 条、S6 追加增殖指数、新增 S8a/S8b（只读输入，输出新目录）
# 登记：DEVIATIONS.md 2026-11-09「补充表新增 Table S8；Table S6 补增殖指数」（运行之前）
suppressPackageStartupMessages({ library(openxlsx); library(data.table) })
setwd(path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"))
DD <- "data_derived"
FI <- file.path(DD, "supp_tables_20261025", "SCI-2_Supplementary_Tables_S1-S7_EN.xlsx")
OD <- file.path(DD, "supp_tables_20261109"); if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
FO <- file.path(OD, "SCI-2_Supplementary_Tables_S1-S8_EN.xlsx")
sink(file.path("logs", "supp_S8_20261109.log"), split = TRUE)
cat("input md5:", unname(tools::md5sum(FI)), "\n")
rd <- function(s) { d <- as.data.table(read.xlsx(FI, s, check.names = FALSE, sep.names = " ")); 
  if (any(vapply(d, function(v) any(grepl("&#[0-9]+;", v)), logical(1)))) stop("HTML entity in sheet ", s); d }
sh <- getSheetNames(FI); IN <- setNames(lapply(sh, rd), sh); OUT <- copy(IN)

# S1：追加 2 条
s1 <- IN$S1_deviation_log; stopifnot(nrow(s1) == 28, max(s1$No.) == 28)
add1 <- data.table(
  No. = 29:30, `Date recorded` = "2026-11-09",
  `Plan section / item` = c("Post hoc summaries (review preparation)", "Supplementary tables"),
  Description = c("Descriptive summaries of existing results: null distribution B for the eight class B passes; 95th percentile of null-A |z|, proportion of random gene sets with |z| >= 1.96 and approximate smallest detectable HR per SD; BH adjustment across all combinations, within class and within signature; overlap of signatures with meta-PCNA genes; correlations with proliferation and fibroblast estimates in TCGA-LUAD; counts of significant fibroblast estimates.",
                  "Table S8 added; proliferation index added to Table S6."),
  `Found before or after results` = "After; registered before running",
  Action = c("Computed from existing result files; no null distribution regenerated; labelled post hoc (Table S8).",
             "Tables S2-S7 otherwise exported unchanged."),
  `Effect on reported results` = c("No reported decision changed. Seven of the eight class B passes passed the co-expression check and six of these also exceeded null B. Passes under alternative FDR families: all combinations 4; within class 0 (A) and 6 (B); within signature 1 (A) and 10 (B).", "—"))
OUT$S1_deviation_log <- rbind(s1, add1, use.names = TRUE)

# S6：追加增殖指数 6 行
s6 <- IN$S6_fibroblast_proliferation; stopifnot(nrow(s6) == 24)
cx <- fread(file.path(DD, "stage3_cox_20261024.csv"))[readout == "Prolif_metaPCNA"]; stopifnot(nrow(cx) == 6)
add6 <- cx[, .(cohort, readout, model = "univariable", n, events, HR, lo, hi, z, p, rho_prolif = NA_real_)]
OUT$S6_fibroblast_proliferation <- rbind(s6, add6, use.names = TRUE)

# S8a / S8b
RS <- file.path(DD, "review_summaries_20261109")
s8a <- fread(file.path(RS, "02_nullA_summary_by_cohort.csv")); stopifnot(nrow(s8a) == 12)
s8b <- fread(file.path(RS, "04_global_bh.csv")); stopifnot(nrow(s8b) == 77)
setnames(s8b, c("p", "pass_family"), c("p_emp_A", "pass_main"))
s8b[, `:=`(pass_global = q_global < 0.05 & abs(z_real) >= 1.96, pass_class = q_class < 0.05 & abs(z_real) >= 1.96,
           pass_signature = q_sig < 0.05 & abs(z_real) >= 1.96)]
stopifnot(sum(s8b$pass_main) == 9, sum(s8b$pass_global) == 4)
OUT$S8a_null_summary_by_cohort <- s8a; OUT$S8b_alternative_FDR <- s8b
fwrite(s8a, file.path(OD, "S8a_null_summary_by_cohort.csv")); fwrite(s8b, file.path(OD, "S8b_alternative_FDR.csv"))

# README
rm <- IN$README; stopifnot("Table S7" %in% rm$Item)
rm[Item == "Table S6", Content := "MCP-counter and xCell fibroblast estimates with and without adjustment for the proliferation index, and the proliferation index itself (post hoc)."]
rm[Item == "Column dictionary", Content := "Definitions of all column names in Tables S2–S8."]
rm[Item == "Source", Content := "Tables S2–S8 are exported unchanged from script-generated result files (see Data Availability Statement); only column definitions and translations were added. Table S8 and the last six rows of Table S6 were added on 2026-11-09 (Table S1, Nos. 29–30)."]
k <- which(rm$Item == "Table S7")
rm <- rbind(rm[1:k], data.table(Item = "Table S8", Content = "Post hoc summaries: (a) by cohort and class, null-A 95th percentile of |z|, proportion of random gene sets with |z| >= 1.96 and approximate smallest detectable HR per SD; (b) passes under alternative false-discovery-rate families."), rm[(k + 1):.N])
OUT$README <- rm

# Column dictionary
cd <- IN$Column_dictionary
add_cd <- data.table(Table = c(rep("S8a", 7), rep("S8b", 5)),
  Column = c("cohort / class / events", "se_med", "p95_med / p95_min / p95_max", "hr_thr_med", "hr_thr_approx",
             "frac_sig_med / frac_sig_min / frac_sig_max", "(all)",
             "z_real / p_emp_A", "pass_main", "q_global / q_class / q_sig", "pass_global / pass_class / pass_signature", "(all)"),
  Definition = c("Cohort, signature class and number of deaths",
                 "Median standard error of the log hazard ratio per SD across readouts in that cohort",
                 "Median, minimum and maximum across signatures of the 95th percentile of |z| in null distribution A",
                 "Median of exp(p95 x se_med): approximate smallest HR per SD that could exceed the null 95th percentile",
                 "Same, approximated as exp(median p95 / sqrt(events))",
                 "Median, minimum and maximum proportion of the 1,000 random gene sets with |z| >= 1.96",
                 "Post hoc summary; see Table S1, No. 29",
                 "Observed z and empirical p against null distribution A",
                 "Pass in the main analysis (BH within cohort and class)",
                 "BH-adjusted p_emp_A across all 77 combinations / within class / within signature",
                 "q < 0.05 and |z| >= 1.96 under each alternative family",
                 "Post hoc sensitivity analysis; see Table S1, No. 29"))
OUT$Column_dictionary <- rbind(cd, add_cd, use.names = TRUE)

ord <- c("README", "S1_deviation_log", "S2_registry", "S3_classB_nullB", "S4_all_combinations", "S5_confounding",
         "S6_fibroblast_proliferation", "S7_SIG012_published_coef", "S8a_null_summary_by_cohort", "S8b_alternative_FDR", "Column_dictionary")
stopifnot(setequal(ord, names(OUT))); OUT <- OUT[ord]
write.xlsx(OUT, FO, colWidths = "auto")

# 核对
cat("\n核对：\n"); RB <- lapply(setNames(ord, ord), function(s) as.data.table(read.xlsx(FO, s, check.names = FALSE, sep.names = " ")))
for (s in c("S2_registry", "S3_classB_nullB", "S4_all_combinations", "S5_confounding", "S7_SIG012_published_coef"))
  cat(sprintf("  %-28s identical to input: %s\n", s, isTRUE(all.equal(RB[[s]], IN[[s]], check.attributes = FALSE))))
cat("  S1 rows 1–28 identical:", isTRUE(all.equal(RB$S1_deviation_log[1:28], IN$S1_deviation_log, check.attributes = FALSE)), "; rows:", nrow(RB$S1_deviation_log), "\n")
cat("  S6 rows 1–24 identical:", isTRUE(all.equal(RB$S6_fibroblast_proliferation[1:24], IN$S6_fibroblast_proliferation, check.attributes = FALSE)), "; rows:", nrow(RB$S6_fibroblast_proliferation), "\n")
han <- sapply(RB, function(d) sum(vapply(d, function(v) sum(grepl("\\p{Han}", as.character(v), perl = TRUE)), numeric(1))))
cat("  含中文字符的单元格数（按 sheet）：\n"); print(han)
cat("\noutput:", FO, "\nmd5:", unname(tools::md5sum(FO)), "\n"); sink()
