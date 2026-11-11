#!/usr/bin/env Rscript
# ==============================================================================
# sci2_make_tables.R —— SCI-2 正文/补充表格（只读已有结果；不写任何数字常量）
#   Table 1  队列基本信息：平台、n、事件、事件率、随访中位（反向 KM）、年龄、性别、分期（I / ≥II / 缺失）
#   Table 2  逐签名汇总：类别、可评估队列数、各信号分类计数、pass_A；A 类另报 pass_B_main / pass_B_sens
#   Table S2 签名注册表（锁定表 + 2026-10-25 角色标注）
#   Table S3 B 类签名的零分布 B 结果（SAP 未规定，仅供参考）
#   Table S4 逐组合结果（A 类、B 类分表）
#   输出：data_derived/tables_<日期>/（已存在即停止）；屏幕打印全部表，供写 Results
#   放置：ch04_benchmark/scripts/   用法：Rscript sci2_make_tables.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD <- file.path(ROOT, "03_signature_registry"); DD <- file.path(CH, "data_derived")
OD <- file.path(DD, sprintf("tables_%s", format(Sys.Date(), "%Y%m%d")))
if (dir.exists(OD)) stop("输出目录已存在，不覆盖：", OD)
dir.create(OD)
suppressMessages(source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")))
newest <- function(pat, dir = DD) { f <- sort(list.files(dir, pat, full.names = TRUE), decreasing = TRUE); stopifnot(length(f) > 0); f[1] }
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
W <- function(x, nm) { fwrite(x, file.path(OD, paste0(nm, ".csv"))); print(x, row.names = FALSE, nrows = 500) }

fE <- newest("^results_extended_.*\\.csv$"); fR <- file.path(DD, "results_by_cohort_role_20261025.csv")
fS <- newest("^m4_stage_binary_2.*\\.csv$"); fA <- file.path(DD, "train_cohort_annotation_20261025.csv")
fL <- newest("^signature_registry_locked_.*\\.csv$", REGD)
for (f in c(fE, fR, fS, fA, fL)) cat(sprintf("  %-46s md5 %s\n", basename(f), unname(tools::md5sum(f))))
E <- fread(fE)[evaluable == TRUE]; L <- fread(fL); RO <- fread(fR); AN <- fread(fA)
STG <- fread(fS, na.strings = c("", "NA"))
stopifnot(nrow(E) == 71L, !anyDuplicated(E[, paste(cohort, sig_id)]), nrow(L) == 12L)
cat("results_extended 列：", paste(names(E), collapse = ", "), "\n")

# ------------------------------------------------------------------ Table 1
hr("Table 1 队列基本信息")
T1 <- rbindlist(lapply(COHORTS, function(a) {
  co <- load_cohort(a, verbose = FALSE)
  stopifnot(!anyDuplicated(colnames(co$expr)), ncol(co$expr) == co$n, co$n == as.integer(co$man$n_final[1]))
  cl <- as.data.table(co$clin); stopifnot(nrow(cl) == co$n)
  age <- suppressWarnings(as.numeric(cl$age)); sx <- tolower(as.character(cl$sex))
  fu <- if (co$has_survival) { k <- survfit(Surv(cl$OS_time, 1 - cl$OS_event) ~ 1); unname(summary(k)$table["median"]) } else NA_real_
  st <- STG[cohort == a]
  if (nrow(st)) stopifnot(!anyDuplicated(st$sample_id), all(st$sample_id %in% colnames(co$expr)))
  data.table(cohort = a, platform = co$man$platform[1], n = co$n,
             events = if (co$has_survival) co$n_events else NA_integer_,
             event_rate = if (co$has_survival) round(co$n_events / co$n, 3) else NA_real_,
             followup_median_yr_reverseKM = round(fu / 365.25, 2),
             age_median = median(age, na.rm = TRUE), age_q1 = quantile(age, .25, na.rm = TRUE), age_q3 = quantile(age, .75, na.rm = TRUE),
             age_missing = sum(is.na(age)),
             male = sum(sx %in% c("male", "m")), female = sum(sx %in% c("female", "f")), sex_missing = sum(!sx %in% c("male", "m", "female", "f")),
             stage_I = if (nrow(st)) sum(st$stage_ge_II == 0, na.rm = TRUE) else NA_integer_,
             stage_geII = if (nrow(st)) sum(st$stage_ge_II == 1, na.rm = TRUE) else NA_integer_,
             stage_missing = if (nrow(st)) co$n - sum(!is.na(st$stage_ge_II)) else NA_integer_)
}))
W(T1, "Table1_cohorts")
cat(sprintf("  生存队列合计：n %d，事件 %d\n", sum(T1$n[!is.na(T1$events)]), sum(T1$events, na.rm = TRUE)))

# ------------------------------------------------------------------ 逐组合 + 分类
E[, sig_class := fifelse(abs(z_real) < 1.96, "no_signal", fifelse(p_emp_A >= 0.05, "not_specific", "stronger_than_random"))]
E <- merge(E, RO[, .(cohort, sig_id, role)], by = c("cohort", "sig_id"), all.x = TRUE)
stopifnot(!anyNA(E$role))

hr("Table 2 逐签名汇总（零分布 B 只对 A 类报告）")
T2 <- E[, .(n_eval = .N, no_signal = sum(sig_class == "no_signal"), not_specific = sum(sig_class == "not_specific"),
            stronger = sum(sig_class == "stronger_than_random"), pass_A = sum(p_emp_A < 0.05),
            pass_B_main = if (class[1] == "A") sum(p_emp_B < 0.05 & coexpr_qc_pass) else NA_integer_,
            pass_B_sens = if (class[1] == "A") sum(p_emp_B < 0.05) else NA_integer_,
            qc_fail = if (class[1] == "A") sum(!coexpr_qc_pass) else NA_integer_,
            median_abs_z = round(median(abs(z_real)), 2)), by = .(sig_id, class)][order(class, sig_id)]
T2 <- merge(L[, .(sig_id, name, n_genes)], T2, by = "sig_id")[order(class, sig_id)]
W(T2, "Table2_by_signature")
cat("\n  合计：\n"); print(E[, .(n = .N, no_signal = sum(sig_class == "no_signal"), not_specific = sum(sig_class == "not_specific"),
                               stronger = sum(sig_class == "stronger_than_random"), below_p95 = sum(abs(z_real) < zA_p95)), by = class])
print(E[, .(n = .N, no_signal = sum(sig_class == "no_signal"), not_specific = sum(sig_class == "not_specific"),
            stronger = sum(sig_class == "stronger_than_random"), below_p95 = sum(abs(z_real) < zA_p95))])

# ------------------------------------------------------------------ Table S2
hr("Table S2 签名注册表")
TS2 <- merge(L[, .(sig_id, name, class, n_genes, source, doi, evidence_level, overlap_cohorts)],
             AN[, .(sig_id, train_cohort)], by = "sig_id", all.x = TRUE)
TS2[sig_id == "SIG-007", `:=`(doi = "10.1016/j.lddd.2026.100333", overlap_cohorts = "TCGA, GSE30219",
                              name = paste0(name, " [first author Li W; corrected 2026-10-25, see Table S1]"))]
W(TS2, "TableS2_registry")

# ------------------------------------------------------------------ Table S3 / S4
keep <- intersect(c("cohort", "sig_id", "class", "role", "n_genes_scored", "n_present", "coverage", "n_events", "z_real", "hr_real", "hr_lo", "hr_hi",
                    "zA_p95", "p_emp_A", "sig_class", "zB_p95", "p_emp_B", "coexpr_real", "coexpr_null_med", "coexpr_ratio", "coexpr_qc_pass"), names(E))
hr("Table S3 B 类签名的零分布 B（SAP 未规定；权重按模块位置分配）")
W(E[class == "B", .SD, .SDcols = intersect(c("cohort", "sig_id", "z_real", "zB_p95", "p_emp_B", "coexpr_ratio", "coexpr_qc_pass"), names(E))][order(sig_id, cohort)],
  "TableS3_classB_nullB")
hr("Table S4a 逐组合结果：A 类")
W(E[class == "A", ..keep][order(sig_id, cohort)], "TableS4a_classA")
hr("Table S4b 逐组合结果：B 类（不含零分布 B 列）")
W(E[class == "B", .SD, .SDcols = setdiff(keep, c("zB_p95", "p_emp_B", "coexpr_real", "coexpr_null_med", "coexpr_ratio", "coexpr_qc_pass"))][order(sig_id, cohort)],
  "TableS4b_classB")

cat("\n写出目录：", OD, "\n"); print(list.files(OD))
