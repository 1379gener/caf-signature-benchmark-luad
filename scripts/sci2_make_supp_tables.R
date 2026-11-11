#!/usr/bin/env Rscript
# ==============================================================================
# sci2_make_supp_tables.R —— SCI-2 稿件 v1 的补充表 S1–S7（只读已有结果；不写任何数字常量）
#   S1 偏差记录（DEVIATIONS.md 原样复制）
#   S2 签名注册表 + 排除的候选（tables_20261025/TableS2_registry.csv、03_signature_registry/excluded_notes.csv）+ SCANI 行
#   S3 B 类签名的零分布 B（tables_20261025/TableS3_classB_nullB.csv）
#   S4 逐组合主结果（A 类含 SCANI：classA_scani_bh_20261025.csv；B 类：sap6_bh_combo_20261025.csv；HR 取自 results_extended / scani 汇总）
#   S5 阶段 2 逐组合衰减（stage2_attenuation_20261017.csv）
#   S6 成纤维细胞读数 ± 增殖校正（stage3_explore_fibro_prolif_20261025.csv）
#   S7 SIG-012 原文权重（sens_SIG012_origweights_20261025.csv）
#   输出：data_derived/supp_tables_<日期>/（已存在即停止）：S1.md、S2–S7.csv；若装有 openxlsx 另出一个 xlsx（每表一页）
#   用法：Rscript sci2_make_supp_tables.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived")
REG <- file.path(ROOT, "03_signature_registry")
OD <- file.path(DD, sprintf("supp_tables_%s", format(Sys.Date(), "%Y%m%d"))); if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
rd <- function(f) { stopifnot(file.exists(f)); cat(sprintf("  %-44s md5 %s\n", basename(f), unname(tools::md5sum(f)))); fread(f, na.strings = c("", "NA")) }
RL <- c("开发队列" = "development", "原文验证队列" = "validation (original study)", "完全独立队列" = "fully independent")
OUT <- list()

# S1
file.copy(file.path(CH, "DEVIATIONS.md"), file.path(OD, "S1_deviation_log.md"))
cat("  S1 DEVIATIONS.md md5 ", unname(tools::md5sum(file.path(CH, "DEVIATIONS.md"))), "\n", sep = "")

# S2
R2 <- rd(file.path(DD, "tables_20261025/TableS2_registry.csv")); X2 <- rd(file.path(REG, "excluded_notes.csv"))
SC <- rd(file.path(DD, "scani_supp_summary_20261025.csv"))
stopifnot(nrow(R2) == 12L, nrow(X2) == 6L, nrow(SC) == 6L)
S2 <- rbind(R2[, .(status = "included", id = sig_id, name, class, n_genes, source, doi, train_cohort, overlap_cohorts, note = NA_character_)],
            data.table(status = "included (added per analysis plan; developed by the authors)", id = "SCANI_myCAF", name = "SCANI myCAF component",
                       class = "A", n_genes = SC$n_genes[1], source = "manuscript submitted", doi = NA_character_,
                       train_cohort = "TCGA (stage II-IV)", overlap_cohorts = "TCGA, GSE72094", note = NA_character_),
            X2[, .(status = paste("excluded:", criterion), id = NA_character_, name, class = NA_character_, n_genes = NA_integer_,
                   source, doi = NA_character_, train_cohort = NA_character_, overlap_cohorts = NA_character_, note = reason)], fill = TRUE)
OUT$S2_registry <- S2

# S3
OUT$S3_classB_nullB <- rd(file.path(DD, "tables_20261025/TableS3_classB_nullB.csv"))
stopifnot(nrow(OUT$S3_classB_nullB) == 41L)

# S4
A <- rd(file.path(DD, "classA_scani_bh_20261025.csv")); B <- rd(file.path(DD, "sap6_bh_combo_20261025.csv"))
E <- rd(sort(list.files(DD, "^results_extended_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1])[evaluable == TRUE]
stopifnot(nrow(A) == 36L, B[class == "B", .N] == 41L)
HR <- rbind(E[, .(cohort, sig_id, hr_per_sd = hr_real)], SC[, .(cohort, sig_id, hr_per_sd = hr_real)])
S4 <- rbind(A[, .(class = "A", cohort, sig_id, role, z = z_real, zA_p95, p_emp_A, q_A, p_emp_B, q_B, coexpr_qc_pass, classification = cls3_bh)],
            B[class == "B", .(class, cohort, sig_id, role, z = z_real, zA_p95, p_emp_A, q_A, p_emp_B = NA_real_, q_B = NA_real_,
                              coexpr_qc_pass = NA, classification = cls3_bh)])
S4 <- merge(S4, HR, by = c("cohort", "sig_id"), all.x = TRUE); stopifnot(nrow(S4) == 77L, !anyNA(S4$hr_per_sd))
S4[, role := RL[role]]; stopifnot(!anyNA(S4$role))
setcolorder(S4, c("class", "sig_id", "cohort", "role", "z", "hr_per_sd")); OUT$S4_all_combinations <- S4[order(class, sig_id, cohort)]

# S5–S7
OUT$S5_confounding <- rd(file.path(DD, "stage2_attenuation_20261017.csv"))
OUT$S6_fibroblast_proliferation <- rd(file.path(DD, "stage3_explore_fibro_prolif_20261025.csv"))
OUT$S7_SIG012_published_coef <- rd(file.path(DD, "sens_SIG012_origweights_20261025.csv")); stopifnot(nrow(OUT$S7_SIG012_published_coef) == 6L)

for (nm in names(OUT)) { fwrite(OUT[[nm]], file.path(OD, paste0(nm, ".csv"))); cat(sprintf("  %-32s %4d 行 %3d 列\n", nm, nrow(OUT[[nm]]), ncol(OUT[[nm]]))) }
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(OUT, file.path(OD, "SCI-2_Supplementary_Tables_S2-S7.xlsx")); cat("  xlsx 已写出\n")
} else cat("  未安装 openxlsx：只输出 csv\n")
cat("\n写出目录 ", OD, "\n", sep = ""); print(list.files(OD))
