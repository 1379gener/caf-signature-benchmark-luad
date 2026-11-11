#!/usr/bin/env Rscript
# ==============================================================================
# sci2_m4_stage_tcga_sens.R —— TCGA 分期敏感性版本（事前登记，阶段 2 任何模型运行之前）
#   主分析：GDC ajcc_pathologic_stage（即底座 stage 列，m4_stage_binary_<日期>.csv）
#   敏感性：GDC 缺失者用 SE colData 的 paper_Tumor.stage 补充；GDC 有值者一律不改
#   不修改底座；输出 data_derived/m4_stage_binary_TCGA_sens_<日期>.csv；日志写 logs/
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_m4_stage_tcga_sens.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(SummarizedExperiment) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
suppressMessages(source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")))
SE_F <- path.expand("~/caf_spatial_project/data/TCGA_LUAD/TCGA_LUAD_SE.rds")
LOG <- file.path(CH, "logs", format(Sys.time(), "m4_stage_tcga_sens_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)
cat("SE md5: ", unname(tools::md5sum(SE_F)), "\n", sep = "")

co <- load_cohort("TCGA_LUAD", verbose = FALSE); cl <- co$clin
stopifnot(!anyDuplicated(cl$patient_id), nrow(cl) == co$n)
fm <- sort(list.files(file.path(CH, "data_derived"), "^m4_stage_binary_[0-9]{8}\\.csv$", full.names = TRUE), decreasing = TRUE)[1]
MAIN <- fread(fm, na.strings = c("", "NA"))[cohort == "TCGA_LUAD"]     # fwrite 把 NA 写成空串，读回须还原
stopifnot(nrow(MAIN) == co$n, setequal(MAIN$sample_id, cl$patient_id))
stopifnot(sum(is.na(MAIN$stage)) == sum(is.na(cl$stage)), identical(is.na(MAIN$stage), is.na(MAIN$stage_ge_II)))
cat("主分析编码文件: ", basename(fm), "\n", sep = "")

se <- readRDS(SE_F); cd <- as.data.table(as.data.frame(colData(se)))
pc <- intersect(c("patient", "bcr_patient_barcode", "submitter_id"), names(cd))[1]; stopifnot(!is.na(pc))
cd[, pid := substr(get(pc), 1, 12)]
P <- unique(cd[pid %in% cl$patient_id, .(pid, paper = paper_Tumor.stage)])[!is.na(paper) & paper != "[Not Available]"]
stopifnot(!anyDuplicated(P$pid))
grp <- function(x) { s <- toupper(sub("^Stage\\s*", "", x, ignore.case = TRUE)); s <- sub("[ABC]$", "", s)
  ifelse(s %in% c("I", "II", "III", "IV"), s, NA_character_) }
P[, paper_grp := grp(paper)]

S <- merge(MAIN[, .(sample_id, stage, stage_ge_II)], P[, .(sample_id = pid, paper_grp)], by = "sample_id", all.x = TRUE)
S[, `:=`(stage_sens = stage, source_sens = fifelse(is.na(stage), NA_character_, "GDC"))]
S[is.na(stage) & !is.na(paper_grp), `:=`(stage_sens = paper_grp, source_sens = "paper_Tumor.stage")]
S[, stage_ge_II_sens := fifelse(is.na(stage_sens), NA_integer_, as.integer(stage_sens != "I"))]
stopifnot(S[!is.na(stage), all(stage_ge_II_sens == stage_ge_II)])    # GDC 有值者未被改动

B <- S[!is.na(stage) & !is.na(paper_grp)]
po <- mean((B$stage != "I") == (B$paper_grp != "I"))
p1 <- mean(B$stage != "I"); p2 <- mean(B$paper_grp != "I"); pe <- p1 * p2 + (1 - p1) * (1 - p2)
cat(sprintf("两源均有值 %d 人：二分类一致 %d（%.1f%%），Cohen kappa %.3f\n", nrow(B),
            sum((B$stage != "I") == (B$paper_grp != "I")), 100 * po, (po - pe) / (1 - pe)))
cat("\n来源分布：\n"); print(S[, .N, by = source_sens])
cat("\n主分析 vs 敏感性 的 stage_ge_II：\n"); print(table(main = S$stage_ge_II, sens = S$stage_ge_II_sens, useNA = "ifany"))

fo <- file.path(CH, "data_derived", sprintf("m4_stage_binary_TCGA_sens_%s.csv", format(Sys.Date(), "%Y%m%d")))
if (file.exists(fo)) stop("已存在，不覆盖：", fo)
fwrite(S[, .(cohort = "TCGA_LUAD", sample_id, stage_main = stage, stage_ge_II_main = stage_ge_II,
             stage_sens, stage_ge_II_sens, source_sens)], fo)
cat("\n写出：", fo, "\n"); sink(); close(zz)
