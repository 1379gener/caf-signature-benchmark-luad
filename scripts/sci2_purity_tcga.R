#!/usr/bin/env Rscript
# ==============================================================================
# sci2_purity_tcga.R —— TCGA 纯度（SAP §7.2：Aran 等 2015 Supplementary Data 1 的 ABSOLUTE 列）
#   1 读入并断言原始表；2 LUAD 原发灶、按病人去重；3 对底座 443 人的覆盖率；
#   4 与 SE 的 paper_Purity.ABSOLUTE.calls（逗号小数）在共同病人上对照；
#   5 M2 / M4 完整病例样本量预估；6 写出 M2/M4 读取的纯度文件
#   只读原始数据与底座；写 data_derived/ 与 logs/
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_purity_tcga.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl); library(SummarizedExperiment) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
suppressMessages(source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")))
LOG <- file.path(CH, "logs", format(Sys.time(), "purity_tcga_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")

hr("1. 原始表")
F <- list.files(file.path(ROOT, "01_bases/A_bulk_OS/raw/purity"), "ncomms9971.*\\.xlsx$", full.names = TRUE)
stopifnot(length(F) == 1)
MD5 <- unname(tools::md5sum(F)); cat(basename(F), "\n文件 md5 ", MD5, "\n", sep = "")
# xlsx 是 zip 包，整文件 md5 会因打包元数据不同而变；内容一致性以内部数据部件的 md5 判定
ux <- file.path(tempdir(), "aran_xlsx"); unzip(F, files = c("xl/worksheets/sheet1.xml", "xl/sharedStrings.xml"), exdir = ux)
PM <- unname(tools::md5sum(file.path(ux, c("xl/worksheets/sheet1.xml", "xl/sharedStrings.xml"))))
cat("sheet1.xml md5 ", PM[1], "\nsharedStrings.xml md5 ", PM[2], "\n", sep = "")
REF <- c("153821725201bbc4f16cf51db7fab760", "85c3c6de97439267d7b1f40dcba4f6b0")   # 2026-10-17 对话中上传版本的内部部件
cat("数据部件与对话中上传版本一致：", identical(PM, REF), "\n", sep = "")
A <- as.data.table(read_excel(F, sheet = "Supp Data 1", skip = 3, na = c("", "NaN")))   # 原表缺失写作字符串 NaN
stopifnot(identical(names(A)[1:7], c("Sample ID", "Cancer type", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")))
setnames(A, c("Sample ID", "Cancer type"), c("sample", "cancer"))
A <- A[!is.na(sample)]
cat("总行数 ", nrow(A), "；重复样本 ID ", anyDuplicated(A$sample), "\n", sep = "")
stopifnot(!anyDuplicated(A$sample), is.numeric(A$ABSOLUTE))
stopifnot(nrow(A) == 9364)                                 # 与本地核对版本同为 9,364 个样本

hr("2. LUAD 原发灶，按病人去重")
L <- A[cancer == "LUAD"]
L[, `:=`(pid = substr(sample, 1, 12), stype = substr(sample, 14, 15))]
cat("LUAD 样本 ", nrow(L), "，病人 ", uniqueN(L$pid), "\n", sep = ""); print(table(sample_type = L$stype))
P <- L[stype == "01" & !is.na(ABSOLUTE)]
k <- P[, .(n_vial = .N, n_val = uniqueN(ABSOLUTE)), by = pid]
cat("原发灶有 ABSOLUTE 的病人 ", nrow(k), "；多管 ", sum(k$n_vial > 1), "；多管取值不一致 ", sum(k$n_val > 1), "\n", sep = "")
stopifnot(all(k$n_val == 1))                               # 同一病人多管取值一致，去重无歧义
P <- unique(P[, .(pid, purity_absolute = ABSOLUTE)])
stopifnot(!anyDuplicated(P$pid), nrow(P) == 354)          # 本地核对：LUAD 原发灶有 ABSOLUTE 者 354 人
print(summary(P$purity_absolute))

hr("3. 对底座 TCGA 443 人的覆盖")
co <- load_cohort("TCGA_LUAD", verbose = FALSE); cl <- co$clin
stopifnot(!anyDuplicated(cl$patient_id), nrow(cl) == 443)
B <- merge(cl[, .(pid = patient_id, OS_time, OS_event, age, sex)], P, by = "pid", all.x = TRUE)
stopifnot(nrow(B) == 443)
cat(sprintf("有 ABSOLUTE：%d / 443（%.1f%%）；其中事件 %d / %d\n", sum(!is.na(B$purity_absolute)),
            100 * mean(!is.na(B$purity_absolute)), sum(B$OS_event[!is.na(B$purity_absolute)]), sum(B$OS_event)))
cat("有无纯度两组的事件率：\n"); print(B[, .(n = .N, event_rate = round(mean(OS_event), 3)), by = .(has_purity = !is.na(purity_absolute))])

hr("4. 与 SE paper_Purity.ABSOLUTE.calls 对照")
se <- readRDS(path.expand("~/caf_spatial_project/data/TCGA_LUAD/TCGA_LUAD_SE.rds"))
cd <- as.data.table(as.data.frame(colData(se)))
pc <- intersect(c("patient", "bcr_patient_barcode", "submitter_id"), names(cd))[1]
cd[, pid := substr(get(pc), 1, 12)]
S <- unique(cd[pid %in% cl$patient_id, .(pid, se_raw = as.character(paper_Purity.ABSOLUTE.calls))])[!is.na(se_raw)]
S[, se_purity := as.numeric(sub(",", ".", se_raw, fixed = TRUE))]
stopifnot(!anyDuplicated(S$pid), !anyNA(S$se_purity))
J <- merge(P, S[, .(pid, se_purity)], by = "pid")
cat(sprintf("SE 列有值 %d 人；两源共同 %d 人；SE 有而 Aran 无 %d 人\n", nrow(S), nrow(J), nrow(S[!pid %in% P$pid])))
if (nrow(J)) cat(sprintf("Spearman %.3f；|差| 中位 %.3f，最大 %.3f；完全相同 %d 人\n",
    cor(J$purity_absolute, J$se_purity, method = "spearman"), median(abs(J$purity_absolute - J$se_purity)),
    max(abs(J$purity_absolute - J$se_purity)), sum(abs(J$purity_absolute - J$se_purity) < 1e-9)))

hr("5. M2 / M4 完整病例预估（不跑模型）")
ST <- fread(sort(list.files(file.path(CH, "data_derived"), "^m4_stage_binary_[0-9]{8}\\.csv$", full.names = TRUE),
                 decreasing = TRUE)[1], na.strings = c("", "NA"))[cohort == "TCGA_LUAD"]
SS <- fread(sort(list.files(file.path(CH, "data_derived"), "^m4_stage_binary_TCGA_sens_[0-9]{8}\\.csv$", full.names = TRUE),
                 decreasing = TRUE)[1], na.strings = c("", "NA"))
B <- merge(B, ST[, .(pid = sample_id, stage_ge_II)], by = "pid")
B <- merge(B, SS[, .(pid = sample_id, stage_ge_II_sens)], by = "pid"); stopifnot(nrow(B) == 443)
cc <- function(ok) sprintf("%d 人 / %d 事件", sum(ok), sum(B$OS_event[ok]))
cat("M2（纯度）                 ", cc(!is.na(B$purity_absolute)), "\n")
cat("M4 主（纯度+GDC 分期+年龄+性别）", cc(!is.na(B$purity_absolute) & !is.na(B$stage_ge_II) & !is.na(B$age) & !is.na(B$sex)), "\n")
cat("M4 敏感性（分期补充版）     ", cc(!is.na(B$purity_absolute) & !is.na(B$stage_ge_II_sens) & !is.na(B$age) & !is.na(B$sex)), "\n")
cat("（增殖评分由表达计算，全覆盖，不影响完整病例数）\n")

hr("6. 写出")
fo <- file.path(CH, "data_derived", sprintf("m2_purity_TCGA_%s.csv", format(Sys.Date(), "%Y%m%d")))
if (file.exists(fo)) stop("已存在，不覆盖：", fo)
fwrite(B[, .(cohort = "TCGA_LUAD", sample_id = pid, purity_absolute)], fo, na = "NA")
cat("写出：", fo, "\n来源：", basename(F), "（文件 md5 ", MD5, "；sheet1.xml md5 ", PM[1], "）\n", sep = "")
sink(); close(zz)
