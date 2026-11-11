# sci2_supp_S2route_20261110.R —— 以 S1–S8 (20261109) 为输入：S2 增 found_via 列与 2 行、S1 追加 No.31、README/Column dictionary 更新
# 登记：DEVIATIONS.md 2026-11-10「签名来源途径重建；Table S2 增 found_via 列」（运行之前）
suppressPackageStartupMessages({ library(openxlsx); library(data.table) })
setwd(path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"))
FI <- "data_derived/supp_tables_20261109/SCI-2_Supplementary_Tables_S1-S8_EN.xlsx"
stopifnot(unname(tools::md5sum(FI)) == "5c5cc692934a1291709e10a8e9c571bb")
OD <- "data_derived/supp_tables_20261110"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
FO <- file.path(OD, "SCI-2_Supplementary_Tables_S1-S8_EN.xlsx")
sink(file.path("logs", "supp_S2route_20261110.log"), split = TRUE)
rd <- function(s) as.data.table(read.xlsx(FI, s, check.names = FALSE, sep.names = " "))
sh <- getSheetNames(FI); IN <- setNames(lapply(sh, rd), sh); OUT <- copy(IN)

# S2：found_via
s2 <- copy(IN$S2_registry); stopifnot(nrow(s2) == 19)
key <- ifelse(is.na(s2$id), s2$name, s2$id)
SC <- "Starting set: cited in SCANI work"; SR <- "Starting set: stromal readout"
AD <- "Added while reviewing sources"; NR <- "Not recorded"
route <- c(`SIG-001` = SC, `SIG-002` = SR, `SIG-003` = SR, `SIG-004` = AD, `SIG-005` = AD, `SIG-006` = SC,
           `SIG-007` = AD, `SIG-008` = AD, `SIG-009` = NR, `SIG-010` = NR, `SIG-011` = AD, `SIG-012` = AD,
           SCANI_myCAF = "Developed by the authors (added per analysis plan)",
           `Cords 2024` = SC, `EPIC CAF (package)` = SR, `xCell Fibroblasts` = SR, `Kong 2024` = NR,
           `Chen 2025` = SC, `Peng 2023` = AD)
stopifnot(setequal(key, names(route)))
s2[, found_via := unname(route[key])]
add2 <- data.table(status = "not assessed (starting set)", id = NA_character_,
  name = c("Liu 2025", "Bagaev TME classification"), class = NA_character_, n_genes = NA,
  source = c("Cited in the SCANI manuscript", "Cancer Cell 2021;39(6):845-865"),
  doi = c(NA_character_, "10.1016/j.ccell.2021.04.014"), train_cohort = NA_character_, overlap_cohorts = NA_character_,
  note = "In the starting set but not assessed as a candidate; reason not recorded (Table S1, No. 31).",
  found_via = c(SC, SR))
OUT$S2_registry <- rbind(s2, add2, use.names = TRUE, fill = TRUE)
cat("found_via 计数（18 个候选）：\n"); print(table(s2$found_via[!is.na(s2$id) & s2$id != "SCANI_myCAF" | grepl("^excluded", s2$status)]))

# S1：No.31
s1 <- IN$S1_deviation_log; stopifnot(nrow(s1) == 30)
OUT$S1_deviation_log <- rbind(s1, data.table(No. = 31L, `Date recorded` = "2026-11-10", `Plan section / item` = "Signature search (routes)",
  Description = "How each candidate was found was reconstructed from the records of the registry work: 4 from the starting literature list, 4 from the starting list of stromal readouts, 7 added while reviewing sources, and not recorded for 3 (SIG-009, SIG-010, Kong 2024). Two items of the starting list (Liu 2025; Bagaev classification) were not assessed as candidates, and the reason was not recorded.",
  `Found before or after results` = "After",
  Action = "found_via column and two rows added to Table S2; Methods section 2.3 describes the search as it was done.",
  `Effect on reported results` = "—"), use.names = TRUE)

# README / Column dictionary
rm <- IN$README
rm[Item == "Table S2", Content := paste0(Content, " Includes how each candidate was found (found_via) and two starting-set items that were not assessed.")]
rm[Item == "Source", Content := paste0(Content, " The found_via column and the last two rows of Table S2 were added on 2026-11-10 (Table S1, No. 31).")]
OUT$README <- rm
OUT$Column_dictionary <- rbind(IN$Column_dictionary, data.table(Table = "S2", Column = "found_via",
  Definition = "How the candidate was found: starting set (signature studies cited in the SCANI work, or widely used stromal readouts), added while reviewing those sources, or not recorded"), use.names = TRUE)

OUT <- OUT[sh]; write.xlsx(OUT, FO, colWidths = "auto")
cat("\n核对：\n"); RB <- lapply(setNames(sh, sh), function(s) as.data.table(read.xlsx(FO, s, check.names = FALSE, sep.names = " ")))
for (s in setdiff(sh, c("README", "S1_deviation_log", "S2_registry", "Column_dictionary")))
  cat(sprintf("  %-28s identical to input: %s\n", s, isTRUE(all.equal(RB[[s]], IN[[s]], check.attributes = FALSE))))
cat("  S1 rows 1–30 identical:", isTRUE(all.equal(RB$S1_deviation_log[1:30], IN$S1_deviation_log, check.attributes = FALSE)), "; rows:", nrow(RB$S1_deviation_log), "\n")
cat("  S2 rows 1–19, original columns identical:", isTRUE(all.equal(RB$S2_registry[1:19, names(IN$S2_registry), with = FALSE], IN$S2_registry, check.attributes = FALSE)), "; rows:", nrow(RB$S2_registry), "\n")
han <- sapply(RB, function(d) sum(vapply(d, function(v) sum(grepl("\\p{Han}", as.character(v), perl = TRUE)), numeric(1))))
cat("  含中文单元格总数:", sum(han), "\n")
print(RB$S2_registry[, .(status, id, name, found_via)])
cat("\noutput:", FO, "\nmd5:", unname(tools::md5sum(FO)), "\n"); sink()
