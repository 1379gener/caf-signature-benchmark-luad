# S1 追加 No.33、No.11 补一句；以 supp_tables_20261110b 为输入，输出 20261110c
suppressPackageStartupMessages({ library(openxlsx); library(data.table) })
setwd(path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"))
FI <- "data_derived/supp_tables_20261110b/SCI-2_Supplementary_Tables_S1-S8_EN.xlsx"
stopifnot(unname(tools::md5sum(FI)) == "e4b9647283a4af0f8cba5767fd547de7")
OD <- "data_derived/supp_tables_20261110c"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
FO <- file.path(OD, "SCI-2_Supplementary_Tables_S1-S8_EN.xlsx")
sh <- getSheetNames(FI); IN <- setNames(lapply(sh, function(s) as.data.table(read.xlsx(FI, s, check.names = FALSE, sep.names = " "))), sh); OUT <- copy(IN)
s1 <- copy(IN$S1_deviation_log); stopifnot(nrow(s1) == 32, grepl("CCNB1 only", s1[No. == 11, Description]))
s1[No. == 11, Description := paste0(Description, " The SAP also stated that 15 meta-PCNA genes had been renamed; the gene table shows 14.")]
OUT$S1_deviation_log <- rbind(s1, data.table(No. = 33L, `Date recorded` = "2026-11-10", `Plan section / item` = "§4.1 common gene pool; §9 circular feature selection",
  Description = "Two items specified in the plan were not carried out: the sensitivity analysis using the common cross-cohort gene pool (11,311 genes; §4.1) and the short section illustrating circular feature selection (§9). The reason was not recorded.",
  `Found before or after results` = "After", Action = "Listed in the note at the front of File S1; no analysis added.", `Effect on reported results` = "—"), use.names = TRUE)
write.xlsx(OUT, FO, colWidths = "auto")
RB <- lapply(setNames(sh, sh), function(s) as.data.table(read.xlsx(FO, s, check.names = FALSE, sep.names = " ")))
for (s in setdiff(sh, "S1_deviation_log")) cat(sprintf("  %-28s identical: %s\n", s, isTRUE(all.equal(RB[[s]], IN[[s]], check.attributes = FALSE))))
cat("  S1 rows:", nrow(RB$S1_deviation_log), "; rows other than 11 unchanged:", isTRUE(all.equal(RB$S1_deviation_log[No. != 11 & No. <= 32], IN$S1_deviation_log[No. != 11], check.attributes = FALSE)), "\n")
cat("output:", FO, "\nmd5:", unname(tools::md5sum(FO)), "\n")
