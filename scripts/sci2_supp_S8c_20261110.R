# sci2_supp_S8c_20261110.R —— 以 S1–S8 (20261110) 为输入：新增 S8c（比例风险检验）、S1 追加 No.32、README/Column dictionary
# 登记：DEVIATIONS.md 2026-11-10「比例风险假定检验」（运行之前）
suppressPackageStartupMessages({ library(openxlsx); library(data.table) })
setwd(path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"))
FI <- "data_derived/supp_tables_20261110/SCI-2_Supplementary_Tables_S1-S8_EN.xlsx"
stopifnot(unname(tools::md5sum(FI)) == "95a9d868a4910b8e3d931bcf30fcc33a")
OD <- "data_derived/supp_tables_20261110b"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
FO <- file.path(OD, "SCI-2_Supplementary_Tables_S1-S8_EN.xlsx")
sink(file.path("logs", "supp_S8c_20261110.log"), split = TRUE)
rd <- function(s) as.data.table(read.xlsx(FI, s, check.names = FALSE, sep.names = " "))
sh <- getSheetNames(FI); IN <- setNames(lapply(sh, rd), sh); OUT <- copy(IN)
ph <- fread("data_derived/ph_check_20261110/ph_results.csv"); stopifnot(nrow(ph) == 77, sum(ph$ph_p < 0.05) == 13)
ph[, ph_violated := ph_p < 0.05]; OUT$S8c_proportional_hazards <- ph
fwrite(ph, file.path(OD, "S8c_proportional_hazards.csv"))
s1 <- IN$S1_deviation_log; stopifnot(nrow(s1) == 31)
OUT$S1_deviation_log <- rbind(s1, data.table(No. = 32L, `Date recorded` = "2026-11-10", `Plan section / item` = "Proportional hazards (post hoc)",
  Description = "The proportional-hazards assumption of the univariable Cox model was checked for all 77 combinations (cox.zph, Kaplan-Meier transform), with hazard ratios for the first 3 years and thereafter.",
  `Found before or after results` = "After; registered before running",
  Action = "Scores recomputed with the main scoring functions; z reproduced the main results for all 77 combinations (difference < 1e-6). Labelled post hoc (Table S8c).",
  `Effect on reported results` = "No decision changed. Assumption rejected (p < 0.05) for 13 of 77 combinations (1 of 36 class A, 12 of 41 class B); in 11 of 13 the hazard ratio was larger in the first 3 years. Four of the nine passes were among them."), use.names = TRUE)
rm <- IN$README; k <- which(rm$Item == "Table S8"); stopifnot(length(k) == 1)
rm[k, Content := paste0(Content, " (c) Proportional-hazards checks for every combination.")]
rm[Item == "Source", Content := paste0(Content, " Table S8c was added on 2026-11-10 (Table S1, No. 32).")]
OUT$README <- rm
OUT$Column_dictionary <- rbind(IN$Column_dictionary, data.table(Table = "S8c",
  Column = c("z / hr", "ph_chisq / ph_p", "hr_le3y / hr_gt3y", "ph_violated"),
  Definition = c("Univariable Cox z and hazard ratio per SD (identical to the main analysis)",
                 "Test statistic and p value of cox.zph (Kaplan-Meier time transform)",
                 "Hazard ratio per SD in the first 3 years of follow-up and thereafter",
                 "ph_p < 0.05")), use.names = TRUE)
ord <- append(sh, "S8c_proportional_hazards", after = which(sh == "S8b_alternative_FDR")); OUT <- OUT[ord]
write.xlsx(OUT, FO, colWidths = "auto")
RB <- lapply(setNames(ord, ord), function(s) as.data.table(read.xlsx(FO, s, check.names = FALSE, sep.names = " ")))
for (s in setdiff(sh, c("README", "S1_deviation_log", "Column_dictionary")))
  cat(sprintf("  %-28s identical to input: %s\n", s, isTRUE(all.equal(RB[[s]], IN[[s]], check.attributes = FALSE))))
cat("  S1 rows 1–31 identical:", isTRUE(all.equal(RB$S1_deviation_log[1:31], IN$S1_deviation_log, check.attributes = FALSE)), "; rows:", nrow(RB$S1_deviation_log), "\n")
cat("  S8c rows:", nrow(RB$S8c_proportional_hazards), " violated:", sum(RB$S8c_proportional_hazards$ph_violated), "\n")
han <- sapply(RB, function(d) sum(vapply(d, function(v) sum(grepl("\\p{Han}", as.character(v), perl = TRUE)), numeric(1)))); cat("  中文单元格:", sum(han), "\n")
cat("output:", FO, "\nmd5:", unname(tools::md5sum(FO)), "\n"); sink()
