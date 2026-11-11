#!/usr/bin/env Rscript
# ==============================================================================
# sci2_m4_stage_coding.R —— M4 的 stage 二分类编码（I vs ≥II）：定义 + 各队列分布登记
#   决定日期 2026-10-17，阶段 2 任何模型运行之前（SAP §14 追加）
#   只读底座；输出 data_derived/m4_stage_binary_<日期>.csv 供 M4 读取；日志写 logs/
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_m4_stage_coding.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
suppressMessages(source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")))
STAMP <- format(Sys.Date(), "%Y%m%d")
LOG <- file.path(CH, "logs", format(Sys.time(), "m4_stage_coding_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)

#' 分期组 -> 二分类。只接受罗马数字分期组（可带 Stage 前缀与 A/B/C 亚分），其余一律缺失
stage_bin <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- sub("^STAGE\\s*", "", s); s <- sub("[ABC][0-9]?$", "", s)
  out <- rep(NA_integer_, length(s))
  out[s == "I"] <- 0L; out[s %in% c("II", "III", "IV")] <- 1L
  out }

OUT <- list(); TAB <- list()
for (a in COHORTS) {
  co <- load_cohort(a, verbose = FALSE)
  if (!co$has_survival) { cat(a, "：无生存数据，不进入 M4\n"); next }
  cl <- co$clin; id <- co$id_col
  stopifnot(!anyDuplicated(cl[[id]]), nrow(cl) == co$n)
  stopifnot("stage" %in% names(cl))
  b <- stage_bin(cl$stage)
  raw_unmapped <- unique(cl$stage[is.na(b) & !is.na(cl$stage) & nzchar(as.character(cl$stage))])
  cat(sprintf("\n%-10s n=%d\n", a, co$n)); print(table(stage = cl$stage, bin = b, useNA = "ifany"))
  if (length(raw_unmapped)) cat("  !! 未映射的非空取值：", paste(raw_unmapped, collapse = " | "), "\n")
  TAB[[a]] <- data.table(cohort = a, n = co$n, stage_I = sum(b == 0L, na.rm = TRUE),
                         stage_ge_II = sum(b == 1L, na.rm = TRUE), stage_missing = sum(is.na(b)),
                         events_ge_II = sum(cl$OS_event[which(b == 1L)], na.rm = TRUE))
  OUT[[a]] <- data.table(cohort = a, sample_id = as.character(cl[[id]]), stage = cl$stage, stage_ge_II = b)
}
T <- rbindlist(TAB); T[, pct_missing := round(100 * stage_missing / n, 1)]
cat("\n汇总：\n"); print(T, row.names = FALSE)
fo <- file.path(CH, "data_derived", sprintf("m4_stage_binary_%s.csv", STAMP))
if (file.exists(fo)) stop("已存在，不覆盖：", fo)
fwrite(rbindlist(OUT), fo); fwrite(T, sub("\\.csv$", "_summary.csv", fo))
cat("\n写出：", fo, "\n")
sink(); close(zz)
