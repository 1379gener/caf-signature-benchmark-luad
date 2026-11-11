#!/usr/bin/env Rscript
# ==============================================================================
# make_train_annotation.R —— 把 analyze_cohort_role.R 里写在代码中的「训练队列人工标注」
#   导出为独立的元数据文件，供 verify 与后续分析读取（标注本身是文献信息，不是数据结果）
#   放置：~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/make_train_annotation.R
#   输出：data_derived/train_cohort_annotation.csv（已存在则不覆盖）
# ==============================================================================
CH <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark")
src <- file.path(CH, "scripts/analyze_cohort_role.R")
out <- file.path(CH, "data_derived/train_cohort_annotation.csv")
if (file.exists(out)) { cat("已存在，不覆盖：", out, "\n"); quit(status = 0) }
x <- readLines(src, warn = FALSE)
s <- grep("^TRAIN <- data.table\\(", x); e <- s
while (!grepl("\\)\\)\\s*(#.*)?$", x[e])) e <- e + 1
env <- new.env(); suppressPackageStartupMessages(library(data.table))
eval(parse(text = x[s:e]), envir = env)
TR <- get("TRAIN", env)
cm <- x[s:e][grepl("SIG-[0-9]{3}", x[s:e]) & grepl("#", x[s:e])]
stopifnot(length(cm) == nrow(TR))
TR[, note := trimws(sub("^[^#]*#", "", cm))]
stopifnot(all(startsWith(TR$note, TR$sig_id)))           # 注释与行一一对应
TR[, source := sprintf("%s L%d–%d", basename(src), s, e)]
stopifnot(nrow(TR) == 12)
fwrite(TR, out, na = "NA"); print(TR, row.names = FALSE); cat("写出：", out, "\n")
