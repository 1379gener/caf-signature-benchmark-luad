#!/usr/bin/env Rscript
# ==============================================================================
# 06_baseA_cptac.R —— 底座 A · CPTAC-LUAD 整备
#
# CPTAC 有四个格式各异的文件，字段名未知。本脚本分三段：
#   A. 打印原脚本 12_cptac_survival_validation.R 的数据加载段（作为参照口径）
#   B. 勘察四个原始文件的结构
#   C. 自动探测字段并整备；探测失败则打印实际列名后停止，不猜。
#
# 输出：~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/CPTAC_LUAD/
# 参考：Paper 1 的 CPTAC 分期独立性用 n=72（软核对，非硬断言）
# ==============================================================================

suppressPackageStartupMessages({ library(data.table) })

SRC  <- path.expand("~/caf_spatial_project/data/CPTAC_LUAD")
OLD  <- path.expand("~/caf_spatial_project/scripts/12_cptac_survival_validation.R")
OUT_ROOT <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived")
D    <- file.path(OUT_ROOT, "CPTAC_LUAD")
dir.create(D, recursive = TRUE, showWarnings = FALSE)
LOG  <- file.path(D, format(Sys.time(), "build_%Y%m%d_%H%M.log"))
sink(LOG, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")
sh <- function(t) cat("\n  --- ", t, "\n", sep="")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k,":"), paste(v, collapse=", ")))

cat("底座 A · CPTAC-LUAD 整备\n时间: ", format(Sys.time()), "\n日志: ", LOG, "\n", sep="")

# ============================================================ A. 原脚本参照
hr("A. 原脚本 12_cptac_survival_validation.R 的数据加载段")
if (file.exists(OLD)) {
  L <- readLines(OLD, warn = FALSE)
  cat("  总行数: ", length(L), "\n\n", sep="")
  for (i in seq_len(min(125, length(L)))) cat(sprintf("  %4d| %s\n", i, L[i]))
} else cat("  [无] ", OLD, "\n", sep="")

# ============================================================ B. 文件勘察
hr("B. 原始文件勘察")
files <- list(
  expr    = file.path(SRC, "cptac_rna.cct"),
  os      = file.path(SRC, "os_data.csv"),
  clin    = file.path(SRC, "clinical_raw.csv"),
  cell20  = file.path(SRC, "cell2020_OS.xlsx")
)
for (n in names(files)) {
  f <- files[[n]]
  if (file.exists(f)) cat(sprintf("  %-8s %8s  %s  %s\n", n,
      paste0(round(file.info(f)$size/1048576, 2), "M"),
      format(file.info(f)$mtime, "%Y-%m-%d"), basename(f)))
  else cat(sprintf("  %-8s !! 不存在  %s\n", n, basename(f)))
}

peek <- function(f, label, nrows = 4) {
  sh(label)
  if (!file.exists(f)) { cat("    !! 不存在\n"); return(NULL) }
  d <- tryCatch(fread(f, nrows = nrows), error = function(e) { cat("    读取失败: ",
                conditionMessage(e), "\n", sep=""); NULL })
  if (is.null(d)) return(NULL)
  cat("    列数: ", ncol(d), "\n", sep="")
  cat("    列名（前 25）: ", paste(head(colnames(d), 25), collapse=", "), "\n", sep="")
  print(d[, seq_len(min(8, ncol(d))), with = FALSE])
  invisible(d)
}
peek(files$expr, "cptac_rna.cct")
peek(files$os,   "os_data.csv")
peek(files$clin, "clinical_raw.csv")

if (requireNamespace("readxl", quietly = TRUE)) {
  sh("cell2020_OS.xlsx")
  sheets <- readxl::excel_sheets(files$cell20)
  kv("  工作表", sheets)
  x <- readxl::read_excel(files$cell20, sheet = sheets[1], n_max = 4)
  cat("    列名: ", paste(colnames(x), collapse=", "), "\n", sep="")
  print(as.data.frame(x)[, seq_len(min(8, ncol(x)))])
} else cat("\n  --- cell2020_OS.xlsx：readxl 未安装，跳过\n")

# ============================================================ C. 整备
hr("C. 整备")

# ---- C1 表达矩阵
sh("C1. 表达矩阵")
ex <- fread(files$expr, sep = "\t", header = TRUE, showProgress = FALSE)
gcol <- colnames(ex)[1]
kv("首列（基因名列）", gcol)
g <- as.character(ex[[gcol]])
M <- as.matrix(ex[, -1, with = FALSE])
if (!is.numeric(M)) storage.mode(M) <- "numeric"
rownames(M) <- g
kv("原始维度", paste(dim(M), collapse=" x "))
kv("重复基因名", sum(duplicated(g)))
if (anyDuplicated(g)) {                       # 与其他队列一致：同名取均表达最大
  rm2 <- rowMeans(M, na.rm = TRUE)
  o <- order(rownames(M), -rm2); M <- M[o, , drop=FALSE]
  M <- M[!duplicated(rownames(M)), , drop=FALSE]
}
M <- M[order(rownames(M)), , drop=FALSE]
kv("去重后维度", paste(dim(M), collapse=" x "))
smp <- M[seq_len(min(2000,nrow(M))), seq_len(min(20,ncol(M))), drop=FALSE]
kv("取值范围（抽样）", sprintf("[%.3g, %.3g]  中位 %.3g",
    min(smp,na.rm=TRUE), max(smp,na.rm=TRUE), median(smp,na.rm=TRUE)))
cat("    -> 中位落在 0-15 通常已是 log2 尺度；此项须写进变量字典\n")
kv("NA 比例", sprintf("%.3f%%", 100*mean(is.na(M))))
for (kg in c("TIMP1","PDPN","PCNA","ACTA2")) 
  cat(sprintf("    %-8s %s\n", kg, ifelse(kg %in% rownames(M), "有", "!! 缺")))
kv("样本 ID 示例", head(colnames(M), 5))

# ---- C2 生存
sh("C2. 生存字段自动探测")
os <- fread(files$os)
pick <- function(dt, pat, what) {
  h <- grep(pat, colnames(dt), ignore.case = TRUE, value = TRUE)
  cat(sprintf("    %-12s 候选: %s\n", what,
      ifelse(length(h), paste(h, collapse=" | "), "（无）")))
  if (length(h)) h[1] else NA_character_
}
id_c   <- pick(os, "patient|case|sample|subject|^id$|idx", "ID 列")
time_c <- pick(os, "os_?time|overall.*surv|surv.*day|surv.*month|days_to|months", "时间列")
ev_c   <- pick(os, "os_?event|vital|status|death|censor", "事件列")

if (any(is.na(c(id_c, time_c, ev_c)))) {
  cat("\n  !! 生存字段探测失败。os_data.csv 的完整列名：\n")
  cat("     ", paste(colnames(os), collapse=", "), "\n", sep="")
  cat("  停止。请把上面的列名贴回，我据此改探测规则。\n")
  quit(status = 0)
}
kv("采用", sprintf("ID=%s  time=%s  event=%s", id_c, time_c, ev_c))

sv <- data.table(sample_id = as.character(os[[id_c]]),
                 t_raw = os[[time_c]], e_raw = os[[ev_c]])
kv("os_data 行数", nrow(sv))
kv("唯一 ID 数", length(unique(sv$sample_id)))
cat("    时间列取值示例: ", paste(head(sv$t_raw, 6), collapse=", "), "\n", sep="")
cat("    事件列取值分布: \n")
tb <- table(sv$e_raw, useNA="ifany")
for (i in seq_along(tb)) cat(sprintf("      %-24s %4d\n", names(tb)[i], tb[i]))

sv[, OS_time_raw := suppressWarnings(as.numeric(as.character(t_raw)))]
sv[, OS_event := as.integer(grepl("dead|decease|^1$|yes", as.character(e_raw), ignore.case=TRUE))]

# 时间单位：中位 < 200 视为月，否则视为天
med <- median(sv$OS_time_raw, na.rm = TRUE)
unit_mult <- if (!is.na(med) && med < 200) 30.4375 else 1
kv("时间中位（原始单位）", round(med, 2))
kv("判定单位", ifelse(unit_mult == 1, "天（不换算）", "月 -> 天（x30.4375）"))
sv[, OS_time := OS_time_raw * unit_mult]

# ---- C3 分期
sh("C3. 分期字段自动探测")
cl <- fread(files$clin)
kv("clinical_raw 维度", paste(dim(cl), collapse=" x "))
# cBioPortal 长表形态：clinicalAttributeId / value
if (all(c("clinicalAttributeId") %in% colnames(cl))) {
  kv("形态", "cBioPortal 长表")
  at <- cl[, .N, by = clinicalAttributeId][order(-N)]
  cat("    属性（前 25）：\n")
  for (i in seq_len(min(25, nrow(at)))) cat(sprintf("      %-42s %4d\n", at$clinicalAttributeId[i], at$N[i]))
  stg_at <- grep("stage", at$clinicalAttributeId, ignore.case=TRUE, value=TRUE)
  kv("分期属性候选", ifelse(length(stg_at), paste(stg_at, collapse=" | "), "（无）"))
  idc <- grep("patient|sample|case", colnames(cl), ignore.case=TRUE, value=TRUE)
  vlc <- grep("^value$|attributeValue", colnames(cl), ignore.case=TRUE, value=TRUE)
  if (length(stg_at) && length(idc) && length(vlc)) {
    st <- cl[clinicalAttributeId == stg_at[1],
             .(sample_id = as.character(get(idc[1])), stage_raw = as.character(get(vlc[1])))]
  } else st <- NULL
} else {
  kv("形态", "宽表")
  cat("    列名: ", paste(colnames(cl), collapse=", "), "\n", sep="")
  idc <- grep("patient|sample|case|^id$", colnames(cl), ignore.case=TRUE, value=TRUE)
  stc <- grep("stage", colnames(cl), ignore.case=TRUE, value=TRUE)
  kv("ID 列候选", ifelse(length(idc), paste(idc, collapse=" | "), "（无）"))
  kv("分期列候选", ifelse(length(stc), paste(stc, collapse=" | "), "（无）"))
  st <- if (length(idc) && length(stc))
    cl[, .(sample_id = as.character(get(idc[1])), stage_raw = as.character(get(stc[1])))] else NULL
}
if (!is.null(st)) {
  st <- unique(st[!is.na(stage_raw) & nzchar(stage_raw)])
  cat("    分期原始取值分布:\n")
  tb <- table(st$stage_raw, useNA="ifany")
  for (i in seq_along(tb)) cat(sprintf("      %-30s %4d\n", names(tb)[i], tb[i]))
  s <- toupper(trimws(st$stage_raw))
  s <- sub("^STAGE\\s*", "", s); s <- sub("[ABC]$", "", s)
  s[s %in% c("1","2","3","4")] <- c("I","II","III","IV")[as.integer(s[s %in% c("1","2","3","4")])]
  s[!s %in% c("I","II","III","IV")] <- NA
  st[, stage := s]
} else cat("    !! 未能提取分期 —— 本队列的 stage 将全部为 NA\n")

# ---- C4 合并
sh("C4. 合并与断言")
ids <- colnames(M)
cat("    ID 匹配情况：\n")
cat(sprintf("      表达矩阵样本数              %d\n", length(ids)))
cat(sprintf("      os_data 中能匹配上的        %d\n", sum(ids %in% sv$sample_id)))
if (sum(ids %in% sv$sample_id) == 0) {
  cat("\n  !! 表达矩阵与 os_data 的 ID 完全不匹配。\n")
  cat("     表达列名示例: ", paste(head(ids, 5), collapse=", "), "\n", sep="")
  cat("     os_data ID 示例: ", paste(head(sv$sample_id, 5), collapse=", "), "\n", sep="")
  cat("  停止。ID 格式需要转换规则，请把上面两行贴回。\n")
  quit(status = 0)
}

d <- data.table(sample_id = ids)
d <- merge(d, sv[!duplicated(sample_id), .(sample_id, OS_time, OS_event)],
           by = "sample_id", all.x = TRUE, sort = FALSE)
if (!is.null(st)) d <- merge(d, st[!duplicated(sample_id), .(sample_id, stage)],
                             by = "sample_id", all.x = TRUE, sort = FALSE)
else d[, stage := NA_character_]
setkey(d, NULL); d <- d[match(ids, sample_id)]
stopifnot(identical(d$sample_id, ids))

d[, valid := !is.na(OS_time) & OS_time > 0 & !is.na(OS_event)]
kv("可用于生存分析", sum(d$valid))
kv("其中事件数", sum(d$OS_event[d$valid]))
kv("有分期者", sum(!is.na(d$stage)))
kv("既有生存又有分期", sum(d$valid & !is.na(d$stage)))
cat("    -> Paper 1 的 CPTAC 分期独立性用 n=72，可与上一行核对（软核对）\n")
cat("    分期分布:\n")
tb <- table(d$stage, useNA="ifany")
for (i in seq_along(tb)) cat(sprintf("      %-6s %4d\n",
    ifelse(is.na(names(tb)[i]),"NA",names(tb)[i]), tb[i]))

# ---- C5 落盘
sh("C5. 输出")
Mo <- M[, d$valid, drop = FALSE]
do <- d[d$valid == TRUE]
stopifnot(identical(colnames(Mo), do$sample_id))
kv("最终分析集", sprintf("%d 基因 x %d 样本，事件 %d", nrow(Mo), ncol(Mo), sum(do$OS_event)))

f_e <- file.path(D, "CPTAC_LUAD_expr_symbol.rds")
f_c <- file.path(D, "CPTAC_LUAD_clinical_clean.csv")
f_m <- file.path(D, "cohort_manifest_CPTAC.csv")
saveRDS(Mo, f_e); write.csv(do, f_c, row.names = FALSE)
write.csv(data.frame(
  cohort = "CPTAC-LUAD", platform = "RNA-seq (LinkedOmics .cct)",
  n_raw = length(ids), n_final = ncol(Mo), n_events = sum(do$OS_event),
  n_genes = nrow(Mo), time_col = time_c, time_mult = unit_mult, time_unit = "days",
  expr_median_sample = round(median(smp, na.rm=TRUE), 3),
  median_followup_days = round(median(do$OS_time, na.rm=TRUE), 1),
  n_with_stage = sum(!is.na(do$stage)),
  source_md5 = tools::md5sum(files$expr),
  built = format(Sys.time(), "%Y-%m-%d %H:%M"), stringsAsFactors = FALSE),
  f_m, row.names = FALSE)
for (f in c(f_e, f_c, f_m)) kv(basename(f), sprintf("%.2f MB", file.info(f)$size/1048576))

hr("sessionInfo"); print(sessionInfo())
cat("\n日志: ", LOG, "\n", sep="")
