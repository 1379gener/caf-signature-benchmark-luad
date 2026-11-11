#!/usr/bin/env Rscript
# ==============================================================================
# 05_baseA_geo.R  (rev2) —— 底座 A · GEO 队列整备
#
# rev2 修正：
#   [1] data.table() 的 key= 是保留参数（设主键用），不能作列名 -> 改名 ckey
#   [2] TNM 正则用大写模式却作用于小写原串（pN0pT1）-> 先 toupper
#   [3] "NA"/"--"/"." 等字符串未被当作缺失 -> 统一 na_str() 处理
#   [4] fread 会把 !series_matrix_table_end 读成一行 -> 显式滤掉 ! 开头行
#   [5] AnnotationDbi::select 传入库中不存在的探针会报错 -> 先取交集
#
# ★ 关键设计：按【键名】解析 characteristics_ch1，不按槽位位置。
#   GSE31210 的槽位跨样本错位，按位置解析会把"距复发天数""距死亡天数"
#   "距复发月数"混进同一列，且不会报错。
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(AnnotationDbi)
})

OUT_ROOT <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT_ROOT, format(Sys.time(), "geo_build_%Y%m%d_%H%M.log"))
sink(LOG, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")
sh <- function(t) cat("\n  --- ", t, "\n", sep="")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k,":"), paste(v, collapse=", ")))

na_str <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("NA","N/A","na","","--","-",".","null","NULL","nd","Unknown")] <- NA
  x
}

# ============================== 队列口径配置 =================================
CFG <- list(
  GSE68465 = list(
    file      = "~/caf_spatial_project/data/GSE68465/GSE68465_series_matrix.txt.gz",
    platform  = "GPL96", anno_pkg = "hgu133a.db", anno_soft = NA,
    tumor_key = "disease_state", tumor_val = "Lung Adenocarcinoma",
    status_key= "vital_status",  dead_pat  = "^Dead$",
    time_key  = "months_to_last_contact_or_death", time_mult = 30.4375,  # 月 -> 天
    stage_key = "disease_stage", stage_type = "tnm",
    age_key   = "age", sex_key = "Sex", extra_keys = NULL,
    assert_n  = NA_integer_, assert_ev = NA_integer_
  ),
  GSE31210 = list(
    file      = "~/caf_spatial_project/results/10_GSE31210_validation/GSE31210_series_matrix.txt.gz",
    platform  = "GPL570", anno_pkg = "hgu133plus2.db", anno_soft = NA,
    tumor_key = "tissue", tumor_val = "primary lung tumor",
    status_key= "death", dead_pat = "^dead$",
    time_key  = "days before death/censor", time_mult = 1,
    stage_key = "pathological stage", stage_type = "roman",
    age_key   = "age (years)", sex_key = "gender",
    extra_keys= "exclude for prognosis analysis due to incomplete resection or adjuvant therapy",
    assert_n  = NA_integer_, assert_ev = NA_integer_
  ),
  GSE72094 = list(
    file      = "~/caf_spatial_project/data/GSE72094_series_matrix.txt.gz",
    platform  = "GPL15048", anno_pkg = NA,
    anno_soft = "~/caf_spatial_project/data/GPL15048.txt",
    tumor_key = NA, tumor_val = NA,                    # 全部为肿瘤
    status_key= "vital_status", dead_pat = "^Dead$",
    time_key  = "survival_time_in_days", time_mult = 1,
    stage_key = "Stage", stage_type = "arabic",
    age_key   = "age_at_diagnosis", sex_key = "gender", extra_keys = NULL,
    assert_n  = 398L, assert_ev = 113L                 # ★ Paper 1 锁定的双数字
  )
)
# =============================================================================

read_header <- function(f) {
  con <- gzfile(f, "rt"); on.exit(close(con))
  out <- character(0)
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (!length(l) || grepl("^!series_matrix_table_begin", l)) break
    out <- c(out, l)
  }
  out
}

# --------------------------------- 按【键名】解析 characteristics（核心）
parse_char_bykey <- function(h) {
  acc <- grep("^!Sample_geo_accession", h, value = TRUE)
  gsm <- gsub('"', '', strsplit(acc[1], "\t", fixed = TRUE)[[1]][-1])
  ch  <- grep("^!Sample_characteristics_ch1", h, value = TRUE)
  recs <- vector("list", length(ch))
  for (i in seq_along(ch)) {
    v <- gsub('"', '', strsplit(ch[i], "\t", fixed = TRUE)[[1]][-1])
    if (length(v) < length(gsm)) v <- c(v, rep("", length(gsm) - length(v)))
    v <- v[seq_along(gsm)]
    hask <- grepl(":", v, fixed = TRUE)
    recs[[i]] <- data.table(
      gsm  = gsm,
      ckey = ifelse(hask, trimws(sub(":.*$", "", v)), NA_character_),
      cval = ifelse(hask, trimws(sub("^[^:]*:\\s*", "", v)), trimws(v))
    )
  }
  long <- rbindlist(recs)
  long <- long[!is.na(ckey) & nzchar(ckey) & nzchar(cval)]
  long <- long[, .(cval = cval[1]), by = .(gsm, ckey)]   # 同键多次取首个非空
  list(gsm = gsm, long = long)
}

get_field <- function(pc, want) {
  if (is.na(want)) return(rep(NA_character_, length(pc$gsm)))
  s <- pc$long[ckey == want]
  if (!nrow(s)) return(rep(NA_character_, length(pc$gsm)))
  s$cval[match(pc$gsm, s$gsm)]
}

# ---------------------------------------------------------- 分期归一
norm_stage <- function(x, type) {
  x <- na_str(x)
  if (type == "roman") {
    s <- sub("^Stage\\s*", "", x, ignore.case = TRUE); s <- sub("[ABC]$", "", s)
  } else if (type == "arabic") {
    s <- sub("[ABC]$", "", x)
    s <- unname(c("1"="I","2"="II","3"="III","4"="IV")[s])
  } else if (type == "tnm") {
    xu <- toupper(x)                                   # rev2 [2]
    tt <- ifelse(grepl("T[0-4X]", xu), sub(".*T([0-4X]).*", "\\1", xu), NA)
    nn <- ifelse(grepl("N[0-3X]", xu), sub(".*?N([0-3X]).*", "\\1", xu), NA)
    s <- rep(NA_character_, length(x))
    s[which(nn=="0" & tt %in% c("1","2"))] <- "I"
    s[which(nn=="0" & tt=="3")]            <- "II"
    s[which(nn=="1" & tt %in% c("1","2"))] <- "II"
    s[which(nn=="1" & tt=="3")]            <- "III"
    s[which(nn %in% c("2","3"))]           <- "III"
    s[which(tt=="4")]                      <- "III"
    bad <- is.na(s) & !is.na(x)
    if (any(bad)) cat("      !! 未能映射的 TNM 取值: ",
                      paste(unique(x[bad]), collapse=" | "), "\n", sep="")
  } else s <- x
  s[!s %in% c("I","II","III","IV")] <- NA
  s
}

# ---------------------------------------------------------- 探针 -> symbol
probe2symbol <- function(cfg, probes) {
  if (!is.na(cfg$anno_pkg)) {
    if (!requireNamespace(cfg$anno_pkg, quietly = TRUE)) {
      cat("      !! 注释包 ", cfg$anno_pkg, " 未安装 —— 跳过本队列\n",
          "         安装: BiocManager::install(\"", cfg$anno_pkg, "\")\n", sep="")
      return(NULL)
    }
    suppressPackageStartupMessages(library(cfg$anno_pkg, character.only = TRUE))
    db <- get(cfg$anno_pkg)
    ok <- intersect(probes, AnnotationDbi::keys(db, keytype = "PROBEID"))  # rev2 [5]
    kv("  库中可查的探针", sprintf("%d / %d", length(ok), length(probes)))
    if (!length(ok)) return(NULL)
    map <- suppressMessages(AnnotationDbi::select(db, keys = ok,
             columns = "SYMBOL", keytype = "PROBEID"))
    map <- map[!is.na(map$SYMBOL), ]
    map <- map[!duplicated(map$PROBEID), ]        # 一探针多基因取第一个
    kv("  注释来源", paste(cfg$anno_pkg, as.character(packageVersion(cfg$anno_pkg))))
    return(setNames(map$SYMBOL, map$PROBEID))
  }
  f <- path.expand(cfg$anno_soft)
  if (!file.exists(f)) { cat("      !! 平台注释缺失: ", f, "\n", sep=""); return(NULL) }
  L <- readLines(f, warn = FALSE)
  b <- grep("^!platform_table_begin", L); e <- grep("^!platform_table_end", L)
  if (!length(b)) { cat("      !! SOFT 中无 platform_table\n"); return(NULL) }
  e <- if (length(e)) e[1] else length(L) + 1L
  tb <- fread(text = L[(b[1]+1L):(e-1L)], sep = "\t", header = TRUE,
              quote = "", showProgress = FALSE)
  kv("  平台表维度", paste(dim(tb), collapse=" x "))
  kv("  平台表列名", head(colnames(tb), 15))
  scol <- grep("symbol", colnames(tb), ignore.case = TRUE, value = TRUE)
  if (!length(scol)) { cat("      !! 未找到 symbol 列，请把上面的列名贴回\n"); return(NULL) }
  kv("  使用的 symbol 列", scol[1])
  idc <- colnames(tb)[1]
  sy  <- trimws(sub("\\s*///.*$", "", as.character(tb[[scol[1]]])))
  keep <- nzchar(sy) & !is.na(sy) & !sy %in% c("---","NA","-")
  kv("  可用 symbol 的探针数", sum(keep))
  setNames(sy[keep], as.character(tb[[idc]])[keep])
}

# ============================================================ 主流程
for (acc in names(CFG)) {
  cfg <- CFG[[acc]]; f <- path.expand(cfg$file)
  hr(paste0(acc, "   ", cfg$platform))
  if (!file.exists(f)) { cat("  !! 文件不存在，跳过\n"); next }

  h  <- read_header(f)
  pc <- parse_char_bykey(h)
  kv("样本数", length(pc$gsm))

  sh("键名解析结果（★ 按键名，非槽位）")
  kt <- pc$long[, .N, by = ckey][order(-N)]
  for (i in seq_len(nrow(kt))) cat(sprintf("    %-62s %4d\n", substr(kt$ckey[i],1,62), kt$N[i]))

  need <- c(cfg$status_key, cfg$time_key, cfg$stage_key)
  miss <- need[!need %in% kt$ckey]
  if (length(miss)) { cat("\n  !! 配置的键不存在: ", paste(miss, collapse=" | "),
                          "\n     按上面打印的实际键名改 CFG 后重跑。本队列跳过。\n", sep=""); next }

  # ---- 临床
  sh("临床字段")
  cl <- data.table(
    gsm       = pc$gsm,
    status    = na_str(get_field(pc, cfg$status_key)),
    tval      = suppressWarnings(as.numeric(na_str(get_field(pc, cfg$time_key)))),
    stage_raw = na_str(get_field(pc, cfg$stage_key)),
    age       = suppressWarnings(as.numeric(na_str(get_field(pc, cfg$age_key)))),
    sex       = na_str(get_field(pc, cfg$sex_key))
  )
  if (!is.na(cfg$tumor_key)) {
    tis <- get_field(pc, cfg$tumor_key)
    cat("    组织类型分布:\n"); tb <- sort(table(tis, useNA="ifany"), decreasing=TRUE)
    for (i in seq_along(tb)) cat(sprintf("      %-46s %4d\n", substr(names(tb)[i],1,46), tb[i]))
    cl$is_tumor <- !is.na(tis) & tis == cfg$tumor_val
  } else cl$is_tumor <- TRUE
  kv("肿瘤样本数", sum(cl$is_tumor))

  cl$OS_event <- as.integer(!is.na(cl$status) & grepl(cfg$dead_pat, cl$status, ignore.case = TRUE))
  cl$OS_time  <- cl$tval * cfg$time_mult
  cat("    分期归一:\n"); cl$stage <- norm_stage(cl$stage_raw, cfg$stage_type)
  tb <- table(cl$stage, useNA="ifany")
  for (i in seq_along(tb)) cat(sprintf("      %-6s %4d\n",
      ifelse(is.na(names(tb)[i]),"NA",names(tb)[i]), tb[i]))
  if (!is.null(cfg$extra_keys)) for (k in cfg$extra_keys) {
    vv <- na_str(get_field(pc, k))
    cl[[paste0("extra_", make.names(substr(k, 1, 12)))]] <- vv
    cat(sprintf("    附加字段 %-30s 非缺失 %d\n", substr(k,1,30), sum(!is.na(vv))))
    if (sum(!is.na(vv))) { t2 <- table(vv); for (i in seq_along(t2))
      cat(sprintf("        %-20s %4d\n", names(t2)[i], t2[i])) }
  }

  cl$valid <- cl$is_tumor & !is.na(cl$OS_time) & cl$OS_time > 0 & !is.na(cl$status)
  kv("时间单位换算", sprintf("原始 x %.4f -> 天", cfg$time_mult))
  kv("可用于生存分析", sum(cl$valid))
  kv("其中事件数", sum(cl$OS_event[cl$valid]))
  kv("中位随访（天）", round(median(cl$OS_time[cl$valid], na.rm=TRUE), 1))

  if (!is.na(cfg$assert_n)) {
    sh("自检")
    okn <- sum(cl$valid) == cfg$assert_n
    oke <- sum(cl$OS_event[cl$valid]) == cfg$assert_ev
    cat(sprintf("    n      : 得到 %d，锁定值 %d  %s\n", sum(cl$valid), cfg$assert_n, ifelse(okn,"OK","!!")))
    cat(sprintf("    事件数 : 得到 %d，锁定值 %d  %s\n", sum(cl$OS_event[cl$valid]), cfg$assert_ev, ifelse(oke,"OK","!!")))
    if (!okn || !oke) cat("    -> 不要改断言迁就结果。先查生存字段的缺失口径。\n")
  }

  # ---- 表达矩阵
  sh("表达矩阵")
  ex <- fread(cmd = paste("zcat", shQuote(f)), skip = "ID_REF", sep = "\t",
              header = TRUE, showProgress = FALSE)
  id1 <- as.character(ex[[1]])
  ex  <- ex[!is.na(id1) & nzchar(id1) & !grepl("^!", id1)]      # rev2 [4]
  pid <- as.character(ex[[1]])
  M   <- suppressWarnings(as.matrix(ex[, -1, with = FALSE]))
  if (!is.numeric(M)) { storage.mode(M) <- "numeric" }
  rownames(M) <- pid; colnames(M) <- gsub('"', '', colnames(M))
  kv("原始维度", paste(dim(M), collapse=" x "))
  smp <- M[seq_len(min(2000,nrow(M))), seq_len(min(20,ncol(M))), drop=FALSE]
  kv("取值范围（抽样）", sprintf("[%.3g, %.3g]  中位 %.3g",
      min(smp,na.rm=TRUE), max(smp,na.rm=TRUE), median(smp,na.rm=TRUE)))
  cat("    -> 中位落在 0-15 通常已是 log2 尺度；此项须写进变量字典\n")

  p2s <- probe2symbol(cfg, rownames(M))
  if (is.null(p2s)) { cat("  -> 注释不可用，本队列跳过输出\n"); rm(ex, M); invisible(gc()); next }

  keep <- rownames(M) %in% names(p2s)
  M2 <- M[keep, , drop = FALSE]; sy <- unname(p2s[rownames(M2)])
  kv("  可映射探针", nrow(M2)); kv("  唯一 symbol", length(unique(sy)))
  rm2 <- rowMeans(M2, na.rm = TRUE)
  o   <- order(sy, -rm2); M2 <- M2[o, , drop=FALSE]; sy <- sy[o]
  k2  <- !duplicated(sy); M2 <- M2[k2, , drop=FALSE]; rownames(M2) <- sy[k2]
  M2  <- M2[order(rownames(M2)), , drop = FALSE]
  kv("  合并后维度（多探针取均表达最大者）", paste(dim(M2), collapse=" x "))

  sel <- cl$gsm[cl$valid]
  sel <- sel[sel %in% colnames(M2)]
  if (length(sel) != sum(cl$valid))
    cat("    !! 有 ", sum(cl$valid)-length(sel), " 个样本在表达矩阵中缺列\n", sep="")
  M2 <- M2[, sel, drop = FALSE]
  clv <- cl[cl$gsm %in% sel]
  setkey(clv, gsm); clv <- clv[match(sel, clv$gsm)]
  stopifnot(identical(colnames(M2), clv$gsm))
  kv("  最终分析集", sprintf("%d 基因 x %d 样本，事件 %d",
                             nrow(M2), ncol(M2), sum(clv$OS_event)))

  d <- file.path(OUT_ROOT, acc); dir.create(d, showWarnings = FALSE)
  saveRDS(M2, file.path(d, paste0(acc, "_expr_symbol.rds")))
  write.csv(clv, file.path(d, paste0(acc, "_clinical_clean.csv")), row.names = FALSE)
  man <- data.frame(cohort = acc, platform = cfg$platform,
    n_raw = length(pc$gsm), n_tumor = sum(cl$is_tumor), n_final = ncol(M2),
    n_events = sum(clv$OS_event), n_genes = nrow(M2),
    time_key = cfg$time_key, time_mult = cfg$time_mult, time_unit = "days",
    stage_type = cfg$stage_type,
    anno = ifelse(is.na(cfg$anno_pkg), "GPL SOFT table", cfg$anno_pkg),
    expr_median_sample = round(median(smp, na.rm=TRUE), 3),
    median_followup_days = round(median(clv$OS_time, na.rm=TRUE), 1),
    source_md5 = tools::md5sum(f), built = format(Sys.time(), "%Y-%m-%d %H:%M"),
    stringsAsFactors = FALSE)
  write.csv(man, file.path(d, paste0("cohort_manifest_", acc, ".csv")), row.names = FALSE)
  cat("    输出 -> ", d, "\n", sep="")
  rm(ex, M, M2); invisible(gc())
}

hr("sessionInfo"); print(sessionInfo())
cat("\n日志: ", LOG, "\n", sep="")
