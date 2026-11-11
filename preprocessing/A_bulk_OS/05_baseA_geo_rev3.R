#!/usr/bin/env Rscript
# ==============================================================================
# 05_baseA_geo.R (rev3) —— 底座 A · GEO 队列整备（5 个队列）
#
# rev3 相对 rev2 的改动：
#   [1] 新增组织学筛选（histo_key / histo_keep）—— GSE30219 与 GSE50081 是
#       NSCLC 混合队列，必须筛腺癌。SAP §2.1 口径：严格筛 ADC。
#   [2] 新增 GSE50081（时间单位为【年】，×365.25）与 GSE30219（月，×30.4375）
#   [3] 组织学与组织类型两级筛选分开报告，便于在论文中给出筛选前后样本数
#
# 保留 rev2 的全部修正：按【键名】解析 characteristics_ch1（非槽位位置）、
# TNM 大写化、"NA"/"--" 字符串归为缺失、滤掉 ! 开头行、探针交集。
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(AnnotationDbi)
})

OUT_ROOT <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived")
RAW      <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/raw")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT_ROOT, format(Sys.time(), "geo_build_%Y%m%d_%H%M.log"))
sink(LOG, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")
sh <- function(t) cat("\n  --- ", t, "\n", sep="")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k,":"), paste(v, collapse=", ")))

na_str <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("NA","N/A","na","","--","-",".","null","NULL","nd","Unknown","NTL")] <- NA
  x
}

# ============================== 队列口径配置 =================================
# histo_key / histo_keep：组织学筛选。NA 表示该队列全为 LUAD，无需筛。
CFG <- list(
  GSE68465 = list(
    file = "~/caf_spatial_project/data/GSE68465/GSE68465_series_matrix.txt.gz",
    platform = "GPL96", anno_pkg = "hgu133a.db", anno_soft = NA,
    tumor_key = "disease_state", tumor_val = "Lung Adenocarcinoma",
    histo_key = NA, histo_keep = NA,
    status_key = "vital_status", dead_pat = "^Dead$",
    time_key = "months_to_last_contact_or_death", time_mult = 30.4375, time_src = "月",
    stage_key = "disease_stage", stage_type = "tnm",
    age_key = "age", sex_key = "Sex", extra_keys = NULL,
    assert_n = NA_integer_, assert_ev = NA_integer_
  ),
  GSE31210 = list(
    file = "~/caf_spatial_project/results/10_GSE31210_validation/GSE31210_series_matrix.txt.gz",
    platform = "GPL570", anno_pkg = "hgu133plus2.db", anno_soft = NA,
    tumor_key = "tissue", tumor_val = "primary lung tumor",
    histo_key = NA, histo_keep = NA,
    status_key = "death", dead_pat = "^dead$",
    time_key = "days before death/censor", time_mult = 1, time_src = "天",
    stage_key = "pathological stage", stage_type = "roman",
    age_key = "age (years)", sex_key = "gender",
    exclude_key = "exclude for prognosis analysis due to incomplete resection or adjuvant therapy",
    exclude_val = "exclude",   # SAP §2.4：遵从原作者的预后分析排除标志
    extra_keys = "exclude for prognosis analysis due to incomplete resection or adjuvant therapy",
    assert_n = NA_integer_, assert_ev = NA_integer_
  ),
  GSE72094 = list(
    file = "~/caf_spatial_project/data/GSE72094_series_matrix.txt.gz",
    platform = "GPL15048", anno_pkg = NA,
    anno_soft = "~/caf_spatial_project/data/GPL15048.txt",
    tumor_key = NA, tumor_val = NA,
    histo_key = NA, histo_keep = NA,
    status_key = "vital_status", dead_pat = "^Dead$",
    time_key = "survival_time_in_days", time_mult = 1, time_src = "天",
    stage_key = "Stage", stage_type = "arabic",
    age_key = "age_at_diagnosis", sex_key = "gender", extra_keys = NULL,
    assert_n = 398L, assert_ev = 113L          # ★ Paper 1 锁定的双数字
  ),
  GSE50081 = list(
    file = file.path(RAW, "GSE50081_series_matrix.txt.gz"),
    platform = "GPL570", anno_pkg = "hgu133plus2.db", anno_soft = NA,
    tumor_key = NA, tumor_val = NA,            # source_name 全部为 Primary NSCLC
    histo_key = "histology", histo_keep = "adenocarcinoma",   # 严格筛：精确匹配
    status_key = "status", dead_pat = "^dead$",
    time_key = "survival time", time_mult = 365.25, time_src = "年",  # ★ 单位是年
    stage_key = "Stage", stage_type = "arabic",
    age_key = "age", sex_key = "Sex", extra_keys = c("t-stage","n-stage","recurrence"),
    assert_n = NA_integer_, assert_ev = NA_integer_
  ),
  GSE30219 = list(
    file = file.path(RAW, "GSE30219_series_matrix.txt.gz"),
    platform = "GPL570", anno_pkg = "hgu133plus2.db", anno_soft = NA,
    tumor_key = NA, tumor_val = NA,            # tissue 字段大小写不一致，改用 histology 筛
    histo_key = "histology", histo_keep = "ADC",
    status_key = "status", dead_pat = "^DEAD$",
    time_key = "follow-up time (months)", time_mult = 30.4375, time_src = "月",
    stage_key = NA, stage_type = "none",       # 只有 pT/pN/pM，无分期组，见下方 §分期
    pt_key = "pt stage", pn_key = "pn stage", pm_key = "pm stage",
    age_key = "age at surgery", sex_key = "gender", extra_keys = NULL,
    assert_n = NA_integer_, assert_ev = NA_integer_
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

# --------------------------------- 按【键名】解析（非槽位位置）
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
      cval = ifelse(hask, trimws(sub("^[^:]*:\\s*", "", v)), trimws(v)))
  }
  long <- rbindlist(recs)
  long <- long[!is.na(ckey) & nzchar(ckey) & nzchar(cval)]
  long <- long[, .(cval = cval[1]), by = .(gsm, ckey)]
  list(gsm = gsm, long = long)
}

get_field <- function(pc, want) {
  if (is.na(want)) return(rep(NA_character_, length(pc$gsm)))
  s <- pc$long[ckey == want]
  if (!nrow(s)) return(rep(NA_character_, length(pc$gsm)))
  s$cval[match(pc$gsm, s$gsm)]
}

tnm2stage <- function(tt, nn, mm = NULL) {
  tt <- toupper(na_str(tt)); nn <- toupper(na_str(nn))
  tt <- ifelse(grepl("T[0-4X]", tt), sub(".*T([0-4X]).*", "\\1", tt), NA)
  nn <- ifelse(grepl("N[0-3X]", nn), sub(".*?N([0-3X]).*", "\\1", nn), NA)
  s <- rep(NA_character_, length(tt))
  s[which(nn=="0" & tt %in% c("1","2"))] <- "I"
  s[which(nn=="0" & tt=="3")]            <- "II"
  s[which(nn=="1" & tt %in% c("1","2"))] <- "II"
  s[which(nn=="1" & tt=="3")]            <- "III"
  s[which(nn %in% c("2","3"))]           <- "III"
  s[which(tt=="4")]                      <- "III"
  if (!is.null(mm)) {
    m1 <- grepl("M1", toupper(na_str(mm)))
    s[which(m1)] <- "IV"
  }
  s
}

norm_stage <- function(x, type) {
  x <- na_str(x)
  if (type == "roman") {
    s <- sub("^Stage\\s*", "", x, ignore.case = TRUE); s <- sub("[ABC]$", "", s)
  } else if (type == "arabic") {
    s <- sub("[ABC]$", "", x)
    s <- unname(c("1"="I","2"="II","3"="III","4"="IV")[s])
  } else if (type == "tnm") {
    s <- tnm2stage(x, x)
    bad <- is.na(s) & !is.na(x)
    if (any(bad)) cat("      !! 未能映射的 TNM 取值: ",
                      paste(unique(x[bad]), collapse=" | "), "\n", sep="")
  } else {
    s <- rep(NA_character_, length(x))
  }
  s[!s %in% c("I","II","III","IV")] <- NA
  s
}

probe2symbol <- function(cfg, probes) {
  if (!is.na(cfg$anno_pkg)) {
    if (!requireNamespace(cfg$anno_pkg, quietly = TRUE)) {
      cat("      !! 注释包 ", cfg$anno_pkg, " 未安装\n", sep=""); return(NULL) }
    suppressPackageStartupMessages(library(cfg$anno_pkg, character.only = TRUE))
    db <- get(cfg$anno_pkg)
    ok <- intersect(probes, AnnotationDbi::keys(db, keytype = "PROBEID"))
    kv("  库中可查的探针", sprintf("%d / %d", length(ok), length(probes)))
    if (!length(ok)) return(NULL)
    map <- suppressMessages(AnnotationDbi::select(db, keys = ok,
             columns = "SYMBOL", keytype = "PROBEID"))
    map <- map[!is.na(map$SYMBOL), ]; map <- map[!duplicated(map$PROBEID), ]
    kv("  注释来源", paste(cfg$anno_pkg, as.character(packageVersion(cfg$anno_pkg))))
    return(setNames(map$SYMBOL, map$PROBEID))
  }
  f <- path.expand(cfg$anno_soft)
  if (!file.exists(f)) { cat("      !! 平台注释缺失\n"); return(NULL) }
  L <- readLines(f, warn = FALSE)
  b <- grep("^!platform_table_begin", L); e <- grep("^!platform_table_end", L)
  if (!length(b)) { cat("      !! SOFT 中无 platform_table\n"); return(NULL) }
  e <- if (length(e)) e[1] else length(L) + 1L
  tb <- fread(text = L[(b[1]+1L):(e-1L)], sep = "\t", header = TRUE,
              quote = "", showProgress = FALSE)
  scol <- grep("symbol", colnames(tb), ignore.case = TRUE, value = TRUE)
  if (!length(scol)) { cat("      !! 未找到 symbol 列: ",
                          paste(colnames(tb), collapse=", "), "\n", sep=""); return(NULL) }
  kv("  使用的 symbol 列", scol[1])
  sy <- trimws(sub("\\s*///.*$", "", as.character(tb[[scol[1]]])))
  keep <- nzchar(sy) & !is.na(sy) & !sy %in% c("---","NA","-")
  kv("  可用 symbol 的探针数", sum(keep))
  setNames(sy[keep], as.character(tb[[colnames(tb)[1]]])[keep])
}

# ============================================================ 主流程
# ---- 可选：只重建指定队列（避免波及其余队列）
.args <- commandArgs(trailingOnly = TRUE)
.h <- grep("^--cohort=", .args, value = TRUE)
ONLY <- if (length(.h)) strsplit(sub("^--cohort=", "", .h[1]), ",")[[1]] else NULL
if (!is.null(ONLY)) cat("\n[限定队列] ", paste(ONLY, collapse = ", "), "\n", sep = "")

for (acc in names(CFG)) {
  if (!is.null(ONLY) && !acc %in% ONLY) next
  cfg <- CFG[[acc]]; f <- path.expand(cfg$file)
  hr(paste0(acc, "   ", cfg$platform))
  if (!file.exists(f)) { cat("  !! 文件不存在，跳过: ", f, "\n", sep=""); next }

  h  <- read_header(f)
  pc <- parse_char_bykey(h)
  kv("样本数（原始）", length(pc$gsm))

  kt <- pc$long[, .N, by = ckey][order(-N)]
  need <- na.omit(c(cfg$status_key, cfg$time_key, cfg$stage_key, cfg$histo_key))
  miss <- need[!need %in% kt$ckey]
  if (length(miss)) {
    cat("\n  !! 配置的键不存在: ", paste(miss, collapse=" | "), "\n", sep="")
    cat("     实际键名: ", paste(kt$ckey, collapse=" | "), "\n", sep="")
    cat("     本队列跳过。\n"); next
  }

  # ---- 两级筛选
  sh("筛选")
  cl <- data.table(gsm = pc$gsm)
  n0 <- nrow(cl)

  if (!is.na(cfg$tumor_key)) {
    tis <- get_field(pc, cfg$tumor_key)
    cl$is_tumor <- !is.na(tis) & tis == cfg$tumor_val
  } else cl$is_tumor <- TRUE
  cat(sprintf("    组织类型筛选后        %d / %d\n", sum(cl$is_tumor), n0))

  if (!is.na(cfg$histo_key)) {
    ht <- na_str(get_field(pc, cfg$histo_key))
    cat("    histology 取值分布:\n")
    tb <- sort(table(ht, useNA="ifany"), decreasing=TRUE)
    for (i in seq_along(tb)) cat(sprintf("      %-46s %4d\n",
        substr(ifelse(is.na(names(tb)[i]),"NA",names(tb)[i]), 1, 46), tb[i]))
    cl$is_histo <- !is.na(ht) & ht == cfg$histo_keep     # 精确匹配，不做模糊
    cat(sprintf("    ★ 严格筛 '%s' 后      %d / %d\n", cfg$histo_keep, sum(cl$is_histo), n0))
    cat("      注：模糊写法（如 'NSCLC-favor adenocarcinoma'）一律排除\n")
  } else cl$is_histo <- TRUE
  cl$keep <- cl$is_tumor & cl$is_histo
  kv("两级筛选后", sum(cl$keep))

  # ---- 临床
  sh("临床字段")
  cl[, status    := na_str(get_field(pc, cfg$status_key))]
  cl[, tval      := suppressWarnings(as.numeric(na_str(get_field(pc, cfg$time_key))))]
  cl[, age       := suppressWarnings(as.numeric(na_str(get_field(pc, cfg$age_key))))]
  cl[, sex       := tolower(na_str(get_field(pc, cfg$sex_key)))]
  cl[, OS_event  := as.integer(!is.na(status) & grepl(cfg$dead_pat, status, ignore.case = TRUE))]
  cl[, OS_time   := tval * cfg$time_mult]
  kv("时间单位换算", sprintf("%s x %.4f -> 天", cfg$time_src, cfg$time_mult))

  # ---- 分期
  cat("    分期:\n")
  if (cfg$stage_type == "none" && !is.null(cfg$pt_key)) {
    cat("      本队列无分期组，由 pT/pN/pM 推导（AJCC 6 近似，标记为 derived）\n")
    cl[, stage_raw := paste(na_str(get_field(pc, cfg$pt_key)),
                            na_str(get_field(pc, cfg$pn_key)),
                            na_str(get_field(pc, cfg$pm_key)), sep="/")]
    cl[, stage := tnm2stage(get_field(pc, cfg$pt_key),
                            get_field(pc, cfg$pn_key),
                            get_field(pc, cfg$pm_key))]
  } else {
    cl[, stage_raw := na_str(get_field(pc, cfg$stage_key))]
    cl[, stage := norm_stage(stage_raw, cfg$stage_type)]
  }
  tb <- table(cl$stage[cl$keep], useNA="ifany")
  for (i in seq_along(tb)) cat(sprintf("      %-6s %4d\n",
      ifelse(is.na(names(tb)[i]),"NA",names(tb)[i]), tb[i]))

  if (!is.null(cfg$extra_keys)) for (k in cfg$extra_keys) {
    vv <- na_str(get_field(pc, k))
    cl[[paste0("extra_", make.names(substr(k, 1, 14)))]] <- vv
    cat(sprintf("    附加字段 %-32s 非缺失 %d\n", substr(k,1,32), sum(!is.na(vv))))
  }

  # ---- SAP §2.4：遵从数据集自带的 exclude 标志（未设该字段的队列不受影响）
  excl <- rep(FALSE, nrow(cl))
  if (!is.null(cfg$exclude_key) && !is.na(cfg$exclude_key)) {
    ev <- na_str(get_field(pc, cfg$exclude_key))
    excl <- !is.na(ev) & ev == cfg$exclude_val
    cat(sprintf("    ★ SAP 2.4 exclude 标志: %d 例标记为 '%s'，予以排除\n",
                sum(excl), cfg$exclude_val))
    cat(sprintf("       该字段非缺失 %d 例，缺失 %d 例（缺失视为不排除）\n",
                sum(!is.na(ev)), sum(is.na(ev))))
  }
  cl[, excluded_by_flag := excl]
  cl[, valid := keep & !excl & !is.na(OS_time) & OS_time > 0 & !is.na(status)]
  kv("可用于生存分析", sum(cl$valid))
  kv("其中事件数", sum(cl$OS_event[cl$valid]))
  kv("事件率", sprintf("%.1f%%", 100*mean(cl$OS_event[cl$valid])))
  kv("中位随访（天）", round(median(cl$OS_time[cl$valid], na.rm=TRUE), 1))

  if (!is.na(cfg$assert_n)) {
    sh("自检")
    cat(sprintf("    n      : 得到 %d，锁定值 %d  %s\n", sum(cl$valid), cfg$assert_n,
                ifelse(sum(cl$valid)==cfg$assert_n, "OK", "!!")))
    cat(sprintf("    事件数 : 得到 %d，锁定值 %d  %s\n", sum(cl$OS_event[cl$valid]),
                cfg$assert_ev, ifelse(sum(cl$OS_event[cl$valid])==cfg$assert_ev, "OK", "!!")))
  }

  # ---- 表达矩阵
  sh("表达矩阵")
  ex <- fread(cmd = paste("zcat", shQuote(f)), skip = "ID_REF", sep = "\t",
              header = TRUE, showProgress = FALSE)
  id1 <- as.character(ex[[1]])
  ex  <- ex[!is.na(id1) & nzchar(id1) & !grepl("^!", id1)]
  pid <- as.character(ex[[1]])
  M   <- suppressWarnings(as.matrix(ex[, -1, with = FALSE]))
  if (!is.numeric(M)) storage.mode(M) <- "numeric"
  rownames(M) <- pid; colnames(M) <- gsub('"', '', colnames(M))
  kv("原始维度", paste(dim(M), collapse=" x "))
  smp <- M[seq_len(min(2000,nrow(M))), seq_len(min(20,ncol(M))), drop=FALSE]
  expr_med <- median(smp, na.rm=TRUE); expr_max <- max(smp, na.rm=TRUE)
  kv("取值范围（抽样）", sprintf("[%.3g, %.3g]  中位 %.3g",
      min(smp,na.rm=TRUE), expr_max, expr_med))
  kv("推定尺度", ifelse(expr_max < 25, "log2（不再变换）", "线性（下游需 log2(x+1)）"))

  p2s <- probe2symbol(cfg, rownames(M))
  if (is.null(p2s)) { cat("  -> 注释不可用，跳过\n"); rm(ex, M); invisible(gc()); next }

  keep <- rownames(M) %in% names(p2s)
  M2 <- M[keep, , drop=FALSE]; sy <- unname(p2s[rownames(M2)])
  rm2 <- rowMeans(M2, na.rm = TRUE)
  o <- order(sy, -rm2); M2 <- M2[o,,drop=FALSE]; sy <- sy[o]
  k2 <- !duplicated(sy); M2 <- M2[k2,,drop=FALSE]; rownames(M2) <- sy[k2]
  M2 <- M2[order(rownames(M2)), , drop=FALSE]
  kv("  合并后维度（同名取均表达最大）", paste(dim(M2), collapse=" x "))

  sel <- cl$gsm[cl$valid]; sel <- sel[sel %in% colnames(M2)]
  if (length(sel) != sum(cl$valid))
    cat("    !! 有 ", sum(cl$valid)-length(sel), " 个样本在表达矩阵中缺列\n", sep="")
  M2 <- M2[, sel, drop=FALSE]
  clv <- cl[match(sel, cl$gsm)]
  stopifnot(identical(colnames(M2), clv$gsm))
  kv("  最终分析集", sprintf("%d 基因 x %d 样本，事件 %d",
                             nrow(M2), ncol(M2), sum(clv$OS_event)))

  d <- file.path(OUT_ROOT, acc); dir.create(d, showWarnings = FALSE)
  saveRDS(M2, file.path(d, paste0(acc, "_expr_symbol.rds")))
  write.csv(clv, file.path(d, paste0(acc, "_clinical_clean.csv")), row.names = FALSE)
  write.csv(data.frame(
    cohort = acc, platform = cfg$platform,
    n_raw = length(pc$gsm), n_after_tumor = sum(cl$is_tumor),
    n_after_histology = sum(cl$keep), n_excluded_by_flag = sum(cl$excluded_by_flag), n_final = ncol(M2),
    histology_filter = ifelse(is.na(cfg$histo_key), "无（全为 LUAD）",
                              paste0(cfg$histo_key, " == '", cfg$histo_keep, "' 精确匹配")),
    n_events = sum(clv$OS_event), n_genes = nrow(M2),
    time_key = cfg$time_key, time_src_unit = cfg$time_src,
    time_mult = cfg$time_mult, time_unit = "days",
    stage_type = cfg$stage_type, n_with_stage = sum(!is.na(clv$stage)),
    has_survival = TRUE,
    anno = ifelse(is.na(cfg$anno_pkg), "GPL SOFT table", cfg$anno_pkg),
    expr_median_sample = round(expr_med, 3),
    expr_scale = ifelse(expr_max < 25, "log2", "linear"),
    need_log = expr_max >= 25,
    gene_collapse = "同名取均表达最大",
    median_followup_days = round(median(clv$OS_time, na.rm=TRUE), 1),
    source_md5 = tools::md5sum(f), built = format(Sys.time(), "%Y-%m-%d %H:%M"),
    stringsAsFactors = FALSE),
    file.path(d, paste0("cohort_manifest_", acc, ".csv")), row.names = FALSE)
  cat("    输出 -> ", d, "\n", sep="")
  rm(ex, M, M2); invisible(gc())
}

hr("sessionInfo"); print(sessionInfo())
cat("\n日志: ", LOG, "\n", sep="")
