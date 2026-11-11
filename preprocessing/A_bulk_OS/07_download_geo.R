#!/usr/bin/env Rscript
# ==============================================================================
# 07_download_geo.R —— 下载底座 A 缺失的两个队列
#
# GSE50081 / GSE30219，均为 GPL570（与 GSE31210 同平台，hgu133plus2.db 已装）。
#
# 下载后立即做两项校验：
#   1) gzip -t —— 本项目有过 7.6G 截断文件的先例，GEO 下载中断很常见
#   2) 样本数与平台 —— 从 series matrix 头部直接读，与 GEO 网页记录核对
#
# 输出目录：~/thesis_LUAD_CAF/01_bases/A_bulk_OS/raw/
# ==============================================================================

suppressPackageStartupMessages(library(GEOquery))
options(timeout = 3600)

ACCS <- c("GSE50081", "GSE30219")
DEST <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/raw")
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(DEST, format(Sys.time(), "download_%Y%m%d_%H%M.log"))
sink(LOG, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 66), "\n", t, "\n", strrep("=", 66), "\n", sep = "")

cat("GEO 队列下载\n时间: ", format(Sys.time()), "\n目标: ", DEST, "\n", sep = "")

# ------------------------------------------------------------------ 下载
hr("1. 下载")
for (acc in ACCS) {
  cat("\n  === ", acc, "  ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
  existing <- list.files(DEST, pattern = paste0(acc, ".*series_matrix.*gz$"),
                         full.names = TRUE, recursive = TRUE)
  if (length(existing)) {
    cat("    已存在，跳过下载：", basename(existing[1]), "\n", sep = "")
    next
  }
  t0 <- Sys.time()
  r <- tryCatch(
    getGEO(acc, destdir = DEST, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) { cat("    !! 失败: ", conditionMessage(e), "\n", sep = ""); NULL })
  if (!is.null(r)) {
    cat(sprintf("    完成，耗时 %.1f 分钟，返回 %d 个 ExpressionSet\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), length(r)))
  }
}

# ------------------------------------------------------------------ 校验
hr("2. 完整性校验")

read_head <- function(f) {
  con <- gzfile(f, "rt"); on.exit(close(con))
  out <- character(0)
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (!length(l) || grepl("^!series_matrix_table_begin", l)) break
    out <- c(out, l)
  }
  out
}

files <- list.files(DEST, pattern = "series_matrix.*\\.gz$",
                    full.names = TRUE, recursive = TRUE)
if (!length(files)) {
  cat("  !! 未找到任何 series matrix 文件\n")
} else {
  for (f in files) {
    cat("\n  --- ", basename(f), "\n", sep = "")
    cat(sprintf("      体量: %.1f MB\n", file.info(f)$size / 1048576))

    ok <- system2("gzip", c("-t", shQuote(f)), stdout = FALSE, stderr = FALSE) == 0
    cat("      gzip -t : ", ifelse(ok, "PASS（未截断）", "!! FAIL —— 文件损坏，需重下"), "\n", sep = "")
    if (!ok) next

    h <- read_head(f)
    one <- function(tag) {
      x <- grep(paste0("^!", tag), h, value = TRUE)
      if (!length(x)) return(NA_character_)
      gsub('"', '', sub("^[^\t]*\t", "", x[1]))
    }
    cat("      标题  : ", substr(one("Series_title"), 1, 80), "\n", sep = "")
    cat("      平台  : ", one("Series_platform_id"), "\n", sep = "")

    a <- grep("^!Sample_geo_accession", h, value = TRUE)
    n <- if (length(a)) length(strsplit(a[1], "\t", fixed = TRUE)[[1]]) - 1 else NA
    cat("      样本数: ", n, "   ★ 须与 GEO 网页记录核对\n", sep = "")
    cat("      md5   : ", tools::md5sum(f), "\n", sep = "")

    ch <- grep("^!Sample_characteristics_ch1", h)
    cat("      characteristics_ch1 槽位数: ", length(ch), "\n", sep = "")
  }
}

hr("3. 下一步")
cat("  1) 核对上面的样本数与 GEO 网页是否一致（不一致即下载不全，须重下）\n")
cat("  2) 跑字段勘察，它会自动扫描本目录下新出现的 series matrix：\n")
cat("     Rscript ~/thesis_LUAD_CAF/00_admin/forensics_20260829/scripts/04_survey_geo.R\n")
cat("  3) 把勘察输出贴回，据实际字段名补 05_baseA_geo_rev2.R 的 CFG 后整备\n")

cat("\n日志: ", LOG, "\n", sep = "")
