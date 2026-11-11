#!/usr/bin/env Rscript
# ==============================================================================
# normalize_symbols.R —— signature 基因符号归一化（锁定前的最后一道工序）
#
# 解决两类问题：
#   [1] signature 用了过时符号（SEPP1 -> SELENOP，WISP2 -> CCN5）
#   [2] 队列的平台注释用了过时符号，导致 signature 基因"看起来缺失"
#
# 两者都会系统性低估覆盖率，进而把"符号老化"误判为"该队列不可评估"。
#
# 处理原则（保守）：
#   - 只在别名唯一映射到一个 ENTREZ 基因时才改写，多重映射一律标记不改
#   - 改写前备份原文件，全部改写记入映射日志
#   - 队列侧的别名命中只报告、不修改底座（底座由 scoring_protocol 统一管理）
#
# 用法：
#   Rscript normalize_symbols.R           # 只报告，不改写
#   Rscript normalize_symbols.R --write   # 改写 signature 基因文件（先备份）
# ==============================================================================

WRITE <- "--write" %in% commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
})

REGDIR <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")
GENDIR <- file.path(REGDIR, "genes")
STAMP  <- format(Sys.Date(), "%Y%m%d")
OUT    <- file.path("/tmp", format(Sys.time(), "symnorm_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
sh <- function(t) cat("\n  --- ", t, "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("signature 基因符号归一化\n时间: ", format(Sys.time()),
    "\n模式: ", ifelse(WRITE, "报告 + 改写", "仅报告（加 --write 改写）"), "\n", sep = "")

# ============================================================ 参照表
hr("0. 参照表")
OFFICIAL <- suppressMessages(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"))
kv("org.Hs.eg.db 版本", as.character(packageVersion("org.Hs.eg.db")))
kv("官方 SYMBOL 数", length(OFFICIAL))

#' 把一个可能过时的符号解析为当前官方符号
#' @return list(status, to, n_map)
#'   status: current / renamed / ambiguous / unknown
resolve <- function(sym) {
  if (sym %in% OFFICIAL) return(list(status = "current", to = sym, n_map = 1L))
  eid <- tryCatch(suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = sym, columns = "ENTREZID",
                          keytype = "ALIAS")$ENTREZID), error = function(e) character(0))
  eid <- unique(na.omit(eid))
  if (!length(eid)) return(list(status = "unknown", to = sym, n_map = 0L))
  if (length(eid) > 1L) return(list(status = "ambiguous", to = sym, n_map = length(eid)))
  new <- tryCatch(suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = eid, columns = "SYMBOL",
                          keytype = "ENTREZID")$SYMBOL[1]), error = function(e) NA)
  if (is.na(new) || new == sym) return(list(status = "unknown", to = sym, n_map = 1L))
  list(status = "renamed", to = new, n_map = 1L)
}

#' 某符号的全部别名（含自身）
aliases_of <- function(sym) {
  eid <- tryCatch(suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = sym, columns = "ENTREZID",
                          keytype = "SYMBOL")$ENTREZID[1]), error = function(e) NA)
  if (is.na(eid)) return(sym)
  al <- tryCatch(suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = eid, columns = "ALIAS",
                          keytype = "ENTREZID")$ALIAS), error = function(e) character(0))
  unique(c(sym, al))
}

# ============================================================ 1. 逐 signature 检查
hr("1. signature 符号检查")
fs <- sort(list.files(GENDIR, pattern = "^SIG-[0-9]+_genes\\.txt$", full.names = TRUE))
if (!length(fs)) stop("未找到任何 SIG-*_genes.txt")
kv("signature 文件数", length(fs))

MAP <- list()
for (f in fs) {
  sid <- sub("_genes\\.txt$", "", basename(f))
  g <- trimws(readLines(f)); g <- g[nzchar(g)]
  res <- lapply(g, resolve)
  st  <- vapply(res, `[[`, "", "status")
  to  <- vapply(res, `[[`, "", "to")
  n_ren <- sum(st == "renamed"); n_amb <- sum(st == "ambiguous"); n_unk <- sum(st == "unknown")
  cat(sprintf("  %-10s %3d 基因   当前 %3d  改写 %2d  歧义 %2d  未知 %2d\n",
              sid, length(g), sum(st == "current"), n_ren, n_amb, n_unk))
  if (n_ren) for (i in which(st == "renamed"))
    cat(sprintf("      改写  %-12s -> %s\n", g[i], to[i]))
  if (n_amb) for (i in which(st == "ambiguous"))
    cat(sprintf("      !! 歧义 %-12s 映射到 %d 个基因，保持原样\n", g[i], res[[i]]$n_map))
  if (n_unk) for (i in which(st == "unknown"))
    cat(sprintf("      !! 未知 %-12s org.Hs.eg.db 中无此符号或别名，保持原样\n", g[i]))
  MAP[[sid]] <- data.table(sig_id = sid, from = g, to = to, status = st)
}
M <- rbindlist(MAP)

sh("汇总")
print(M[, .N, by = status], row.names = FALSE)
if (any(M$status == "renamed")) {
  sh("全部改写项")
  print(M[status == "renamed", .(sig_id, from, to)], row.names = FALSE)
}

# ============================================================ 2. 队列侧别名命中
hr("2. 队列侧别名解析（signature 基因缺失时，检查其别名是否在矩阵中）")
sp <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/scripts/scoring_protocol.R")
if (!file.exists(sp)) {
  cat("  未找到 scoring_protocol.R，跳过本节\n")
} else {
  suppressMessages(source(sp))
  RESCUE <- list()
  for (a in COHORTS) {
    co <- tryCatch(suppressMessages(load_cohort(a, verbose = FALSE)), error = function(e) NULL)
    if (is.null(co)) next
    rn <- rownames(co$expr)
    rows <- list()
    for (sid in names(MAP)) {
      g <- MAP[[sid]]$to                     # 用归一化后的符号去查
      miss <- setdiff(g, rn)
      for (mg in miss) {
        al <- setdiff(aliases_of(mg), mg)
        hit <- intersect(al, rn)
        rows[[length(rows) + 1L]] <- data.table(
          cohort = a, sig_id = sid, gene = mg,
          rescued = length(hit) > 0, alias_in_matrix = paste(hit, collapse = ","))
      }
    }
    if (length(rows)) RESCUE[[a]] <- rbindlist(rows)
  }
  if (length(RESCUE)) {
    R <- rbindlist(RESCUE)
    sh("缺失基因及其别名命中情况")
    cat(sprintf("  缺失组合总数 %d，其中可由别名救回 %d\n", nrow(R), sum(R$rescued)))
    if (any(R$rescued)) {
      cat("\n  ★ 可救回的（说明是符号老化，非平台缺失）：\n")
      print(R[rescued == TRUE, .(cohort, sig_id, gene, alias_in_matrix)], row.names = FALSE)
      cat("\n  -> 这些应在 scoring_protocol 的 gene_coverage/打分环节做别名回退，\n")
      cat("     或在底座整备时把矩阵行名归一化。两者择一，须写入 SAP §4.1。\n")
    }
    cat("\n  真实缺失（平台无探针，无别名可救）按 signature 计数：\n")
    print(R[rescued == FALSE, .N, by = .(sig_id, cohort)][order(sig_id, cohort)], row.names = FALSE)
    fR <- file.path(REGDIR, sprintf("alias_rescue_report_%s.csv", STAMP))
    if (WRITE) { fwrite(R, fR); cat("\n  别名报告 -> ", fR, "\n", sep = "") }
  } else cat("  无缺失基因\n")
}

# ============================================================ 3. 改写
hr("3. 改写")
if (!WRITE) {
  cat("  未加 --write，本次只报告不改写。\n")
  cat("  确认上面的改写项无误后，重跑并加 --write。\n")
} else {
  n_files <- 0L
  bakdir <- file.path(GENDIR, paste0("pre_symnorm_", STAMP))
  dir.create(bakdir, showWarnings = FALSE)
  for (sid in names(MAP)) {
    m <- MAP[[sid]]
    if (!any(m$status == "renamed")) next
    fg <- file.path(GENDIR, paste0(sid, "_genes.txt"))
    file.copy(fg, file.path(bakdir, basename(fg)), overwrite = FALSE)
    writeLines(m$to, fg)
    fw <- file.path(GENDIR, paste0(sid, "_weights.csv"))
    if (file.exists(fw)) {
      file.copy(fw, file.path(bakdir, basename(fw)), overwrite = FALSE)
      w <- fread(fw); w$gene <- m$to[match(w$gene, m$from)]
      fwrite(w, fw)
    }
    n_files <- n_files + 1L
    cat(sprintf("  %s 已改写（%d 处），原文件备份至 %s/\n", sid,
                sum(m$status == "renamed"), basename(bakdir)))
  }
  fM <- file.path(REGDIR, sprintf("symbol_mapping_%s.csv", STAMP))
  fwrite(M, fM)
  cat("\n  改写文件数 ", n_files, "\n  映射日志 -> ", fM, "\n", sep = "")
  cat("  备份目录 -> ", bakdir, "\n", sep = "")
}

hr("4. 下一步")
cat("  1) 确认改写项无误 -> 加 --write 执行\n")
cat("  2) 决定队列侧别名回退的实现位置（scoring_protocol 或底座整备），写入 SAP §4.1\n")
cat("  3) 补齐 signature 至 SAP §3.3 的目标区间\n")
cat("  4) 锁定清单：导出 signature_registry_locked_%s.csv，算 md5 写入 SAP §14\n")
cat("  5) 锁定后方可运行 run_benchmark.R\n")
cat("\n输出: ", OUT, "\n", sep = "")
