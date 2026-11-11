#!/usr/bin/env Rscript
# ==============================================================================
# patch_alias_fallback.R —— 给 scoring_protocol.R 加别名回退层
#
# 背景：符号归一化把 signature 基因统一到了当前官方符号（如 WISP1 -> CCN4），
#   但部分队列的平台注释仍用旧符号。归一化之后，这些基因在那些队列里反而
#   从"能找到"变成"找不到"，覆盖率被我们自己弄低了。
#
#   实测受影响 11 处：
#     SIG-003 x GSE72094  CCN4/ADGRA2/PLPPR4/NME8/TENM4  （矩阵里是旧符号）
#     SIG-006 x TCGA/GSE72094/CPTAC  TMT1A  （矩阵里是 METTL7A）
#     SIG-007 x GSE72094  CEMIP   （矩阵里是 KIAA1199）
#     SIG-009 x GSE72094  C5AR2   （矩阵里是 GPR77）
#     SIG-003 x CPTAC     CCN4
#
# 做法（保守）：
#   基因在矩阵中缺失时，检索其在 org.Hs.eg.db 中的全部别名；
#   恰有一个别名存在于该矩阵中 -> 以该别名取值并记录；
#   零个或多个 -> 视为真实缺失，不做任何猜测。
#
# 本脚本修改 scoring_protocol.R（先备份），并跑一次前后对比。
# ==============================================================================

SP <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/scripts/scoring_protocol.R")
STAMP <- format(Sys.Date(), "%Y%m%d")
OUT <- file.path("/tmp", format(Sys.time(), "alias_patch_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")

cat("给 scoring_protocol.R 加别名回退层\n时间: ", format(Sys.time()), "\n", sep = "")

# ============================================================ 1. 备份
hr("1. 备份")
if (!file.exists(SP)) stop("未找到 ", SP)
bak <- paste0(SP, ".bak_", STAMP)
if (!file.exists(bak)) file.copy(SP, bak)
cat("  原文件 -> ", basename(bak), "\n", sep = "")

src <- readLines(SP, warn = FALSE)
if (any(grepl("resolve_aliases", src))) {
  cat("\n  已包含别名回退层，无需重复修改。\n")
} else {

# ============================================================ 2. 注入
hr("2. 注入别名回退层")

PATCH <- c(
'',
'# ============================== 别名回退（SAP §4.1 补充条款）==================',
'#\' signature 基因在矩阵中缺失时，尝试用其别名取值',
'#\' 保守规则：恰有一个别名存在于矩阵中才回退；零个或多个均视为真实缺失。',
'#\' @return list(present, rescued, missing, map)',
'#\'   present 可用于取值的行名（原符号 + 回退到的别名）',
'#\'   map     data.frame(gene, used) 记录全部回退映射',
'resolve_aliases <- function(genes, rownames_expr) {',
'  genes <- unique(as.character(genes))',
'  direct <- intersect(genes, rownames_expr)',
'  miss   <- setdiff(genes, rownames_expr)',
'  if (!length(miss) || !requireNamespace("org.Hs.eg.db", quietly = TRUE))',
'    return(list(present = direct, rescued = character(0), missing = miss,',
'                map = data.frame(gene = character(0), used = character(0),',
'                                 stringsAsFactors = FALSE)))',
'  suppressPackageStartupMessages({',
'    requireNamespace("AnnotationDbi", quietly = TRUE) })',
'  rescued <- character(0); used <- character(0); still <- character(0)',
'  for (g in miss) {',
'    eid <- tryCatch(suppressMessages(AnnotationDbi::select(',
'             org.Hs.eg.db::org.Hs.eg.db, keys = g, columns = "ENTREZID",',
'             keytype = "SYMBOL")$ENTREZID[1]), error = function(e) NA)',
'    al <- if (is.na(eid)) character(0) else tryCatch(suppressMessages(',
'            AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = eid,',
'              columns = "ALIAS", keytype = "ENTREZID")$ALIAS),',
'            error = function(e) character(0))',
'    hit <- setdiff(intersect(unique(al), rownames_expr), c(direct, used))',
'    if (length(hit) == 1L) {                 # 恰有一个才回退',
'      rescued <- c(rescued, g); used <- c(used, hit)',
'    } else {',
'      still <- c(still, g)',
'    }',
'  }',
'  list(present = c(direct, used), rescued = rescued, missing = still,',
'       map = data.frame(gene = rescued, used = used, stringsAsFactors = FALSE))',
'}',
'')

# gene_coverage 改为带回退版本
NEW_COV <- c(
'#\' signature 基因在该队列的覆盖率（含别名回退，SAP §4.1）',
'#\' 低于 PROTOCOL$min_coverage 的组合标记为不可评估，不给经验 p 值',
'gene_coverage <- function(genes, expr, use_alias = TRUE) {',
'  genes <- unique(as.character(genes))',
'  if (use_alias) {',
'    r <- resolve_aliases(genes, rownames(expr))',
'    present <- r$present; missing <- r$missing',
'    rescued <- r$rescued; map <- r$map',
'  } else {',
'    present <- intersect(genes, rownames(expr)); missing <- setdiff(genes, present)',
'    rescued <- character(0)',
'    map <- data.frame(gene = character(0), used = character(0), stringsAsFactors = FALSE)',
'  }',
'  frac <- length(present) / length(genes)',
'  list(n_total = length(genes), n_present = length(present),',
'       frac = frac, evaluable = frac >= PROTOCOL$min_coverage,',
'       present = present, missing = missing,',
'       n_rescued = length(rescued), alias_map = map)',
'}')

i0 <- grep("^gene_coverage <- function", src)
if (!length(i0)) stop("未找到 gene_coverage 定义")
# 找到该函数的结束（下一个顶层 } ）
j <- i0; depth <- 0L
repeat {
  depth <- depth + lengths(regmatches(src[j], gregexpr("\\{", src[j]))) -
                   lengths(regmatches(src[j], gregexpr("\\}", src[j])))
  if (depth <= 0L && j > i0) break
  j <- j + 1L
  if (j > length(src)) stop("未能定位 gene_coverage 的结束")
}
# 注释行往上找函数的文档块起点
i_doc <- i0
while (i_doc > 1L && grepl("^#'", src[i_doc - 1L])) i_doc <- i_doc - 1L

cat("  gene_coverage 原定义: 第 ", i_doc, " - ", j, " 行\n", sep = "")
src2 <- c(src[seq_len(i_doc - 1L)], PATCH, NEW_COV, src[(j + 1L):length(src)])

# 打分函数：也用回退后的基因
src2 <- sub('gs <- lapply\\(gene_sets, function\\(g\\) intersect\\(unique\\(as.character\\(g\\)\\), rownames\\(expr\\)\\)\\)',
            'gs <- lapply(gene_sets, function(g) resolve_aliases(g, rownames(expr))$present)',
            src2)

# PROTOCOL 加一条记录
src2 <- sub('  drop_constant  = FALSE',
            '  drop_constant  = FALSE,           # 不预先剔除常数行，交由 GSVA 处理\n  alias_fallback = TRUE',
            src2, fixed = TRUE)

writeLines(src2, SP)
cat("  已注入 resolve_aliases()，并改写 gene_coverage() 与 score_ssgsea()\n")
}

# ============================================================ 3. 前后对比
hr("3. 覆盖率前后对比")
suppressMessages(source(SP))
d <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")
sigs <- sort(list.files(file.path(d, "genes"), pattern = "^SIG-[0-9]+_genes\\.txt$"))
sigs <- sub("_genes\\.txt$", "", sigs)

rows <- list(); RESC <- list()
for (a in COHORTS) {
  co <- tryCatch(suppressMessages(load_cohort(a, verbose = FALSE)), error = function(e) NULL)
  if (is.null(co)) next
  for (s in sigs) {
    g <- readLines(file.path(d, "genes", paste0(s, "_genes.txt")))
    g <- trimws(g); g <- g[nzchar(g)]
    old <- gene_coverage(g, co$expr, use_alias = FALSE)
    new <- gene_coverage(g, co$expr, use_alias = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = a, sig_id = s, n = length(g),
      before = old$n_present, after = new$n_present, rescued = new$n_rescued,
      frac_before = round(old$frac, 3), frac_after = round(new$frac, 3),
      eval_before = old$evaluable, eval_after = new$evaluable,
      stringsAsFactors = FALSE)
    if (new$n_rescued > 0L) RESC[[length(RESC) + 1L]] <-
      cbind(cohort = a, sig_id = s, new$alias_map)
  }
}
R <- do.call(rbind, rows)

chg <- R[R$rescued > 0L, ]
cat("  发生回退的组合: ", nrow(chg), " / ", nrow(R), "\n\n", sep = "")
if (nrow(chg)) print(chg[, c("cohort","sig_id","n","before","after","rescued",
                             "frac_before","frac_after")], row.names = FALSE)

flip <- R[!R$eval_before & R$eval_after, ]
if (nrow(flip)) {
  cat("\n  ★ 由『不可评估』翻转为『可评估』的组合:\n")
  print(flip[, c("cohort","sig_id","frac_before","frac_after")], row.names = FALSE)
} else cat("\n  无组合发生可评估性翻转\n")

still <- R[!R$eval_after, ]
if (nrow(still)) {
  cat("\n  回退后仍不可评估（真实平台缺失）:\n")
  print(still[, c("cohort","sig_id","n","after","frac_after")], row.names = FALSE)
} else cat("\n  回退后全部组合均可评估\n")

if (length(RESC)) {
  M <- do.call(rbind, RESC)
  f <- file.path(d, sprintf("alias_fallback_map_%s.csv", STAMP))
  write.csv(M, f, row.names = FALSE)
  cat("\n  全部回退映射 (", nrow(M), " 条) -> ", f, "\n", sep = "")
  cat("  该文件须作为补充材料随论文提交。\n")
}

hr("4. 待办")
cat("  1) 在 SAP §14 追加修订记录：\n\n")
cat("     §4.1 补充：signature 基因在某队列表达矩阵中缺失时，检索其在\n")
cat("     org.Hs.eg.db 中的全部别名；恰有一个别名存在于该矩阵中则以该别名\n")
cat("     取值并记录映射；零个或多个均视为真实缺失。全部回退映射见补充材料。\n")
cat("     修订理由：符号归一化把 signature 统一到当前官方符号后，部分队列的\n")
cat("     平台注释仍用旧符号，导致覆盖率被系统性低估。\n")
cat("     发生于结果产出之前。\n\n")
cat("  2) 重跑 lock_registry.R 确认，然后 --lock\n")

cat("\n输出: ", OUT, "\n", sep = "")
