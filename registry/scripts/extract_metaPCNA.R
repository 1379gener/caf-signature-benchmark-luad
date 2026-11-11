#!/usr/bin/env Rscript
# ==============================================================================
# extract_metaPCNA.R —— 从 Venet 2011 Table S1 提取 meta-PCNA 基因集
#
# 来源：Venet D, Dumont JE, Detours V. Most Random Gene Expression Signatures Are
#       Significantly Associated with Breast Cancer Outcome.
#       PLoS Comput Biol. 2011;7(10):e1002240. doi:10.1371/journal.pcbi.1002240
#       Table S1（补充材料 PDF）
#
# 原始定义（原文 Methods）：
#   在 Ge 等 2005（GDS1096，36 个正常组织样本 / 27 个器官）上计算 PCNA 与全部基因的
#   Pearson 相关，取最高的 1%，即 13,077 个中的 131 个。
#   meta-PCNA 指数 = 这 131 个基因表达值的【中位数】（非均值，非 ssGSEA）。
#
# ★ 该基因集在【正常组织】上定义，故用它校正肿瘤数据不构成循环。
#
# 本脚本：解析 PDF 文本 -> 提取符号与相关系数 -> 归一化到当前 HGNC -> 落盘
#
# 用法：
#   Rscript extract_metaPCNA.R                # 只报告
#   Rscript extract_metaPCNA.R --write        # 落盘
# ==============================================================================

WRITE <- "--write" %in% commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
})

PDF <- path.expand("~/thesis_LUAD_CAF/03_signature_registry/pdfs/metaPCNA/TableS1_metaPCNA.pdf")
OUTD <- path.expand("~/thesis_LUAD_CAF/03_signature_registry/confounders")
STAMP <- format(Sys.Date(), "%Y%m%d")
OUT <- file.path("/tmp", format(Sys.time(), "metaPCNA_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-32s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("meta-PCNA 基因集提取\n时间: ", format(Sys.time()),
    "\n模式: ", ifelse(WRITE, "报告 + 落盘", "仅报告"), "\n", sep = "")

# ============================================================ 1. 解析 PDF
hr("1. 解析 Table S1")
if (!file.exists(PDF)) stop("未找到 ", PDF)
kv("源文件", basename(PDF))
kv("md5", unname(tools::md5sum(PDF)))

txt <- system2("pdftotext", c("-layout", shQuote(PDF), "-"), stdout = TRUE)
kv("提取行数", length(txt))
# PDF 分页处的首行会带换页符 \f (\x0c)，会让 ^ 锚定失败 —— 先剥掉
n_ff <- sum(grepl("\f", txt, fixed = TRUE))
txt <- gsub("\f", "", txt, fixed = TRUE)
kv("剥除换页符的行数", n_ff)

# 数据行形态：SYMBOL  GENENAME...  ENTREZID  CORRELATION
# 基因名可能换行，故只认「行首是符号、行尾是两个数字」的行
pat <- "^([A-Za-z0-9._-]+)\\s+.*?\\s+([0-9]+)\\s+(-?[0-9.]+)\\s*$"
hit <- grep(pat, txt, value = TRUE)
kv("匹配到的数据行", length(hit))

d <- data.table(
  symbol_2011 = sub(pat, "\\1", hit),
  entrez      = as.integer(sub(pat, "\\2", hit)),
  pcna_cor    = as.numeric(sub(pat, "\\3", hit)))
d <- d[!symbol_2011 %in% c("SYMBOL")]
d <- unique(d, by = "symbol_2011")
kv("解析出的基因数", nrow(d))
cat(sprintf("  相关系数范围 [%.2f, %.2f]\n", min(d$pcna_cor), max(d$pcna_cor)))

if (nrow(d) != 131L) {
  cat("\n  !! 解析出 ", nrow(d), " 个，原文为 131 个。解析规则可能漏行。\n", sep = "")
  cat("     未匹配上的非空行（供人工核对）：\n")
  miss <- setdiff(trimws(txt), c("", trimws(hit)))
  miss <- miss[nzchar(miss) & !grepl("^(The meta-PCNA|SYMBOL)", miss)]
  for (l in head(miss, 20)) cat("       ", substr(l, 1, 100), "\n", sep = "")
} else cat("  ✔ 与原文的 131 个一致\n")

cat("\n  相关性最高的 10 个：\n")
print(head(d[order(-pcna_cor)], 10), row.names = FALSE)

# ============================================================ 2. 符号归一化
hr("2. 符号归一化（2011 年 HGNC -> 当前）")
OFFICIAL <- suppressMessages(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"))
kv("org.Hs.eg.db 版本", as.character(packageVersion("org.Hs.eg.db")))

# 优先用 Entrez ID 映射（比符号可靠，Entrez ID 极少变动）
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
         keys = as.character(d$entrez), columns = "SYMBOL", keytype = "ENTREZID"))
setDT(map); setnames(map, c("entrez_chr", "symbol_now"))
d[, entrez_chr := as.character(entrez)]
d <- merge(d, map, by = "entrez_chr", all.x = TRUE, sort = FALSE)

d[, status := fifelse(is.na(symbol_now), "entrez 已退役",
              fifelse(symbol_now == symbol_2011, "未变", "已改名"))]
print(d[, .N, by = status][order(-N)], row.names = FALSE)

chg <- d[status == "已改名"]
if (nrow(chg)) {
  cat("\n  改名明细：\n")
  print(chg[order(-pcna_cor), .(symbol_2011, symbol_now, entrez, pcna_cor)], row.names = FALSE)
}
dead <- d[status == "entrez 已退役"]
if (nrow(dead)) {
  cat("\n  !! Entrez ID 已退役（保留 2011 符号，标记待人工核对）：\n")
  print(dead[, .(symbol_2011, entrez, pcna_cor)], row.names = FALSE)
  d[status == "entrez 已退役", symbol_now := symbol_2011]
}

d[, symbol := symbol_now]
kv("最终唯一符号数", length(unique(d$symbol)))
if (anyDuplicated(d$symbol)) {
  cat("  !! 归一化后出现重复符号：",
      paste(unique(d$symbol[duplicated(d$symbol)]), collapse = ", "), "\n")
}

# ============================================================ 3. 逐队列覆盖率
hr("3. 逐队列覆盖率")
sp <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/scripts/scoring_protocol.R")
if (file.exists(sp)) {
  suppressMessages(source(sp))
  for (a in COHORTS) {
    co <- tryCatch(suppressMessages(load_cohort(a, verbose = FALSE)), error = function(e) NULL)
    if (is.null(co)) next
    cv <- gene_coverage(d$symbol, co$expr)
    cat(sprintf("  %-12s %3d/%3d  %5.1f%%  回退 %d  %s\n", a, cv$n_present, cv$n_total,
        100 * cv$frac, cv$n_rescued,
        ifelse(cv$frac >= 0.70, "可用", "!! 覆盖率不足，该队列的 M3 不可算")))
  }
  cat("\n  注：meta-PCNA 是混杂读数，不是被测 signature，故不适用 SAP §4.2 的\n")
  cat("      70% 阈值判定；但覆盖率过低会使增殖评分不可靠，须逐队列报告。\n")
} else cat("  未找到 scoring_protocol.R，跳过\n")

# ============================================================ 4. 与被测 signature 的重叠
hr("4. 与被测 signature 的基因重叠（SAP §7.2：>30% 须标注结构性重叠）")
gd <- path.expand("~/thesis_LUAD_CAF/03_signature_registry/genes")
fs <- sort(list.files(gd, pattern = "^SIG-[0-9]+_genes\\.txt$", full.names = TRUE))
if (length(fs)) {
  ov <- rbindlist(lapply(fs, function(f) {
    sid <- sub("_genes\\.txt$", "", basename(f))
    g <- trimws(readLines(f)); g <- g[nzchar(g)]
    i <- intersect(g, d$symbol)
    data.table(sig_id = sid, n_genes = length(g), overlap = length(i),
               frac = round(length(i) / length(g), 3),
               genes = paste(i, collapse = ","))
  }))
  print(ov[order(-frac)], row.names = FALSE)
  bad <- ov[frac > 0.30]
  if (nrow(bad)) {
    cat("\n  ★ 以下组合存在结构性重叠（>30%），M3 结果须单独讨论、不计入主汇总：\n")
    print(bad[, .(sig_id, n_genes, overlap, frac, genes)], row.names = FALSE)
  } else cat("\n  无 signature 与 meta-PCNA 的重叠超过 30%\n")
} else cat("  未找到 signature 基因文件\n")

# ============================================================ 5. 落盘
hr("5. 输出")
if (!WRITE) {
  cat("  未加 --write，本次只报告。确认无误后重跑并加 --write。\n")
} else {
  dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
  fg <- file.path(OUTD, "metaPCNA_genes.txt")
  writeLines(sort(unique(d$symbol)), fg)
  ft <- file.path(OUTD, "metaPCNA_table.csv")
  fwrite(d[order(-pcna_cor), .(symbol, symbol_2011, entrez, pcna_cor, status)], ft)
  fm <- file.path(OUTD, "metaPCNA_provenance.txt")
  writeLines(c(
    "meta-PCNA gene set",
    "",
    "Source: Venet D, Dumont JE, Detours V. Most Random Gene Expression Signatures",
    "        Are Significantly Associated with Breast Cancer Outcome.",
    "        PLoS Comput Biol. 2011;7(10):e1002240. Table S1.",
    "        doi:10.1371/journal.pcbi.1002240",
    paste0("Source PDF md5: ", unname(tools::md5sum(PDF))),
    paste0("Extracted:      ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    "Original definition (paper Methods):",
    "  Pearson correlation between PCNA and all genes in Ge et al. 2005 (GDS1096,",
    "  36 normal tissue samples spanning 27 organs); top 1% most positively",
    "  correlated retained, i.e. 131 of 13,077 genes.",
    "  meta-PCNA index = MEDIAN expression of these genes (not mean, not ssGSEA).",
    "",
    "  Defined on NORMAL tissue, therefore adjusting tumour data for it",
    "  does not introduce circularity.",
    "",
    paste0("Genes parsed:      ", nrow(d)),
    paste0("Renamed to current HGNC: ", sum(d$status == "已改名")),
    paste0("Retired Entrez IDs:      ", sum(d$status == "entrez 已退役")),
    paste0("Final unique symbols:    ", length(unique(d$symbol))),
    "",
    "NOTE ON ADJUSTMENT METHOD:",
    "  Venet et al. adjusted the expression matrix gene-by-gene (lm on the",
    "  meta-PCNA index, replacing each value by residual + cohort mean) and",
    "  then re-scored the signature.",
    "  This study follows SAP 7.1 instead: the proliferation index enters the",
    "  Cox model as a covariate (model M3). The SAP approach is more",
    "  conservative; attenuation reported here is therefore not directly",
    "  comparable to the near-complete abrogation reported by Venet et al."), fm)
  for (f in c(fg, ft, fm)) kv(basename(f), sprintf("%.1f KB", file.info(f)$size / 1024))
  kv("目录", OUTD)
}

hr("6. 待办")
cat("  1) SAP 附录 A：把 meta-PCNA 的【待核实】改为已核实，填入完整引用\n")
cat("     Venet D, Dumont JE, Detours V. PLoS Comput Biol. 2011;7(10):e1002240.\n")
cat("  2) SAP 附录 A：Venet 2011 的著者更正为 Venet, Dumont, Detours\n")
cat("     （此前记录的 'Venet, Dutoit & Michiels' 有误，勿使用）\n")
cat("  3) 增殖评分的计算方式按原文：131 个基因表达的【中位数】\n")
cat("     不用 ssGSEA —— 这与被测 signature 的打分方式不同，须在 Methods 中写明\n")
cat("  4) 阶段 2 的另一前置：TCGA 的 ABSOLUTE 纯度仍未获取\n")

cat("\n输出: ", OUT, "\n", sep = "")
