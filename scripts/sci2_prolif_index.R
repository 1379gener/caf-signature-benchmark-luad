#!/usr/bin/env Rscript
# ==============================================================================
# sci2_prolif_index.R —— M3/M4 的增殖指数（SAP v1.2 §7.2）
#   定义：Venet 2011 meta-PCNA 131 基因在各样本中表达值（log2，按 load_cohort 协议）的中位数
#   基因匹配：先按现行 HGNC 符号精确匹配；未命中者用 resolve_aliases() 的别名回退（与覆盖率口径同一函数）
#   另报：逐队列覆盖、与 MKI67/TOP2A 的 Spearman（描述性效度检查）、与 12 个被测签名的基因重叠
#   只读底座与注册表；写 data_derived/ 与 logs/
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_prolif_index.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD <- file.path(ROOT, "03_signature_registry")
suppressMessages(source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")))
STAMP <- format(Sys.Date(), "%Y%m%d")
LOG <- file.path(CH, "logs", format(Sys.time(), "prolif_index_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")

hr("1. meta-PCNA 基因集")
fg <- file.path(REGD, "confounders/metaPCNA_genes.txt")
cat("文件 ", fg, "\nmd5 ", unname(tools::md5sum(fg)), "\n", sep = "")
G <- trimws(readLines(fg, warn = FALSE)); G <- G[nzchar(G) & !grepl("^#", G)]
cat("条目 ", length(G), "；去重后 ", length(unique(G)), "\n", sep = "")
stopifnot(length(G) == 131, !anyDuplicated(G))
cat("\nprovenance：\n"); cat(readLines(file.path(REGD, "confounders/metaPCNA_provenance.txt"), warn = FALSE), sep = "\n")

hr("2. 与 12 个被测签名的基因重叠（SAP §7.2 第 5 条：>30% 标注结构性重叠）")
ov <- rbindlist(lapply(sprintf("SIG-%03d", 1:12), function(s) {
  g <- trimws(readLines(file.path(REGD, "genes", paste0(s, "_genes.txt")), warn = FALSE)); g <- unique(g[nzchar(g)])
  i <- intersect(g, G)
  data.table(sig_id = s, n_sig = length(g), n_overlap = length(i), frac = round(length(i) / length(g), 3),
             structural_overlap = length(i) / length(g) > 0.30, genes = paste(i, collapse = ";")) }))
print(ov, row.names = FALSE)

hr("3. 逐队列计算")
OUT <- list(); COV <- list()
for (a in COHORTS) {
  co <- load_cohort(a, verbose = FALSE); X <- co$expr
  if (!co$has_survival) { cat(a, "：无生存数据，不进入 M3/M4，跳过\n"); next }
  id <- co$id_col; stopifnot(!anyDuplicated(colnames(X)), ncol(X) == co$n)
  r <- resolve_aliases(G, rownames(X))
  exact <- intersect(G, rownames(X))
  used <- r$present                                        # 精确命中 + 唯一别名命中（矩阵中的名字）
  stopifnot(!anyDuplicated(used))
  px <- apply(X[used, , drop = FALSE], 2, median)
  stopifnot(all(is.finite(px)))
  rho <- function(g) if (g %in% rownames(X)) round(cor(px, X[g, ], method = "spearman"), 3) else NA_real_
  COV[[a]] <- data.table(cohort = a, n = co$n, n_exact = length(exact), n_alias = length(r$rescued),
                         n_used = length(used), coverage = round(length(used) / 131, 3),
                         missing = paste(r$missing, collapse = ";"),
                         alias_map = paste(sprintf("%s->%s", r$map$gene, r$map$used), collapse = ";"),
                         rho_MKI67 = rho("MKI67"), rho_TOP2A = rho("TOP2A"),
                         prolif_median = round(median(px), 3), prolif_iqr = round(IQR(px), 3))
  OUT[[a]] <- data.table(cohort = a, sample_id = colnames(X), prolif_metaPCNA = px)
}
C <- rbindlist(COV)
print(C[, !c("missing", "alias_map"), with = FALSE], row.names = FALSE)
cat("\n缺失基因与别名映射：\n")
for (i in seq_len(nrow(C))) cat(sprintf("  %-10s 缺: %s\n             别名: %s\n", C$cohort[i],
                                        ifelse(nzchar(C$missing[i]), C$missing[i], "无"),
                                        ifelse(nzchar(C$alias_map[i]), C$alias_map[i], "无")))

hr("4. 写出")
fo <- file.path(CH, "data_derived", sprintf("m3_prolif_metaPCNA_%s.csv", STAMP))
if (file.exists(fo)) stop("已存在，不覆盖：", fo)
fwrite(rbindlist(OUT), fo)
fwrite(C, sub("\\.csv$", "_coverage.csv", fo)); fwrite(ov, sub("\\.csv$", "_sig_overlap.csv", fo))
cat("写出：", fo, "（及 _coverage / _sig_overlap）\n", sep = "")
sink(); close(zz)
