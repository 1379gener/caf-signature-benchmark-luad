#!/usr/bin/env Rscript
# ==============================================================================
# sci2_stage3_redundancy.R —— 阶段 3：冗余与相关性分析（SAP v1.2 §8）
#
#   读数（2026-10-17 登记）
#   - 12 个被测签名：与阶段 1/2 同一打分路径（候选池 + score_sets_*），以阶段 1 z_real 为闸门
#   - 外部读数（6 队列）：MCP-counter 全部 10 类；xCell 的 Fibroblasts / ImmuneScore / StromaScore /
#     MicroenvironmentScore；ESTIMATE ImmuneSignature（estimate 包 SI_geneset，走同一 ssGSEA 路径）；
#     meta-PCNA 增殖指数
#   - 仅 TCGA：EPIC 的 CAFs / Endothelial（log2(TPM+1) 还原为 TPM 后输入）；ABSOLUTE 纯度
#   分析
#   1 各队列内全部读数两两 Spearman
#   2 层次聚类：距离 1 − |rho|，ward.D2
#   3 跨队列一致性：两两队列在共有读数上（a）重新聚类后的 cophenetic 相关（dendextend::cor_cophenetic），
#     （b）相关矩阵上三角的 Spearman；报中位与范围。TCGA 独有读数不进入跨队列比较
#   4 预后横向比较：每个读数在每个队列的单变量 Cox（per SD，同一终点、同一有效样本定义）
#   输出：data_derived/stage3_*_<日期>.*；日志写 logs/
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_stage3_redundancy.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA)
  library(MCPcounter); library(xCell); library(EPIC); library(estimate); library(dendextend) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD <- file.path(ROOT, "03_signature_registry"); DD <- file.path(CH, "data_derived")
suppressMessages({ source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
                   source(file.path(CH, "scripts/null_engine.R")) })
STAMP <- format(Sys.Date(), "%Y%m%d")
LOG <- file.path(CH, "logs", format(Sys.time(), "stage3_redundancy_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
OUTF <- function(stem, ext = "csv") { f <- file.path(DD, sprintf("stage3_%s_%s.%s", stem, STAMP, ext))
  if (file.exists(f)) stop("已存在，不覆盖：", f); f }
latest <- function(pat) { f <- sort(list.files(DD, pat, full.names = TRUE), decreasing = TRUE); stopifnot(length(f) >= 1); f[1] }
for (p in c("MCPcounter", "xCell", "EPIC", "estimate", "dendextend", "GSVA"))
  cat(sprintf("  %-11s %s\n", p, as.character(packageVersion(p))))

LOCK <- require_locked_registry(REGD); stopifnot(LOCK$md5 == "4ad3218931d95dfb6ee67a7ba07b9a75"); reg <- LOCK$table
PRO <- fread(latest("^m3_prolif_metaPCNA_[0-9]{8}\\.csv$"))
PUR <- fread(latest("^m2_purity_TCGA_[0-9]{8}\\.csv$"), na.strings = c("", "NA"))
read_sig <- function(sid, cls) {
  g <- trimws(readLines(file.path(REGD, "genes", paste0(sid, "_genes.txt")))); g <- g[nzchar(g)]
  w <- NULL
  if (cls == "B") { d <- fread(file.path(REGD, "genes", paste0(sid, "_weights.csv")))
                    w <- setNames(as.numeric(d$weight), d$gene)[g]; stopifnot(!anyNA(w)) }
  list(genes = g, weights = w) }

hr("0. ESTIMATE 基因集")
data(SI_geneset, package = "estimate", envir = environment())
gs_row <- function(r) { v <- as.character(unlist(SI_geneset[r, ])); v <- v[nzchar(v) & !is.na(v)]; v[grepl("^[A-Za-z0-9.-]+$", v) & v != r] }
EST_IMM <- unique(gs_row("ImmuneSignature")); EST_STR <- unique(gs_row("StromalSignature"))
s3 <- trimws(readLines(file.path(REGD, "genes", "SIG-003_genes.txt"))); s3 <- s3[nzchar(s3)]
cat(sprintf("  ImmuneSignature %d 基因；StromalSignature %d 基因，与 SIG-003 共有 %d / %d\n",
            length(EST_IMM), length(EST_STR), length(intersect(EST_STR, s3)), length(s3)))
stopifnot(length(EST_IMM) >= 130)

XC_USE  <- c("Fibroblasts", "ImmuneScore", "StromaScore", "MicroenvironmentScore")
EPIC_USE <- c("CAFs", "Endothelial")

hr("1. 逐队列计算全部读数")
rds <- list.files(file.path(DD, "null_benchmark"), "\\.rds$", full.names = TRUE)
EV <- rbindlist(lapply(rds, function(f) { x <- readRDS(f)
  data.table(cohort = x$cohort, sig_id = x$sig_id, class = x$class, evaluable = isTRUE(x$evaluable),
             z_real = if (is.null(x$z_real)) NA_real_ else x$z_real) }))
LONG <- list(); SURV <- list(); NOTE <- list()
for (a in setdiff(COHORTS, "CPTAC_LUAD")) {
  co <- load_cohort(a, verbose = FALSE); X <- co$expr
  stopifnot(!anyDuplicated(colnames(X)), ncol(X) == co$n)
  id <- colnames(X); pool <- build_pool(X)$genes
  R <- list()
  for (i in which(EV$cohort == a & EV$evaluable)) {                           # 签名
    s <- EV$sig_id[i]; cls <- EV$class[i]; sg <- read_sig(s, cls); gp <- intersect(sg$genes, pool)
    sc <- suppressWarnings(if (cls == "A") score_sets_ssgsea(X, list(gp), prefix = "real")[1, ] else
                                           score_sets_weighted(X, list(gp), sg$weights[gp])[1, ])
    z <- cox_one(sc, co$clin$OS_time, co$clin$OS_event)$z; stopifnot(abs(z - EV$z_real[i]) < 1e-6)
    R[[s]] <- as.numeric(sc) }
  gi <- intersect(EST_IMM, pool)                                              # ESTIMATE 免疫
  R[["ESTIMATE_Immune"]] <- as.numeric(suppressWarnings(score_sets_ssgsea(X, list(gi), prefix = "est")[1, ]))
  NOTE[[length(NOTE) + 1]] <- data.table(cohort = a, readout = "ESTIMATE_Immune", n_genes = length(gi), note = "")
  mc <- MCPcounter.estimate(X, featuresType = "HUGO_symbols")                 # MCP-counter
  for (r in rownames(mc)) R[[paste0("MCP_", gsub("[^A-Za-z0-9]+", "_", r))]] <- as.numeric(mc[r, id])
  xc <- suppressWarnings(xCellAnalysis(X, rnaseq = (a == "TCGA_LUAD"), parallel.sz = 1))   # xCell
  stopifnot(all(XC_USE %in% rownames(xc)))
  for (r in XC_USE) R[[paste0("xCell_", r)]] <- as.numeric(xc[r, id])
  R[["Prolif_metaPCNA"]] <- PRO[cohort == a][match(id, sample_id), prolif_metaPCNA]
  if (a == "TCGA_LUAD") {
    w <- character(0)
    ep <- withCallingHandlers(EPIC(bulk = 2^X - 1), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    for (r in EPIC_USE) R[[paste0("EPIC_", r)]] <- as.numeric(ep$cellFractions[id, r])
    NOTE[[length(NOTE) + 1]] <- data.table(cohort = a, readout = "EPIC", n_genes = NA_integer_, note = paste(unique(w), collapse = " | "))
    R[["Purity_ABSOLUTE"]] <- PUR[match(id, sample_id), purity_absolute]
  }
  M <- as.data.table(R); M[, sample_id := id]
  LONG[[a]] <- melt(M, id.vars = "sample_id", variable.name = "readout", value.name = "value")[, cohort := a]
  SURV[[a]] <- data.table(cohort = a, sample_id = id, t = co$clin$OS_time, e = co$clin$OS_event)
  cat(sprintf("  %-10s 样本 %d，读数 %d\n", a, length(id), length(R)))
}
L <- rbindlist(LONG); S <- rbindlist(SURV)
fwrite(L, OUTF("readouts_long")); fwrite(rbindlist(NOTE), OUTF("readout_notes"))

hr("2. 各队列相关矩阵与层次聚类")
COR <- list(); HC <- list()
for (a in names(LONG)) {
  W <- dcast(L[cohort == a], sample_id ~ readout, value.var = "value"); m <- as.matrix(W[, -1])
  v <- apply(m, 2, function(x) sd(x, na.rm = TRUE)); drop <- names(v)[!is.finite(v) | v == 0]
  if (length(drop)) cat("  ", a, " 无变异、剔除：", paste(drop, collapse = ", "), "\n", sep = "")
  m <- m[, setdiff(colnames(m), drop), drop = FALSE]
  rho <- cor(m, method = "spearman", use = "pairwise.complete.obs"); COR[[a]] <- rho
  HC[[a]] <- hclust(as.dist(1 - abs(rho)), method = "ward.D2")
  fwrite(as.data.table(rho, keep.rownames = "readout"), OUTF(paste0("cor_", a)))
  cat(sprintf("\n  --- %s：聚类叶序\n  %s\n", a, paste(HC[[a]]$labels[HC[[a]]$order], collapse = " | ")))
}
saveRDS(list(cor = COR, hclust = HC), OUTF("cor_hclust", "rds"))

hr("3. 跨队列一致性（共有读数，TCGA 独有读数除外）")
ch <- names(COR); P <- list()
for (i in seq_along(ch)) for (j in seq_along(ch)) if (i < j) {
  com <- setdiff(intersect(colnames(COR[[ch[i]]]), colnames(COR[[ch[j]]])), c(paste0("EPIC_", EPIC_USE), "Purity_ABSOLUTE"))
  r1 <- COR[[ch[i]]][com, com]; r2 <- COR[[ch[j]]][com, com]
  h1 <- hclust(as.dist(1 - abs(r1)), "ward.D2"); h2 <- hclust(as.dist(1 - abs(r2)), "ward.D2")
  ut <- upper.tri(r1)
  P[[length(P) + 1]] <- data.table(cohort_1 = ch[i], cohort_2 = ch[j], n_readouts = length(com),
    cophenetic = cor_cophenetic(as.dendrogram(h1), as.dendrogram(h2)),
    cor_matrix_spearman = cor(r1[ut], r2[ut], method = "spearman"))
}
P <- rbindlist(P); fwrite(P, OUTF("consistency"))
print(P[, .(cohort_1, cohort_2, n_readouts, cophenetic = round(cophenetic, 3), cor_matrix = round(cor_matrix_spearman, 3))], row.names = FALSE)
cat(sprintf("\n  cophenetic 中位 %.3f（%.3f 至 %.3f）；相关矩阵 Spearman 中位 %.3f（%.3f 至 %.3f）；队列对 %d\n",
    median(P$cophenetic), min(P$cophenetic), max(P$cophenetic),
    median(P$cor_matrix_spearman), min(P$cor_matrix_spearman), max(P$cor_matrix_spearman), nrow(P)))

hr("4. 预后横向比较（单变量 Cox，per SD）")
CX <- L[, { d <- merge(.SD, S[cohort == .BY$cohort], by = "sample_id")
            d <- d[!is.na(t) & t > 0 & !is.na(e) & is.finite(value)]
            if (nrow(d) < 20 || sd(d$value) == 0) list(n = nrow(d), events = sum(d$e), HR = NA_real_, lo = NA_real_, hi = NA_real_, z = NA_real_, p = NA_real_) else {
              r <- cox_one(d$value, d$t, d$e); list(n = r$n, events = r$events, HR = r$hr, lo = r$lo, hi = r$hi, z = r$z, p = r$p) } },
        by = .(cohort, readout)]
CX[, type := fifelse(grepl("^SIG-", readout), "被测签名", "外部读数")]
fwrite(CX, OUTF("cox"))
options(width = 200)
print(dcast(CX, readout + type ~ cohort, value.var = "z", fun.aggregate = function(x) if (length(x)) round(x[1], 2) else NA_real_, fill = NA)[order(type, readout)], row.names = FALSE)
print(sessionInfo())
sink(); close(zz)
