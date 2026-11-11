#!/usr/bin/env Rscript
# ==============================================================================
# sci2_stage3_facts.R —— 阶段 3 事实清单（只读已写出的 stage3_* 文件，全部数字现算）
#   输出：data_derived/stage3_facts_<日期>.md
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_stage3_facts.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
CH <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived")
STAMP <- format(Sys.Date(), "%Y%m%d")
latest <- function(pat) { f <- sort(list.files(DD, pat, full.names = TRUE), decreasing = TRUE); stopifnot(length(f) >= 1); f[1] }
FL <- latest("^stage3_readouts_long_[0-9]{8}\\.csv$"); FC <- latest("^stage3_cox_[0-9]{8}\\.csv$")
FK <- latest("^stage3_consistency_[0-9]{8}\\.csv$"); FN <- latest("^stage3_readout_notes_[0-9]{8}\\.csv$")
L <- fread(FL); CX <- fread(FC); K <- fread(FK)
stopifnot(!anyDuplicated(L[, .(cohort, sample_id, readout)]), !anyDuplicated(CX[, .(cohort, readout)]))
COH <- sort(unique(L$cohort)); stopifnot(length(COH) == 6)

f2 <- function(x) formatC(x, format = "f", digits = 2)
rng <- function(x) { x <- x[is.finite(x)]; paste0(f2(min(x)), " 至 ", f2(max(x))) }
rho <- function(a, x, y) { W <- dcast(L[cohort == a & readout %in% c(x, y)], sample_id ~ readout, value.var = "value")
  if (!all(c(x, y) %in% names(W))) return(NA_real_); cor(W[[x]], W[[y]], method = "spearman", use = "pairwise.complete.obs") }
rr <- function(x, y, coh = COH) vapply(coh, function(a) rho(a, x, y), 0)
zq <- function(r) CX[readout == r][match(COH, cohort), z]

SIGNAL <- sprintf("SIG-%03d", c(1, 4, 6, 7, 8, 11, 12)); STROMA <- sprintf("SIG-%03d", c(2, 3, 9, 10))
r_mx <- rr("MCP_Fibroblasts", "xCell_Fibroblasts"); r_m2 <- rr("MCP_Fibroblasts", "SIG-002")
r_mp <- rr("MCP_Fibroblasts", "Prolif_metaPCNA"); r_xp <- rr("xCell_Fibroblasts", "Prolif_metaPCNA")
r_xs <- rr("xCell_Fibroblasts", "xCell_StromaScore")
z_m <- zq("MCP_Fibroblasts"); z_x <- zq("xCell_Fibroblasts"); z_p <- zq("Prolif_metaPCNA")
st_mcp <- unlist(lapply(STROMA, function(s) rr(s, "MCP_Fibroblasts")))
st_pro <- unlist(lapply(STROMA, function(s) rr(s, "Prolif_metaPCNA")))
sg_pro <- unlist(lapply(SIGNAL, function(s) rr(s, "Prolif_metaPCNA")))
st_pur <- vapply(STROMA, function(s) rho("TCGA_LUAD", s, "Purity_ABSOLUTE"), 0)
nsig <- function(r) CX[readout == r, sum(abs(z) >= 1.96, na.rm = TRUE)]
cmp <- CX[grepl("^SIG-", readout) & readout %in% SIGNAL][
  CX[readout == "Prolif_metaPCNA", .(cohort, z_prolif = z)], on = "cohort"][, above := abs(z) > abs(z_prolif)]
geo <- cmp[cohort != "TCGA_LUAD" & !is.na(z)]; tc <- cmp[cohort == "TCGA_LUAD"]

T <- c(
  sprintf("# 阶段 3 事实清单（%s，sci2_stage3_facts.R 现算）", STAMP), "",
  "来源文件：", paste0("- ", basename(c(FL, FC, FK, FN)), "（md5 ", unname(tools::md5sum(c(FL, FC, FK, FN))), "）"), "",
  "## G1 相关结构的跨队列一致性（共有读数，15 对队列）",
  sprintf("- cophenetic 相关中位 %s（%s）；相关矩阵上三角 Spearman 中位 %s（%s）。最低的 3 对：%s。",
          f2(median(K$cophenetic)), rng(K$cophenetic), f2(median(K$cor_matrix_spearman)), rng(K$cor_matrix_spearman),
          paste(K[order(cophenetic)][1:3, sprintf("%s–%s %s", cohort_1, cohort_2, f2(cophenetic))], collapse = "；")),
  "- 不可写：『签名分为 X 类』（未设类数 k）。叶序相邻不等于同类。", "",
  "## G2 合理性检查",
  sprintf("- SIG-002（MCP-counter Fibroblasts 基因）与 MCP_Fibroblasts 的 Spearman：%s。", rng(r_m2)),
  sprintf("- TCGA：MCP_Fibroblasts 与 EPIC_CAFs %s；xCell_Fibroblasts 与 EPIC_CAFs %s。",
          f2(rho("TCGA_LUAD", "MCP_Fibroblasts", "EPIC_CAFs")), f2(rho("TCGA_LUAD", "xCell_Fibroblasts", "EPIC_CAFs"))), "",
  "## G3 两个『成纤维细胞』读数",
  sprintf("- MCP_Fibroblasts 与 xCell_Fibroblasts 的 Spearman：%s（%s）。", rng(r_mx), paste(sprintf("%s %s", COH, f2(r_mx)), collapse = "，")),
  sprintf("- 与增殖指数：MCP_Fibroblasts %s；xCell_Fibroblasts %s（为负的队列 %d/6）。xCell_Fibroblasts 与 xCell_StromaScore %s。",
          rng(r_mp), rng(r_xp), sum(r_xp < 0), rng(r_xs)),
  sprintf("- 单变量 Cox z：MCP_Fibroblasts %s（>0 的队列 %d/6，|z|≥1.96 的 %d）；xCell_Fibroblasts %s（<0 的队列 %d/6，|z|≥1.96 的 %d）。",
          rng(z_m), sum(z_m > 0), nsig("MCP_Fibroblasts"), rng(z_x), sum(z_x < 0), nsig("xCell_Fibroblasts")),
  "- 可写：两种成纤维细胞读数中度相关，但与增殖的相关方向不同，预后关联方向相反。",
  "- 不可写：『成纤维细胞具有保护 / 促进作用』；『xCell 的负向关联由增殖导致』（未做校正检验）。", "",
  "## G4 基质含量类签名（SIG-002/003/009/010）",
  sprintf("- 与 MCP_Fibroblasts 的 Spearman %s；与增殖指数 %s；TCGA 中与 ABSOLUTE 纯度 %s。", rng(st_mcp), rng(st_pro), rng(st_pur)),
  sprintf("- 6 队列中 |z|≥1.96 的队列数：%s。", paste(sprintf("%s %d", STROMA, vapply(STROMA, nsig, 0L)), collapse = "，")),
  "- SIG-003 与 estimate 包 StromalSignature 共有 136/141 基因（阶段 3 日志第 0 节）。", "",
  "## G5 有信号签名与增殖指数（阶段 2 F3 的阶段 3 复核）",
  sprintf("- SIG-001/004/006/007/008/011/012 与增殖指数的 Spearman（6 队列）：%s。", rng(sg_pro)), "",
  "## G6 预后横向比较：增殖指数作为参照读数",
  sprintf("- 增殖指数单变量 z：%s，|z|≥1.96 的队列 %d/6。", rng(z_p), nsig("Prolif_metaPCNA")),
  sprintf("- GEO 5 队列中，上述 7 个签名的可评估组合 %d 个，|z| 超过同队列增殖指数的 %d 个。", nrow(geo), sum(geo$above)),
  sprintf("- TCGA 中 7 个签名 |z| 超过增殖指数（%s）的 %d 个：%s。TCGA 是 SIG-001/004/006/011/012 的开发队列。",
          f2(tc$z_prolif[1]), sum(tc$above), paste(tc[above == TRUE, sprintf("%s %s", readout, f2(z))], collapse = "，")),
  "- 不可写：『签名不优于 / 优于增殖指数』（SAP 未预设 z 的比较检验，仅描述）。"
)
fo <- file.path(DD, sprintf("stage3_facts_%s.md", STAMP)); if (file.exists(fo)) stop("已存在，不覆盖：", fo)
writeLines(T, fo); cat(T, sep = "\n"); cat("\n写出 ", fo, "\n", sep = "")
