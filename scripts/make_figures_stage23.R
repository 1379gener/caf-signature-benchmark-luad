#!/usr/bin/env Rscript
# ==============================================================================
# make_figures_stage23.R —— SCI-2 阶段 2 / 阶段 3 图（只读 data_derived/ 下已写出的文件，无字面量数据）
#   Fig5  阶段 2：A 有信号组合 M1 → M3 的 HR（同样本，per SD）；B TCGA 的 M1 / M2 / M4 / M4 敏感性
#   Fig6  阶段 3：A TCGA 全部读数相关矩阵（按 ward.D2 叶序）；B 两个成纤维细胞读数与增殖指数的 Spearman；
#         C 两读数单变量 vs 加增殖指数后的 z（事后探索性）
#   FigS  阶段 3：其余 5 个队列的相关矩阵，各一页
#   输出：figures/Fig5_*_<日期>.pdf/png、Fig6_*、FigS_cor_<队列>_*；同名不覆盖
#   拼版用 cowplot::plot_grid（服务器 ggplot2 4.0 与已装 patchwork 不兼容）
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/make_figures_stage23.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
CH <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived"); FIG <- file.path(CH, "figures")
STAMP <- format(Sys.Date(), "%Y%m%d"); dir.create(FIG, showWarnings = FALSE)
latest <- function(pat) { f <- sort(list.files(DD, pat, full.names = TRUE), decreasing = TRUE); stopifnot(length(f) >= 1)
  cat(sprintf("  %-45s md5 %s\n", basename(f[1]), unname(tools::md5sum(f[1])))); f[1] }
save2 <- function(p, stem, w, h) { for (ext in c("pdf", "png")) {
  f <- file.path(FIG, sprintf("%s_%s.%s", stem, STAMP, ext)); if (file.exists(f)) stop("已存在，不覆盖：", f)
  ggsave(f, p, width = w, height = h, dpi = 300, device = if (ext == "pdf") cairo_pdf else "png") }
  cat("  写出 ", stem, "\n", sep = "") }
TH <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
                                     plot.tag = element_text(face = "bold"))
cat("输入：\n")
MOD <- fread(latest("^stage2_models_[0-9]{8}\\.csv$"))
DIA <- fread(latest("^stage2_diagnostics_[0-9]{8}\\.csv$"))
CH3 <- readRDS(latest("^stage3_cor_hclust_[0-9]{8}\\.rds$"))
EXP <- fread(latest("^stage3_explore_fibro_prolif_[0-9]{8}\\.csv$"))
LNG <- fread(latest("^stage3_readouts_long_[0-9]{8}\\.csv$"))
stopifnot(!anyDuplicated(MOD[, .(cohort, sig_id, model, sample)]), !anyDuplicated(DIA[, .(cohort, sig_id, comparison)]),
          !anyDuplicated(EXP[, .(cohort, readout, model)]))
COH <- c("TCGA_LUAD", "GSE68465", "GSE72094", "GSE31210", "GSE50081", "GSE30219")
coh_lab <- function(x) factor(sub("_LUAD", "", x), levels = sub("_LUAD", "", COH))

# ------------------------------------------------------------------ Fig5A
sig3 <- DIA[comparison == "M3@S3" & M1_signal == TRUE, .(cohort, sig_id, VIF_M3)]
A <- merge(MOD[sample == "S3" & model %in% c("M1", "M3")], sig3, by = c("cohort", "sig_id"))
stopifnot(nrow(A) == 2 * nrow(sig3))
A[, `:=`(coh = coh_lab(cohort), model = factor(model, c("M1", "M3"), c("M1  score", "M3  score + proliferation")),
         hi_vif = VIF_M3 >= 2.5)]
A[, row := factor(paste(sig_id, coh), levels = rev(unique(A[order(sig_id, coh), paste(sig_id, coh)])))]
pA <- ggplot(A, aes(x = HR, y = row, colour = model)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(aes(group = row), colour = "grey70", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4, orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  geom_point(data = A[hi_vif == TRUE & model == levels(model)[1]], aes(x = 0.3, y = row), shape = 8, colour = "black", size = 1.2, inherit.aes = FALSE) +
  scale_x_log10() + scale_colour_manual(values = c("grey20", "#C2410C"), name = NULL) +
  labs(x = "Hazard ratio per SD (log scale)", y = NULL,
       subtitle = sprintf("%d pairs with |z| >= 1.96 under M1 (same sample)\n* VIF >= 2.5 (score vs proliferation index)", nrow(sig3))) +
  TH + theme(legend.position = "top", axis.text.y = element_text(size = 6))

# ------------------------------------------------------------------ Fig5B
sigT <- DIA[cohort == "TCGA_LUAD" & comparison == "M2@S2" & M1_signal == TRUE, sig_id]
B <- MOD[cohort == "TCGA_LUAD" & sig_id %in% sigT &
         ((model == "M1" & sample == "full") | (model == "M2" & sample == "S2") | (model == "M4" & sample %in% c("S4", "S4_sens")))]
B[, lab := factor(fcase(model == "M1", "M1 (all samples)", model == "M2", "M2 + purity",
                        sample == "S4", "M4 full adj.", default = "M4 full adj. (stage-supplemented)"),
                  levels = c("M1 (all samples)", "M2 + purity", "M4 full adj.", "M4 full adj. (stage-supplemented)"))]
pB <- ggplot(B, aes(x = HR, y = lab)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4, orientation = "y") + geom_point(size = 1.6) +
  facet_wrap(~ sig_id, ncol = 1, strip.position = "left") + scale_x_log10() + scale_y_discrete(limits = rev) +
  labs(x = "Hazard ratio per SD (log scale)", y = NULL, subtitle = "TCGA-LUAD\npurity: ABSOLUTE (Aran 2015)") +
  TH + theme(strip.text.y.left = element_text(angle = 0), axis.text.y = element_text(size = 6))
nB <- unique(B[, .(lab, n, events)])[order(lab)]
pB <- pB + labs(caption = paste(sprintf("%s: n = %d, events = %d", nB$lab, nB$n, nB$events), collapse = "\n"))
save2(plot_grid(pA, pB, labels = c("A", "B"), nrow = 1, rel_widths = c(1.35, 1)), "Fig5_stage2_confounding", 10, 9)

# ------------------------------------------------------------------ Fig6 / FigS
heat <- function(a) {
  r <- CH3$cor[[a]]; h <- CH3$hclust[[a]]; o <- h$labels[h$order]; stopifnot(setequal(o, rownames(r)))
  D <- as.data.table(as.table(r[o, o])); setnames(D, c("x", "y", "rho"))
  D[, `:=`(x = factor(x, o), y = factor(y, rev(o)))]
  ggplot(D, aes(x, y, fill = rho)) + geom_tile() +
    scale_fill_gradient2(low = "#2563EB", mid = "white", high = "#C2410C", limits = c(-1, 1), name = "Spearman") +
    coord_fixed() + labs(x = NULL, y = NULL, subtitle = sprintf("%s (n = %d); order: ward.D2 on 1 - |rho|", sub("_LUAD", "", a),
                                                              uniqueN(LNG[cohort == a, sample_id]))) +
    TH + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5), axis.text.y = element_text(size = 5.5))
}
rp <- LNG[readout %in% c("MCP_Fibroblasts", "xCell_Fibroblasts", "Prolif_metaPCNA")]
rp <- dcast(rp, cohort + sample_id ~ readout, value.var = "value")
R6 <- melt(rp[, .(MCP_Fibroblasts = cor(MCP_Fibroblasts, Prolif_metaPCNA, method = "spearman"),
                  xCell_Fibroblasts = cor(xCell_Fibroblasts, Prolif_metaPCNA, method = "spearman")), by = cohort],
           id.vars = "cohort", variable.name = "readout", value.name = "rho")
R6[, coh := coh_lab(cohort)]
pBf <- ggplot(R6, aes(x = rho, y = coh, colour = readout, shape = readout)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2) + scale_y_discrete(limits = rev) + scale_colour_manual(values = c("grey20", "#0F766E"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(x = "Spearman with meta-PCNA proliferation index", y = NULL) + TH + theme(legend.position = "top")
E <- copy(EXP)[, `:=`(coh = coh_lab(cohort), model = factor(model, c("univariable", "+prolif"), c("univariable", "+ proliferation index")))]
pC <- ggplot(E, aes(x = z, y = coh, colour = model)) +
  annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, colour = "grey50") + geom_line(aes(group = coh), colour = "grey70") + geom_point(size = 1.8) +
  facet_wrap(~ readout, nrow = 1) + scale_y_discrete(limits = rev) + scale_colour_manual(values = c("grey20", "#C2410C"), name = NULL) +
  labs(x = "Cox z (per SD); shaded |z| < 1.96", y = NULL, subtitle = "Post hoc exploratory analysis (registered 2026-10-25)") +
  TH + theme(legend.position = "top")
right <- plot_grid(pBf, pC, labels = c("B", "C"), ncol = 1, rel_heights = c(1, 1.1))
save2(plot_grid(heat("TCGA_LUAD"), right, labels = c("A", ""), nrow = 1, rel_widths = c(1.5, 1)), "Fig6_stage3_readouts", 13, 7.5)
for (a in setdiff(COH, "TCGA_LUAD")) save2(heat(a), sprintf("FigS_cor_%s", a), 7.5, 7)
print(sessionInfo())
