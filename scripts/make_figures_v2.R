#!/usr/bin/env Rscript
# ==============================================================================
# make_figures_v1.R —— SCI-2 稿件 v1 的 Figure 1–4（按 SAP 主口径：BH、A/B 分开、含 SCANI）
#   原 make_figures.R 保留不动。本脚本不内嵌任何数值，全部从 data_derived/ 读取。
#   输入：classA_scani_bh_20261025.csv（A 类 36 组合，含 SCANI，BH 家族 = 6）
#         sap6_bh_combo_20261025.csv（取 B 类 41 组合）
#         posthoc_checks_20261025.csv（过拟合比 z/HR 版；类别 × 角色通过数与 95% CI）
#         results_extended_*.csv（队列事件数，用于排序）；锁定表（基因数）；scani_supp_summary_20261025.csv
#   输出：figures/Fig1–4_*_v1_<日期>.pdf/png（同名不覆盖）
#   用法：Rscript make_figures_v1.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
DD <- file.path(CH, "data_derived"); FIG <- file.path(CH, "figures"); REG <- file.path(ROOT, "03_signature_registry")
ST <- format(Sys.Date(), "%Y%m%d")
NAMES <- c("Fig1_signal_heatmap", "Fig2_real_vs_null", "Fig3_pass_by_role", "Fig4_overfit_gradient")
for (nm in NAMES) for (ext in c("pdf", "png")) {
  f <- file.path(FIG, sprintf("%s_v2_%s.%s", nm, ST, ext)); if (file.exists(f)) stop("已存在，不覆盖：", f) }
save2 <- function(p, nm, w, h) for (ext in c("pdf", "png"))
  ggsave(file.path(FIG, sprintf("%s_v2_%s.%s", nm, ST, ext)), p, width = w, height = h, dpi = 300, device = if (ext == "pdf") cairo_pdf else ext)
md <- function(f) { cat(sprintf("  %-40s md5 %s\n", basename(f), unname(tools::md5sum(f)))); f }

A  <- fread(md(file.path(DD, "classA_scani_bh_20261025.csv")))
BB <- fread(md(file.path(DD, "sap6_bh_combo_20261025.csv")))
P  <- fread(md(file.path(DD, "posthoc_checks_20261025.csv")))
S  <- fread(md(file.path(DD, "scani_supp_summary_20261025.csv")))
E  <- fread(md(sort(list.files(DD, "^results_extended_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1]))[evaluable == TRUE]
L  <- fread(md(sort(list.files(REG, "^signature_registry_locked_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1]))
stopifnot(nrow(A) == 36L, !anyDuplicated(A[, paste(cohort, sig_id)]), BB[class == "B", .N] == 41L)

D <- rbind(A[, .(cohort, sig_id, class = "A", role, z_real, zA_p95, cls3_bh)],
           BB[class == "B", .(cohort, sig_id, class, role, z_real, zA_p95, cls3_bh)])
stopifnot(nrow(D) == 77L, !anyNA(D$cls3_bh))
NG <- rbind(L[, .(sig_id, n_genes)], S[, .(sig_id, n_genes)][1])
D <- merge(D, unique(NG), by = "sig_id", all.x = TRUE); stopifnot(!anyNA(D$n_genes))
D[, sig_lab := fifelse(sig_id == "SCANI_myCAF", sprintf("SCANI myCAF (n=%d) *", n_genes), sprintf("%s (n=%d)", sig_id, n_genes))]
D[, sig_lab := factor(sig_lab, levels = unique(D[order(class, sig_id == "SCANI_myCAF", sig_id), sig_lab]))]
EV <- unique(E[, .(cohort, n_events)]); stopifnot(!anyDuplicated(EV$cohort))
D[, cohort := factor(cohort, levels = rev(EV[order(-n_events), cohort]))]
D[, cls := factor(cls3_bh, levels = c("no_signal", "not_specific", "stronger"),
                  labels = c("No prognostic signal", "Signal, not stronger than random", "Stronger than random (q < 0.05)"))]
D[, class_lab := factor(paste("Class", class), levels = c("Class A", "Class B"))]
PAL <- setNames(c("#BDBDBD", "#F6C667", "#C0492B"), levels(D$cls))
TH <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
                                       legend.position = "bottom", plot.title = element_text(face = "bold", size = 11),
                                       plot.subtitle = element_text(size = 8.5, colour = "grey30"))

# ---------------------------------------------------------------- Fig 1
f1 <- ggplot(D, aes(sig_lab, cohort, fill = cls)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_point(data = D[role == "开发队列"], shape = 4, size = 2.6, stroke = 1.1, inherit.aes = FALSE, aes(sig_lab, cohort)) +
  geom_point(data = D[role == "原文验证队列"], shape = 1, size = 2.8, stroke = 1.1, inherit.aes = FALSE, aes(sig_lab, cohort)) +
  facet_grid(. ~ class_lab, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = PAL, name = NULL, drop = FALSE) +
  labs(title = "Prognostic signal and its specificity",
       subtitle = "False-discovery-rate control within cohort and class. \u00d7 development cohort; \u25cb validation cohort of the original study; * developed by the authors.\nEmpty cell = not evaluable (coverage < 70%).",
       x = NULL, y = NULL) + TH + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill = guide_legend(nrow = 1))
save2(f1, "Fig1_signal_heatmap", 10, 4.6)

# ---------------------------------------------------------------- Fig 2
mx <- max(c(abs(D$z_real), D$zA_p95)) * 1.05
f2 <- ggplot(D, aes(zA_p95, abs(z_real), colour = cls)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") + geom_hline(yintercept = 1.96, linetype = 3, colour = "grey55") +
  geom_point(aes(shape = cohort == "TCGA_LUAD"), size = 2, alpha = 0.85) + facet_wrap(~ class_lab) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("Other cohorts", "TCGA-LUAD"), name = NULL) + coord_equal(xlim = c(0, mx), ylim = c(0, mx)) +
  scale_colour_manual(values = PAL, name = NULL, drop = FALSE) +
  labs(title = "Each signature against its own null distribution",
       subtitle = "Diagonal: |z| equal to the 95th percentile of |z| in 1,000 expression-matched random gene sets. Dotted: |z| = 1.96.",
       x = "95th percentile of |z| in null distribution A", y = "|z| of the signature") + TH + guides(colour = guide_legend(nrow = 1))
save2(f2, "Fig2_real_vs_null", 9, 5.2)

# ---------------------------------------------------------------- Fig 3
RC <- P[part == "role_ci"]; stopifnot(nrow(RC) == 6L)
RL <- c("开发队列" = "Development", "原文验证队列" = "Validation (original study)", "完全独立队列" = "Fully independent")
R3 <- rbind(RC[, .(class, role, n, k = raw, lo = raw_lo, hi = raw_hi, adj = "Unadjusted")],
            RC[, .(class, role, n, k = bh,  lo = bh_lo,  hi = bh_hi,  adj = "FDR-adjusted")])
R3[, `:=`(rate = k / n, role_en = factor(RL[role], levels = rev(unname(RL))), class_lab = paste("Class", class),
          adj = factor(adj, levels = c("Unadjusted", "FDR-adjusted")))]
pd <- position_dodge(width = 0.55)
f3 <- ggplot(R3, aes(rate, role_en, colour = adj)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, position = pd, linewidth = 0.7) +
  geom_point(size = 2.4, position = pd) +
  geom_text(aes(x = hi, label = sprintf("%d/%d", k, n)), position = pd, hjust = -0.6, size = 2.9, show.legend = FALSE) +
  facet_wrap(~ class_lab) + scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.25), labels = scales::percent) +
  scale_colour_manual(values = setNames(c("grey55", "#C0492B"), levels(R3$adj)), name = NULL) +
  labs(title = "Passes against null distribution A by cohort role (post hoc)",
       subtitle = "Points: proportion of combinations passing; bars: exact binomial 95% CI.", x = "Proportion passing", y = NULL) + TH
save2(f3, "Fig3_pass_by_role", 9, 3.6)

# ---------------------------------------------------------------- Fig 4
G <- P[part == "overfit" & sig_id != "SCANI_myCAF"]; stopifnot(nrow(G) == 6L)
S4x <- fread(md(file.path(DD, "supp_tables_20261025", "S4_all_combinations.csv")))
mk <- function(dev, ext, ratio, ylab, ttl) {
  Gl <- rbind(G[, .(sig_id, where = "Development (TCGA-LUAD)", v = get(dev))], G[, .(sig_id, where = "External (median)", v = get(ext))])
  Gl[, where := factor(where, levels = c("Development (TCGA-LUAD)", "External (median)"))]
  lab <- G[, .(sig_id, where = factor("Development (TCGA-LUAD)", levels = levels(Gl$where)), v = get(dev), r = get(ratio))]
  X <- S4x[sig_id %in% G$sig_id & cohort != "TCGA_LUAD"]
  X[, v := if (dev == "z_dev") abs(z) else abs(log(hr_per_sd))]
  chk <- merge(X[, .(m = median(v)), by = sig_id], G[, .(sig_id, e = get(ext))], by = "sig_id"); stopifnot(max(abs(chk$m - chk$e)) < 1e-6)
  X[, where := factor("External (median)", levels = levels(Gl$where))]
  ggplot(Gl, aes(where, v, group = sig_id, colour = sig_id)) +
    geom_point(data = X, aes(where, v, colour = sig_id), inherit.aes = FALSE, size = 1.1, alpha = 0.45,
               position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
    geom_text(data = lab, aes(label = sprintf("%s (\u00d7%.2f)", sig_id, r)), hjust = 1.12, size = 3, show.legend = FALSE) +
    scale_x_discrete(expand = expansion(add = c(1.4, 0.3))) + labs(title = ttl, x = NULL, y = ylab) + TH + theme(legend.position = "none")
}
f4 <- plot_grid(mk("z_dev", "z_ext_med", "ratio_z", "|z|", "A  Wald z"),
                mk("lhr_dev", "lhr_ext_med", "ratio_hr", "|log hazard ratio| per SD", "B  Hazard ratio"), nrow = 1)
save2(f4, "Fig4_overfit_gradient", 11, 5)

# ---------------------------------------------------------------- 核对
cat("\nFig 1 计数（按类别）：\n"); print(dcast(D[, .N, by = .(class, cls3_bh)], class ~ cls3_bh, value.var = "N"))
cat("\nFig 2 低于对角线：\n"); print(D[, .(n = .N, below = sum(abs(z_real) < zA_p95)), by = class])
cat("\nFig 3：\n"); print(RC[, .(class, role, n, raw, bh)])
cat("\nFig 4：\n"); print(G[, .(sig_id, ratio_z = round(ratio_z, 2), ratio_hr = round(ratio_hr, 2))])
cat("\n输出：\n"); print(list.files(FIG, sprintf("_v2_%s", ST)))
