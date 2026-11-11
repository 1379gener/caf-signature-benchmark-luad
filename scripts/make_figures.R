#!/usr/bin/env Rscript
# ==============================================================================
# make_figures.R —— Paper 2 主图
#
# ★ SAP §11.1 第 7 条：出图脚本不得出现字面量数据向量，所有图必须从数据文件读。
#   本脚本只从以下四个文件读数，不内嵌任何数值：
#     results_extended_*.csv        逐组合结果（含 z_real、zA_p95、p_emp）
#     results_by_cohort_role_*.csv  队列角色标注
#     overfit_gradient_*.csv        过拟合梯度
#     signature_registry_locked_*.csv  signature 元信息
#
# 四个面板：
#   Fig 1  signature × 队列的信号分类热图（无信号 / 有信号不特异 / 强于随机）
#   Fig 2  |z_real| vs 零分布 |z| 的 95% 分位 —— 区分「无信号」与「不特异」
#   Fig 3  通过率：开发队列 vs 原文验证队列 vs 完全独立队列
#   Fig 4  过拟合梯度：开发队列 z 与外部队列 z 的配对比较
#
# 图内文字用英文（论文为英文），控制台输出用中文。
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
if (!requireNamespace("scales", quietly = TRUE))
  stop("需要 scales 包：install.packages(\"scales\")")

CH    <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark")
DD    <- file.path(CH, "data_derived")
FIG   <- file.path(CH, "figures")
REG   <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")
STAMP <- format(Sys.Date(), "%Y%m%d")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path("/tmp", format(Sys.time(), "figures_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

newest <- function(dir, pat) {
  f <- sort(list.files(dir, pattern = pat, full.names = TRUE), decreasing = TRUE)
  if (!length(f)) stop("未找到 ", pat, " in ", dir)
  f[1]
}
save2 <- function(p, name, w, h) {
  for (ext in c("pdf", "png")) {
    f <- file.path(FIG, sprintf("%s_%s.%s", name, STAMP, ext))
    ggsave(f, p, width = w, height = h, dpi = 300,
           device = ext, limitsize = FALSE)
  }
  cat(sprintf("  %-28s %.0f x %.0f in\n", name, w, h))
}

cat("Paper 2 主图生成\n时间: ", format(Sys.time()), "\n", sep = "")

# ============================================================ 载入
hr("1. 载入数据（不内嵌任何数值）")
fE <- newest(DD, "^results_extended_.*\\.csv$");        kv("逐组合结果", basename(fE))
fR <- newest(DD, "^results_by_cohort_role_.*\\.csv$");   kv("队列角色", basename(fR))
fG <- newest(DD, "^overfit_gradient_.*\\.csv$");         kv("过拟合梯度", basename(fG))
fL <- newest(REG, "^signature_registry_locked_.*\\.csv$"); kv("锁定表", basename(fL))
kv("锁定表 md5", unname(tools::md5sum(fL)))

E <- fread(fE)[evaluable == TRUE]
R <- fread(fR)
G <- fread(fG)
L <- fread(fL)

D <- merge(E, R[, .(cohort, sig_id, role, is_train)], by = c("cohort", "sig_id"), all.x = TRUE)
D <- merge(D, L[, .(sig_id, class_L = class, n_genes)], by = "sig_id", all.x = TRUE)
kv("可评估组合", nrow(D))

# 队列按事件数排序；signature 按类别 + 编号排序
co_ord <- D[, .(ev = n_events[1]), by = cohort][order(-ev), cohort]
D[, cohort := factor(cohort, levels = rev(co_ord))]
sg_ord <- unique(D[order(class_L, sig_id), sig_id])
D[, sig_lab := factor(sprintf("%s (%s, n=%d)", sig_id, class_L, n_genes),
                      levels = unique(sprintf("%s (%s, n=%d)",
                        D[order(class_L, sig_id), sig_id],
                        D[order(class_L, sig_id), class_L],
                        D[order(class_L, sig_id), n_genes])))]

# 信号分类（与 diagnose_results.R 同一判据）
D[, sig_class := fifelse(abs(z_real) < 1.96, "No prognostic signal",
                  fifelse(p_emp_A >= 0.05, "Signal, not specific",
                          "Signal, stronger than random"))]
D[, sig_class := factor(sig_class, levels = c(
  "No prognostic signal", "Signal, not specific", "Signal, stronger than random"))]

PAL <- c("No prognostic signal" = "#BDBDBD",
         "Signal, not specific" = "#F6C667",
         "Signal, stronger than random" = "#C0492B")

TH <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey96", colour = NA),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom", legend.title = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))

# ============================================================ Fig 1
hr("2. Fig 1 — 信号分类热图")
f1 <- ggplot(D, aes(x = sig_lab, y = cohort, fill = sig_class)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_point(data = D[is_train == TRUE], aes(x = sig_lab, y = cohort),
             shape = 4, size = 2.6, stroke = 1.1, colour = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = PAL, name = NULL, drop = FALSE) +
  labs(title = "Prognostic signal and its specificity across cohorts",
       subtitle = paste("x = development cohort of that signature.",
                        "Empty cells = not evaluable (coverage < 70%)."),
       x = NULL, y = NULL) +
  TH + theme(axis.text.x = element_text(angle = 45, hjust = 1),
             legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 1))
save2(f1, "Fig1_signal_heatmap", 9.5, 4.2)

# ============================================================ Fig 2
hr("3. Fig 2 — |z_real| vs 零分布 95% 分位")
mx <- max(c(abs(D$z_real), D$zA_p95), na.rm = TRUE) * 1.05
f2 <- ggplot(D, aes(x = zA_p95, y = abs(z_real), colour = sig_class, shape = class_L)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_hline(yintercept = 1.96, linetype = 3, colour = "grey55") +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = PAL, name = NULL, drop = FALSE) +
  scale_shape_manual(values = c(A = 16, B = 17), name = "Class") +
  coord_equal(xlim = c(0, mx), ylim = c(0, mx)) +
  labs(title = "Real signature vs its own null distribution",
       subtitle = paste0("Diagonal: |z_real| equals the 95th percentile of |z| from 1000\n",
                         "expression-matched random gene sets. Points below the diagonal\n",
                         "are weaker than most random sets. Dotted line: |z| = 1.96."),
       x = "95th percentile of |z| in null distribution A",
       y = "|z| of the real signature") +
  TH + guides(colour = guide_legend(nrow = 1, order = 1),
              shape = guide_legend(nrow = 1, order = 2))
save2(f2, "Fig2_real_vs_null", 6.4, 6.8)

# ============================================================ Fig 3
hr("4. Fig 3 — 按队列角色的通过率")
RL <- c("开发队列" = "Development cohort",
        "原文验证队列" = "Validation cohort (original paper)",
        "完全独立队列" = "Fully independent cohort")
D[, role_en := factor(RL[role], levels = rev(unname(RL)))]
A3 <- D[!is.na(role_en), .(n = .N, pass = sum(p_emp_A < 0.05, na.rm = TRUE)),
        by = role_en][, rate := pass / n][]
f3 <- ggplot(A3, aes(x = role_en, y = rate, fill = role_en)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d/%d  (%.0f%%)", pass, n, 100 * rate)),
            hjust = -0.08, size = 3.2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.42),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c("#7F8C8D", "#95A5A6", "#C0492B")) +
  coord_flip() +
  labs(title = "Pass rate against null distribution A, by cohort role",
       subtitle = "A signature passes when its empirical p < 0.05 in that cohort.",
       x = NULL, y = "Proportion of combinations passing") + TH
save2(f3, "Fig3_pass_by_role", 7.8, 2.9)

# ============================================================ Fig 4
hr("5. Fig 4 — 过拟合梯度")
GG <- G[!is.na(z_train)]
if (nrow(GG)) {
  Gl <- melt(GG[, .(sig_id, train_cohort,
                    Development = abs(z_train), External = abs(z_ext_median))],
             id.vars = c("sig_id", "train_cohort"),
             variable.name = "where", value.name = "abs_z")
  Gl[, where := factor(where, levels = c("Development", "External"))]
  Gl[, sig_id := factor(sig_id, levels = GG[order(-ratio), sig_id])]
  # 标签放在 Development 侧：该侧数值分散，不会像 External 侧那样重叠
  LAB <- Gl[where == "Development"]
  LAB <- merge(LAB, GG[, .(sig_id, ratio)], by = "sig_id", all.x = TRUE)
  LAB[, lbl := sprintf("%s (x%.2f)", sig_id, ratio)]
  f4 <- ggplot(Gl, aes(x = where, y = abs_z, group = sig_id)) +
    geom_line(aes(colour = sig_id), linewidth = 0.9, alpha = 0.9) +
    geom_point(aes(colour = sig_id), size = 2.6) +
    geom_text(data = LAB, aes(x = where, y = abs_z, label = lbl),
              hjust = 1.12, size = 3.1, inherit.aes = FALSE) +
    scale_x_discrete(expand = expansion(add = c(1.35, 0.25))) +
    labs(title = "Performance drop from development to external cohorts",
         subtitle = paste0("Each line is one signature. External value is the median |z|\n",
                           "across its non-development cohorts.\n",
                           "Ratio = |z_development| / |z_external|."),
         x = NULL, y = "|z| of the signature") +
    TH + theme(legend.position = "none")
  save2(f4, "Fig4_overfit_gradient", 6.6, 5.2)
} else cat("  无可用的过拟合梯度数据，跳过\n")

# ============================================================ 数据核对
hr("6. 图中数值的来源核对（供 verify 脚本引用）")
cat("\n  Fig 1 各类别的组合数：\n")
print(D[, .N, by = sig_class][order(-N)], row.names = FALSE)
cat("\n  Fig 2 落在对角线下方（|z_real| < zA_p95）的组合数：\n")
cat(sprintf("    %d / %d  (%.0f%%)\n", sum(abs(D$z_real) < D$zA_p95, na.rm = TRUE),
            nrow(D), 100 * mean(abs(D$z_real) < D$zA_p95, na.rm = TRUE)))
cat("\n  Fig 3 数据：\n"); print(A3, row.names = FALSE)
cat("\n  Fig 4 数据：\n")
if (nrow(GG)) print(GG[order(-ratio), .(sig_id, train_cohort, z_train,
                                        z_ext_median, ratio)], row.names = FALSE)

hr("7. 输出")
for (f in sort(list.files(FIG, pattern = STAMP, full.names = TRUE)))
  cat(sprintf("  %-46s %.0f KB\n", basename(f), file.info(f)$size / 1024))
cat("\n  目录: ", FIG, "\n", sep = "")
cat("\n  ★ 本脚本未内嵌任何数值，全部从 data_derived/ 的 csv 读取（SAP §11.1-7）。\n")
cat("  ★ 论文中出现的每个数字，verify 脚本须有对应 check（SAP §12）。\n")

cat("\n输出: ", OUT, "\n", sep = "")
