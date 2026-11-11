# 由 v1 脚本生成 v2（不改原脚本）；每处替换必须精确匹配预期次数，否则停止
setwd(path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts"))
patch <- function(src, dst, reps) {
  if (file.exists(dst)) stop("已存在，不覆盖：", dst)
  x <- paste(readLines(src, encoding = "UTF-8"), collapse = "\n")
  for (r in reps) {
    n <- lengths(regmatches(x, gregexpr(r[[1]], x, fixed = TRUE)))
    if (n != r[[3]]) stop(sprintf("%s：预期匹配 %d 次，实际 %d 次：\n%s", src, r[[3]], n, substr(r[[1]], 1, 120)))
    x <- gsub(r[[1]], r[[2]], x, fixed = TRUE) }
  writeLines(x, dst, useBytes = TRUE); cat("写出", dst, "\n") }

patch("make_figures_v1.R", "make_figures_v2.R", list(
  list(r"---(_v1_%s)---", r"---(_v2_%s)---", 3),
  list(r"---(dpi = 300, device = ext))---", r"---(dpi = 300, device = if (ext == "pdf") cairo_pdf else ext))---", 1),
  list(r"---(  geom_point(data = D[role == "开发队列"], shape = 4, size = 2.6, stroke = 1.1, inherit.aes = FALSE, aes(sig_lab, cohort)) +)---",
       r"---(  geom_point(data = D[role == "开发队列"], shape = 4, size = 2.6, stroke = 1.1, inherit.aes = FALSE, aes(sig_lab, cohort)) +
  geom_point(data = D[role == "原文验证队列"], shape = 1, size = 2.8, stroke = 1.1, inherit.aes = FALSE, aes(sig_lab, cohort)) +)---", 1),
  list(r"---(x = development cohort. * developed by the authors.)---",
       r"---(\u00d7 development cohort; \u25cb validation cohort of the original study; * developed by the authors.)---", 1),
  list(r"---(  geom_point(size = 2, alpha = 0.85) + facet_wrap(~ class_lab))---",
       r"---(  geom_point(aes(shape = cohort == "TCGA_LUAD"), size = 2, alpha = 0.85) + facet_wrap(~ class_lab) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("Other cohorts", "TCGA-LUAD"), name = NULL))---", 1),
  list(r"---(  ggplot(Gl, aes(where, v, group = sig_id, colour = sig_id)) + geom_line(linewidth = 0.9) + geom_point(size = 2.4) +)---",
       r"---(  X <- S4x[sig_id %in% G$sig_id & cohort != "TCGA_LUAD"]
  X[, v := if (dev == "z_dev") abs(z) else abs(log(hr_per_sd))]
  chk <- merge(X[, .(m = median(v)), by = sig_id], G[, .(sig_id, e = get(ext))], by = "sig_id"); stopifnot(max(abs(chk$m - chk$e)) < 1e-6)
  X[, where := factor("External (median)", levels = levels(Gl$where))]
  ggplot(Gl, aes(where, v, group = sig_id, colour = sig_id)) +
    geom_point(data = X, aes(where, v, colour = sig_id), inherit.aes = FALSE, size = 1.1, alpha = 0.45,
               position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) +)---", 1),
  list(r"---(sprintf("%s (x%.2f)", sig_id, r))---", r"---(sprintf("%s (\u00d7%.2f)", sig_id, r))---", 1),
  list(r"---(G <- P[part == "overfit" & sig_id != "SCANI_myCAF"]; stopifnot(nrow(G) == 6L))---",
       r"---(G <- P[part == "overfit" & sig_id != "SCANI_myCAF"]; stopifnot(nrow(G) == 6L)
S4x <- fread(md(file.path(DD, "supp_tables_20261025", "S4_all_combinations.csv"))))---", 1)
))

patch("make_figures_stage23.R", "make_figures_stage23_v2.R", list(
  list(r"---(sprintf("%d pairs with |z| >= 1.96 under M1 (same sample)\n* VIF >= 2.5 (score vs proliferation index)", nrow(sig3)))---",
       r"---(sprintf("%d combinations with |z| \u2265 1.96 under M1 (same sample)\n* VIF \u2265 2.5 (score vs proliferation index)", nrow(sig3)))---", 1),
  list(r"---(         ((model == "M1" & sample == "full") | (model == "M2" & sample == "S2") | (model == "M4" & sample %in% c("S4", "S4_sens")))]
B[, lab := factor(fcase(model == "M1", "M1 (all samples)", model == "M2", "M2 + purity",
                        sample == "S4", "M4 full adj.", default = "M4 full adj. (stage-supplemented)"),
                  levels = c("M1 (all samples)", "M2 + purity", "M4 full adj.", "M4 full adj. (stage-supplemented)"))])---",
       r"---(         ((model == "M1" & sample %in% c("full", "S2", "S4", "S4_sens")) | (model == "M2" & sample == "S2") | (model == "M4" & sample %in% c("S4", "S4_sens")))]
LV <- c("M1, all samples", "M1, purity sample", "M2 + purity", "M1, M4 sample", "M4 full adj.",
        "M1, stage-supplemented sample", "M4 full adj. (stage-supplemented)")
B[, lab := factor(fcase(model == "M1" & sample == "full", LV[1], model == "M1" & sample == "S2", LV[2], model == "M2", LV[3],
                        model == "M1" & sample == "S4", LV[4], model == "M4" & sample == "S4", LV[5],
                        model == "M1" & sample == "S4_sens", LV[6], default = LV[7]), levels = LV)]
print(B[, .N, by = sig_id]); stopifnot(all(B[, .N, by = sig_id]$N == 7)))---", 1),
  list(r"---("Fig5_stage2_confounding", 10, 9))---", r"---("Fig5_stage2_confounding_v2", 10, 12))---", 1),
  list(r"---(subtitle = "Post hoc exploratory analysis (registered 2026-10-25)")---", r"---(subtitle = "Post hoc exploratory analysis")---", 1),
  list(r"---(heat <- function(a) {)---",
       r"---(pretty_lab <- function(x) { x <- sub("^MCP_", "MCP-counter ", x); x <- sub("^xCell_", "xCell ", x); x <- sub("^EPIC_", "EPIC ", x)
  x <- sub("^ESTIMATE_", "ESTIMATE ", x); x <- sub("Prolif_metaPCNA", "Proliferation index", x)
  x <- sub("Purity_ABSOLUTE", "Purity (ABSOLUTE)", x); gsub("_", " ", x) }
heat <- function(a) {)---", 1),
  list(r"---(  D[, `:=`(x = factor(x, o), y = factor(y, rev(o)))])---",
       r"---(  D[, `:=`(x = factor(pretty_lab(as.character(x)), pretty_lab(o)), y = factor(pretty_lab(as.character(y)), rev(pretty_lab(o))))])---", 1),
  list(r"---(R6[, coh := coh_lab(cohort)])---", r"---(R6[, `:=`(coh = coh_lab(cohort), readout = pretty_lab(as.character(readout)))])---", 1),
  list(r"---(E <- copy(EXP)[, `:=`(coh = coh_lab(cohort), model = )---", r"---(E <- copy(EXP)[, `:=`(coh = coh_lab(cohort), readout = pretty_lab(readout), model = )---", 1),
  list(r"---("Fig6_stage3_readouts", 13, 7.5))---", r"---("Fig6_stage3_readouts_v2", 13, 7.5))---", 1),
  list(r"---(for (a in setdiff(COH, "TCGA_LUAD")) save2(heat(a), sprintf("FigS_cor_%s", a), 7.5, 7))---",
       r"---(# v2：Figures S1–S5 不重出)---", 1)
))
