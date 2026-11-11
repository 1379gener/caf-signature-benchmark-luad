# 由 v2 脚本生成 v3（不改 v2）；每处替换必须精确匹配预期次数，否则停止
setwd(path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts"))
patch <- function(src, dst, reps) {
  if (file.exists(dst)) stop("已存在，不覆盖：", dst)
  x <- paste(readLines(src, encoding = "UTF-8"), collapse = "\n")
  for (r in reps) {
    n <- lengths(regmatches(x, gregexpr(r[[1]], x, fixed = TRUE)))
    if (n != r[[3]]) stop(sprintf("%s：预期匹配 %d 次，实际 %d 次：\n%s", src, r[[3]], n, substr(r[[1]], 1, 120)))
    x <- gsub(r[[1]], r[[2]], x, fixed = TRUE) }
  writeLines(x, dst, useBytes = TRUE); cat("写出", dst, "\n") }

patch("make_figures_v2.R", "make_figures_v3.R", list(
  list(r"---(_v2_%s)---", r"---(_v3_%s)---", 3),
  list(r"---(scale_x_discrete(expand = expansion(add = c(1.4, 0.3))))---",
       r"---(scale_x_discrete(limits = levels(Gl$where), expand = expansion(add = c(1.4, 0.3))))---", 1),
  list(r"---(facet_grid(. ~ class_lab, scales = "free_x", space = "free_x") +)---",
       r"---(facet_grid(. ~ class_lab, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = function(x) sub("_LUAD", "-LUAD", x)) +)---", 1)
))

patch("make_figures_stage23_v2.R", "make_figures_stage23_v3.R", list(
  list(r"---("Fig5_stage2_confounding_v2")---", r"---("Fig5_stage2_confounding_v3")---", 1),
  list(r"---("Fig6_stage3_readouts_v2")---", r"---("Fig6_stage3_readouts_v3")---", 1),
  list(r"---(labs(x = "Spearman with meta-PCNA proliferation index", y = NULL) + TH + theme(legend.position = "top"))---",
       r"---(labs(x = "Spearman with meta-PCNA proliferation index", y = NULL) + TH +
  theme(legend.position = "top", legend.text = element_text(margin = margin(l = 2, r = 14))))---", 1),
  list(r"---("%s (n = %d); order: ward.D2 on 1 - |rho|")---", r"---("%s (n = %d); order: Ward (ward.D2) on 1 \u2212 |\u03c1|")---", 1),
  list(r"---(axis.text.y = element_text(size = 5.5)))---",
       r"---(axis.text.y = element_text(size = 5.5), plot.background = element_rect(fill = "white", colour = NA)))---", 1)
))
