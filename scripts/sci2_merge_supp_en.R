#!/usr/bin/env Rscript
# sci2_merge_supp_en.R — 把英文 S1 / S2 排除理由 / 列说明并入 S2–S7 工作簿（只读两个输入，输出新文件，不覆盖）
#   输入：data_derived/supp_tables_20261025/SCI-2_Supplementary_Tables_S2-S7.xlsx
#         data_derived/supp_tables_20261025/SCI-2_Supplementary_S1_and_notes_EN.xlsx
#   输出：data_derived/supp_tables_20261025/SCI-2_Supplementary_Tables_S1-S7_EN.xlsx
#   v5：SIG-009/010 来源作者更正为 Koeck（S1 第 28 条 → 共 28 条）
#   v4：S1 按 SAP v1.2 修订表与 DEVIATIONS.md 定稿（27 条），无【核对】
#   v3：读入后解码 HTML 数字实体（v2 输出含 &#167; 等乱码，已作废）
#   v2：S7 的 class_main / class_w 中文类别译为英文；残留检查改为逐表报告（perl 正则 + 固定串「【」）
suppressPackageStartupMessages({ library(openxlsx); library(data.table) })
OD <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/data_derived/supp_tables_20261025")
F1 <- file.path(OD, "SCI-2_Supplementary_Tables_S2-S7.xlsx"); F2 <- file.path(OD, "SCI-2_Supplementary_S1_and_notes_EN.xlsx")
FO <- file.path(OD, "SCI-2_Supplementary_Tables_S1-S7_EN.xlsx")
stopifnot(file.exists(F1), file.exists(F2)); if (file.exists(FO)) stop("已存在，不覆盖：", FO)
for (f in c(F1, F2)) cat(sprintf("  %-48s md5 %s\n", basename(f), unname(tools::md5sum(f))))
# openxlsx 读 openpyxl 写出的 inlineStr 时，非 ASCII 字符以 &#NNN; 实体原样返回 → 读入后统一解码
dec <- function(x) {
  if (!is.character(x)) return(x)
  m <- gregexpr("&#[0-9]+;", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(e) vapply(as.integer(gsub("\\D", "", e)), intToUtf8, ""))
  for (k in list(c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&apos;", "'"), c("&amp;", "&"))) x <- gsub(k[1], k[2], x, fixed = TRUE)
  x
}
rx <- function(f, s, ...) { d <- as.data.table(read.xlsx(f, s, check.names = FALSE, sep.names = " ", ...)); d[] <- lapply(d, dec); setnames(d, dec(names(d))); d }

S1 <- rx(F2, "Table S1", startRow = 3); X <- rx(F2, "S2 exclusion reasons", startRow = 3)
CD <- rx(F2, "Column dictionary"); RM <- rx(F2, "README"); S2 <- rx(F1, "S2_registry")
stopifnot(nrow(S1) == 28L, nrow(X) == 6L, nrow(S2) == 19L, sum(grepl("^excluded", S2$status)) == 6L, "Reason (English)" %in% names(X))

# S2：排除行名称与理由改为英文（按 name 匹配；匹配不上即停）
S2[name == "EPIC CAF (包内)", name := "EPIC CAF (package)"]
# SIG-009/010：注册表作者名记为 Ruhle，DOI 10.1080/2162402X.2023.2274130 的第一作者为 Koeck（锁定表不改，只改补充表显示名）
stopifnot(S2[id %in% c("SIG-009", "SIG-010"), all(grepl("^Ruhle 2023", name))])
S2[id %in% c("SIG-009", "SIG-010"), name := sub("^Ruhle 2023", "Koeck 2023", name)]
ex <- grepl("^excluded", S2$status); m <- match(S2$name[ex], X$Candidate)
if (anyNA(m)) { print(S2$name[ex]); print(X$Candidate); stop("S2 排除行 name 与英文表 Candidate 不一致") }
S2[ex, note := X$`Reason (English)`[m]]

OUT <- list(README = RM, S1_deviation_log = S1, S2_registry = S2)
for (nm in setdiff(getSheetNames(F1), "S2_registry")) OUT[[nm]] <- rx(F1, nm)

# S7：类别标签译为英文（与 S4 的 classification 取值一致）；出现未知标签即停
CL <- c("无信号" = "no_signal", "不特异" = "not_specific", "强于随机" = "stronger")
for (cc in c("class_main", "class_w")) {
  v <- OUT$S7_SIG012_published_coef[[cc]]; stopifnot(!is.null(v))
  tr <- ifelse(v %in% names(CL), CL[v], v); bad <- setdiff(unique(tr), unname(CL))
  if (length(bad)) stop("S7 ", cc, " 出现未知标签：", paste(bad, collapse = ", "))
  set(OUT$S7_SIG012_published_coef, j = cc, value = unname(tr))
}
cat("S7 类别（unadjusted）：\n"); print(OUT$S7_SIG012_published_coef[, .(cohort, class_main, class_w)])
OUT$Column_dictionary <- CD

write.xlsx(OUT, FO, colWidths = "auto")
cat("\n写出 ", FO, "\n", sep = ""); print(sapply(OUT, nrow))

# 投稿前检查：逐表列出仍含中文或【 的单元格
cat("\n仍含中文 / 【核对】 / 未解码实体的单元格（按表）：\n"); tot <- 0L
for (nm in names(OUT)) {
  v <- as.character(unlist(OUT[[nm]], use.names = FALSE)); v <- v[!is.na(v)]
  hit <- v[grepl("\u3010", v, fixed = TRUE) | grepl("\\p{Han}", v, perl = TRUE) | grepl("&#[0-9]+;", v)]
  if (length(hit)) { cat(sprintf("\n[%s] %d 处\n", nm, length(hit))); print(substr(hit, 1, 110)) }
  tot <- tot + length(hit)
}
cat("\n合计：", tot, "\n")
