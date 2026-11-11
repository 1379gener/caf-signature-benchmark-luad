#!/usr/bin/env Rscript
# ==============================================================================
# recover_ren2023_coef.R —— 从已算好的 riskscore 反解 Ren 2023 的模型系数
#
# 背景：Ren 2023 (Front Immunol 14:1201573) 是 4 基因加权模型，按 SAP §3.3
#   归入 B 类。B 类必须有显式权重，否则按 §3.2 排除。
#   原文 Methods 与 Results 给出的 CLEC3B 系数【符号相反】，须以数据为准。
#
# 方法：riskscore 是四个基因表达的线性组合。对 500 个样本做最小二乘回归，
#   若 R² = 1 且残差为机器精度量级，则系数被精确恢复。
#
# 只读，不写任何文件。
# ==============================================================================

f <- path.expand("~/thesis_LUAD_CAF/03_signature_registry/pdfs/Ren2023_code/results/tcga_cox_datas.txt")
if (!file.exists(f)) {
  f2 <- path.expand("~/thesis_LUAD_CAF/03_signature_registry/pdfs/Ren2023_code/07.CAF.module.cli/results/tcga_cox_datas.txt")
  if (file.exists(f2)) f <- f2 else stop("未找到 tcga_cox_datas.txt")
}

hr <- function(t) cat("\n", strrep("=", 66), "\n", t, "\n", strrep("=", 66), "\n", sep="")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k,":"), paste(v, collapse=", ")))

cat("Ren 2023 模型系数反解\n文件: ", f, "\n", sep="")

d <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
kv("维度", paste(dim(d), collapse=" x "))
G <- c("CLEC3B","EXO1","CCNB1","CD302")
stopifnot(all(c(G, "riskscore") %in% colnames(d)))
for (g in c(G, "riskscore")) d[[g]] <- as.numeric(d[[g]])
d <- d[complete.cases(d[, c(G, "riskscore")]), ]
kv("完整样本", nrow(d))

hr("1. 四个基因表达值的描述")
for (g in G) cat(sprintf("  %-8s 范围 [%7.4f, %7.4f]  中位 %7.4f\n",
    g, min(d[[g]]), max(d[[g]]), median(d[[g]])))
kv("riskscore 范围", sprintf("[%.4f, %.4f]", min(d$riskscore), max(d$riskscore)))

hr("2. 最小二乘反解")
fit <- lm(riskscore ~ CLEC3B + EXO1 + CCNB1 + CD302, data = d)
s <- summary(fit)
kv("R²", sprintf("%.10f", s$r.squared))
kv("残差标准误", sprintf("%.3e", s$sigma))
kv("残差绝对值最大", sprintf("%.3e", max(abs(resid(fit)))))
exact <- s$r.squared > 0.999999 && max(abs(resid(fit))) < 1e-6
cat(sprintf("\n  判定: %s\n", ifelse(exact,
  "✅ 系数被精确恢复（R²=1，残差为机器精度量级）",
  "⚠️ 非精确线性关系 —— riskscore 可能还含其他项或做过变换")))

cat("\n  恢复的系数：\n")
co <- coef(fit)
for (n in names(co)) cat(sprintf("    %-14s %+.6f\n", n, co[n]))

hr("3. 与原文报告的系数比对")
cat("  原文 Methods §2.3:  +0.123*CLEC3B + 0.114*EXO1 + 0.103*CCNB1 - 0.177*CD302\n")
cat("  原文 Results §3.3:  -0.123*CLEC3B + 0.114*EXO1 + 0.103*CCNB1 - 0.177*CD302\n")
cat("                      ↑ 两处 CLEC3B 符号相反\n\n")
pub <- c(CLEC3B = 0.123, EXO1 = 0.114, CCNB1 = 0.103, CD302 = -0.177)
cat(sprintf("  %-8s %12s %12s %12s %s\n", "基因", "恢复值", "Methods", "Results", "判定"))
for (g in G) {
  r <- unname(co[g]); m <- unname(pub[g]); rs <- if (g == "CLEC3B") -pub[g] else pub[g]
  verdict <- if (abs(r - m) < 0.002) "与 Methods 一致" else
             if (abs(r - rs) < 0.002) "与 Results 一致" else "两者皆不符"
  cat(sprintf("  %-8s %+12.6f %+12.3f %+12.3f  %s\n", g, r, m, rs, verdict))
}

hr("4. 方向自洽性检查（系数符号 vs 单变量效应）")
cat("  原文 Discussion 将 CLEC3B 与 CD302 描述为保护性基因。\n")
cat("  保护性基因在【风险】评分中应为负系数。\n\n")
if (requireNamespace("survival", quietly = TRUE)) {
  suppressPackageStartupMessages(library(survival))
  tcol <- if ("OS.time" %in% colnames(d)) "OS.time" else "time"
  ecol <- if ("OS" %in% colnames(d)) "OS" else "status"
  d[[tcol]] <- as.numeric(d[[tcol]]); d[[ecol]] <- as.numeric(d[[ecol]])
  ok <- !is.na(d[[tcol]]) & d[[tcol]] > 0 & !is.na(d[[ecol]])
  cat(sprintf("  %-8s %10s %10s %10s  %s\n", "基因", "系数", "单变量HR", "p", "方向自洽"))
  for (g in G) {
    cx <- summary(coxph(Surv(d[[tcol]][ok], d[[ecol]][ok]) ~ d[[g]][ok]))
    hrv <- cx$conf.int[1]; pv <- cx$coefficients[5]
    consist <- (co[g] > 0 && hrv > 1) || (co[g] < 0 && hrv < 1)
    cat(sprintf("  %-8s %+10.4f %10.3f %10.4f  %s\n", g, co[g], hrv, pv,
                ifelse(consist, "是", "!! 否")))
  }
  cat("\n  系数为正而 HR<1（或反之）= 方向矛盾，须核对原文。\n")
}

hr("5. 供收录表填写的结论")
cat("  SIG-00X   Ren 等, Front Immunol, 2023;14:1201573\n")
cat("  DOI       10.3389/fimmu.2023.1201573\n")
cat("  类别      B（3-9 基因加权模型）\n")
cat("  基因数    4\n")
cat("  基因与权重（以本次反解为准，非以原文文本为准）：\n")
for (g in G) cat(sprintf("    %-8s %+.6f\n", g, co[g]))
cat("\n  备注：原文 Methods 与 Results 的 CLEC3B 系数符号自相矛盾；\n")
cat("        本研究采用从作者提供的 riskscore 反解所得的系数，\n")
cat("        并在论文中注明该不一致。\n")
cat("\n  可写出为 SIG-00X_weights.csv 的内容：\n")
cat("    gene,weight\n")
for (g in G) cat(sprintf("    %s,%.6f\n", g, co[g]))
