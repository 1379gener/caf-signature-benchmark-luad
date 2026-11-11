#!/usr/bin/env Rscript
# ==============================================================================
# recover_weights.R —— 通用权重反解工具
#
# 用途：多数 signature 论文的补充数据会给出「各基因表达值 + 算好的 risk score」。
#   若 risk score 是表达值的线性组合，可用最小二乘精确反解出系数 ——
#   这比抄论文正文可靠：Ren 2023 的 Methods 与 Results 中 CLEC3B 符号相反，
#   反解结果证明 Results 版正确。
#
# 用法：
#   Rscript recover_weights.R <数据文件> [选项]
#
#   必需：
#     <数据文件>            .txt / .csv / .tsv / .xlsx，宽表，一行一个样本
#   可选：
#     --score=<列名>        风险评分列名。省略则自动探测
#     --genes=A,B,C         基因列名，逗号分隔。省略则自动探测
#     --time=<列名> --event=<列名>   提供则做方向自洽性检查
#     --sheet=<n|名>        xlsx 的工作表
#     --id=SIG-00X          写出文件的编号；提供则落盘到 registry
#     --out=<目录>          落盘目录，默认 ~/thesis_LUAD_CAF/03_signature_registry/genes
#
# 示例：
#   Rscript recover_weights.R data.txt
#   Rscript recover_weights.R data.txt --score=riskScore --time=OS.time --event=OS --id=SIG-002
#
# 不带 --id 时只读不写。
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] %in% c("-h", "--help")) {
  cat(paste(readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))[3:32],
            collapse = "\n"), "\n")
  quit(status = 0)
}

opt <- function(k, d = NA) {
  h <- grep(paste0("^--", k, "="), args, value = TRUE)
  if (length(h)) sub(paste0("^--", k, "="), "", h[1]) else d
}
FILE   <- args[1]
SCORE  <- opt("score"); GENES <- opt("genes")
TCOL   <- opt("time");  ECOL  <- opt("event")
SHEET  <- opt("sheet"); SIGID <- opt("id")
OUTDIR <- opt("out", "~/thesis_LUAD_CAF/03_signature_registry/genes")

hr <- function(t) cat("\n", strrep("=", 66), "\n", t, "\n", strrep("=", 66), "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("通用权重反解\n文件: ", FILE, "\n", sep = "")

# ---------------------------------------------------------------- 读入
FILE <- path.expand(FILE)
if (!file.exists(FILE)) stop("文件不存在: ", FILE)
d <- if (grepl("\\.xlsx?$", FILE, ignore.case = TRUE)) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("需要 readxl")
  sh <- if (is.na(SHEET)) 1 else if (grepl("^[0-9]+$", SHEET)) as.integer(SHEET) else SHEET
  as.data.frame(suppressMessages(readxl::read_excel(FILE, sheet = sh)))
} else {
  sep <- if (grepl("\\.csv$", FILE, ignore.case = TRUE)) "," else "\t"
  read.delim(FILE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
}
kv("维度", paste(dim(d), collapse = " x "))

# ---------------------------------------------------------------- 探测列
num_ok <- vapply(d, function(x) {
  v <- suppressWarnings(as.numeric(as.character(x))); mean(!is.na(v)) > 0.9 }, TRUE)
num_cols <- names(num_ok)[num_ok]
for (c0 in num_cols) d[[c0]] <- suppressWarnings(as.numeric(as.character(d[[c0]])))

hr("1. 列探测")
if (is.na(SCORE)) {
  cand <- grep("risk[_.]?score|riskscore|^score$|^risk$|signature[_.]?score|CAF[_.]?Index",
               num_cols, ignore.case = TRUE, value = TRUE)
  cand <- cand[!grepl("z$|_z$|group|cutoff|norm", cand, ignore.case = TRUE)]  # 排除 Z 化 / 分组列
  if (!length(cand)) {
    cat("  !! 未能自动识别风险评分列。数值型列名如下，请用 --score= 指定：\n")
    cat("     ", paste(num_cols, collapse = ", "), "\n", sep = "")
    quit(status = 0)
  }
  SCORE <- cand[1]
  if (length(cand) > 1) kv("评分列候选（取首个）", cand)
}
kv("风险评分列", SCORE)

if (is.na(GENES)) {
  # 基因列的启发式：全大写字母数字、2-12 字符、非临床字段
  ex <- "^(age|sex|gender|time|status|os|dfs|rfs|pfs|stage|grade|risk|score|sample|id|cutoff|t|n|m)"
  g <- num_cols[!grepl(ex, num_cols, ignore.case = TRUE)]
  g <- setdiff(g, SCORE)
  g <- g[grepl("^[A-Za-z][A-Za-z0-9._-]{1,14}$", g)]
  if (!length(g)) {
    cat("  !! 未能自动识别基因列，请用 --genes= 指定。数值型列：\n")
    cat("     ", paste(num_cols, collapse = ", "), "\n", sep = "")
    quit(status = 0)
  }
  GENES <- g
} else GENES <- trimws(strsplit(GENES, ",")[[1]])
miss <- setdiff(GENES, colnames(d))
if (length(miss)) stop("以下基因列不存在: ", paste(miss, collapse = ", "))
kv("候选基因列", sprintf("%d 个：%s", length(GENES), paste(GENES, collapse = ", ")))

use <- d[complete.cases(d[, c(GENES, SCORE)]), ]
kv("完整样本", nrow(use))
if (nrow(use) < length(GENES) + 5) stop("样本数不足以反解 ", length(GENES), " 个系数")

# ---------------------------------------------------------------- 反解
hr("2. 最小二乘反解")
f  <- as.formula(paste(SCORE, "~", paste(sprintf("`%s`", GENES), collapse = " + ")))
fit <- lm(f, data = use); s <- summary(fit)
kv("R²", sprintf("%.10f", s$r.squared))
kv("残差绝对值最大", sprintf("%.3e", max(abs(resid(fit)))))
exact <- s$r.squared > 0.999999 && max(abs(resid(fit))) < 1e-6
cat(sprintf("\n  判定: %s\n", ifelse(exact,
  "✅ 精确恢复（R²=1，残差为机器精度量级）",
  "⚠️ 非精确线性关系")))

co <- coef(fit); icpt <- co[1]; co <- co[-1]; names(co) <- GENES
cat("\n  截距: ", sprintf("%+.6f", icpt),
    ifelse(abs(icpt) < 1e-8, "（为 0，符合无截距的加权和）",
           "（非 0 —— 原模型可能先做过中心化）"), "\n", sep = "")
cat("\n  恢复的系数：\n")
for (n in names(co)) cat(sprintf("    %-14s %+.6f\n", n, co[n]))

if (!exact) {
  cat("\n  非精确的常见原因（按概率）：\n")
  cat("    1) 基因列不全 —— 模型还含未被识别为基因列的项\n")
  cat("    2) 原模型用的是 Z 标准化后的表达值，而本表给的是原始值\n")
  cat("    3) 评分列选错（如选到了 Z 化后的版本）\n")
  cat("    4) 表中的表达值经过舍入，精度不足\n")
  cat("\n  逐基因贡献（|系数| x 该基因标准差）：\n")
  ctb <- sort(abs(co) * vapply(use[GENES], sd, 0), decreasing = TRUE)
  for (n in names(ctb)) cat(sprintf("    %-14s %.4f\n", n, ctb[n]))
  cat("\n  尝试 Z 标准化后重解：\n")
  z <- use; for (g in GENES) z[[g]] <- as.numeric(scale(z[[g]]))
  f2 <- summary(lm(f, data = z))
  cat(sprintf("    Z 化后 R² = %.10f  %s\n", f2$r.squared,
      ifelse(f2$r.squared > 0.999999, "<<< 原模型使用 Z 标准化后的表达值", "")))
  if (f2$r.squared > 0.999999) {
    co <- coef(lm(f, data = z))[-1]; names(co) <- GENES
    cat("    Z 尺度下的系数：\n")
    for (n in names(co)) cat(sprintf("      %-14s %+.6f\n", n, co[n]))
  }
}

# ---------------------------------------------------------------- 方向自洽
if (!is.na(TCOL) && !is.na(ECOL) && requireNamespace("survival", quietly = TRUE)) {
  hr("3. 方向自洽性检查")
  if (!all(c(TCOL, ECOL) %in% colnames(use))) {
    cat("  !! 时间或事件列不存在，跳过\n")
  } else {
    suppressPackageStartupMessages(library(survival))
    tt <- as.numeric(use[[TCOL]]); ee <- as.numeric(use[[ECOL]])
    ok <- !is.na(tt) & tt > 0 & !is.na(ee)
    cat(sprintf("  可用样本 %d，事件 %d\n\n", sum(ok), sum(ee[ok])))
    cat(sprintf("  %-14s %10s %10s %10s  %s\n", "基因", "系数", "单变量HR", "p", "方向自洽"))
    nbad <- 0
    for (g in GENES) {
      cx <- summary(coxph(Surv(tt[ok], ee[ok]) ~ use[[g]][ok]))
      h <- cx$conf.int[1]; p <- cx$coefficients[5]
      cons <- (co[g] > 0 && h > 1) || (co[g] < 0 && h < 1)
      if (!cons) nbad <- nbad + 1
      cat(sprintf("  %-14s %+10.4f %10.3f %10.4f  %s\n", g, co[g], h, p,
                  ifelse(cons, "是", "!! 否")))
    }
    cat(sprintf("\n  方向矛盾的基因数: %d %s\n", nbad,
        ifelse(nbad == 0, "（全部自洽）", "—— 须核对原文，或说明该基因在多变量中翻转")))
    cxr <- summary(coxph(Surv(tt[ok], ee[ok]) ~ use[[SCORE]][ok]))
    cat(sprintf("  评分整体：HR %.3f (%.3f-%.3f)，p = %.4g\n",
        cxr$conf.int[1], cxr$conf.int[3], cxr$conf.int[4], cxr$coefficients[5]))
  }
}

# ---------------------------------------------------------------- 落盘
hr("4. 供收录表填写")
cat("  类别      B（3-9 基因加权模型）\n")
cat(sprintf("  基因数    %d %s\n", length(GENES),
    ifelse(length(GENES) >= 10, "  ← ≥10，按 SAP §3.3 应归入 A 类，改用等权 ssGSEA",
    ifelse(length(GENES) < 3, "  ← <3，按 SAP §3.2 排除", ""))))
cat("\n  gene,weight\n")
for (n in names(co)) cat(sprintf("  %s,%.6f\n", n, co[n]))

if (!is.na(SIGID)) {
  OUTDIR <- path.expand(OUTDIR)
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  fg <- file.path(OUTDIR, paste0(SIGID, "_genes.txt"))
  fw <- file.path(OUTDIR, paste0(SIGID, "_weights.csv"))
  writeLines(names(co), fg)
  write.csv(data.frame(gene = names(co), weight = round(unname(co), 6)),
            fw, row.names = FALSE, quote = FALSE)
  cat("\n  已写出:\n    ", fg, "\n    ", fw, "\n", sep = "")
  cat("\n  收录表备注建议：权重由作者提供的风险评分反解所得",
      sprintf("（R²=%.6f）", s$r.squared), "，非取自原文文本。\n", sep = "")
} else {
  cat("\n  （未提供 --id，本次只读不写。加 --id=SIG-00X 可落盘）\n")
}
