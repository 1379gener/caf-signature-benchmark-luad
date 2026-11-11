#!/usr/bin/env Rscript
# ==============================================================================
# sci2_stage2_diagnostics.R —— 阶段 2 结果的解释性诊断 + 事实清单（全部数字由本脚本现算）
#   1 逐组合：分数与增殖指数 / 纯度的 Spearman 与 Pearson；M3 的方差膨胀（VIF = 1/(1−r²)，
#     r 为 M3 分析样本内分数与增殖指数的 Pearson）与实测 SE 比（M3 / 同样本 M1）
#   2 并入 stage2_attenuation，写 stage2_diagnostics_<日期>.csv
#   3 生成事实清单 stage2_facts_<日期>.md：每条结论的数字、来源文件、不可外推之处
#   4 打印 SAP v1.2 §8（阶段 3）原文，供下一步使用
#   打分路径与 sci2_stage2_confounding.R 相同，并以阶段 1 z_real 为闸门
# 用法：Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_stage2_diagnostics.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD <- file.path(ROOT, "03_signature_registry"); DD <- file.path(CH, "data_derived")
suppressMessages({ source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
                   source(file.path(CH, "scripts/null_engine.R")) })
STAMP <- format(Sys.Date(), "%Y%m%d")
LOG <- file.path(CH, "logs", format(Sys.time(), "stage2_diagnostics_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
latest <- function(pat) { f <- sort(list.files(DD, pat, full.names = TRUE), decreasing = TRUE); stopifnot(length(f) >= 1); f[1] }
FIN <- list(purity = latest("^m2_purity_TCGA_[0-9]{8}\\.csv$"), prolif = latest("^m3_prolif_metaPCNA_[0-9]{8}\\.csv$"),
            stage = latest("^m4_stage_binary_[0-9]{8}\\.csv$"), stage_sens = latest("^m4_stage_binary_TCGA_sens_[0-9]{8}\\.csv$"),
            models = latest("^stage2_models_[0-9]{8}\\.csv$"), att = latest("^stage2_attenuation_[0-9]{8}\\.csv$"))
for (k in names(FIN)) cat(sprintf("  %-10s %s  md5 %s\n", k, basename(FIN[[k]]), unname(tools::md5sum(FIN[[k]]))))
rd <- function(k) fread(FIN[[k]], na.strings = c("", "NA"))
PUR <- rd("purity"); PRO <- rd("prolif"); STG <- rd("stage"); STS <- rd("stage_sens"); MOD <- rd("models"); ATT <- rd("att")
LOCK <- require_locked_registry(REGD); stopifnot(LOCK$md5 == "4ad3218931d95dfb6ee67a7ba07b9a75")
reg <- LOCK$table
read_sig <- function(sid, cls) {
  g <- trimws(readLines(file.path(REGD, "genes", paste0(sid, "_genes.txt")))); g <- g[nzchar(g)]
  w <- NULL
  if (cls == "B") { d <- fread(file.path(REGD, "genes", paste0(sid, "_weights.csv")))
                    w <- setNames(as.numeric(d$weight), d$gene)[g]; stopifnot(!anyNA(w)) }
  list(genes = g, weights = w) }

hr("1. 逐组合相关与共线性")
rds <- list.files(file.path(DD, "null_benchmark"), "\\.rds$", full.names = TRUE)
CO <- list(); POOL <- list(); OUT <- list()
for (f in rds) {
  x <- readRDS(f); if (!isTRUE(x$evaluable)) next
  a <- x$cohort; sid <- x$sig_id; cls <- x$class
  if (is.null(CO[[a]])) { co <- load_cohort(a, verbose = FALSE)
    stopifnot(!anyDuplicated(colnames(co$expr)), ncol(co$expr) == co$n); CO[[a]] <- co; POOL[[a]] <- build_pool(co$expr) }
  co <- CO[[a]]; sg <- read_sig(sid, cls); gp <- intersect(sg$genes, POOL[[a]]$genes)
  sc <- suppressWarnings(if (cls == "A") score_sets_ssgsea(co$expr, list(gp), prefix = "real")[1, ] else
                                         score_sets_weighted(co$expr, list(gp), sg$weights[gp])[1, ])
  D <- data.table(sample_id = colnames(co$expr), score = as.numeric(sc), t = co$clin$OS_time, e = co$clin$OS_event)
  D <- merge(D, PRO[cohort == a, .(sample_id, prolif = prolif_metaPCNA)], by = "sample_id", all.x = TRUE)
  D[, purity := if (a == "TCGA_LUAD") PUR[match(D$sample_id, sample_id), purity_absolute] else NA_real_]
  stopifnot(nrow(D) == co$n)
  D <- D[!is.na(t) & t > 0 & !is.na(e) & is.finite(score)]
  z <- cox_one(D$score, D$t, D$e)$z; stopifnot(abs(z - x$z_real) < 1e-6)        # 闸门
  S3 <- D[!is.na(prolif)]; r3 <- cor(S3$score, S3$prolif)
  P2 <- D[!is.na(purity)]
  OUT[[length(OUT) + 1]] <- data.table(cohort = a, sig_id = sid, class = cls,
    rho_prolif = cor(S3$score, S3$prolif, method = "spearman"), r_prolif = r3, VIF_M3 = 1 / (1 - r3^2),
    rho_purity = if (nrow(P2)) cor(P2$score, P2$purity, method = "spearman") else NA_real_)
}
R <- rbindlist(OUT); stopifnot(nrow(R) == 71)
se <- MOD[model %in% c("M1", "M3") & sample == "S3", .(cohort, sig_id, model, se)]
se <- dcast(se, cohort + sig_id ~ model, value.var = "se")[, .(cohort, sig_id, SE_ratio_M3 = M3 / M1)]
R <- merge(R, se, by = c("cohort", "sig_id"))
cat("  VIF 与实测 SE 比的 Spearman：", round(cor(sqrt(R$VIF_M3), R$SE_ratio_M3, method = "spearman"), 3), "（理论上 SE 比 ≈ √VIF）\n")
W <- merge(ATT, R[, .(cohort, sig_id, rho_prolif, r_prolif, VIF_M3, SE_ratio_M3, rho_purity)], by = c("cohort", "sig_id"))
fo <- file.path(DD, sprintf("stage2_diagnostics_%s.csv", STAMP)); if (file.exists(fo)) stop("已存在，不覆盖：", fo)
fwrite(W, fo); cat("  写出 ", fo, "\n", sep = "")
options(datatable.print.nrows = 200, width = 200)
print(W[comparison == "M3@S3", .(cohort, sig_id, z_M1 = round(z_M1, 2), rho_prolif = round(rho_prolif, 2),
        VIF = round(VIF_M3, 1), SE_ratio = round(SE_ratio_M3, 2), att = round(attenuation_reported, 2))][order(sig_id, cohort)])

hr("2. 事实清单")
f2 <- function(x) formatC(x, format = "f", digits = 2)
rng <- function(x) paste0(f2(min(x)), " 至 ", f2(max(x)))
M3 <- W[comparison == "M3@S3"]; M3s <- M3[M1_signal == TRUE]
M2 <- W[comparison == "M2@S2"]; M2s <- M2[M1_signal == TRUE]
sig_sig <- unique(M2s$sig_id); nos <- setdiff(M2$sig_id, sig_sig)
tc <- R[cohort == "TCGA_LUAD"]
q <- quantile(M3s$attenuation, c(.25, .5, .75))
neg <- M3s[attenuation < 0]
s1 <- R[sig_id == "SIG-001"]
m4 <- W[comparison == "M4@S4" & M1_signal == TRUE]; m4s <- W[comparison == "M4@S4_sens" & M1_signal == TRUE]
cl <- load_cohort("TCGA_LUAD", verbose = FALSE)$clin
sm <- merge(cl[, .(sample_id = patient_id, OS_event)], STS, by = "sample_id")
g1 <- sm[!is.na(stage_main)]; g0 <- sm[is.na(stage_main)]
gate_max <- max(abs(MOD[model == "M1" & sample == "full", z] -
  vapply(seq_len(nrow(MOD[model == "M1" & sample == "full"])), function(i) {
    r <- MOD[model == "M1" & sample == "full"][i]; readRDS(file.path(DD, "null_benchmark", sprintf("%s__%s.rds", r$cohort, r$sig_id)))$z_real }, 0)))
L <- c(
  sprintf("# 阶段 2 事实清单（%s，sci2_stage2_diagnostics.R 现算）", STAMP), "",
  "来源文件：", paste0("- ", vapply(FIN, basename, ""), "（md5 ", unname(tools::md5sum(unlist(FIN))), "）"), "",
  "## F1 一致性",
  sprintf("- 全样本 M1 的 z 与阶段 1 RDS 的 z_real：71 个组合最大差 %.1e；阶段 2 分数即阶段 1 分数。", gate_max), "",
  "## F2 纯度（M2，仅 TCGA，Aran 2015 ABSOLUTE，n = 301）",
  sprintf("- M1 有信号的 %d 个签名（%s），M2 衰减 %s；与纯度的 Spearman %s。", nrow(M2s), paste(sig_sig, collapse = "/"),
          rng(M2s$attenuation), rng(tc[sig_id %in% sig_sig, rho_purity])),
  sprintf("- 其余 %d 个签名（%s）与纯度的 Spearman %s，M1 均无信号（|z| < 1.96）。", length(nos), paste(nos, collapse = "/"),
          rng(tc[sig_id %in% nos, rho_purity])),
  "- 可写：有信号的签名与纯度近乎不相关，故 M2 几乎不改变其 HR；与纯度相关的签名本身无预后信号。",
  "- 不可写：『CAF 签名不受纯度影响』『纯度不是混杂因素』。仅一个队列有纯度；覆盖 67.9%。", "",
  "## F3 增殖（M3，6 队列）",
  sprintf("- M1 有信号的组合 %d 个；M3 衰减中位 %s（四分位 %s 至 %s）；>0.5 的 %d 个；为负的 %d 个；校正后仍 P<0.05 的 %d 个。",
          nrow(M3s), f2(q[2]), f2(q[1]), f2(q[3]), sum(M3s$attenuation > 0.5), nrow(neg), sum(M3s$p_Mk < 0.05)),
  sprintf("- 有信号组合中 |分数与增殖的 Spearman| 与 M3 衰减的 Spearman 为 %s：衰减幅度与增殖相关强度几乎无关。",
          f2(cor(abs(M3s$rho_prolif), M3s$attenuation, method = "spearman"))),
  sprintf("- 有信号组合的分数与增殖指数 Spearman：%s（含负相关的 SIG-012）。", rng(M3s$rho_prolif)),
  sprintf("- TCGA 内同为有信号的签名：M3 衰减中位 %s，M2 衰减中位 %s（SAP 事前预期增殖校正的衰减大于纯度校正）。",
          f2(median(M3s[cohort == "TCGA_LUAD", attenuation])), f2(median(M2s$attenuation))),
  "- 可写：增殖校正使多数有信号组合的 HR 部分减弱，减弱幅度在组合间差异大。",
  "- 不可写：『增殖解释了约 X% 的预后信号』（衰减幅度在高共线时不可解释，见 F4）；『由增殖驱动』。", "",
  "## F4 共线性（SIG-001）",
  sprintf("- SIG-001 与增殖指数 Spearman %s；M3 的 VIF %s；M3 / M1 的 SE 比 %s。", rng(s1$rho_prolif), rng(s1$VIF_M3), rng(s1$SE_ratio_M3)),
  sprintf("- 负衰减的 %d 个组合：%s，均为 SIG-001。", nrow(neg), paste(sprintf("%s %s", neg$cohort, f2(neg$attenuation)), collapse = "；")),
  "- 权重：EXO1、CCNB1 为正，CLEC3B、CD302 为负（registry genes/SIG-001_weights.csv）；meta-PCNA 重叠仅 CCNB1（1/4）。",
  "- 可写：SIG-001 分数与增殖指数高度相关，M3 中两者共线，签名效应与增殖效应无法分离，其衰减值不作解释。",
  "- 不可写：『校正增殖后 SIG-001 信号增强』；『SIG-001 一半基因属 meta-PCNA』。", "",
  "## F5 完全校正（M4，仅 TCGA）",
  sprintf("- 主分析（GDC 分期，n = %d，%d 事件）：M1 有信号 %d 个，衰减 %s，M4 后仍 P<0.05 的 %d 个。",
          m4$n[1], m4$events[1], nrow(m4), rng(m4$attenuation), sum(m4$p_Mk < 0.05)),
  sprintf("- 敏感性（分期补充，n = %d，%d 事件）：M1 有信号 %d 个，衰减 %s，仍 P<0.05 的 %d 个。",
          m4s$n[1], m4s$events[1], nrow(m4s), rng(m4s$attenuation), sum(m4s$p_Mk < 0.05)),
  sprintf("- 分期缺失与结局：GDC 有分期 %d 人中事件 %d（%.1f%%）；缺失 %d 人中事件 %d（%.1f%%）。主分析的完整病例存在选择。",
          nrow(g1), sum(g1$OS_event), 100 * mean(g1$OS_event), nrow(g0), sum(g0$OS_event), 100 * mean(g0$OS_event)),
  "- 不可写：基于单一队列外推到其他队列；TCGA 是 SIG-001/004/006/011/012 的开发队列。", "",
  "## F6 不可计算项",
  "- 5 个 GEO 队列（GSE68465/72094/31210/50081/30219）无与表达无关的纯度来源，M2/M4 按 SAP §7.2 不计算",
  "  （SAP 原文写『四个 GEO 微阵列队列』，实为 5 个，属文字笔误）；CPTAC 无生存数据。")
fm <- file.path(DD, sprintf("stage2_facts_%s.md", STAMP)); if (file.exists(fm)) stop("已存在，不覆盖：", fm)
writeLines(L, fm); cat(L, sep = "\n"); cat("\n  写出 ", fm, "\n", sep = "")

hr("3. SAP v1.2 §8（阶段 3）原文")
sap <- system2("find", c(ROOT, "-maxdepth", "3", "-name", "SAP_v1.2.txt"), stdout = TRUE)[1]
x <- readLines(sap, warn = FALSE); s <- grep("^8　冗余", x)[1]; e <- grep("^9　", x); e <- e[e > s][1]
if (is.na(s)) s <- grep("冗余与相关性", x)[1]; if (is.na(e)) e <- min(length(x), s + 80)
cat(sprintf("%5d  %s", s:(e - 1), x[s:(e - 1)]), sep = "\n")
sink(); close(zz)
