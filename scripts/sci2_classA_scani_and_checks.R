#!/usr/bin/env Rscript
# ==============================================================================
# sci2_classA_scani_and_checks.R —— 只读已有结果，不重算任何 z / p_emp
#   登记：DEVIATIONS.md 2026-10-25「SCANI 并入 A 类 BH；事后检查：方向、HR 版过拟合比、通过率 CI」（运行之前）
#
#   第 1 部分（SAP §3.4 + §6）：SCANI_myCAF 并入 A 类，每队列 BH 家族 = 6；q_A、q_B 重算；
#            零分布 B 主计数要求 coexpr_qc_pass。主结果 = 含 SCANI；不含 SCANI 的结果来自 sap6_bh_combo_20261025.csv。
#            SCANI 的零分布 B「第 51 次直接采用」成员，按 relaxed 口径做最坏情形上界（同 sci2_nullB_relax_bound.R）。
#   第 2 部分（事后，描述性）：
#     (a) 方向：B 类为风险分（原文高分 = 高风险），预期 z > 0；SCANI 预期 z > 0（SCI-1）；
#         其余 A 类原文未给等权分数的方向，记为「未定义」，不判定。
#     (b) HR 版过拟合比 = |log HR 开发队列| / 外部队列 |log HR| 中位数（per SD），与 z 版并列；
#         对象：开发队列为 TCGA 的签名（SIG-001/004/006/007/011/012；SCANI 另列并注明其开发样本为 II–IV 期子集）。
#     (c) 按「类别 × 队列角色」的通过数给出精确二项 95% CI（未校正与 BH 两套）。
#   输出（data_derived/，已存在即停止）：classA_scani_bh_<日期>.csv、posthoc_checks_<日期>.csv
#   放置：ch04_benchmark/scripts/   用法：Rscript sci2_classA_scani_and_checks.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
CH <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived")
ST <- format(Sys.Date(), "%Y%m%d")
F1 <- file.path(DD, sprintf("classA_scani_bh_%s.csv", ST)); F2 <- file.path(DD, sprintf("posthoc_checks_%s.csv", ST))
for (f in c(F1, F2)) if (file.exists(f)) stop("输出已存在，不覆盖：", f)
md <- function(f) { cat(sprintf("  %-38s md5 %s\n", basename(f), unname(tools::md5sum(f)))); f }
fB <- md(file.path(DD, "sap6_bh_combo_20261025.csv")); fS <- md(file.path(DD, "scani_supp_summary_20261025.csv"))
fE <- md(sort(list.files(DD, "^results_extended_.*\\.csv$", full.names = TRUE), decreasing = TRUE)[1])
B <- fread(fB); S <- fread(fS); E <- fread(fE)[evaluable == TRUE]
stopifnot(nrow(B) == 71L, nrow(S) == 6L, all(S$evaluable), !anyDuplicated(S$cohort))

# ------------------------------------------------------------------ 1. SCANI 并入 A 类
cat("\n== 1. A 类 BH（含 SCANI，家族 = 6）\n")
role_scani <- function(a) fifelse(a == "TCGA_LUAD", "开发队列", fifelse(a == "GSE72094", "原文验证队列", "完全独立队列"))
A <- rbind(B[class == "A", .(cohort, sig_id, role, z_real, zA_p95, p_emp_A, p_emp_B, coexpr_qc_pass, q_A_noSCANI = q_A, q_B_noSCANI = q_B)],
           S[, .(cohort, sig_id, role = role_scani(cohort), z_real, zA_p95, p_emp_A, p_emp_B, coexpr_qc_pass, q_A_noSCANI = NA_real_, q_B_noSCANI = NA_real_)])
stopifnot(nrow(A) == 36L, A[, .N, by = cohort][, all(N == 6L)])
A[, `:=`(q_A = p.adjust(p_emp_A, "BH"), q_B = p.adjust(p_emp_B, "BH")), by = cohort]
A[, `:=`(signal = abs(z_real) >= 1.96, pass_A_raw = p_emp_A < 0.05, pass_A_bh = q_A < 0.05,
         pass_B_raw = p_emp_B < 0.05 & coexpr_qc_pass, pass_B_bh = q_B < 0.05 & coexpr_qc_pass)]
A[, cls3_bh := fifelse(!signal, "no_signal", fifelse(pass_A_bh, "stronger", "not_specific"))]
print(A[sig_id == "SCANI_myCAF", .(cohort, role, z = round(z_real, 2), p_emp_A, q_A = round(q_A, 4), p_emp_B, q_B = round(q_B, 4), cls3_bh)], row.names = FALSE)
cat("\n  含 SCANI 后 BH 判定与不含 SCANI 时不同的原 A 类组合：\n")
print(A[sig_id != "SCANI_myCAF" & ((q_A < 0.05) != (q_A_noSCANI < 0.05) | (q_B < 0.05) != (q_B_noSCANI < 0.05)),
        .(cohort, sig_id, q_A_noSCANI, q_A, q_B_noSCANI, q_B)], row.names = FALSE)
SA <- A[, .(n = .N, no_signal = sum(!signal), not_specific_raw = sum(signal & !pass_A_raw), stronger_raw = sum(pass_A_raw),
            not_specific_bh = sum(signal & !pass_A_bh), stronger_bh = sum(pass_A_bh), below_p95 = sum(abs(z_real) < zA_p95),
            passB_raw = sum(pass_B_raw), passB_bh = sum(pass_B_bh))]
cat("\n  A 类合计（含 SCANI，36 个组合）：\n"); print(SA, row.names = FALSE)

# SCANI 零分布 B 放宽成员的最坏情形上界（relaxed 口径）
cat("\n  SCANI 零分布 B「第 51 次采用」上界：\n")
RB <- rbindlist(lapply(S$cohort, function(a) {
  x <- readRDS(file.path(DD, "null_benchmark_scani_supp", sprintf("%s__SCANI_myCAF.rds", a)))
  zn <- x$z_null_B[is.finite(x$z_null_B)]; N <- length(zn); c0 <- sum(abs(zn) >= abs(x$z_real)); m <- as.integer(x$n_relax_B)
  jlo <- max(0, m - (N - c0)); jhi <- min(m, c0); p <- (1 + c0) / (1 + N)
  lo <- (1 + c0 - jhi) / (1 + N - m); hi <- (1 + c0 - jlo) / (1 + N - m)
  data.table(cohort = a, m_relaxed = m, p_emp_B = p, p_lo = lo, p_hi = hi, flip_possible = (p < .05) != (lo < .05) | (p < .05) != (hi < .05)) }))
print(RB, row.names = FALSE); stopifnot(all(abs(RB$p_emp_B - S[match(RB$cohort, cohort), p_emp_B]) < 1e-12))
fwrite(A[order(sig_id, cohort)], F1)

# ------------------------------------------------------------------ 2. 事后检查
cat("\n== 2a. 方向（B 类与 SCANI 预期 z > 0）\n")
B[, expected_dir := fifelse(class == "B", "positive", "undefined")]
D <- rbind(B[, .(cohort, sig_id, class, z_real, expected_dir, signal = abs(z_real) >= 1.96, pass_A_raw, pass_A_bh)],
           A[sig_id == "SCANI_myCAF", .(cohort, sig_id, class = "A", z_real, expected_dir = "positive", signal, pass_A_raw, pass_A_bh)])
D[, dir_ok := fifelse(expected_dir == "undefined", NA, z_real > 0)]
print(D[expected_dir != "undefined", .(n = .N, z_pos = sum(dir_ok), signal = sum(signal), signal_wrong_dir = sum(signal & !dir_ok),
                                       pass_raw_wrong_dir = sum(pass_A_raw & !dir_ok), pass_bh_wrong_dir = sum(pass_A_bh & !dir_ok)), by = class])
print(D[expected_dir != "undefined" & !dir_ok, .(cohort, sig_id, z = round(z_real, 2), signal)], row.names = FALSE)

cat("\n== 2b. 过拟合比：z 版与 HR 版（per SD）\n")
H <- rbind(E[, .(cohort, sig_id, z_real, loghr = log(hr_real))],
           S[, .(cohort, sig_id, z_real, loghr = log(hr_real))])
dev <- c("SIG-001", "SIG-004", "SIG-006", "SIG-007", "SIG-011", "SIG-012", "SCANI_myCAF")
G <- H[sig_id %in% dev, .(z_dev = abs(z_real[cohort == "TCGA_LUAD"]), z_ext_med = median(abs(z_real[cohort != "TCGA_LUAD"])),
                          lhr_dev = abs(loghr[cohort == "TCGA_LUAD"]), lhr_ext_med = median(abs(loghr[cohort != "TCGA_LUAD"])),
                          n_ext = sum(cohort != "TCGA_LUAD")), by = sig_id]
G[, `:=`(ratio_z = z_dev / z_ext_med, ratio_hr = lhr_dev / lhr_ext_med)]
print(G[order(-ratio_z), lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)], row.names = FALSE)
cat("  注：SCANI 的开发样本为 TCGA II–IV 期（166 例），此处 TCGA 为全分期 443 例。\n")

cat("\n== 2c. 按类别 × 队列角色的通过数与精确二项 95% CI\n")
RA <- rbind(B[class == "B", .(class, role, pass_A_raw, pass_A_bh)], A[, .(class = "A", role, pass_A_raw, pass_A_bh)])
CI <- RA[, .(n = .N, raw = sum(pass_A_raw), bh = sum(pass_A_bh)), by = .(class, role)][order(class, role)]
CI[, c("raw_lo", "raw_hi") := as.list(binom.test(raw, n)$conf.int), by = .(class, role)]
CI[, c("bh_lo", "bh_hi") := as.list(binom.test(bh, n)$conf.int), by = .(class, role)]
print(CI[, lapply(.SD, function(v) if (is.numeric(v) && !is.integer(v)) round(v, 3) else v)], row.names = FALSE)

fwrite(rbindlist(list(D[, .(part = "direction", cohort, sig_id, class, z_real, expected_dir, dir_ok, signal)],
                      G[, .(part = "overfit", sig_id, z_dev, z_ext_med, ratio_z, lhr_dev, lhr_ext_med, ratio_hr, n_ext)],
                      CI[, .(part = "role_ci", class, role, n, raw, raw_lo, raw_hi, bh, bh_lo, bh_hi)]), fill = TRUE), F2)
cat("\n写出：\n  ", F1, "\n  ", F2, "\n", sep = "")
