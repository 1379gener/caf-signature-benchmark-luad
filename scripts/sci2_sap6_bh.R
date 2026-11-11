#!/usr/bin/env Rscript
# ==============================================================================
# sci2_sap6_bh.R —— 补执行 SAP §6 的 BH FDR（登记：DEVIATIONS.md 2026-10-25「SAP §6 / §3.3 / §3.4」）
#
#   SAP §6：p_emp < 0.05 判定通过；跨 signature（同一队列、同一类别内）做 BH FDR，q < 0.05；
#           跨队列不校正；零分布 A 与 B 之间不校正。
#   SAP §3.3：A、B 两类结果不得合并计算通过比例 → 本脚本所有汇总均按类别分开。
#   零分布 B 的 BH 只对 A 类做（DEVIATIONS 2026-10-25 G）；主计数仍要求 coexpr_qc_pass。
#   只读已有结果，不重算任何 z 或 p_emp。
#   输出（data_derived/，已存在即停止）：sap6_bh_combo_<日期>.csv、sap6_bh_summary_<日期>.csv
#   放置：ch04_benchmark/scripts/   用法：Rscript sci2_sap6_bh.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
CH <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived")
ST <- format(Sys.Date(), "%Y%m%d")
F1 <- file.path(DD, sprintf("sap6_bh_combo_%s.csv", ST)); F2 <- file.path(DD, sprintf("sap6_bh_summary_%s.csv", ST))
for (f in c(F1, F2)) if (file.exists(f)) stop("输出已存在，不覆盖：", f)
newest <- function(p) { f <- sort(list.files(DD, p, full.names = TRUE), decreasing = TRUE); stopifnot(length(f) > 0); f[1] }
fE <- newest("^results_extended_.*\\.csv$"); fR <- file.path(DD, "results_by_cohort_role_20261025.csv")
for (f in c(fE, fR)) cat(sprintf("  %-40s md5 %s\n", basename(f), unname(tools::md5sum(f))))
E <- fread(fE)[evaluable == TRUE]; R <- fread(fR)
stopifnot(nrow(E) == 71L, !anyDuplicated(E[, paste(cohort, sig_id)]))
E <- merge(E, R[, .(cohort, sig_id, role)], by = c("cohort", "sig_id")); stopifnot(nrow(E) == 71L, !anyNA(E$role))

E[, q_A := p.adjust(p_emp_A, "BH"), by = .(cohort, class)]
E[, n_family_A := .N, by = .(cohort, class)]
E[class == "A", q_B := p.adjust(p_emp_B, "BH"), by = cohort]
E[, `:=`(signal = abs(z_real) >= 1.96,
         pass_A_raw = p_emp_A < 0.05, pass_A_bh = q_A < 0.05,
         pass_B_raw = class == "A" & p_emp_B < 0.05 & coexpr_qc_pass,
         pass_B_bh  = class == "A" & !is.na(q_B) & q_B < 0.05 & coexpr_qc_pass)]
E[, cls3_raw := fifelse(!signal, "no_signal", fifelse(pass_A_raw, "stronger", "not_specific"))]
E[, cls3_bh  := fifelse(!signal, "no_signal", fifelse(pass_A_bh,  "stronger", "not_specific"))]

cat("\n== 每个 (队列, 类别) 的 BH 家族大小\n"); print(dcast(unique(E[, .(cohort, class, n_family_A)]), cohort ~ class, value.var = "n_family_A"))

cat("\n== 未校正通过、BH 后不再通过的组合\n")
print(E[pass_A_raw & !pass_A_bh, .(cohort, sig_id, class, role, z = round(z_real, 2), p_emp_A, q_A = round(q_A, 4))][order(class, sig_id)], row.names = FALSE)
cat("\n== BH 后仍通过的组合\n")
print(E[pass_A_bh == TRUE, .(cohort, sig_id, class, role, z = round(z_real, 2), p_emp_A, q_A = round(q_A, 4))][order(class, sig_id)], row.names = FALSE)

S1 <- E[, .(n = .N, no_signal = sum(cls3_raw == "no_signal"),
            not_specific_raw = sum(cls3_raw == "not_specific"), stronger_raw = sum(cls3_raw == "stronger"),
            not_specific_bh = sum(cls3_bh == "not_specific"), stronger_bh = sum(cls3_bh == "stronger"),
            below_p95 = sum(abs(z_real) < zA_p95),
            passB_raw = if (class[1] == "A") sum(pass_B_raw) else NA_integer_,
            passB_bh  = if (class[1] == "A") sum(pass_B_bh)  else NA_integer_), by = class][order(class)]
S2 <- E[, .(n = .N, pass_A_raw = sum(pass_A_raw), pass_A_bh = sum(pass_A_bh)), by = .(class, role)][order(class, role)]
S3 <- E[, .(pass_A_raw = sum(pass_A_raw), pass_A_bh = sum(pass_A_bh)), by = .(class, sig_id)][order(class, sig_id)]
cat("\n== 按类别汇总（SAP §3.3：不合并）\n"); print(S1, row.names = FALSE)
cat("\n== 按类别 × 队列角色\n"); print(S2, row.names = FALSE)
cat("\n== 按签名\n"); print(S3, row.names = FALSE)

fwrite(E[, .(cohort, sig_id, class, role, z_real, zA_p95, p_emp_A, q_A, n_family_A, p_emp_B, q_B, coexpr_qc_pass,
             pass_A_raw, pass_A_bh, pass_B_raw, pass_B_bh, cls3_raw, cls3_bh)][order(class, sig_id, cohort)], F1)
fwrite(rbindlist(list(S1[, .(table = "by_class", class, role = NA_character_, sig_id = NA_character_, n, no_signal, not_specific_raw, stronger_raw,
                             not_specific_bh, stronger_bh, below_p95, passB_raw, passB_bh)],
                      S2[, .(table = "by_class_role", class, role, n, pass_A_raw, pass_A_bh)],
                      S3[, .(table = "by_signature", class, sig_id, pass_A_raw, pass_A_bh)]), fill = TRUE), F2)
cat("\n写出：\n  ", F1, "\n  ", F2, "\n", sep = "")
