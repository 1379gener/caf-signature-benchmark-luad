#!/usr/bin/env Rscript
# sci2_nullB_relax_bound.R —— 零分布 B「第 51 次直接采用」成员对 p_emp_B 判定的影响上界（只读已有结果）
#   每个组合：N = 有限 z_null_B 个数，c = #{|z_null_B| ≥ |z_real|}，m = 受影响成员数；
#   j ∈ [max(0, m−(N−c)), min(m, c)] 个属于「超过」，剔除后 p' = (1+c−j)/(1+N−m)；
#   p' 整个区间与原 p 在 0.05 同侧 → 判定不可能改变。m 取 out35 与 relaxed 两种口径。
#   登记：DEVIATIONS.md 2026-10-25「零分布引擎复核」(2)。输出 data_derived/nullB_relax_bound_<日期>.csv
suppressPackageStartupMessages(library(data.table))
CH <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived")
OUT <- file.path(DD, sprintf("nullB_relax_bound_%s.csv", format(Sys.Date(), "%Y%m%d")))
if (file.exists(OUT)) stop("输出已存在，不覆盖：", OUT)
V <- fread(file.path(DD, "verify_nullengine_combo_20261025.csv"))
stopifnot(nrow(V) == 71L, !anyDuplicated(V[, paste(cohort, sig_id)]))
rds <- list.files(file.path(DD, "null_benchmark"), "\\.rds$", full.names = TRUE)
R <- rbindlist(lapply(rds, function(f) { x <- readRDS(f)
  if (!isTRUE(x$evaluable) || is.null(x$z_null_B)) return(NULL)
  zn <- x$z_null_B[is.finite(x$z_null_B)]
  data.table(cohort = x$cohort, sig_id = x$sig_id, z_real = x$z_real, N = length(zn),
             c = sum(abs(zn) >= abs(x$z_real)), p_rds = x$p_emp_B, qc = isTRUE(x$coexpr_qc_pass)) }))
stopifnot(nrow(R) == 71L, !anyDuplicated(R[, paste(cohort, sig_id)]))
R[, p := (1 + c) / (1 + N)]
stopifnot(max(abs(R$p - R$p_rds)) < 1e-12)
X <- merge(R, V[, .(cohort, sig_id, class, k, out35 = B_relaxed_out35, relaxed = B_n_relaxed)], by = c("cohort", "sig_id"))
stopifnot(nrow(X) == 71L)
bound <- function(N, c, m) { jlo <- max(0, m - (N - c)); jhi <- min(m, c)
  c(lo = (1 + c - jhi) / (1 + N - m), hi = (1 + c - jlo) / (1 + N - m)) }
for (tag in c("out35", "relaxed")) {
  b <- t(mapply(bound, X$N, X$c, X[[tag]]))
  X[, paste0("p_lo_", tag) := b[, "lo"]][, paste0("p_hi_", tag) := b[, "hi"]]
  X[, paste0("flip_possible_", tag) := ((p < 0.05) != (b[, "lo"] < 0.05)) | ((p < 0.05) != (b[, "hi"] < 0.05))]
}
cat("受影响的组合（relaxed > 0）：\n")
print(X[relaxed > 0, .(cohort, sig_id, class, k, z_real = round(z_real, 2), c, p = round(p, 4), qc,
                       out35, p_lo_out35 = round(p_lo_out35, 4), p_hi_out35 = round(p_hi_out35, 4), flip_possible_out35,
                       relaxed, flip_possible_relaxed)][order(cohort, sig_id)], row.names = FALSE)
cat(sprintf("\n可能改变 p_emp_B < 0.05 判定的组合数：out35 口径 %d；relaxed 口径 %d（共 %d 个组合）\n",
            sum(X$flip_possible_out35), sum(X$flip_possible_relaxed), nrow(X)))
cat(sprintf("其中计入 pass_B_main（QC 通过）的：out35 %d；relaxed %d\n",
            sum(X$flip_possible_out35 & X$qc), sum(X$flip_possible_relaxed & X$qc)))
fwrite(X, OUT); cat("写出 ", OUT, "\n", sep = "")
