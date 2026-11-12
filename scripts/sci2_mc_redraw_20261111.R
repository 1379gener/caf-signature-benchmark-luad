# sci2_mc_redraw_20261111.R —— 蒙特卡洛误差：边界组合（0.02 < p_emp < 0.10）用新种子加抽 10,000 次；全部组合给出 p_emp 的蒙特卡洛标准误
# 登记：DEVIATIONS.md 2026-11-11「可重复性复核：阳性对照、蒙特卡洛加抽、GSE68465 按中心分层」(2)（运行之前）
source("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_common_20261111.R")
OD <- "data_derived/mc_redraw_20261111"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
sink("logs/mc_redraw_20261111.log", split = TRUE); cat("run:", format(Sys.time()), "\n")
M <- read_main(); print(M[, .N, by = class]); stopifnot(nrow(M) == 77, all(M$sig_id %in% names(SIGS)))
M[, se_A := sqrt(p_emp_A * (1 - p_emp_A) / 1000)][, se_B := sqrt(p_emp_B * (1 - p_emp_B) / 1000)]
BL <- rbind(M[p_emp_A > 0.02 & p_emp_A < 0.10, .(cohort, sig_id, class, null = "A", p_1000 = p_emp_A)],
            M[class == "A" & qc & p_emp_B > 0.02 & p_emp_B < 0.10, .(cohort, sig_id, class, null = "B", p_1000 = p_emp_B)])
cat("边界组合数：", nrow(BL), "\n"); print(BL)
N <- 10000L; res <- list()
for (a in unique(BL$cohort)) {
  co <- get_cohort(a); tt <- co$clin$OS_time; ee <- co$clin$OS_event
  pool <- build_pool(co$expr); rz <- build_rank_z(co$expr, pool)
  for (i in which(BL$cohort == a)) {
    sid <- BL$sig_id[i]; nt <- BL$null[i]; gp <- intersect(SIGS[[sid]]$genes, pool$genes); t0 <- Sys.time()
    rc <- cox_one(score_one(co$expr, gp, sid), tt, ee)
    zr <- M[cohort == a & sig_id == sid, z_real]; stopifnot(abs(rc$z - zr) < 1e-6)
    seed <- make_seed(a, sid, paste0(nt, "10k"))
    S <- if (nt == "A") sample_null_A(pool, gp, n_null = N, seed = seed) else sample_null_B(pool, rz, gp, n_null = N, seed = seed)
    zn <- cox_z_batch(score_sets(co$expr, S$sets, sid, gp, BP), tt, ee)
    p <- empirical_p(rc$z, zn); nf <- sum(is.finite(zn))
    res[[length(res) + 1]] <- data.table(BL[i], p_10000 = p, n_finite = nf, se_10000 = sqrt(p * (1 - p) / nf), seed_10000 = seed,
      qc_10000 = if (nt == "B") S$qc_pass else NA, minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
    print(res[[length(res)]])
  }
}
R <- rbindlist(res)
# 用 10,000 次的 p 替换边界组合后，重算主分析的 BH（队列 × 类别；零分布 B 只对 A 类），比较判定
M2 <- copy(M)
for (i in seq_len(nrow(R))) { col <- if (R$null[i] == "A") "p_emp_A" else "p_emp_B"
  M2[cohort == R$cohort[i] & sig_id == R$sig_id[i], (col) := R$p_10000[i]] }
q <- function(D) { D <- copy(D); D[, q_A := p.adjust(p_emp_A, "BH"), by = .(cohort, class)]
  D[class == "A", q_B := p.adjust(p_emp_B, "BH"), by = cohort]; D }
Q1 <- q(M); Q2 <- q(M2)
CMP <- merge(Q1[, .(cohort, sig_id, class, z_real, q_A_1000 = q_A, q_B_1000 = q_B, qc)], Q2[, .(cohort, sig_id, q_A_mix = q_A, q_B_mix = q_B)], by = c("cohort", "sig_id"))
CMP[, pass_A_1000 := q_A_1000 < 0.05][, pass_A_mix := q_A_mix < 0.05]
CMP[, pass_B_1000 := class == "A" & qc & q_B_1000 < 0.05][, pass_B_mix := class == "A" & qc & q_B_mix < 0.05]
cat(sprintf("\n主分析通过数（零分布 A）：1,000 次 A 类 %d、B 类 %d；替换后 A 类 %d、B 类 %d\n",
  sum(CMP[class == "A", pass_A_1000]), sum(CMP[class == "B", pass_A_1000]), sum(CMP[class == "A", pass_A_mix]), sum(CMP[class == "B", pass_A_mix])))
cat("判定改变的组合：\n"); print(CMP[pass_A_1000 != pass_A_mix | pass_B_1000 != pass_B_mix])
cat(sprintf("全部 77 个组合 p_emp_A 的蒙特卡洛标准误：中位 %.4f，最大 %.4f\n", median(M$se_A), max(M$se_A)))
fwrite(R, file.path(OD, "mc_redraw_borderline_20261111.csv")); fwrite(CMP, file.path(OD, "mc_redraw_decisions_20261111.csv"))
fwrite(M[, .(cohort, sig_id, class, p_emp_A, se_A, p_emp_B, se_B)], file.path(OD, "mc_se_all_20261111.csv"))
cat("done:", format(Sys.time()), "\n"); sink()
