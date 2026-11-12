# sci2_posctrl_20261111.R —— 阳性对照：meta-PCNA 增殖基因集（131 基因）作为 A 类基因集，走与主分析相同的零分布 A、B
# 登记：DEVIATIONS.md 2026-11-11「可重复性复核：阳性对照、蒙特卡洛加抽、GSE68465 按中心分层」(1)（运行之前）
source("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_common_20261111.R")
OD <- "data_derived/posctrl_20261111"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
sink("logs/posctrl_20261111.log", split = TRUE); cat("run:", format(Sys.time()), "\n")
g <- trimws(readLines(file.path(REGDIR, "confounders", "metaPCNA_genes.txt"))); g <- unique(g[nzchar(g)])
cat("meta-PCNA genes:", length(g), "\n"); stopifnot(length(g) == 131)
SID <- "POS-METAPCNA"; out <- list()
for (a in names(NEXP)) {
  co <- get_cohort(a); tt <- co$clin$OS_time; ee <- co$clin$OS_event
  pool <- build_pool(co$expr); rz <- build_rank_z(co$expr, pool); gp <- intersect(g, pool$genes)
  rs <- score_sets_ssgsea(co$expr, list(gp), prefix = "real")[1, ]; rc <- cox_one(rs, tt, ee)
  sA <- sample_null_A(pool, gp, n_null = 1000L, seed = make_seed(a, SID, "A"))
  zA <- cox_z_batch(score_sets_ssgsea(co$expr, sA$sets, prefix = "A", bp = BP), tt, ee)
  sB <- sample_null_B(pool, rz, gp, n_null = 1000L, seed = make_seed(a, SID, "B"))
  zB <- cox_z_batch(score_sets_ssgsea(co$expr, sB$sets, prefix = "B", bp = BP), tt, ee)
  out[[a]] <- data.table(cohort = a, n_genes = length(g), n_used = length(gp), z_real = rc$z, hr_per_sd = rc$hr,
    zA_p95 = unname(quantile(abs(zA), .95, na.rm = TRUE)), p_emp_A = empirical_p(rc$z, zA), n_null_A = sum(is.finite(zA)),
    zB_p95 = unname(quantile(abs(zB), .95, na.rm = TRUE)), p_emp_B = empirical_p(rc$z, zB), n_null_B = sum(is.finite(zB)),
    coexpr_real = sB$real_coexpr, coexpr_ratio = sB$coexpr_ratio, qc_pass = sB$qc_pass, n_relax_B = sB$n_relax,
    seed_A = make_seed(a, SID, "A"), seed_B = make_seed(a, SID, "B"))
  print(out[[a]][, .(cohort, n_used, z = round(z_real, 2), zA_p95 = round(zA_p95, 2), p_emp_A, zB_p95 = round(zB_p95, 2), p_emp_B, qc_pass)])
}
R <- rbindlist(out)
R[, q_A := p.adjust(p_emp_A, "BH")]   # 仅作描述：6 个队列内一族
fwrite(R, file.path(OD, "posctrl_metaPCNA_20261111.csv"))
cat(sprintf("\n通过零分布 A（p_emp < 0.05）：%d / %d；零分布 B（QC 通过且 p_emp < 0.05）：%d / %d\n",
            sum(R$p_emp_A < 0.05), nrow(R), sum(R$qc_pass & R$p_emp_B < 0.05), sum(R$qc_pass)))
cat("done:", format(Sys.time()), "\n"); sink()
