# sci2_gse68465_site_20261111.R —— GSE68465 按中心分层的 Cox（真实分数与原零分布 A 的 1,000 个成员同样分层），及分层后的比例风险检验
# 登记：DEVIATIONS.md 2026-11-11「可重复性复核：阳性对照、蒙特卡洛加抽、GSE68465 按中心分层」(3)（运行之前）
source("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_common_20261111.R")
OD <- "data_derived/gse68465_site_20261111"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
sink("logs/gse68465_site_20261111.log", split = TRUE); cat("run:", format(Sys.time()), "\n")
A <- "GSE68465"; co <- get_cohort(A); tt <- co$clin$OS_time; ee <- co$clin$OS_event
L <- readLines(gzfile(path.expand("~/caf_spatial_project/data/GSE68465/GSE68465_series_matrix.txt.gz")), n = 400)
fld <- function(tag) gsub('"', "", strsplit(L[startsWith(L, tag)][1], "\t")[[1]][-1])
map <- setNames(fld("!Sample_source_name_ch1"), fld("!Sample_geo_accession"))
ids <- if ("gsm" %in% names(co$clin)) co$clin$gsm else colnames(co$expr)
stopifnot(identical(as.character(ids), colnames(co$expr)), all(ids %in% names(map)))
site <- factor(map[ids]); print(table(site, event = ee)); stopifnot(!"Stratagene" %in% levels(site))
z_strat_batch <- function(S) apply(S, 1, function(s) {
  if (!is.finite(sd(s)) || sd(s) == 0) return(NA_real_)
  f <- tryCatch(coxph(Surv(tt, ee) ~ as.numeric(scale(s)) + strata(site)), error = function(e) NULL)
  if (is.null(f)) NA_real_ else unname(summary(f)$coefficients[1, "z"]) })
M <- read_main()[cohort == A]; pool <- build_pool(co$expr); out <- list()
for (sid in M$sig_id) {
  gp <- intersect(SIGS[[sid]]$genes, pool$genes); s <- score_one(co$expr, gp, sid)
  z0 <- cox_one(s, tt, ee)$z; stopifnot(abs(z0 - M[sig_id == sid, z_real]) < 1e-6)
  fs <- coxph(Surv(tt, ee) ~ as.numeric(scale(s)) + strata(site)); zs <- unname(summary(fs)$coefficients[1, "z"])
  ph <- cox.zph(fs, transform = "km")$table[1, "p"]
  sA <- sample_null_A(pool, gp, n_null = 1000L, seed = make_seed(A, sid, "A")); Sm <- score_sets(co$expr, sA$sets, sid, gp, BP)
  zn0 <- cox_z_batch(Sm, tt, ee); orig <- readRDS(file.path("data_derived", if (sid == "SCANI_myCAF") "null_benchmark_scani_supp" else "null_benchmark", M[sig_id == sid, file]))$z_null_A
  same <- isTRUE(max(abs(zn0 - orig), na.rm = TRUE) < 1e-6)
  zns <- z_strat_batch(Sm)
  out[[sid]] <- data.table(cohort = A, sig_id = sid, class = CLS[[sid]], z_unstrat = z0, p_emp_A_unstrat = M[sig_id == sid, p_emp_A],
    z_site_strat = zs, p_emp_A_site_strat = empirical_p(zs, zns), zA_p95_strat = unname(quantile(abs(zns), .95, na.rm = TRUE)),
    ph_p_unstrat = cox.zph(coxph(Surv(tt, ee) ~ as.numeric(scale(s))), transform = "km")$table[1, "p"], ph_p_site_strat = ph,
    null_reproduced = same)
  print(out[[sid]])
}
R <- rbindlist(out); fwrite(R, file.path(OD, "gse68465_site_strat_20261111.csv"))
cat(sprintf("\n|z| ≥ 1.96：未分层 %d，分层 %d；p_emp_A < 0.05：未分层 %d，分层 %d；PH p < 0.05：未分层 %d，分层 %d；零分布复现 %d/%d\n",
  sum(abs(R$z_unstrat) >= 1.96), sum(abs(R$z_site_strat) >= 1.96), sum(R$p_emp_A_unstrat < 0.05), sum(R$p_emp_A_site_strat < 0.05),
  sum(R$ph_p_unstrat < 0.05), sum(R$ph_p_site_strat < 0.05), sum(R$null_reproduced), nrow(R)))
cat("done:", format(Sys.time()), "\n"); sink()
