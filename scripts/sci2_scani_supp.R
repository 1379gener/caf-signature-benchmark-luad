#!/usr/bin/env Rscript
# ==============================================================================
# sci2_scani_supp.R —— 补执行 SAP §3.4：SCANI 的 myCAF ssGSEA 组分作为 A 类基因集过零分布
#   登记：DEVIATIONS.md 2026-10-25「SAP §3.4 补执行：SCANI myCAF」（运行之前）
#
#   基因集：~/caf_spatial_project/results/01_scRNA/CAF_subtype_signatures.rds 的 $myCAF
#           （SCI-1 step6e_final_model.R L14、L28 所用；只读）
#   做法与 run_benchmark.R 的 A 类完全相同：覆盖率（含别名回退）≥ 70%；打分基因 = 基因集 ∩ 候选池（精确符号）；
#     ssGSEA（ssgsea.norm = FALSE）；单变量 Cox per SD；零分布 A、B 各 1000 次；p_emp 双侧；共表达 QC。
#   seed = make_seed(队列, "SCANI_myCAF", "A"/"B")。已锁定的 71 个组合不改动；结果另存。
#   资源：GSVA 2 进程（BiocParallel MulticoreParam(2)）；6 个队列顺序执行。
#   输出：data_derived/null_benchmark_scani_supp/<队列>__SCANI_myCAF.rds；scani_supp_summary_<日期>.csv
#   放置：ch04_benchmark/scripts/   用法：nohup Rscript sci2_scani_supp.R > /tmp/scani_supp.out 2>&1 &
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA); library(BiocParallel) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark"); DD <- file.path(CH, "data_derived")
OD <- file.path(DD, "null_benchmark_scani_supp"); dir.create(OD, showWarnings = FALSE)
FS <- file.path(DD, sprintf("scani_supp_summary_%s.csv", format(Sys.Date(), "%Y%m%d")))
if (file.exists(FS)) stop("输出已存在，不覆盖：", FS)
LOG <- file.path(CH, "logs", format(Sys.time(), "scani_supp_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
suppressMessages({
  source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
  source(file.path(CH, "scripts/null_engine.R"))
})
BP <- MulticoreParam(workers = 2)
SID <- "SCANI_myCAF"; N_NULL <- 1000L

fg <- path.expand("~/caf_spatial_project/results/01_scRNA/CAF_subtype_signatures.rds")
cat("基因集文件 ", fg, "\nmd5 ", unname(tools::md5sum(fg)), "\n", sep = "")
stopifnot(unname(tools::md5sum(fg)) == "17d3458f1a9f55c0e90316aed443cd3e")
SIGS <- readRDS(fg); cat("元素：", paste(names(SIGS), lengths(SIGS), sep = "=", collapse = ", "), "\n")
stopifnot("myCAF" %in% names(SIGS))
g <- trimws(as.character(SIGS$myCAF)); g <- g[nzchar(g)]
stopifnot(!anyDuplicated(g), length(g) >= 10L)                       # A 类：≥ 10 基因
cat("myCAF 基因数 ", length(g), "\n", paste(g, collapse = " "), "\n", sep = "")

coxz <- function(sc, t, e) { ok <- !is.na(t) & t > 0 & !is.na(e) & is.finite(sc)
  f <- tryCatch(coxph(Surv(t[ok], e[ok]) ~ as.numeric(scale(sc[ok]))), error = function(x) NULL)
  if (is.null(f)) return(c(z = NA, hr = NA, lo = NA, hi = NA))
  s <- summary(f); c(z = s$coefficients[1, "z"], hr = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4]) }

OUT <- list(); t0 <- Sys.time()
for (a in COHORTS) {
  co <- load_cohort(a, verbose = FALSE)
  stopifnot(!anyDuplicated(colnames(co$expr)), ncol(co$expr) == co$n, co$n == as.integer(co$man$n_final[1]))
  if (!co$has_survival) next
  stopifnot(nrow(co$clin) == co$n, co$n_events == as.integer(co$man$n_events[1]))
  fo <- file.path(OD, sprintf("%s__%s.rds", a, SID)); if (file.exists(fo)) stop("已存在，不覆盖：", fo)
  cv <- gene_coverage(g, co$expr)
  row <- list(cohort = a, sig_id = SID, class = "A", n_genes = cv$n_total, n_present = cv$n_present, coverage = cv$frac,
              evaluable = cv$evaluable, n = co$n, n_events = co$n_events, missing_genes = paste(cv$missing, collapse = ";"))
  if (!cv$evaluable) { saveRDS(row, fo); OUT[[a]] <- as.data.table(row); cat(a, "不可评估\n"); next }
  pool <- build_pool(co$expr); gp <- intersect(g, pool$genes); rz <- build_rank_z(co$expr, pool)
  tt <- co$clin$OS_time; ee <- co$clin$OS_event
  real <- score_sets_ssgsea(co$expr, list(gp), prefix = "real")[1, ]
  cr <- coxz(real, tt, ee)
  run_null <- function(sets, type) {
    S <- score_sets_ssgsea(co$expr, sets, prefix = type, bp = BP)
    zn <- apply(S, 1, function(s) coxz(s, tt, ee)[["z"]])
    list(z = zn, p = empirical_p(cr[["z"]], zn))
  }
  sA <- make_seed(a, SID, "A"); nA <- sample_null_A(pool, gp, n_null = N_NULL, seed = sA); rA <- run_null(nA$sets, "A")
  sB <- make_seed(a, SID, "B"); nB <- sample_null_B(pool, rz, gp, n_null = N_NULL, seed = sB); rB <- run_null(nB$sets, "B")
  row <- c(row, list(n_in_pool = length(gp), z_real = cr[["z"]], hr_real = cr[["hr"]], hr_lo = cr[["lo"]], hr_hi = cr[["hi"]],
                     p_emp_A = rA$p, seed_A = sA, n_fallback_A = nA$n_fallback, zA_p95 = unname(quantile(abs(rA$z), .95, na.rm = TRUE)),
                     p_emp_B = rB$p, seed_B = sB, n_relax_B = nB$n_relax, zB_p95 = unname(quantile(abs(rB$z), .95, na.rm = TRUE)),
                     coexpr_real = nB$real_coexpr, coexpr_null_med = nB$null_coexpr_median, coexpr_ratio = nB$coexpr_ratio,
                     coexpr_qc_pass = nB$qc_pass, n_null = N_NULL, engine_version = ENGINE$version))
  saveRDS(c(row, list(z_null_A = rA$z, z_null_B = rB$z)), fo)
  OUT[[a]] <- as.data.table(row)
  cat(sprintf("  %-10s 基因 %d/%d  z %+.2f  p_emp_A %.4f (p95 %.2f)  p_emp_B %.4f  QC %s  第51次采用 %d  %.1f 分钟\n",
              a, length(gp), length(g), cr[["z"]], rA$p, row$zA_p95, rB$p, nB$qc_pass, nB$n_relax,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
X <- rbindlist(OUT, fill = TRUE); fwrite(X, FS)
cat("\n"); print(X[, .(cohort, coverage = round(coverage, 3), n_in_pool, z_real = round(z_real, 2), hr_real = round(hr_real, 2),
                       zA_p95 = round(zA_p95, 2), p_emp_A, p_emp_B, coexpr_qc_pass)], row.names = FALSE)
cat("写出 ", FS, "\n", sep = ""); print(sessionInfo())
