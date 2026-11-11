# sci2_ph_check_20261110.R —— 77 个组合的比例风险假定检验（cox.zph, KM 变换）；只读输入，输出新目录
# 登记：DEVIATIONS.md 2026-11-10「比例风险假定检验」（运行之前）
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT <- path.expand("~/thesis_LUAD_CAF"); CHDIR <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGDIR <- file.path(ROOT, "03_signature_registry"); setwd(CHDIR)
OD <- "data_derived/ph_check_20261110"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
sink(file.path("logs", "ph_check_20261110.log"), split = TRUE); cat("run:", format(Sys.time()), "\n")
suppressMessages({ source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")); source("scripts/null_engine.R") })
norm <- function(x) gsub("-", "_", x)
ref <- rbind(fread("data_derived/sap6_bh_combo_20261025.csv")[, .(cohort, sig_id, class, z_real)],
             fread("data_derived/classA_scani_bh_20261025.csv")[sig_id == "SCANI_myCAF", .(cohort, sig_id, class = "A", z_real)])
stopifnot(nrow(ref) == 77)
reg <- fread(file.path(REGDIR, "signature_registry_locked_20260903.csv")); stopifnot(unname(tools::md5sum(file.path(REGDIR, "signature_registry_locked_20260903.csv"))) == "4ad3218931d95dfb6ee67a7ba07b9a75")
read_sig <- function(sid, cls) {
  g <- trimws(readLines(file.path(REGDIR, "genes", paste0(sid, "_genes.txt")))); g <- g[nzchar(g)]; w <- NULL
  if (cls == "B") { d <- fread(file.path(REGDIR, "genes", paste0(sid, "_weights.csv"))); w <- setNames(as.numeric(d$weight), d$gene)[g]; stopifnot(!anyNA(w)) }
  list(genes = g, weights = w) }
SIGS <- setNames(lapply(seq_len(nrow(reg)), function(i) read_sig(reg$sig_id[i], reg$class[i])), reg$sig_id)
cand <- list.files(ROOT, pattern = "(?i)mycaf.*\\.(txt|csv)$", recursive = TRUE, full.names = TRUE)
cand <- cand[!grepl("ARCHIVE|superseded|bak", cand)]
nG <- vapply(cand, function(f) { x <- if (grepl("csv$", f)) fread(f)[[1]] else readLines(f); length(unique(trimws(x[nzchar(trimws(x))]))) }, numeric(1))
cat("SCANI 候选文件：\n"); print(data.table(file = cand, n = nG))
sc_file <- cand[nG == 50]
if (length(sc_file) == 1) { x <- if (grepl("csv$", sc_file)) fread(sc_file)[[1]] else readLines(sc_file)
  SIGS$SCANI_myCAF <- list(genes = unique(trimws(x[nzchar(trimws(x))])), weights = NULL); cat("SCANI genes from:", sc_file, "\n")
} else cat("!! SCANI 基因文件不唯一或未找到，SCANI 跳过\n")
CLS <- c(setNames(reg$class, reg$sig_id), SCANI_myCAF = "A")
NEXP <- c(TCGA_LUAD = 443, GSE68465 = 442, GSE72094 = 398, GSE31210 = 204, GSE50081 = 127, GSE30219 = 83)
out <- list()
for (a in COHORTS) {
  co <- tryCatch(suppressMessages(load_cohort(a, verbose = FALSE)), error = function(e) NULL)
  if (is.null(co) || !isTRUE(co$has_survival)) next
  an <- norm(a); if (!an %in% names(NEXP)) next
  stopifnot(co$n == NEXP[[an]], nrow(co$clin) == co$n)
  pool <- build_pool(co$expr); tt <- co$clin$OS_time; ee <- co$clin$OS_event
  for (sid in intersect(names(SIGS), ref[cohort == an, sig_id])) {
    sg <- SIGS[[sid]]; gp <- intersect(sg$genes, pool$genes)
    s <- if (CLS[[sid]] == "A") score_sets_ssgsea(co$expr, list(gp), prefix = "real")[1, ] else score_sets_weighted(co$expr, list(gp), sg$weights[gp])[1, ]
    ok <- !is.na(tt) & tt > 0 & !is.na(ee) & is.finite(s)
    d <- data.frame(time = tt[ok], event = ee[ok], zs = as.numeric(scale(s[ok])))
    f <- coxph(Surv(time, event) ~ zs, data = d); z <- unname(coef(f) / sqrt(diag(vcov(f))))
    zref <- ref[cohort == an & sig_id == sid, z_real]
    if (abs(z - zref) > 1e-6) stop(sprintf("z 不一致：%s %s 重算 %.6f vs 主结果 %.6f", an, sid, z, zref))
    ph <- cox.zph(f, transform = "km"); p <- ph$table["zs", "p"]
    sp <- survSplit(Surv(time, event) ~ zs, data = d, cut = 3 * 365.25, episode = "period")
    f2 <- tryCatch(coxph(Surv(tstart, time, event) ~ zs:strata(period), data = sp), error = function(e) NULL)
    hr2 <- if (is.null(f2)) c(NA, NA) else unname(exp(coef(f2)))
    out[[length(out) + 1]] <- data.table(cohort = an, sig_id = sid, class = CLS[[sid]], n = nrow(d), events = sum(d$event),
      z = z, hr = exp(unname(coef(f))), ph_chisq = ph$table["zs", "chisq"], ph_p = p, hr_le3y = hr2[1], hr_gt3y = hr2[2])
  }
  cat(sprintf("  %-10s done (%d 组合)\n", an, sum(vapply(out, function(r) r$cohort == an, logical(1)))))
}
R <- rbindlist(out); setorder(R, cohort, sig_id); fwrite(R, file.path(OD, "ph_results.csv"))
cat("\n检验的组合数:", nrow(R), "（预期 77；若跳过 SCANI 则 71）\nz 值与主结果全部一致（差 < 1e-6）\n")
cat("cox.zph p < 0.05 的组合数:", sum(R$ph_p < 0.05), "\n"); print(R[ph_p < 0.05][order(ph_p)], digits = 3)
cat("\n按类别：\n"); print(R[, .(n = .N, ph_viol = sum(ph_p < 0.05)), by = class])
cat("\n按队列：\n"); print(R[, .(n = .N, ph_viol = sum(ph_p < 0.05)), by = cohort])
cat("\ndone:", format(Sys.time()), "\n"); sink()
