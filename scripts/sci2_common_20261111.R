# sci2_common_20261111.R —— 11-11 三项事后分析的共用设置（只读输入）
suppressPackageStartupMessages({ library(data.table); library(survival) })
ROOT   <- path.expand("~/thesis_LUAD_CAF")
CHDIR  <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGDIR <- file.path(ROOT, "03_signature_registry")
setwd(CHDIR)
suppressMessages({ source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R")); source("scripts/null_engine.R") })
REG_F <- file.path(REGDIR, "signature_registry_locked_20260903.csv")
stopifnot(unname(tools::md5sum(REG_F)) == "4ad3218931d95dfb6ee67a7ba07b9a75")
reg <- fread(REG_F)
read_sig <- function(sid, cls) {
  g <- trimws(readLines(file.path(REGDIR, "genes", paste0(sid, "_genes.txt")))); g <- g[nzchar(g)]; w <- NULL
  if (cls == "B") { d <- fread(file.path(REGDIR, "genes", paste0(sid, "_weights.csv")))
    w <- setNames(as.numeric(d$weight), d$gene)[g]; stopifnot(!anyNA(w)) }
  list(genes = g, weights = w) }
SIGS <- setNames(lapply(seq_len(nrow(reg)), function(i) read_sig(reg$sig_id[i], reg$class[i])), reg$sig_id)
sc_f <- file.path(ROOT, "05_manuscripts/SCI-1_IJMS/tables/S1_Table_myCAF_signature_50genes.csv")
sc_d <- fread(sc_f); sc_col <- which(vapply(sc_d, function(v) is.character(v) && all(grepl("^[A-Za-z0-9.-]+$", v)), TRUE))[1]
scani <- unique(trimws(sc_d[[sc_col]])); stopifnot(length(scani) == 50)
SIGS$SCANI_myCAF <- list(genes = scani, weights = NULL)
CLS  <- c(setNames(reg$class, reg$sig_id), SCANI_myCAF = "A")
NEXP <- c(TCGA_LUAD = 443, GSE68465 = 442, GSE72094 = 398, GSE31210 = 204, GSE50081 = 127, GSE30219 = 83)
get_cohort <- function(a) { co <- suppressMessages(load_cohort(a, verbose = FALSE))
  stopifnot(isTRUE(co$has_survival), co$n == NEXP[[a]], nrow(co$clin) == co$n, !anyDuplicated(colnames(co$expr))); co }
score_one <- function(expr, genes, sid) if (CLS[[sid]] == "A") score_sets_ssgsea(expr, list(genes), prefix = "real")[1, ] else
  score_sets_weighted(expr, list(genes), SIGS[[sid]]$weights[genes])[1, ]
score_sets <- function(expr, sets, sid, gp, bp) if (CLS[[sid]] == "A") score_sets_ssgsea(expr, sets, prefix = "n", bp = bp) else
  score_sets_weighted(expr, sets, SIGS[[sid]]$weights[gp])
BP <- if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::MulticoreParam(workers = 4) else NULL
# 主分析结果（只读）：71 个注册签名组合 + 6 个 SCANI 组合
read_main <- function() {
  f1 <- list.files("data_derived/null_benchmark", "\\.rds$", full.names = TRUE)
  f2 <- list.files("data_derived/null_benchmark_scani_supp", "\\.rds$", full.names = TRUE)
  rbindlist(lapply(c(f1, f2), function(f) { x <- readRDS(f)
    if (!isTRUE(x$evaluable)) return(NULL)
    data.table(cohort = x$cohort, sig_id = x$sig_id, class = if (is.null(x$class)) "A" else x$class,
               z_real = x$z_real, p_emp_A = x$p_emp_A, p_emp_B = if (is.null(x$p_emp_B)) NA_real_ else x$p_emp_B,
               qc = isTRUE(x$coexpr_qc_pass), file = basename(f)) }), fill = TRUE)
}
