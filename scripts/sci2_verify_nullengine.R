#!/usr/bin/env Rscript
# ==============================================================================
# sci2_verify_nullengine.R —— 零分布引擎的独立复核（对照表 ⚠️ E / F / G / H）
#   用法：Rscript sci2_verify_nullengine.R [--workers=6]
#   H  SAP §12：随机抽 20 个零分布成员，从表达矩阵独立打分、独立拟合 Cox，与 RDS 的 z_null 比对（1e-6）
#   R  重现性闸门：重抽的零分布 A 与引擎 sample_null_A 相同；零分布 B 前 200 模块共表达中位数与 RDS 相同
#   E  零分布 A 退让分层计数（相邻箱 / 全池兜底）；零分布 B 第 51 次直接采用数及其中超出 ±35% 数
#   F  共表达 QC 用全部 1000 个模块重算，与前 200 个的判定比对
#   G  B 类：sample_null_weighted 与 sample_null_A 在同一 seed 下基因集与权重是否相同
#   输出（data_derived/，已存在即停止）：verify_nullengine_combo_<日期>.csv、verify_nullengine_members_<日期>.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(GSVA); library(parallel) })
arg <- commandArgs(TRUE)
WORKERS <- suppressWarnings(as.integer(sub("^--workers=", "", grep("^--workers=", arg, value = TRUE)[1])))
if (is.na(WORKERS)) WORKERS <- 6L
ROOT <- path.expand("~/thesis_LUAD_CAF"); CH <- file.path(ROOT, "02_chapters/ch04_benchmark")
REGD <- file.path(ROOT, "03_signature_registry"); DD <- file.path(CH, "data_derived")
STAMP <- format(Sys.Date(), "%Y%m%d")
F_COMBO <- file.path(DD, sprintf("verify_nullengine_combo_%s.csv", STAMP))
F_MEM   <- file.path(DD, sprintf("verify_nullengine_members_%s.csv", STAMP))
for (f in c(F_COMBO, F_MEM)) if (file.exists(f)) stop("输出已存在，不覆盖：", f)
dir.create(file.path(CH, "logs"), showWarnings = FALSE)
LOG <- file.path(CH, "logs", format(Sys.time(), "verify_nullengine_%Y%m%d_%H%M.log"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
CHK <- list()
check <- function(id, what, expected, observed, tol = 0) {
  ok <- isTRUE(all(abs(expected - observed) <= tol))
  CHK[[length(CHK) + 1]] <<- data.table(id, what, expected, observed, status = ifelse(ok, "PASS", "FAIL"))
  cat(sprintf("  [%s] %-5s %-58s 期望 %-8s 实得 %s\n", ifelse(ok, "PASS", "FAIL"), id, what, expected, observed))
}
cat("verify_nullengine  ", format(Sys.time()), "  workers ", WORKERS, "\n日志 ", LOG, "\n", sep = "")
suppressMessages({
  source(file.path(ROOT, "01_bases/A_bulk_OS/scripts/scoring_protocol.R"))
  source(file.path(CH, "scripts/null_engine.R"))
})
cat("null_engine.R md5 ", unname(tools::md5sum(file.path(CH, "scripts/null_engine.R"))), "\n", sep = "")
cat("ENGINE: match_tol ", ENGINE$match_tol, " / relax ", ENGINE$match_tol_relax, " / max_resample ",
    ENGINE$max_resample, " / coexpr_min_ratio ", ENGINE$coexpr_min_ratio, "\n", sep = "")

hr("0. 队列")
CO <- list(); POOL <- list(); RZ <- list()
for (a in COHORTS) {
  co <- load_cohort(a, verbose = FALSE)
  stopifnot(!anyDuplicated(colnames(co$expr)), ncol(co$expr) == co$n, co$n == as.integer(co$man$n_final[1]))
  if (!co$has_survival) next
  stopifnot(nrow(co$clin) == co$n, co$n_events == as.integer(co$man$n_events[1]))
  CO[[a]] <- co; POOL[[a]] <- build_pool(co$expr); RZ[[a]] <- build_rank_z(co$expr, POOL[[a]])
  cat(sprintf("  %-10s n %4d  事件 %3d  候选池 %5d\n", a, co$n, co$n_events, length(POOL[[a]]$genes)))
}
read_sig <- function(sid, cls) {
  g <- trimws(readLines(file.path(REGD, "genes", paste0(sid, "_genes.txt")))); g <- g[nzchar(g)]
  w <- NULL
  if (cls == "B") { d <- fread(file.path(REGD, "genes", paste0(sid, "_weights.csv")))
                    w <- setNames(as.numeric(d$weight), d$gene)[g]; stopifnot(!anyNA(w)) }
  list(genes = g, weights = w)
}
rds <- list.files(file.path(DD, "null_benchmark"), "\\.rds$", full.names = TRUE)
META <- rbindlist(lapply(rds, function(f) { x <- readRDS(f)
  data.table(file = f, cohort = x$cohort, sig_id = x$sig_id, class = x$class, evaluable = isTRUE(x$evaluable),
             has_null = !is.null(x$z_null_A) && !is.null(x$z_null_B)) }))
E <- META[evaluable == TRUE & has_null == TRUE]
check("R0", "可评估且含零分布的组合数", 71, nrow(E))

null_A_trace <- function(pool, real_genes, n_null, seed) {
  set.seed(seed)
  real <- intersect(real_genes, pool$genes); k <- length(real)
  idx_real <- match(real, pool$genes); bins_real <- pool$bin[idx_real]
  keep <- setdiff(seq_along(pool$genes), idx_real); avail <- split(keep, pool$bin[keep])
  nb_bins <- function(b) { bm <- b %/% 100L; bv <- b %% 100L
    cand <- c((bm - 1L) * 100L + bv, (bm + 1L) * 100L + bv, (bm - 2L) * 100L + bv, (bm + 2L) * 100L + bv)
    cand[cand > 100L & cand < 1100L] }
  fb1 <- 0L; fb2 <- 0L; sets <- vector("list", n_null)
  for (i in seq_len(n_null)) {
    used <- integer(0); picked <- integer(k)
    for (j in seq_len(k)) {
      b <- bins_real[j]; cand <- setdiff(avail[[as.character(b)]], used)
      if (!length(cand)) for (nb in nb_bins(b)) { cand <- setdiff(avail[[as.character(nb)]], used)
                                                  if (length(cand)) { fb1 <- fb1 + 1L; break } }
      if (!length(cand)) { cand <- setdiff(seq_along(pool$genes), c(idx_real, used)); fb2 <- fb2 + 1L }
      picked[j] <- if (length(cand) == 1L) cand else sample(cand, 1L)
      used <- c(used, picked[j])
    }
    sets[[i]] <- pool$genes[picked]
  }
  list(sets = sets, fb1 = fb1, fb2 = fb2)
}
null_B_trace <- function(pool, rz, real_genes, n_null, seed) {
  set.seed(seed)
  real <- intersect(real_genes, pool$genes); k <- length(real)
  idx_real <- match(real, pool$genes); bins_real <- pool$bin[idx_real]
  tm <- mean(pool$mean[idx_real]); tv <- mean(pool$var[idx_real]); n <- ncol(rz)
  keep <- setdiff(seq_along(pool$genes), idx_real); avail <- split(keep, pool$bin[keep])
  sets <- vector("list", n_null)
  tries_v <- integer(n_null); relaxed <- logical(n_null); extreme <- logical(n_null)
  dm <- numeric(n_null); dv <- numeric(n_null)
  for (i in seq_len(n_null)) {
    tol <- ENGINE$match_tol; got <- NULL
    for (tries in seq_len(ENGINE$max_resample + 1L)) {
      tries_v[i] <- tries
      if (tries > ENGINE$max_resample) { tol <- ENGINE$match_tol_relax; relaxed[i] <- TRUE }
      b <- bins_real[sample.int(k, 1L)]
      cand <- avail[[as.character(b)]]
      if (!length(cand)) next
      s <- if (length(cand) == 1L) cand else sample(cand, 1L)
      cr <- as.numeric(rz %*% rz[s, ]) / (n - 1); cr[s] <- -Inf; cr[idx_real] <- -Inf
      mod <- c(s, order(cr, decreasing = TRUE)[seq_len(k - 1L)])
      mm <- mean(pool$mean[mod]); mv <- mean(pool$var[mod])
      if (abs(mm - tm) <= tol * abs(tm) && abs(mv - tv) <= tol * abs(tv)) { got <- mod; break }
      if (tries > ENGINE$max_resample) { got <- mod; break }
    }
    if (is.null(got)) { got <- sample(keep, k); extreme[i] <- TRUE }
    sets[[i]] <- pool$genes[got]
    dm[i] <- abs(mean(pool$mean[got]) - tm) / abs(tm); dv[i] <- abs(mean(pool$var[got]) - tv) / abs(tv)
  }
  list(sets = sets, tries = tries_v, relaxed = relaxed, extreme = extreme, dev_mean = dm, dev_var = dv)
}

set.seed(20261025)
GRID <- CJ(row = seq_len(nrow(E)), type = c("A", "B"))
HS <- GRID[sample.int(nrow(GRID), 20L)][, member := sample.int(1000L, 20L, replace = TRUE)]
HS[, `:=`(cohort = E$cohort[row], sig_id = E$sig_id[row], class = E$class[row])]
cat("\nH 抽中的 20 个成员（抽样 seed 20261025）：\n"); print(HS[, .(cohort, sig_id, class, type, member)], row.names = FALSE)

hr(sprintf("3. 逐组合重抽（%d 个组合，%d 进程）", nrow(E), WORKERS))
t0 <- Sys.time()
one <- function(r) {
  x <- readRDS(E$file[r]); a <- x$cohort; sid <- x$sig_id; cls <- x$class
  pool <- POOL[[a]]; rz <- RZ[[a]]; sg <- read_sig(sid, cls)
  gp <- intersect(sg$genes, pool$genes); k <- length(gp)
  sA <- make_seed(a, sid, "A"); sB <- make_seed(a, sid, "B")
  A  <- null_A_trace(pool, gp, 1000L, sA)
  eng_A <- sample_null_A(pool, gp, n_null = 1000L, seed = sA)
  g_sets <- g_w <- NA
  if (cls == "B") {
    W <- sample_null_weighted(pool, gp, sg$weights[gp], n_null = 1000L, seed = sA)
    g_sets <- identical(W$sets, A$sets); g_w <- isTRUE(all.equal(unname(W$weights), unname(sg$weights[gp])))
  }
  B  <- null_B_trace(pool, rz, gp, 1000L, sB)
  rc <- mean_pairwise_cor(rz, gp)
  nc <- vapply(B$sets, function(g) mean_pairwise_cor(rz, g), 0)
  m200 <- median(nc[1:200], na.rm = TRUE); m1000 <- median(nc, na.rm = TRUE)
  q <- function(m) { rt <- m / rc; is.finite(rt) && rt >= ENGINE$coexpr_min_ratio }
  row <- data.table(cohort = a, sig_id = sid, class = cls, k = k,
    seedA_ok = identical(as.numeric(sA), as.numeric(x$seed_A)), seedB_ok = identical(as.numeric(sB), as.numeric(x$seed_B)),
    A_same_as_engine = identical(A$sets, eng_A$sets),
    A_fb_neighbor = A$fb1, A_fb_neighbor_rds = as.integer(x$n_fallback_A), A_fb_wholepool = A$fb2,
    G_weighted_same_sets = g_sets, G_weighted_same_weights = g_w,
    coexpr_real = rc, coexpr_real_rds = x$coexpr_real, coexpr_med200 = m200, coexpr_med200_rds = x$coexpr_null_med,
    coexpr_med1000 = m1000, ratio200 = m200 / rc, ratio1000 = m1000 / rc,
    qc200 = q(m200), qc1000 = q(m1000), qc_rds = isTRUE(x$coexpr_qc_pass),
    B_n_relax_rds = if (is.null(x$n_relax_B)) NA_integer_ else as.integer(x$n_relax_B),
    B_mean_tries = mean(B$tries), B_n_relaxed = sum(B$relaxed),
    B_relaxed_out20 = sum(B$relaxed & (B$dev_mean > ENGINE$match_tol | B$dev_var > ENGINE$match_tol)),
    B_relaxed_out35 = sum(B$relaxed & (B$dev_mean > ENGINE$match_tol_relax | B$dev_var > ENGINE$match_tol_relax)),
    B_n_extreme = sum(B$extreme), B_max_dev_mean = max(B$dev_mean), B_max_dev_var = max(B$dev_var))
  keepH <- HS[cohort == a & sig_id == sid]
  sets <- NULL
  if (nrow(keepH)) {
    ia <- unique(keepH[type == "A", member]); ib <- unique(keepH[type == "B", member])
    sets <- list(A = setNames(A$sets[ia], ia), B = setNames(B$sets[ib], ib))
  }
  message(sprintf("  完成 %-10s %-8s  %.1f 分钟", a, sid, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(row = row, sets = sets)
}
RES <- mclapply(seq_len(nrow(E)), one, mc.cores = WORKERS, mc.preschedule = FALSE)
bad <- which(vapply(RES, function(z) inherits(z, "try-error") || is.null(z$row), TRUE))
if (length(bad)) { print(RES[bad]); stop("有组合出错，见上") }
C <- rbindlist(lapply(RES, `[[`, "row"))
cat(sprintf("用时 %.1f 分钟\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

hr("4. H：20 个零分布成员的 Cox z 独立复算")
wsum_pos <- function(expr, genes, w) {
  X <- expr[genes, , drop = FALSE]
  Z <- t(apply(X, 1, function(v) { s <- sd(v); if (!is.finite(s) || s == 0) v * 0 else (v - mean(v)) / s }))
  as.numeric(crossprod(unname(w), Z))
}
coxz <- function(sc, t, e) { ok <- !is.na(t) & t > 0 & !is.na(e) & is.finite(sc)
  summary(coxph(Surv(t[ok], e[ok]) ~ as.numeric(scale(sc[ok]))))$coefficients[1, "z"] }
names(RES) <- paste(E$cohort, E$sig_id)
M <- rbindlist(lapply(seq_len(nrow(HS)), function(i) {
  h <- HS[i]; co <- CO[[h$cohort]]; x <- readRDS(E$file[h$row])
  g <- RES[[paste(h$cohort, h$sig_id)]]$sets[[h$type]][[as.character(h$member)]]
  sc <- if (h$class == "A") as.numeric(GSVA::gsva(co$expr, list(m = g), method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)[1, ]) else {
    sg <- read_sig(h$sig_id, h$class); wsum_pos(co$expr, g, sg$weights[intersect(sg$genes, POOL[[h$cohort]]$genes)]) }
  z  <- coxz(sc, co$clin$OS_time, co$clin$OS_event)
  zr <- x[[paste0("z_null_", h$type)]][h$member]
  data.table(h[, .(cohort, sig_id, class, type, member)], k = length(g), z_recomputed = z, z_rds = zr, abs_diff = abs(z - zr))
}))
print(M, row.names = FALSE)

hr("5. 检查")
check("H1", "20 个成员 |z 复算 − z_RDS| > 1e-6 的个数", 0, sum(!(M$abs_diff <= 1e-6)))
check("R1", "seed 与 RDS 记录不一致的组合数", 0, sum(!C$seedA_ok | !C$seedB_ok))
check("R2", "零分布 A 重抽与引擎不同的组合数", 0, sum(!C$A_same_as_engine))
check("R3", "零分布 A 相邻箱退让计数与 RDS 不一致的组合数", 0, sum(C$A_fb_neighbor != C$A_fb_neighbor_rds))
check("R4", "共表达（真实签名）与 RDS 不一致的组合数", 0, sum(abs(C$coexpr_real - C$coexpr_real_rds) > 1e-10, na.rm = TRUE))
check("R5", "前 200 模块共表达中位数与 RDS 不一致的组合数", 0, sum(abs(C$coexpr_med200 - C$coexpr_med200_rds) > 1e-10, na.rm = TRUE))
check("R6", "前 200 模块的 QC 判定与 RDS 不一致的组合数", 0, sum(C$qc200 != C$qc_rds))
check("R7", "第 51 次采用数与 RDS 的 n_relax_B 不一致的组合数", 0, sum(C$B_n_relaxed != C$B_n_relax_rds, na.rm = TRUE))
cat("\n  —— 以下为计数（0 只是「无影响」的参照，不是断言）\n")
check("E1", "零分布 A 全池兜底（引擎不记录）总次数", 0, sum(C$A_fb_wholepool))
check("E2", "零分布 B 第 51 次直接采用的成员总数", 0, sum(C$B_n_relaxed))
check("E3", "  其中偏离超过 ±35% 的成员总数", 0, sum(C$B_relaxed_out35))
check("E4", "零分布 B 极端兜底（随机 k 个）成员总数", 0, sum(C$B_n_extreme))
check("F1", "QC 用 1000 个模块与用前 200 个判定不同的组合数", 0, sum(C$qc1000 != C$qc200))
GB <- C[class == "B"]
check("G1", "B 类：sample_null_weighted 与 sample_null_A 基因集不同的组合数", 0, sum(!GB$G_weighted_same_sets))
check("G2", "B 类：sample_null_weighted 的权重与打分权重不同的组合数", 0, sum(!GB$G_weighted_same_weights))
cat("\n零分布 B 放宽情况（有放宽的组合）：\n")
print(C[B_n_relaxed > 0, .(cohort, sig_id, class, k, B_mean_tries, B_n_relaxed, B_relaxed_out20, B_relaxed_out35,
                           B_max_dev_mean = round(B_max_dev_mean, 3), B_max_dev_var = round(B_max_dev_var, 3))], row.names = FALSE)
cat("\nQC 判定 200 vs 1000 不同的组合：\n")
print(C[qc1000 != qc200, .(cohort, sig_id, class, ratio200 = round(ratio200, 3), ratio1000 = round(ratio1000, 3), qc200, qc1000)], row.names = FALSE)
cat("\n零分布 A 有退让的组合：\n")
print(C[A_fb_neighbor > 0 | A_fb_wholepool > 0, .(cohort, sig_id, class, k, A_fb_neighbor, A_fb_wholepool)], row.names = FALSE)
fwrite(C, F_COMBO); fwrite(M, F_MEM)
K <- rbindlist(CHK); cat(sprintf("\n闸门（H1、R1–R7）：%d / 8 PASS\n", K[id %in% c("H1", paste0("R", 1:7)), sum(status == "PASS")]))
cat("写出：\n  ", F_COMBO, "\n  ", F_MEM, "\n", sep = "")
cat("\n"); print(sessionInfo())
