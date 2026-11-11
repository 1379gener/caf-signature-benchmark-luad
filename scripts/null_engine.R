#!/usr/bin/env Rscript
# ==============================================================================
# null_engine.R —— 零分布引擎（函数库，供 run_benchmark.R source）
#
# 实现 SAP v1.1 §5 的三种零分布：
#   零分布 A（§5.1）朴素三重匹配随机基因集    —— 尺寸/平均表达/方差 匹配
#   零分布 B（§5.2）共表达结构匹配随机模块    —— 另加模块内共表达强度匹配
#   加权复合零分布（§5.3）                    —— B 类小型加权模型专用
#
# 关键实现约束：
#   - 分箱统计量在 log2 尺度上计算（§5.1 v1.1 补充条款）
#   - 每「队列 × signature × 零分布类型」只调用一次 GSVA，1000 个基因集打包
#   - 全部抽样使用显式 seed，seed 规则在 make_seed() 中定死
#   - 零分布 B 有强制质量验证，不通过则该 signature 的 B 类结论作废
#
# 本文件只定义函数，不执行任何分析。
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(survival); library(GSVA)
})

# ============================== 常量 ==========================================
ENGINE <- list(
  version      = "1.0",
  n_null       = 1000L,      # §5.1/5.2 重抽次数
  n_bins       = 10L,        # 平均表达与方差各分 10 个等频箱
  detect_frac  = 0.50,       # 候选池：在 >= 50% 样本中表达非零
  match_tol    = 0.20,       # §5.2 模块整体匹配容差 ±20%
  match_tol_relax = 0.35,    # 重抽超限后放宽至 ±35%
  max_resample = 50L,        # §5.2 重抽上限
  coexpr_min_ratio = 0.70,   # §5.2 强制验证：模块内平均相关须 >= 真实值的 70%
  seed_base    = 20260830L
)

#' seed 规则：由队列名、signature id、零分布类型唯一确定，不得事后更换
make_seed <- function(cohort, sig_id, null_type) {
  s <- paste(cohort, sig_id, null_type, sep = "|")
  ENGINE$seed_base + sum(utf8ToInt(s) * seq_along(utf8ToInt(s))) %% 100000L
}

# ============================== 候选基因池与分箱 ==============================
#' 构建候选基因池并做联合分箱
#' @param expr 已按协议处理的表达矩阵（log2 尺度，基因 x 样本）
#' @return list(genes, mean, var, bin, bin_mean, bin_var)
build_pool <- function(expr) {
  detected <- rowMeans(expr > 0, na.rm = TRUE) >= ENGINE$detect_frac
  m <- expr[detected, , drop = FALSE]
  mu <- rowMeans(m, na.rm = TRUE)
  vv <- apply(m, 1, var, na.rm = TRUE)
  ok <- is.finite(mu) & is.finite(vv) & vv > 0
  mu <- mu[ok]; vv <- vv[ok]
  q <- function(x) {
    b <- cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = ENGINE$n_bins + 1),
                                  na.rm = TRUE), include.lowest = TRUE, labels = FALSE)
    b[is.na(b)] <- 1L; as.integer(b)
  }
  bm <- q(mu); bv <- q(vv)
  list(genes = names(mu), mean = mu, var = vv,
       bin = bm * 100L + bv,      # 联合箱 id
       bin_mean = bm, bin_var = bv)
}

#' 预先计算秩标准化矩阵（Spearman 相关用）
#' cor_spearman(a, b) = crossprod(Zr[a,], Zr[b,]) / (n-1)
build_rank_z <- function(expr, pool) {
  m <- expr[pool$genes, , drop = FALSE]
  r <- t(apply(m, 1, rank, na.last = "keep"))
  r <- t(scale(t(r)))                       # 行标准化
  r[!is.finite(r)] <- 0
  r
}

# ============================== 零分布 A ======================================
#' §5.1 朴素三重匹配随机基因集
#' @return list(sets = list of 1000 character vectors, n_fallback = int)
sample_null_A <- function(pool, real_genes, n_null = ENGINE$n_null, seed) {
  set.seed(seed)
  real <- intersect(real_genes, pool$genes)
  k <- length(real)
  if (k < 1L) stop("真实 signature 在候选池中无基因")
  idx_real <- match(real, pool$genes)
  bins_real <- pool$bin[idx_real]

  # 预建：每个联合箱内的可用基因下标（排除真实 signature 自身）
  keep  <- setdiff(seq_along(pool$genes), idx_real)
  avail <- split(keep, pool$bin[keep])

  # 相邻平均表达箱的退让映射
  neighbor_bins <- function(b) {
    bm <- b %/% 100L; bv <- b %% 100L
    cand <- c((bm - 1L) * 100L + bv, (bm + 1L) * 100L + bv,
              (bm - 2L) * 100L + bv, (bm + 2L) * 100L + bv)
    cand[cand > 100L & cand < 1100L]
  }

  n_fb <- 0L
  sets <- vector("list", n_null)
  for (i in seq_len(n_null)) {
    used <- integer(0)
    picked <- integer(k)
    for (j in seq_len(k)) {
      b <- bins_real[j]
      cand <- setdiff(avail[[as.character(b)]], used)
      if (!length(cand)) {                       # 退让至相邻平均表达箱
        for (nb in neighbor_bins(b)) {
          cand <- setdiff(avail[[as.character(nb)]], used)
          if (length(cand)) { n_fb <- n_fb + 1L; break }
        }
      }
      if (!length(cand)) cand <- setdiff(seq_along(pool$genes), c(idx_real, used))
      picked[j] <- if (length(cand) == 1L) cand else sample(cand, 1L)
      used <- c(used, picked[j])
    }
    sets[[i]] <- pool$genes[picked]
  }
  list(sets = sets, n_fallback = n_fb, k = k, real_used = real)
}

# ============================== 零分布 B ======================================
#' 模块内平均两两 Spearman 相关
mean_pairwise_cor <- function(rz, genes) {
  i <- match(intersect(genes, rownames(rz)), rownames(rz))
  if (length(i) < 2L) return(NA_real_)
  n <- ncol(rz)
  C <- tcrossprod(rz[i, , drop = FALSE]) / (n - 1)
  mean(C[upper.tri(C)], na.rm = TRUE)
}

#' §5.2 共表达结构匹配随机模块
sample_null_B <- function(pool, rz, real_genes, n_null = ENGINE$n_null, seed) {
  set.seed(seed)
  real <- intersect(real_genes, pool$genes)
  k <- length(real)
  idx_real <- match(real, pool$genes)
  bins_real <- pool$bin[idx_real]
  target_mean <- mean(pool$mean[idx_real])
  target_var  <- mean(pool$var[idx_real])
  n <- ncol(rz)

  keep <- setdiff(seq_along(pool$genes), idx_real)
  avail <- split(keep, pool$bin[keep])

  sets <- vector("list", n_null)
  n_relax <- 0L; n_try_tot <- 0L
  for (i in seq_len(n_null)) {
    tol <- ENGINE$match_tol
    got <- NULL
    for (tries in seq_len(ENGINE$max_resample + 1L)) {
      n_try_tot <- n_try_tot + 1L
      if (tries > ENGINE$max_resample) { tol <- ENGINE$match_tol_relax; n_relax <- n_relax + 1L }
      # 种子基因：与真实 signature 中随机一个基因同联合箱
      b <- bins_real[sample.int(k, 1L)]
      cand <- avail[[as.character(b)]]
      if (!length(cand)) next
      s <- if (length(cand) == 1L) cand else sample(cand, 1L)
      # 与候选池全部基因的 Spearman 相关，取 top (k-1)
      cr <- as.numeric(rz %*% rz[s, ]) / (n - 1)
      cr[s] <- -Inf
      cr[idx_real] <- -Inf                        # 不取真实 signature 自身基因
      top <- order(cr, decreasing = TRUE)[seq_len(k - 1L)]
      mod <- c(s, top)
      # 模块整体的平均表达与平均方差是否落在真实值 ±tol 内
      mm <- mean(pool$mean[mod]); mv <- mean(pool$var[mod])
      if (abs(mm - target_mean) <= tol * abs(target_mean) &&
          abs(mv - target_var)  <= tol * abs(target_var)) { got <- mod; break }
      if (tries > ENGINE$max_resample) { got <- mod; break }   # 放宽后仍取该次
    }
    if (is.null(got)) got <- sample(keep, k)      # 极端兜底，会在日志中体现
    sets[[i]] <- pool$genes[got]
  }

  # ---- 强制质量验证（§5.2）
  real_coexpr <- mean_pairwise_cor(rz, real)
  null_coexpr <- vapply(sets[seq_len(min(200L, n_null))],
                        function(g) mean_pairwise_cor(rz, g), 0)
  med_null <- median(null_coexpr, na.rm = TRUE)
  ratio <- med_null / real_coexpr
  pass <- is.finite(ratio) && ratio >= ENGINE$coexpr_min_ratio

  list(sets = sets, k = k, n_relax = n_relax, mean_tries = n_try_tot / n_null,
       real_coexpr = real_coexpr, null_coexpr_median = med_null,
       coexpr_ratio = ratio, qc_pass = pass)
}

# ============================== 加权复合零分布 ================================
#' §5.3 B 类小型加权模型：随机抽同数量、联合箱匹配的基因，套用同一组权重
sample_null_weighted <- function(pool, real_genes, weights, n_null = ENGINE$n_null,
                                 seed, coexpr = FALSE, rz = NULL) {
  set.seed(seed)
  real <- intersect(real_genes, pool$genes)
  w <- weights[match(real, real_genes)]
  k <- length(real)
  idx_real <- match(real, pool$genes)
  bins_real <- pool$bin[idx_real]
  keep <- setdiff(seq_along(pool$genes), idx_real)
  avail <- split(keep, pool$bin[keep])

  sets <- vector("list", n_null)
  for (i in seq_len(n_null)) {
    used <- integer(0); picked <- integer(k)
    for (j in seq_len(k)) {
      cand <- setdiff(avail[[as.character(bins_real[j])]], used)
      if (!length(cand)) cand <- setdiff(keep, used)
      picked[j] <- if (length(cand) == 1L) cand else sample(cand, 1L)
      used <- c(used, picked[j])
    }
    sets[[i]] <- pool$genes[picked]
  }
  list(sets = sets, weights = w, k = k, real_used = real)
}

# ============================== 打分 ==========================================
ENGINE_SSNORM <- function() FALSE     # 与 scoring_protocol.R 的 PROTOCOL 一致

#' A 类：一次 GSVA 调用，1000 个基因集打包
#' @return 矩阵 n_sets x n_samples
score_sets_ssgsea <- function(expr, sets, prefix = "null", bp = NULL) {
  gsl <- setNames(sets, sprintf("%s_%04d", prefix, seq_along(sets)))
  gsl <- lapply(gsl, function(g) intersect(g, rownames(expr)))
  gsl <- gsl[vapply(gsl, length, 0L) > 0L]
  args <- list(expr, gsl, method = "ssgsea",
               ssgsea.norm = ENGINE_SSNORM(), verbose = FALSE)
  if (!is.null(bp)) args$BPPARAM <- bp
  do.call(GSVA::gsva, args)
}

#' B 类：加权和（各基因先在队列内 Z 标准化）
score_sets_weighted <- function(expr, sets, weights) {
  z <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s }
  t(vapply(sets, function(g) {
    i <- match(g, rownames(expr))
    ok <- !is.na(i)
    if (!any(ok)) return(rep(NA_real_, ncol(expr)))
    as.numeric(crossprod(weights[ok], t(apply(expr[i[ok], , drop = FALSE], 1, z))))
  }, numeric(ncol(expr))))
}

# ============================== Cox 与经验 p ==================================
#' 对一批分数逐个拟合单变量 Cox，返回 Wald z
cox_z_batch <- function(score_mat, time, event) {
  ok <- !is.na(time) & time > 0 & !is.na(event)
  tt <- time[ok]; ee <- event[ok]
  sv <- Surv(tt, ee)
  apply(score_mat[, ok, drop = FALSE], 1, function(s) {
    if (!is.finite(sd(s, na.rm = TRUE)) || sd(s, na.rm = TRUE) == 0) return(NA_real_)
    f <- tryCatch(coxph(sv ~ scale(s)), error = function(e) NULL)
    if (is.null(f)) NA_real_ else unname(summary(f)$coefficients[1, "z"])
  })
}

cox_one <- function(score, time, event) {
  ok <- !is.na(time) & time > 0 & !is.na(event) & is.finite(score)
  f <- coxph(Surv(time[ok], event[ok]) ~ scale(score[ok]))
  s <- summary(f)
  list(z = unname(s$coefficients[1, "z"]),
       hr = unname(s$conf.int[1, 1]), lo = unname(s$conf.int[1, 3]),
       hi = unname(s$conf.int[1, 4]), p = unname(s$coefficients[1, 5]),
       n = sum(ok), events = sum(event[ok]))
}

#' §4.3 经验 p 值（双侧）
empirical_p <- function(z_real, z_null) {
  z_null <- z_null[is.finite(z_null)]
  if (!length(z_null) || !is.finite(z_real)) return(NA_real_)
  (1 + sum(abs(z_null) >= abs(z_real))) / (1 + length(z_null))
}


# ============================== 清单锁定闸（SAP §3.6 / §11.1-2）===============
#' 引擎拒绝在 signature 清单锁定之前运行。
#' SAP §3 规定：清单必须在任何零分布运行之前完全锁定，锁定后不得增删。
#' 该规则由本函数强制，不依赖人工遵守。
require_locked_registry <- function(
    dir = "~/thesis_LUAD_CAF/03_signature_registry") {
  d <- path.expand(dir)
  f <- list.files(d, pattern = "^signature_registry_locked_[0-9]{8}\\.csv$",
                  full.names = TRUE)
  if (!length(f)) {
    stop("\n  signature 清单尚未锁定，引擎拒绝运行。\n",
         "  SAP §3.6：清单填写完毕后须导出为 signature_registry_locked_YYYYMMDD.csv，\n",
         "  计算 md5 并记入 SAP §14；此后分析只能读取该锁定文件。\n",
         "  当前目录: ", d, "\n", call. = FALSE)
  }
  f <- sort(f, decreasing = TRUE)[1]
  list(file = f, md5 = unname(tools::md5sum(f)), table = data.table::fread(f))
}

# ============================== 自检 ==========================================
#' 用模拟数据快速验证三种零分布的行为是否符合预期。
#' 不读取任何真实数据，不受锁定闸限制。
engine_selftest <- function(n_gene = 3000, n_samp = 200, k = 20, n_null = 100, seed = 1) {
  set.seed(seed)
  hr <- function(t) cat("\n", strrep("-", 62), "\n", t, "\n", strrep("-", 62), "\n", sep = "")
  cat("null_engine 自检  v", ENGINE$version, "\n", sep = "")

  # --- 造一个含共表达模块结构的表达矩阵
  base <- matrix(rnorm(n_gene * n_samp, 6, 1.5), n_gene, n_samp)
  for (b in 1:30) {                     # 30 个共表达模块，每个 40 基因
    i <- ((b - 1) * 40 + 1):(b * 40)
    f <- rnorm(n_samp, 0, 1.2)
    base[i, ] <- base[i, ] + matrix(rep(f, each = length(i)), length(i))
  }
  base[base < 0] <- 0
  rownames(base) <- sprintf("G%04d", seq_len(n_gene))
  colnames(base) <- sprintf("S%03d", seq_len(n_samp))

  hr("1. 候选池与分箱")
  pool <- build_pool(base)
  cat("  候选基因数        ", length(pool$genes), "\n")
  cat("  联合箱数          ", length(unique(pool$bin)), " （上限 100）\n")
  tb <- table(pool$bin)
  cat("  每箱基因数 中位/最小 ", median(tb), " / ", min(tb), "\n", sep = "")

  rz <- build_rank_z(base, pool)
  cat("  秩标准化矩阵      ", paste(dim(rz), collapse = " x "), "\n")

  # 真实 signature：取自一个共表达模块（模拟真实 signature 的共表达特性）
  real <- rownames(base)[1:k]

  hr("2. 零分布 A（朴素三重匹配）")
  sd <- make_seed("TEST", "SIG-000", "A")
  a <- sample_null_A(pool, real, n_null = n_null, seed = sd)
  cat("  seed              ", sd, "\n")
  cat("  生成集合数        ", length(a$sets), "  每集尺寸 ", length(a$sets[[1]]), "\n", sep = "")
  cat("  退让次数          ", a$n_fallback, "\n")
  cat("  集合内有无重复基因 ", ifelse(any(vapply(a$sets, function(s) anyDuplicated(s) > 0, TRUE)),
      "!! 有", "无"), "\n")
  cat("  与真实 signature 的交集（应为 0）: ",
      sum(vapply(a$sets, function(s) length(intersect(s, real)), 0L)), "\n")
  mm <- vapply(a$sets, function(s) mean(pool$mean[match(s, pool$genes)]), 0)
  cat(sprintf("  平均表达匹配: 真实 %.3f  零分布中位 %.3f  偏差 %.1f%%\n",
      mean(pool$mean[match(real, pool$genes)]), median(mm),
      100 * (median(mm) - mean(pool$mean[match(real, pool$genes)])) /
        mean(pool$mean[match(real, pool$genes)])))

  hr("3. 零分布 B（共表达结构匹配）")
  sdb <- make_seed("TEST", "SIG-000", "B")
  b <- sample_null_B(pool, rz, real, n_null = n_null, seed = sdb)
  cat("  seed              ", sdb, "\n")
  cat("  放宽容差次数      ", b$n_relax, "  平均尝试次数 ", round(b$mean_tries, 2), "\n", sep = "")
  cat(sprintf("  真实 signature 模块内平均相关  %.4f\n", b$real_coexpr))
  cat(sprintf("  零分布 B 模块内平均相关 中位   %.4f\n", b$null_coexpr_median))
  cat(sprintf("  比值 %.3f （阈值 %.2f）  ->  %s\n", b$coexpr_ratio,
      ENGINE$coexpr_min_ratio, ifelse(b$qc_pass, "质量验证通过", "!! 未通过，B 类结论作废")))
  ca <- median(vapply(a$sets[1:min(50, n_null)], function(g) mean_pairwise_cor(rz, g), 0))
  cat(sprintf("  对照：零分布 A 模块内平均相关 中位 %.4f\n", ca))
  cat("  -> B 的共表达强度应显著高于 A，这正是设立零分布 B 的理由\n")

  hr("4. 加权复合零分布（B 类）")
  w <- round(rnorm(5, 0, 0.3), 3)
  rg <- rownames(base)[c(1, 45, 90, 200, 500)]
  sdw <- make_seed("TEST", "SIG-000", "W")
  wn <- sample_null_weighted(pool, rg, w, n_null = n_null, seed = sdw)
  cat("  基因/权重         ", paste(sprintf("%s=%.3f", rg, w), collapse = "  "), "\n")
  cat("  生成集合数        ", length(wn$sets), "  每集尺寸 ", wn$k, "\n", sep = "")
  sc <- score_sets_weighted(base, wn$sets, wn$weights)
  cat("  加权分数矩阵      ", paste(dim(sc), collapse = " x "), "\n")

  hr("5. Cox 与经验 p")
  time <- rexp(n_samp, 1 / 800); event <- rbinom(n_samp, 1, 0.4)
  real_sc <- score_sets_weighted(base, list(rg), w)[1, ]
  r <- cox_one(real_sc, time, event)
  zn <- cox_z_batch(sc, time, event)
  p <- empirical_p(r$z, zn)
  cat(sprintf("  真实 z = %.3f   HR %.3f (%.3f-%.3f)  n=%d 事件=%d\n",
      r$z, r$hr, r$lo, r$hi, r$n, r$events))
  cat(sprintf("  零分布 z: 中位 %.3f  |z| 90%%分位 %.3f  有效 %d/%d\n",
      median(zn, na.rm = TRUE), quantile(abs(zn), .9, na.rm = TRUE),
      sum(is.finite(zn)), length(zn)))
  cat(sprintf("  经验 p = %.4f   （随机数据下应接近均匀分布，即 p 无系统性偏小）\n", p))
  cat(sprintf("  p 的理论下界 = 1/(1+%d) = %.4f\n", n_null, 1 / (1 + n_null)))

  hr("6. 锁定闸")
  g <- tryCatch(require_locked_registry(), error = function(e) conditionMessage(e))
  if (is.character(g)) cat("  按预期拒绝运行：\n", g, "\n", sep = "") else
    cat("  已找到锁定文件: ", basename(g$file), "  md5 ", g$md5, "\n", sep = "")

  hr("自检结束")
  cat("  以上仅为模拟数据下的行为验证，不构成任何研究结论。\n")
  invisible(list(pool = pool, A = a, B = b, W = wn))
}

if (sys.nframe() == 0L) engine_selftest()
