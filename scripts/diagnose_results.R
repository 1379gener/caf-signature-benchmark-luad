#!/usr/bin/env Rscript
# ==============================================================================
# diagnose_results.R —— 零分布基准结果的诊断与扩展汇总
#
# 解决三件事：
#
# [1] 区分「无预后信号」与「信号不特异」
#     经验 p 大有两种截然不同的原因：
#       (a) 真实 z 接近 0 —— 该 signature 本身就没有预后信号
#       (b) 真实 z 不小，但随机集同样强 —— 信号存在但不特异
#     只报 p_emp 无法区分。本脚本补上零分布 |z| 的分位数与真实 z 并列。
#
# [2] 问题二的敏感性分析
#     coexpr_ratio 指标在真实 signature 内部共表达接近零时不稳定
#     （分母趋零或为负，比值出现负值与荒谬值）。
#     主分析照 SAP 字面执行；此处并列报告【不做 QC 过滤】的 p_emp_B。
#
# [3] 生成可直接用于论文的扩展汇总表
#
# 只读 rds，不修改任何已有结果。
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# R 4.3.1 尚无 base 版 %||%（4.4 才加入），自行定义
`%||%` <- function(a, b) if (is.null(a)) b else a

CH   <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark")
RDIR <- file.path(CH, "data_derived/null_benchmark")
REG  <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")
STAMP <- format(Sys.Date(), "%Y%m%d")
OUT  <- file.path("/tmp", format(Sys.time(), "diag_results_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 74), "\n", t, "\n", strrep("=", 74), "\n", sep = "")
sh <- function(t) cat("\n  --- ", t, "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("零分布基准结果诊断\n时间: ", format(Sys.time()), "\n", sep = "")

# ============================================================ 1. 载入
hr("1. 载入")
fs <- sort(list.files(RDIR, pattern = "\\.rds$", full.names = TRUE))
kv("结果文件数", length(fs))
if (!length(fs)) stop("无结果文件")

rows <- list()
for (f in fs) {
  x <- readRDS(f)
  zA <- x$z_null_A; zB <- x$z_null_B
  q <- function(z, p) if (is.null(z) || !any(is.finite(z))) NA_real_ else
                      unname(quantile(abs(z[is.finite(z)]), p))
  rows[[length(rows) + 1L]] <- data.table(
    cohort = x$cohort, sig_id = x$sig_id, class = x$class,
    evaluable = x$evaluable %||% NA, low_power = x$low_power %||% NA,
    coverage = x$coverage %||% NA_real_,
    n = x$n %||% NA_integer_, n_events = x$n_events %||% NA_integer_,
    z_real = x$z_real %||% NA_real_, hr_real = x$hr_real %||% NA_real_,
    p_cox = x$p_cox %||% NA_real_,
    p_emp_A = x$p_emp_A %||% NA_real_, p_emp_B = x$p_emp_B %||% NA_real_,
    zA_med = q(zA, .50), zA_p95 = q(zA, .95), zA_max = q(zA, 1),
    zB_med = q(zB, .50), zB_p95 = q(zB, .95), zB_max = q(zB, 1),
    coexpr_real = x$coexpr_real %||% NA_real_,
    coexpr_null = x$coexpr_null_med %||% NA_real_,
    coexpr_ratio = x$coexpr_ratio %||% NA_real_,
    coexpr_qc_pass = x$coexpr_qc_pass %||% NA,
    n_relax_B = x$n_relax_B %||% NA_integer_,
    n_null = x$n_null %||% NA_integer_)
}
S <- rbindlist(rows, fill = TRUE)
kv("组合数", nrow(S))
kv("可评估", sum(S$evaluable, na.rm = TRUE))
kv("低功效队列组合", sum(S$low_power, na.rm = TRUE))

# 附上 signature 元信息
fr <- list.files(REG, pattern = "^signature_registry_locked_.*\\.csv$", full.names = TRUE)
if (length(fr)) {
  L <- fread(sort(fr, decreasing = TRUE)[1])
  keep <- intersect(c("sig_id","name","caf_definition","evidence_level"), names(L))
  S <- merge(S, L[, ..keep], by = "sig_id", all.x = TRUE)
  kv("已并入锁定表元信息", basename(sort(fr, decreasing = TRUE)[1]))
}

E <- S[evaluable == TRUE]

# ============================================================ 2. 无信号 vs 不特异
hr("2. 区分「无预后信号」与「信号不特异」")
cat("  判据（描述性，非检验）：\n")
cat("    |z_real| < 1.96                     -> 无预后信号（单变量 Cox 本身就不显著）\n")
cat("    |z_real| >= 1.96 且 p_emp_A >= 0.05  -> 有信号但不特异（随机集同样强）\n")
cat("    |z_real| >= 1.96 且 p_emp_A <  0.05  -> 有信号且强于随机集\n")

E[, sig_class := fifelse(abs(z_real) < 1.96, "无信号",
                  fifelse(p_emp_A >= 0.05, "有信号·不特异", "有信号·强于随机"))]

sh("按 signature 汇总（零分布 A）")
tb <- dcast(E, sig_id + class ~ sig_class, fun.aggregate = length, value.var = "cohort")
setnames(tb, old = names(tb), new = gsub("·", "_", names(tb)))
print(tb, row.names = FALSE)

sh("逐组合明细：真实 z 与零分布 |z| 的对比")
show <- c("cohort","sig_id","class","z_real","zA_med","zA_p95","p_emp_A","sig_class")
D <- E[order(sig_id, cohort), ..show]
D[, `:=`(z_real = round(z_real, 2), zA_med = round(zA_med, 2), zA_p95 = round(zA_p95, 2),
         p_emp_A = round(p_emp_A, 4))]
print(D, row.names = FALSE, nrows = 100)

sh("★ 关键读法")
cat("  zA_p95 是零分布 |z| 的 95% 分位。\n")
cat("  若 |z_real| 明显小于 zA_p95，说明真实 signature 弱于多数随机基因集 ——\n")
cat("  这不是『不特异』，是『本身没有信号』。两者的论文表述完全不同。\n")

sh("零分布本身的强度（各队列的 zA_p95 中位）")
print(E[, .(zA_p95_median = round(median(zA_p95, na.rm = TRUE), 2),
            zB_p95_median = round(median(zB_p95, na.rm = TRUE), 2),
            n_events = n_events[1]), by = cohort][order(-n_events)], row.names = FALSE)
cat("\n  零分布强度随事件数增加而增加，属预期（z 统计量含样本量信息）。\n")

# ============================================================ 3. QC 敏感性
hr("3. 共表达质量验证（coexpr_ratio）的敏感性分析")
cat("  该指标 = 零分布 B 模块内平均两两 Spearman ÷ 真实 signature 的同一指标。\n")
cat("  当真实 signature 内部不共表达（甚至反相关）时，分母趋零或为负，比值失效。\n")

sh("coexpr_real 的分布（分母）")
print(E[, .(coexpr_real_min = round(min(coexpr_real, na.rm=TRUE), 4),
            coexpr_real_med = round(median(coexpr_real, na.rm=TRUE), 4),
            coexpr_real_max = round(max(coexpr_real, na.rm=TRUE), 4),
            n_negative = sum(coexpr_real < 0, na.rm = TRUE),
            n_near_zero = sum(abs(coexpr_real) < 0.05, na.rm = TRUE)),
        by = sig_id][order(coexpr_real_med)], row.names = FALSE)
cat("\n  n_negative > 0 或 n_near_zero > 0 的 signature，其 coexpr_ratio 不可解释。\n")

sh("受影响的组合")
bad <- E[coexpr_real < 0 | abs(coexpr_real) < 0.05 | coexpr_ratio < 0 | coexpr_ratio > 5]
if (nrow(bad)) {
  print(bad[order(sig_id, cohort),
            .(cohort, sig_id, coexpr_real = round(coexpr_real, 4),
              coexpr_null = round(coexpr_null, 4),
              coexpr_ratio = round(coexpr_ratio, 2), coexpr_qc_pass)],
        row.names = FALSE, nrows = 60)
  cat(sprintf("\n  共 %d / %d 个组合的 coexpr_ratio 不可解释\n", nrow(bad), nrow(E)))
} else cat("  无受影响组合\n")

sh("通过队列数：主分析 vs 敏感性分析")
ok <- E[low_power == FALSE]
agg <- ok[, .(n_cohort = .N,
              pass_A = sum(p_emp_A < 0.05, na.rm = TRUE),
              pass_B_main = sum(p_emp_B < 0.05 & coexpr_qc_pass, na.rm = TRUE),
              pass_B_sens = sum(p_emp_B < 0.05, na.rm = TRUE)),
          by = .(sig_id, class)][order(sig_id)]
agg[, diff_B := pass_B_sens - pass_B_main]
print(agg, row.names = FALSE)
cat("\n  pass_B_main  按 SAP §5.2 字面执行（仅计入 coexpr_qc_pass = TRUE 的组合）\n")
cat("  pass_B_sens  不做 QC 过滤的敏感性分析\n")
cat("  diff_B       两者之差 —— 该值不为 0 的 signature，其 B 类结论依赖 QC 判定\n")

# ============================================================ 4. 分组
hr("4. 按 CAF 定义方式分组（SAP §10 路径 b 的探索性分析）")
if ("caf_definition" %in% names(E)) {
  E[, grp := fifelse(grepl("^统计筛选", caf_definition), "统计筛选",
              fifelse(grepl("解卷积", caf_definition), "解卷积衍生",
              fifelse(grepl("原代", caf_definition), "原代细胞对比",
              fifelse(grepl("综述", caf_definition), "教科书标志", "泛基质"))))]
  g <- E[low_power == FALSE, .(n_sig = uniqueN(sig_id), n_combo = .N,
          pass_A = sum(p_emp_A < 0.05, na.rm = TRUE),
          pass_rate_A = sprintf("%.0f%%", 100*mean(p_emp_A < 0.05, na.rm = TRUE)),
          median_abs_z = round(median(abs(z_real), na.rm = TRUE), 2)),
        by = grp][order(-pass_A)]
  print(g, row.names = FALSE)
  cat("\n  ★ 该分析为探索性（SAP §10 路径 b 涵盖），须在论文中明确标注。\n")
  cat("    组数少、每组 signature 数不均，不做统计检验，仅描述。\n")
}

# ============================================================ 5. 落盘
hr("5. 输出")
f1 <- file.path(CH, "data_derived", sprintf("results_extended_%s.csv", STAMP))
fwrite(S, f1); kv("扩展汇总表", f1)
f2 <- file.path(CH, "data_derived", sprintf("results_pass_counts_%s.csv", STAMP))
fwrite(agg, f2); kv("通过队列数（含敏感性）", f2)

hr("6. 供论文使用的三点")
cat("  1) 报 p_emp 时必须同时报 z_real 与 zA_p95，否则读者无法区分\n")
cat("     『无信号』与『信号不特异』。这两者的结论完全不同。\n\n")
cat("  2) coexpr_ratio 的失效须在 Limitations 中如实说明：该指标在真实\n")
cat("     signature 内部共表达接近零时不稳定，主分析与敏感性分析并列报告。\n")
cat("     须注明该问题发现于结果产出之后（SAP §14 记录）。\n\n")
cat("  3) 按 CAF 定义方式分组的通过率对比为探索性分析，不做统计检验。\n")

cat("\n输出: ", OUT, "\n", sep = "")
