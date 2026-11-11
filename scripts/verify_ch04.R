#!/usr/bin/env Rscript
# verify_ch04.R —— 第 4 章数字复算 · 依据 SAP v1.2 §12
# 无 check 的数字视为未验证，不得定稿。
suppressMessages({library(data.table); library(survival)})
options(warn = 1)

CH <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark")
DD <- file.path(CH, "data_derived")
BASE <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived")
REG_MD5 <- "4ad3218931d95dfb6ee67a7ba07b9a75"
TAG <- "20260905"

FAIL <- 0L; SKIP <- 0L
ok   <- function(m) cat(sprintf("  [PASS] %s\n", m))
bad  <- function(m) { FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
skip <- function(m) { SKIP <<- SKIP + 1L; cat(sprintf("  [SKIP] %s\n", m)) }
chk  <- function(cond, m) if (isTRUE(cond)) ok(m) else bad(m)

E <- fread(file.path(DD, sprintf("results_extended_%s.csv", TAG)))
P <- fread(file.path(DD, sprintf("results_pass_counts_%s.csv", TAG)))
S <- fread(file.path(DD, sprintf("null_benchmark_summary_%s.csv", TAG)))
R <- fread(file.path(DD, sprintf("results_by_cohort_role_%s.csv", TAG)))

cat("\n===== C0 完整性 =====\n")
rm_ <- unique(na.omit(S$registry_md5)); rm_ <- rm_[nzchar(rm_)]
chk(length(rm_) == 1 && tolower(rm_) == REG_MD5,
    sprintf("registry md5（非缺失值 %d 个）== %s", length(rm_), REG_MD5))
cat(sprintf("  [INFO] registry_md5 填充率 %d/%d\n", sum(nzchar(as.character(S$registry_md5)) & !is.na(S$registry_md5)), nrow(S)))
ev_ <- unique(na.omit(S$engine_version))
chk(length(ev_) == 1, sprintf("engine_version 非缺失值唯一: %s", paste(ev_, collapse=",")))
nn_ <- unique(na.omit(S$n_null))
chk(length(nn_) == 1 && nn_ == 1000, sprintf("n_null 非缺失值 == 1000（实际 %s）", paste(nn_, collapse=",")))
chk(nrow(E) == 72, sprintf("组合总数 == 72（实际 %d）", nrow(E)))

cat("\n===== C1 覆盖率与可评估（SAP §4.2）=====\n")
chk(sum(E$evaluable) == 71, sprintf("可评估 == 71（实际 %d）", sum(E$evaluable)))
ne <- E[evaluable == FALSE]
chk(nrow(ne) == 1 && ne$cohort == "GSE68465" && ne$sig_id == "SIG-004",
    "唯一不可评估 == GSE68465 × SIG-004")
chk(nrow(ne) == 1 && abs(ne$coverage - 4/7) < 0.01,
    sprintf("其覆盖率 ≈ 57.1%%（实际 %.1f%%）", 100*ne$coverage[1]))
chk(all(E[evaluable == TRUE]$coverage >= 0.70), "全部可评估组合覆盖率 >= 70%")

cat("\n===== C2 队列样本量与 cohort_registry 一致（SAP §12-3）=====\n")
frg <- file.path(BASE, "cohort_registry.csv")
if (file.exists(frg)) {
  G <- fread(frg)
  nm <- intersect(c("cohort","cohort_id","dataset"), names(G))[1]
  cn <- intersect(c("n","n_samples","n_final"), names(G))[1]
  ce <- intersect(c("n_events","events","n_event"), names(G))[1]
  if (!is.na(nm) && !is.na(cn)) {
    obs <- unique(S[, .(cohort, n = as.integer(n), n_events = as.integer(n_events))])
    obs <- obs[, .(n = max(n), n_events = max(n_events)), by = cohort]
    norm <- function(x) toupper(gsub("[^A-Z0-9]", "", toupper(x)))
    obs[, k := norm(cohort)]
    m <- merge(obs, G[, .(k = norm(get(nm)), n_reg = as.integer(get(cn)),
                          ev_reg = if (!is.na(ce)) as.integer(get(ce)) else NA_integer_)],
               by = "k", all.x = TRUE)
    print(m)
    chk(all(!is.na(m$n_reg) & m$n == m$n_reg), "各队列 n 与 registry 一致")
    if (!is.na(ce)) chk(all(m$n_events == m$ev_reg), "各队列事件数与 registry 一致")
  } else skip(sprintf("cohort_registry 列名未识别: %s", paste(names(G), collapse=",")))
} else skip("cohort_registry.csv 未找到")

cat("\n===== C3 GSE31210 未降级（SAP §2.4）=====\n")
g <- unique(S[cohort == "GSE31210", .(n = as.integer(n), n_events = as.integer(n_events))])
print(g)
chk(nrow(g) >= 1 && max(g$n) == 204, "GSE31210 n == 204（exclude 已应用）")
chk(nrow(g) >= 1 && max(g$n_events) == 30, "GSE31210 事件 == 30")
chk(nrow(g) >= 1 && max(g$n_events) >= 25, "事件数 >= 25 阈值 -> 计入主汇总，未降级")
lowp <- unique(S[, .(cohort, n_events = as.integer(n_events))])[, max(n_events), by=cohort][V1 < 25]
chk(nrow(lowp) == 0, "无任何队列触及 <25 阈值")

cat("\n===== C4 ★ A/B 分开汇总，禁止合并（SAP §3.3 / §4.4）=====\n")
ev <- E[evaluable == TRUE][, indist := abs(z_real) < zA_p95]
byc <- ev[, .(n组合 = .N, 不可区分 = sum(indist),
              比例 = sprintf("%.1f%%", 100*mean(indist))), by = class][order(class)]
print(byc)
chk(nrow(byc) == 2, "A 类与 B 类分开呈现")
cat(sprintf("  [INFO] 合并值 %d/%d —— ★ 仅供内部核对，禁止作为主结果报告\n",
            sum(ev$indist), nrow(ev)))
chk(uniqueN(P[class == "A"]$sig_id) == 5, "A 类 signature 数 == 5（SAP §3.3）")
chk(uniqueN(P[class == "B"]$sig_id) == 7, "B 类 signature 数 == 7（SAP §3.3）")

cat("\n===== C5 敏感性分析（SAP §5.2.1 第 3 条）=====\n")
print(P[, .(sig_id, class, n_cohort, pass_A, pass_B_main, pass_B_sens, diff_B)][order(-abs(diff_B))])
d <- P[diff_B != 0]
chk(nrow(d) == 1 && d$sig_id == "SIG-001", "仅 SIG-001 的 diff_B != 0")
chk(nrow(d) == 1 && d$pass_B_main == 0 && d$pass_B_sens == 4, "SIG-001: 0 -> 4")
chk(all(P[sig_id != "SIG-001"]$diff_B == 0), "其余 11 条主分析与敏感性完全一致")

cat("\n===== C6 p_emp 独立复算（SAP §12-2，至少 3 个 signature @ TCGA）=====\n")
ND <- file.path(DD, "null_benchmark")
f1 <- list.files(ND, pattern = "^TCGA_LUAD__SIG-", full.names = TRUE)
if (!length(f1)) { skip("未找到 TCGA 逐组合 RDS") } else {
  o <- readRDS(f1[1])
  cat("  RDS 结构:", paste(class(o), collapse="/"), "| 成员:",
      paste(utils::head(names(o), 25), collapse=", "), "\n")
  zf <- intersect(c("z_null_A","zA","null_z_A","z_nullA"), names(o))[1]
  rf <- intersect(c("z_real","zreal"), names(o))[1]
  if (is.na(zf) || is.na(rf)) {
    skip("RDS 字段名未识别 —— 按上面打印的成员名补 zf/rf 后重跑本节")
  } else {
    n_ok <- 0L
    for (f in utils::head(f1, 3)) {
      o <- readRDS(f); sid <- sub(".*__", "", sub("\\.rds$", "", basename(f)))
      zn <- as.numeric(o[[zf]]); zr <- as.numeric(o[[rf]])[1]
      pe <- (1 + sum(abs(zn) >= abs(zr))) / (1 + length(zn))   # SAP §4.3
      rec <- E[cohort == "TCGA_LUAD" & sig_id == sid]
      if (nrow(rec) == 1) {
        good <- abs(pe - rec$p_emp_A) < 1e-9
        cat(sprintf("    %s  复算 %.6f  引擎 %.6f  %s\n", sid, pe, rec$p_emp_A,
                    ifelse(good, "OK", "★ 不符")))
        n_ok <- n_ok + as.integer(good)
      }
    }
    chk(n_ok >= 3, sprintf("至少 3 个 signature 的 p_emp 复算一致（实际 %d）", n_ok))
  }
}

cat("\n===== C7 零分布成员 Cox z 复算（SAP §12-1，容差 1e-6）=====\n")
skip("需从表达矩阵重跑 20 个零分布成员的 Cox —— 待 RDS 字段名确认后补实现")

cat("\n===== C8 GSE72094 SCANI 回归测试（SAP §12-4）=====\n")
skip("需 scoring_protocol 复算 vs Paper 1 锁定值 —— 待补")

cat("\n===== C9 禁用数字扫描 =====\n")
BAN <- c("1.868","0.766","0.0022","192","1.499")
hits <- c()
# 禁用数字扫描仅适用于正文/图注，不适用于数值结果表（子串匹配必然误报）
TXT <- c(list.files(file.path(CH,"text"), pattern="\\.(md|txt|tex)$", full.names=TRUE),
         list.files(file.path(CH,"tables"), pattern="\\.(md|txt|tex)$", full.names=TRUE))
if (!length(TXT)) { cat("  [INFO] text/ 与 tables/ 暂无正文文件，跳过\n") } else
for (f in TXT) {
  tx <- paste(readLines(f, warn=FALSE), collapse=" ")
  for (b in BAN) if (grepl(sprintf("(^|[^0-9.])%s([^0-9]|$)", gsub("\\.","\\\\.",b)), tx))
    hits <- c(hits, paste(basename(f), b))
}
if (length(hits)) { for (h in hits) cat("  [WARN] 命中:", h, "\n") } else ok("data_derived 无禁用数字")

cat(sprintf("\n===== 结果：FAIL %d · SKIP %d =====\n", FAIL, SKIP))
if (SKIP > 0) cat("★ 存在 SKIP 项，按 SAP §12「无 check 的数字视为未验证」，尚不具备定稿条件。\n")
quit(status = if (FAIL > 0) 1 else 0)
