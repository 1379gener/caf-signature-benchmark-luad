#!/usr/bin/env Rscript
# ==============================================================================
# lock_registry.R —— signature 清单锁定（SAP §3.6）
#
# 这是分析开工前的最后一道闸。锁定后按 SAP §11.1 第 2 条不得增删任何 signature。
# 因此锁定前必须确认清单是完整、自洽、可复算的。
#
# 检查项（BLOCK 不通过则拒绝锁定，WARN 仅提示）：
#   B1 每个 signature 有基因文件，B 类另有权重文件
#   B2 基因数与登记表一致，且无重复、无空行
#   B3 类别与基因数规则一致（A: >=10；B: 3-9 且有权重）
#   B4 B 类权重与基因列表一一对应
#   B5 无任何一对 signature 重叠率 > 90%（SAP §3.2）
#   B6 符号归一化已执行（存在 symbol_mapping_*.csv）
#   W1 数量是否落在 SAP §3.3 的目标区间
#   W2 是否存在证据等级 4（原文有文本矛盾）的条目
#
# 用法：
#   Rscript lock_registry.R                    # 只检查，不锁定
#   Rscript lock_registry.R --lock             # 全部 BLOCK 通过则锁定
#   Rscript lock_registry.R --lock --waive=B6 --reason="..."   # 豁免某项（须给理由）
# ==============================================================================

args   <- commandArgs(trailingOnly = TRUE)
opt    <- function(k, d = NA) { h <- grep(paste0("^--", k, "="), args, value = TRUE)
                                if (length(h)) sub(paste0("^--", k, "="), "", h[1]) else d }
DO_LOCK <- "--lock" %in% args
WAIVE   <- if (is.na(opt("waive"))) character(0) else strsplit(opt("waive"), ",")[[1]]
REASON  <- opt("reason")

suppressPackageStartupMessages(library(data.table))

REGDIR <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")
GENDIR <- file.path(REGDIR, "genes")
STAMP  <- format(Sys.Date(), "%Y%m%d")
OUT    <- file.path("/tmp", format(Sys.time(), "lock_check_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

RESULTS <- list()
check <- function(id, desc, ok, detail = "") {
  waived <- id %in% WAIVE
  status <- if (ok) "PASS" else if (waived) "WAIVED" else if (substr(id,1,1) == "B") "BLOCK" else "WARN"
  RESULTS[[id]] <<- list(id = id, desc = desc, status = status, detail = detail)
  cat(sprintf("  [%-6s] %-3s %s\n", status, id, desc))
  if (nzchar(detail)) for (l in strsplit(detail, "\n")[[1]]) cat("            ", l, "\n", sep = "")
  invisible(ok)
}

cat("signature 清单锁定检查\n时间: ", format(Sys.time()),
    "\n模式: ", ifelse(DO_LOCK, "检查 + 锁定", "仅检查"), "\n", sep = "")
if (length(WAIVE)) {
  kv("豁免项", WAIVE)
  if (is.na(REASON)) stop("使用 --waive 时必须同时给出 --reason=\"...\"，理由会写入锁定文件")
  kv("豁免理由", REASON)
}

# ============================================================ 载入
hr("0. 载入登记表")
fr <- file.path(REGDIR, "registry_notes.csv")
if (!file.exists(fr)) stop("未找到 registry_notes.csv")
reg <- fread(fr)
kv("登记条目", nrow(reg))
kv("列", names(reg))
fe <- file.path(REGDIR, "excluded_notes.csv")
n_exc <- if (file.exists(fe)) nrow(fread(fe)) else 0L
kv("排除条目", n_exc)

GEN <- list(); WT <- list()
for (i in seq_len(nrow(reg))) {
  sid <- reg$sig_id[i]
  fg <- file.path(GENDIR, paste0(sid, "_genes.txt"))
  if (file.exists(fg)) { g <- trimws(readLines(fg)); GEN[[sid]] <- g[nzchar(g)] }
  fw <- file.path(GENDIR, paste0(sid, "_weights.csv"))
  if (file.exists(fw)) WT[[sid]] <- fread(fw)
}

# ============================================================ 检查
hr("1. 检查")

# ---- B1 文件齐备
miss <- character(0)
for (i in seq_len(nrow(reg))) {
  sid <- reg$sig_id[i]
  if (is.null(GEN[[sid]])) miss <- c(miss, paste0(sid, "_genes.txt"))
  if (reg$class[i] == "B" && is.null(WT[[sid]])) miss <- c(miss, paste0(sid, "_weights.csv"))
}
check("B1", "每个 signature 有基因文件，B 类另有权重文件",
      length(miss) == 0, if (length(miss)) paste("缺:", paste(miss, collapse = ", ")) else "")

# ---- B2 基因数一致、无重复
bad <- character(0)
for (i in seq_len(nrow(reg))) {
  sid <- reg$sig_id[i]; g <- GEN[[sid]]; if (is.null(g)) next
  if (length(g) != reg$n_genes[i])
    bad <- c(bad, sprintf("%s 文件 %d 个 vs 登记 %d 个", sid, length(g), reg$n_genes[i]))
  if (anyDuplicated(g)) bad <- c(bad, sprintf("%s 有重复基因: %s", sid,
      paste(unique(g[duplicated(g)]), collapse = ",")))
}
check("B2", "基因数与登记表一致，无重复基因",
      length(bad) == 0, paste(bad, collapse = "\n"))

# ---- B3 类别规则
bad <- character(0)
for (i in seq_len(nrow(reg))) {
  n <- reg$n_genes[i]; cl <- reg$class[i]; sid <- reg$sig_id[i]
  if (cl == "A" && n < 10L) bad <- c(bad, sprintf("%s 归 A 类但只有 %d 个基因（应 >=10）", sid, n))
  if (cl == "B" && (n < 3L || n > 9L)) bad <- c(bad, sprintf("%s 归 B 类但有 %d 个基因（应 3-9）", sid, n))
  if (!cl %in% c("A","B")) bad <- c(bad, sprintf("%s 类别为 '%s'，只允许 A 或 B", sid, cl))
}
check("B3", "类别与基因数规则一致（SAP §3.3）",
      length(bad) == 0, paste(bad, collapse = "\n"))

# ---- B4 权重与基因一一对应
bad <- character(0)
for (sid in names(WT)) {
  w <- WT[[sid]]; g <- GEN[[sid]]
  if (!all(c("gene","weight") %in% names(w))) { bad <- c(bad, paste(sid, "权重文件缺 gene/weight 列")); next }
  if (!setequal(w$gene, g)) bad <- c(bad, sprintf("%s 权重基因与基因列表不符（仅权重侧: %s；仅列表侧: %s）",
      sid, paste(setdiff(w$gene, g), collapse=","), paste(setdiff(g, w$gene), collapse=",")))
  if (any(!is.finite(w$weight))) bad <- c(bad, paste(sid, "权重含非有限值"))
}
check("B4", "B 类权重与基因列表一一对应且为有限值",
      length(bad) == 0, paste(bad, collapse = "\n"))

# ---- B5 重叠率
ids <- names(GEN); pairs <- character(0); ovtab <- list()
for (i in seq_along(ids)) for (j in seq_along(ids)) {
  if (j <= i) next
  a <- GEN[[ids[i]]]; b <- GEN[[ids[j]]]; ov <- length(intersect(a, b))
  if (ov == 0L) next
  fa <- ov/length(a); fb <- ov/length(b)
  ovtab[[length(ovtab)+1L]] <- data.table(A = ids[i], nA = length(a), B = ids[j], nB = length(b),
                                          overlap = ov, fracA = round(fa,3), fracB = round(fb,3))
  if (max(fa, fb) > 0.90) pairs <- c(pairs, sprintf("%s vs %s: %d 个交集（%.0f%% / %.0f%%）",
      ids[i], ids[j], ov, 100*fa, 100*fb))
}
check("B5", "无任何一对重叠率 > 90%（SAP §3.2）",
      length(pairs) == 0, paste(pairs, collapse = "\n"))
if (length(ovtab)) { cat("\n  重叠明细：\n"); print(rbindlist(ovtab)[order(-overlap)], row.names = FALSE) }

# ---- B6 符号归一化
fm <- list.files(REGDIR, pattern = "^symbol_mapping_[0-9]{8}\\.csv$", full.names = TRUE)
b6_detail <- if (length(fm)) {
  paste("找到:", basename(sort(fm, decreasing = TRUE)[1]))
} else {
  "未找到。先跑 normalize_symbols.R --write"
}
check("B6", "符号归一化已执行（存在 symbol_mapping_*.csv）", length(fm) > 0, b6_detail)

# ---- W1 数量
n <- nrow(reg)
check("W1", "数量落在 SAP §3.3 目标区间 12-20",
      n >= 12L && n <= 20L, sprintf("当前 %d 个（A 类 %d，B 类 %d）", n,
      sum(reg$class=="A"), sum(reg$class=="B")))

# ---- W2 证据等级
if ("evidence_level" %in% names(reg)) {
  lv4 <- reg$sig_id[reg$evidence_level == 4L]
  check("W2", "无证据等级 4（原文存在文本矛盾）的条目",
        length(lv4) == 0,
        if (length(lv4)) paste("等级 4:", paste(lv4, collapse = ", "),
                               "\n结果解读时须打折，并在论文中如实报告") else "")
}

# ============================================================ 判定
hr("2. 判定")
st <- vapply(RESULTS, `[[`, "", "status")
n_block <- sum(st == "BLOCK"); n_warn <- sum(st == "WARN"); n_waive <- sum(st == "WAIVED")
cat(sprintf("  PASS %d   BLOCK %d   WARN %d   WAIVED %d\n",
            sum(st == "PASS"), n_block, n_warn, n_waive))

if (n_block > 0L) {
  cat("\n  存在 BLOCK 项，拒绝锁定。逐项修正后重跑。\n")
  cat("  确有正当理由需豁免的，用 --waive=B# --reason=\"...\"，理由会写入锁定文件。\n")
  quit(status = 1)
}
if (n_warn > 0L) cat("\n  存在 WARN 项，不阻止锁定，但会记入锁定文件。\n")

if (!DO_LOCK) {
  cat("\n  全部 BLOCK 项通过。未加 --lock，本次不锁定。\n")
  cat("  确认无误后重跑并加 --lock。\n")
  quit(status = 0)
}

# ============================================================ 锁定
hr("3. 锁定")
md5 <- function(p) if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
L <- copy(reg)
L[, genes_file   := paste0(sig_id, "_genes.txt")]
L[, genes_md5    := vapply(sig_id, function(s) md5(file.path(GENDIR, paste0(s, "_genes.txt"))), "")]
L[, weights_file := ifelse(class == "B", paste0(sig_id, "_weights.csv"), NA_character_)]
L[, weights_md5  := vapply(sig_id, function(s) md5(file.path(GENDIR, paste0(s, "_weights.csv"))), "")]
L[, locked_date  := STAMP]
L[, n_excluded   := n_exc]
L[, waived       := paste(WAIVE, collapse = ";")]
L[, waive_reason := ifelse(is.na(REASON), "", REASON)]
L[, warn_items   := paste(names(st)[st == "WARN"], collapse = ";")]

fL <- file.path(REGDIR, sprintf("signature_registry_locked_%s.csv", STAMP))
if (file.exists(fL)) {
  cat("  !! 锁定文件已存在: ", basename(fL), "\n", sep = "")
  cat("     SAP §11.1 第 2 条：清单锁定后不得增删。\n")
  cat("     如确需重新锁定，须先在 SAP §14 记录修订，再手工移走旧文件。\n")
  quit(status = 1)
}
fwrite(L, fL)
lock_md5 <- md5(fL)

kv("锁定文件", basename(fL))
kv("md5", lock_md5)
kv("signature 数", nrow(L))
cat("\n  各文件 md5：\n")
for (i in seq_len(nrow(L)))
  cat(sprintf("    %-10s genes %s%s\n", L$sig_id[i], substr(L$genes_md5[i], 1, 12),
      ifelse(is.na(L$weights_md5[i]) || !nzchar(L$weights_md5[i]), "",
             paste0("   weights ", substr(L$weights_md5[i], 1, 12)))))

hr("4. 请立即完成的两件事")
cat("  1) 把下面这行写入 SAP §14「signature 清单锁定文件 md5」：\n\n")
cat("       ", lock_md5, "\n\n", sep = "")
cat("  2) 把锁定日期 ", STAMP, " 一并记入 §14 修订记录。\n\n", sep = "")
cat("  此后：\n")
cat("    - 不得增删任何 signature（SAP §11.1 第 2 条）\n")
cat("    - run_benchmark.R 的锁定闸将读取该文件并校验\n")
cat("    - 如发现某 signature 基因列表有误，按 §11.2 处理：移出主分析、记录、全部结果重跑\n")

cat("\n输出: ", OUT, "\n", sep = "")
