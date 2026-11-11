#!/usr/bin/env Rscript
# ==============================================================================
# fix_gse31210_exclude.R —— 修复 SAP §2.4 的实际偏离
#
# 问题：SAP v1.1 §2.4 规定遵从 GSE31210 原作者的 exclude 标志、排除 22 例，
#   但 05_baseA_geo_rev3.R 只把该字段提取存进 extra_ 列，未用于过滤。
#   cl$valid 的定义里没有这一项，导致该队列仍是 226 例 / 35 事件。
#
# 本脚本对 05_baseA_geo_rev3.R 做两处修改：
#   [1] CFG 中为 GSE31210 增加 exclude_key / exclude_val
#   [2] cl$valid 的定义中加入排除条件（对未设该字段的队列无影响）
#   [3] 顺带增加 --cohort= 过滤，避免重建时波及其余四个队列
#
# 修改前备份。不自动重建，只打印后续命令供人工确认后执行。
# ==============================================================================

SP <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/scripts/05_baseA_geo_rev3.R")
STAMP <- format(Sys.Date(), "%Y%m%d")
OUT <- file.path("/tmp", format(Sys.time(), "fix31210_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")

cat("修复 GSE31210 的 exclude 标志（SAP §2.4）\n时间: ", format(Sys.time()), "\n", sep = "")

# ============================================================ 0. 现状核对
hr("0. 现状核对")
DER <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/GSE31210")
fc <- file.path(DER, "GSE31210_clinical_clean.csv")
if (file.exists(fc)) {
  cl0 <- read.csv(fc, stringsAsFactors = FALSE, check.names = FALSE)
  cat("  当前 GSE31210 分析集: ", nrow(cl0), " 例 / ",
      sum(cl0$OS_event, na.rm = TRUE), " 事件\n", sep = "")
  ec <- grep("^extra_", colnames(cl0), value = TRUE)
  cat("  extra 列: ", ifelse(length(ec), paste(ec, collapse = ", "), "（无）"), "\n", sep = "")
  for (e in ec) {
    tb <- table(cl0[[e]], useNA = "ifany")
    cat("    ", e, " 取值: ", paste(sprintf("%s=%d",
        ifelse(is.na(names(tb)), "NA", names(tb)), tb), collapse = "  "), "\n", sep = "")
  }
  cat("\n  -> 上表中标记为 'exclude' 的例数即 SAP §2.4 应排除但未排除的部分\n")
} else cat("  !! 未找到 ", fc, "\n", sep = "")

# ============================================================ 1. 备份
hr("1. 备份")
if (!file.exists(SP)) stop("未找到 ", SP)
bak <- paste0(SP, ".bak_", STAMP)
if (!file.exists(bak)) file.copy(SP, bak)
cat("  ", basename(SP), " -> ", basename(bak), "\n", sep = "")
src <- readLines(SP, warn = FALSE)

if (any(grepl("exclude_key", src))) {
  cat("\n  已包含 exclude 过滤逻辑，跳过修改。\n")
} else {

# ============================================================ 2. 修改 CFG
hr("2. 修改 GSE31210 的 CFG")
i <- grep('extra_keys = "exclude for prognosis analysis', src)
if (!length(i)) {
  i <- grep('^\\s*extra_keys\\s*=\\s*"exclude', src)
}
if (!length(i)) stop("未在 CFG 中定位 GSE31210 的 extra_keys 行")
cat("  定位到第 ", i[1], " 行:\n    ", trimws(src[i[1]]), "\n", sep = "")

EXKEY <- "exclude for prognosis analysis due to incomplete resection or adjuvant therapy"
add <- c(
  sprintf('    exclude_key = "%s",', EXKEY),
  '    exclude_val = "exclude",   # SAP §2.4：遵从原作者的预后分析排除标志')
src <- append(src, add, after = i[1] - 1L)
cat("  已在其前插入 exclude_key / exclude_val\n")

# ============================================================ 3. 修改 valid 定义
hr("3. 修改 cl$valid 的定义")
j <- grep('valid\\s*:=\\s*keep', src)
if (!length(j)) j <- grep('cl\\$valid\\s*<-\\s*keep', src)
if (!length(j)) stop("未定位到 cl$valid 的定义行")
cat("  定位到第 ", j[1], " 行:\n    ", trimws(src[j[1]]), "\n", sep = "")

NEWBLK <- c(
'  # ---- SAP §2.4：遵从数据集自带的 exclude 标志（未设该字段的队列不受影响）',
'  excl <- rep(FALSE, nrow(cl))',
'  if (!is.null(cfg$exclude_key) && !is.na(cfg$exclude_key)) {',
'    ev <- na_str(get_field(pc, cfg$exclude_key))',
'    excl <- !is.na(ev) & ev == cfg$exclude_val',
'    cat(sprintf("    ★ SAP 2.4 exclude 标志: %d 例标记为 \'%s\'，予以排除\\n",',
'                sum(excl), cfg$exclude_val))',
'    cat(sprintf("       该字段非缺失 %d 例，缺失 %d 例（缺失视为不排除）\\n",',
'                sum(!is.na(ev)), sum(is.na(ev))))',
'  }',
'  cl[, excluded_by_flag := excl]',
'  cl[, valid := keep & !excl & !is.na(OS_time) & OS_time > 0 & !is.na(status)]')

src <- c(src[seq_len(j[1] - 1L)], NEWBLK, src[(j[1] + 1L):length(src)])
cat("  已替换为含排除条件的版本\n")

# ============================================================ 4. 增加队列过滤
hr("4. 增加 --cohort= 过滤")
k <- grep('^for \\(acc in names\\(CFG\\)\\) \\{', src)
if (!length(k)) stop("未定位到主循环")
FLT <- c(
'# ---- 可选：只重建指定队列（避免波及其余队列）',
'.args <- commandArgs(trailingOnly = TRUE)',
'.h <- grep("^--cohort=", .args, value = TRUE)',
'ONLY <- if (length(.h)) strsplit(sub("^--cohort=", "", .h[1]), ",")[[1]] else NULL',
'if (!is.null(ONLY)) cat("\\n[限定队列] ", paste(ONLY, collapse = ", "), "\\n", sep = "")',
'')
src <- append(src, FLT, after = k[1] - 1L)
k <- grep('^for \\(acc in names\\(CFG\\)\\) \\{', src)
src <- append(src, '  if (!is.null(ONLY) && !acc %in% ONLY) next', after = k[1])
cat("  已加入 --cohort= 支持\n")

# ============================================================ 5. manifest 记录
hr("5. manifest 增加排除计数")
m <- grep('n_after_histology = sum\\(cl\\$keep\\)', src)
if (length(m)) {
  src[m[1]] <- sub('n_after_histology = sum\\(cl\\$keep\\),',
    'n_after_histology = sum(cl$keep), n_excluded_by_flag = sum(cl$excluded_by_flag),',
    src[m[1]])
  cat("  已在 manifest 中加入 n_excluded_by_flag\n")
} else cat("  未定位到 manifest 的 n_after_histology 字段，跳过（不影响主逻辑）\n")

writeLines(src, SP)
cat("\n  修改完成，共 ", length(src), " 行\n", sep = "")
}

# ============================================================ 6. 校验
hr("6. 语法校验")
p <- tryCatch({ parse(SP); TRUE }, error = function(e) { cat("  !! 解析失败: ",
     conditionMessage(e), "\n", sep = ""); FALSE })
cat("  parse(): ", ifelse(p, "通过", "失败 —— 请从备份恢复"), "\n", sep = "")
if (!p) { cat("  恢复命令: cp ", bak, " ", SP, "\n", sep = ""); quit(status = 1) }

cat("\n  改动确认：\n")
for (pat in c("exclude_key", "exclude_val", "excluded_by_flag", "ONLY")) {
  n <- sum(grepl(pat, readLines(SP, warn = FALSE), fixed = TRUE))
  cat(sprintf("    %-20s 出现 %d 次\n", pat, n))
}

# ============================================================ 7. 后续命令
hr("7. 请依次执行")
cat("  # 1) 只重建 GSE31210（其余四个 GEO 队列不动）\n")
cat("  Rscript ", SP, " --cohort=GSE31210\n\n", sep = "")
cat("  #    看输出中这两行：\n")
cat("  #      ★ SAP 2.4 exclude 标志: N 例标记为 'exclude'，予以排除\n")
cat("  #      可用于生存分析 / 其中事件数\n")
cat("  #    预期 226 -> 约 204 例，35 -> 约 30 事件\n")
cat("  #    ★ 若事件数 < 25，按 SAP §2.4 该队列降为敏感性分析队列\n\n")
cat("  # 2) 重新统一 manifest（新增字段 + 样本量变化）\n")
cat("  Rscript ~/thesis_LUAD_CAF/01_bases/A_bulk_OS/scripts/08_harmonize_manifests.R\n\n")
cat("  # 3) 删除该队列的 12 个旧结果\n")
cat("  rm ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/data_derived/null_benchmark/GSE31210__*.rds\n\n")
cat("  # 4) 断点续跑（只会重算 GSE31210 的 12 个组合）\n")
cat("  nohup Rscript ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/run_benchmark.R \\\n")
cat("        --workers=8 > /tmp/bm_fix31210.log 2>&1 &\n\n")
cat("  # 5) 在 SAP §14 追加修订记录（见下）\n")

hr("8. SAP §14 修订记录（请抄入）")
cat("  # 3 | 2026-09-04 | 修正实现与 §2.4 的偏离：GSE31210 的 exclude 标志\n")
cat("  #                 此前仅被提取存档，未用于过滤；现已应用。\n")
cat("  #                 理由：实现未遵循条款，非条款变更。\n")
cat("  #                 结果产出前/后：★ 之后（首轮结果已产出，该队列 12 个组合作废重算）\n")
cat("\n  注：这是【实现修正】不是【条款修订】。SAP 条款本身未变，\n")
cat("      变的是实现终于符合了条款。仍须记录，因为它导致部分结果作废重算。\n")

cat("\n输出: ", OUT, "\n", sep = "")
