#!/usr/bin/env Rscript
# ==============================================================================
# 03_baseA_tcga.R —— 底座 A · TCGA-LUAD 整备
#
# 依据 02_verify_tcga_se.R 的核验结果编写。五条口径写在下面 CONFIG 中，
# 每一条都是显式的、可复算的，不接受"默认行为"。
#
# 内置自检：Stage II-IV 去重后唯一病人数必须 == 166（Paper 1 锁定值）。
# 该断言失败即停止 —— 说明去重规则与重建时不一致，须查明后再继续。
#
# 输入（只读）：~/caf_spatial_project/data/TCGA_LUAD/TCGA_LUAD_SE.rds
# 输出：        ~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/TCGA_LUAD/
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

# ============================== CONFIG（口径，改这里，不改下面） ==============
CFG <- list(
  se_path      = "~/caf_spatial_project/data/TCGA_LUAD/TCGA_LUAD_SE.rds",
  out_dir      = "~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/TCGA_LUAD",

  assay_use    = "tpm_unstrand",   # 口径 1：与 Paper 1 的 TPM 一致
  dedup_rule   = "max_libsize",    # 口径 2：同病人多 aliquot 取文库最大者
                                   #         备选 "first_barcode"（字典序首位）
  gene_id      = "symbol",         # 口径 3：Ensembl -> HGNC symbol
  symbol_dedup = "max_mean",       # 同 symbol 多探针/多 ID 取均表达最大者
  biotype_keep = "protein_coding", # 口径 4：NA 表示不过滤
  sample_set   = "all_stages",     # 口径 5：全分期（"stage_II_IV" 为 Paper 1 口径）

  # 自检靶子：Stage II-IV 去重后应为 166（Paper 1 锁定值）
  assert_stage24_n = 166L,

  time_unit    = "days"            # 全底座统一为天
)
# =============================================================================

CFG$se_path <- path.expand(CFG$se_path)
CFG$out_dir <- path.expand(CFG$out_dir)
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(CFG$out_dir, format(Sys.time(), "build_%Y%m%d_%H%M.log"))
sink(LOG, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("-", 66), "\n", t, "\n", strrep("-", 66), "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-32s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("底座 A · TCGA-LUAD 整备\n时间: ", format(Sys.time()), "\n", sep = "")
hr("口径（本次运行使用的全部约定）")
for (n in names(CFG)) kv(n, as.character(CFG[[n]]))

# ------------------------------------------------------------------ 1. 载入
hr("1. 载入")
se <- readRDS(CFG$se_path)
kv("维度", paste(dim(se), collapse = " x "))
stopifnot(CFG$assay_use %in% assayNames(se))

bc    <- colnames(se)
patid <- substr(bc, 1, 12)
stype <- vapply(strsplit(bc, "-"), function(x) substr(x[4], 1, 2), "")
stopifnot(all(stype == "01"))          # 核验已确认：全部为原发肿瘤
kv("列数 / 唯一病人数", paste(length(bc), length(unique(patid)), sep = " / "))

# ------------------------------------------------------ 2. 同病人多 aliquot 去重
hr("2. 去重（口径 2）")
libsize <- colSums(assay(se, "unstranded"))
dup_pat <- names(which(table(patid) > 1))
kv("有多列的病人数", length(dup_pat))

cat("\n  多列病人的完整 barcode 与文库大小（★ 标记为被保留者）：\n")
for (p in dup_pat) {
  idx <- which(patid == p)
  idx <- idx[order(libsize[idx], decreasing = TRUE)]
  keep_i <- if (CFG$dedup_rule == "max_libsize") idx[1] else which(patid == p)[order(bc[which(patid == p)])][1]
  for (i in idx) {
    cat(sprintf("    %s %-30s libsize=%12.0f\n",
                if (i == keep_i) "★" else " ", bc[i], libsize[i]))
  }
}

keep <- vapply(split(seq_along(patid), patid), function(idx) {
  if (length(idx) == 1L) return(idx)
  if (CFG$dedup_rule == "max_libsize") idx[which.max(libsize[idx])] else idx[order(bc[idx])][1]
}, integer(1))
keep <- sort(unname(keep))

se2 <- se[, keep]
pat2 <- substr(colnames(se2), 1, 12)
kv("去重后列数", ncol(se2))
kv("去重后唯一病人数", length(unique(pat2)))
stopifnot(!anyDuplicated(pat2))
cat("  去重后无重复病人 —— 通过\n")

# ---------------------------------------------------------- 3. 基因 ID 与类型
hr("3. 基因映射（口径 3、4）")
rd <- as.data.frame(rowData(se2))
kv("rowData 可用列", colnames(rd))
stopifnot(all(c("gene_name", "gene_type") %in% colnames(rd)))

cat("\n  gene_type 分布（前 8）：\n")
tt <- sort(table(rd$gene_type), decreasing = TRUE)
for (i in seq_len(min(8, length(tt)))) cat(sprintf("    %-26s %6d\n", names(tt)[i], tt[i]))

sel <- rep(TRUE, nrow(se2))
if (!is.na(CFG$biotype_keep)) {
  sel <- rd$gene_type == CFG$biotype_keep
  kv(paste0("保留 ", CFG$biotype_keep), sum(sel))
}
sel <- sel & !is.na(rd$gene_name) & nzchar(rd$gene_name)

mat <- assay(se2, CFG$assay_use)[sel, , drop = FALSE]
sym <- rd$gene_name[sel]
kv("过滤后行数", nrow(mat))
kv("唯一 symbol 数", length(unique(sym)))
kv("需合并的重复 symbol", sum(duplicated(sym)))

if (anyDuplicated(sym)) {
  rmean <- rowMeans(mat)
  ord   <- order(sym, -rmean)                 # 同 symbol 内按均表达降序
  mat   <- mat[ord, , drop = FALSE]; sym <- sym[ord]
  keepr <- !duplicated(sym)                   # 取每个 symbol 的第一行 = 均表达最大
  mat   <- mat[keepr, , drop = FALSE]; sym <- sym[keepr]
}
rownames(mat) <- sym
mat <- mat[order(rownames(mat)), , drop = FALSE]
colnames(mat) <- pat2                          # 列名统一为病人 ID
kv("最终表达矩阵", paste(dim(mat), collapse = " x "))
kv("取值范围（抽样）", sprintf("[%.3g, %.3g]", min(mat[1:min(2000,nrow(mat)),1:min(20,ncol(mat))]),
                                max(mat[1:min(2000,nrow(mat)),1:min(20,ncol(mat))])))

# ------------------------------------------------------------------ 4. 临床表
hr("4. 临床与生存（统一为天）")
cd <- as.data.frame(colData(se2))
need <- c("vital_status", "days_to_death", "days_to_last_follow_up", "ajcc_pathologic_stage")
miss <- setdiff(need, colnames(cd))
if (length(miss)) { cat("  !! colData 缺少字段: ", paste(miss, collapse=", "), "\n"); quit(status=1) }

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
cl <- data.frame(
  patient_id = pat2,
  OS_event   = as.integer(grepl("^dead$", cd$vital_status, ignore.case = TRUE)),
  OS_time    = ifelse(grepl("^dead$", cd$vital_status, ignore.case = TRUE),
                      num(cd$days_to_death), num(cd$days_to_last_follow_up)),
  stage_raw  = as.character(cd$ajcc_pathologic_stage),
  age        = if ("age_at_index" %in% colnames(cd)) num(cd$age_at_index) else NA_real_,
  sex        = if ("gender" %in% colnames(cd)) as.character(cd$gender) else NA_character_,
  stringsAsFactors = FALSE
)

# 分期归一：Stage IA/IB -> I，IIA/IIB -> II，IIIA/IIIB -> III
cl$stage <- sub("^Stage\\s+", "", cl$stage_raw)
cl$stage <- sub("[ABC]$", "", cl$stage)
cl$stage[!cl$stage %in% c("I","II","III","IV")] <- NA
cat("\n  分期归一后分布：\n")
tb <- table(cl$stage, useNA = "ifany")
for (i in seq_along(tb)) cat(sprintf("    %-6s %5d\n", ifelse(is.na(names(tb)[i]),"NA",names(tb)[i]), tb[i]))

kv("OS_time 缺失", sum(is.na(cl$OS_time)))
kv("OS_time <= 0", sum(!is.na(cl$OS_time) & cl$OS_time <= 0))
kv("死亡事件数（全分期）", sum(cl$OS_event == 1, na.rm = TRUE))

valid <- !is.na(cl$OS_time) & cl$OS_time > 0
kv("有可用生存信息者", sum(valid))

# --------------------------------------------------- 5. ★ 自检：复现 n=166
hr("5. 自检 —— Stage II-IV 去重后应为 166")
s24 <- cl$stage %in% c("II","III","IV")
kv("Stage II-IV 病人数（去重后）", sum(s24))
kv("其中有可用生存信息者", sum(s24 & valid))
kv("其中死亡事件数", sum(cl$OS_event[s24 & valid] == 1))

got <- sum(s24 & valid)
if (identical(as.integer(got), CFG$assert_stage24_n)) {
  cat("\n  ✅ 通过：复现 n=", got, "，与 Paper 1 锁定值一致。\n", sep = "")
  cat("     → 去重规则正确，本文件可作为底座 A 的 TCGA 来源。\n")
} else {
  cat("\n  ⚠️  未复现：得到 ", got, "，锁定值为 ", CFG$assert_stage24_n, "。\n", sep = "")
  cat("     不要改断言去迁就结果。可能原因（按概率）：\n")
  cat("       1) Paper 1 的生存信息来自 clinical.csv 而非 colData，两者取值口径不同\n")
  cat("       2) Paper 1 的分期筛选用的是原始 ajcc 字符串而非归一后的 I/II/III/IV\n")
  cat("       3) 去重规则不同（Paper 1 可能取 barcode 首位而非文库最大）\n")
  cat("     先查明差异来源，再决定改哪一边。本脚本继续执行，但结果标记为待确认。\n")
}

# ------------------------------------------------------------------ 6. 落盘
hr("6. 输出")
if (CFG$sample_set == "stage_II_IV") {
  sel_s <- s24 & valid
} else {
  sel_s <- valid
}
cl_out  <- cl[sel_s, ]
mat_out <- mat[, sel_s, drop = FALSE]
stopifnot(identical(colnames(mat_out), cl_out$patient_id))
kv("最终分析集", sprintf("%d 病人 x %d 基因，事件 %d",
                        nrow(cl_out), nrow(mat_out), sum(cl_out$OS_event)))

f_expr <- file.path(CFG$out_dir, "TCGA_LUAD_expr_symbol_tpm.rds")
f_clin <- file.path(CFG$out_dir, "TCGA_LUAD_clinical_clean.csv")
saveRDS(mat_out, f_expr)
write.csv(cl_out, f_clin, row.names = FALSE)

man <- data.frame(
  cohort = "TCGA-LUAD", platform = "Illumina RNA-seq (STAR counts, GDC)",
  n_raw_cols = ncol(se), n_dedup = ncol(se2), n_final = nrow(cl_out),
  n_events = sum(cl_out$OS_event), n_genes = nrow(mat_out),
  assay = CFG$assay_use, dedup_rule = CFG$dedup_rule,
  biotype = ifelse(is.na(CFG$biotype_keep), "all", CFG$biotype_keep),
  sample_set = CFG$sample_set, time_unit = CFG$time_unit,
  median_followup_days = round(median(cl_out$OS_time, na.rm = TRUE), 1),
  source_md5 = tools::md5sum(CFG$se_path),
  built = format(Sys.time(), "%Y-%m-%d %H:%M"),
  stringsAsFactors = FALSE
)
f_man <- file.path(CFG$out_dir, "cohort_manifest_TCGA.csv")
write.csv(man, f_man, row.names = FALSE)

for (f in c(f_expr, f_clin, f_man)) kv(basename(f), sprintf("%.1f MB", file.info(f)$size/1048576))

hr("sessionInfo")
print(sessionInfo())
cat("\n日志: ", LOG, "\n", sep = "")
