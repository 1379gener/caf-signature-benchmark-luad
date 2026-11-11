#!/usr/bin/env Rscript
# ==============================================================================
# 08_harmonize_manifests.R —— 统一七个队列的 manifest 字段模式
#
# 问题：TCGA 的 manifest 由最早写的 03_baseA_tcga.R 产出，缺 expr_scale /
#       need_log / time_src_unit / n_with_stage 等字段。协议的 load_cohort()
#       因此对 TCGA 走【取值域推断】而非读声明 —— 这次推断对了，但协议的
#       意义就是不靠推断。CPTAC 用 n_samples、其余用 n_final，字段名也不一致。
#
# 做法：能从数据直接算的字段一律现算（不引用旧 manifest 的表面值），
#       只有出处类信息（平台、时间字段名、注释来源、md5）从旧 manifest 继承。
#       旧 manifest 备份为 .bak_<日期> 后覆盖。
#
# 只写 derived/ 之下。不触碰表达矩阵与临床表。
# ==============================================================================

suppressPackageStartupMessages({ library(data.table) })

DER <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived")
STAMP <- format(Sys.Date(), "%Y%m%d")
OUT <- file.path("/tmp", format(Sys.time(), "harmonize_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")
kv <- function(k, v) cat(sprintf("  %-28s %s\n", paste0(k,":"), paste(v, collapse=", ")))

cat("manifest 字段模式统一\n时间: ", format(Sys.time()), "\n输出: ", OUT, "\n", sep="")

# ---- 各队列的出处类补充信息（旧 manifest 里没有、又无法从数据算出来的）
SUPP <- list(
  TCGA_LUAD = list(
    cohort = "TCGA-LUAD", platform = "Illumina RNA-seq (STAR counts, GDC)",
    n_raw = 540L, n_after_tumor = 517L, n_after_histology = 517L,
    histology_filter = "无（全部为 LUAD）",
    biotype_filter = "protein_coding",
    assay_used = "tpm_unstrand",
    sample_dedup = "同病人多 aliquot 取文库大小最大者",
    gene_collapse = "同名 symbol 取均表达最大",
    anno = "rowData(se)$gene_name（Ensembl -> HGNC）",
    time_key = "days_to_death / days_to_last_follow_up",
    time_src_unit = "天", time_mult = 1, stage_type = "ajcc",
    source_file = "~/caf_spatial_project/data/TCGA_LUAD/TCGA_LUAD_SE.rds"),
  CPTAC_LUAD = list(
    cohort = "CPTAC-LUAD", platform = "RNA-seq (LinkedOmics .cct)",
    n_raw = 110L, n_after_tumor = 110L, n_after_histology = 110L,
    histology_filter = "无（全部为 Tumor）",
    biotype_filter = "无", assay_used = "cptac_rna.cct",
    sample_dedup = "无重复", gene_collapse = "同名取均表达最大",
    anno = "源表已为 gene x sample",
    time_key = NA_character_, time_src_unit = NA_character_,
    time_mult = NA_real_, stage_type = "arabic",
    source_file = "~/caf_spatial_project/data/CPTAC_LUAD/cptac_rna.cct")
)

# ---- 规范字段顺序
SCHEMA <- c(
  "cohort","platform",
  "n_raw","n_after_tumor","n_after_histology","n_final","n_events","event_rate",
  "n_genes",
  "expr_scale","expr_median","expr_max","need_log",
  "has_survival","time_key","time_src_unit","time_mult","time_unit",
  "median_followup_days",
  "stage_type","n_with_stage",
  "histology_filter","biotype_filter","assay_used","sample_dedup",
  "gene_collapse","anno","source_md5","built","harmonized")

COHORTS <- c("TCGA_LUAD","GSE68465","GSE72094","GSE31210",
             "GSE50081","GSE30219","CPTAC_LUAD")

# R 4.3.1 尚无 base 版 %||%（4.4 才加入），自行定义
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

pick <- function(m, nm, default = NA) {
  if (nrow(m) && nm %in% colnames(m) && !is.na(m[[nm]][1])) m[[nm]][1] else default
}

rows <- list()
for (acc in COHORTS) {
  hr(acc)
  d <- file.path(DER, acc)
  if (!dir.exists(d)) { cat("  !! 目录不存在，跳过\n"); next }

  fe <- list.files(d, pattern = "_expr_symbol.*\\.rds$", full.names = TRUE)
  fc <- list.files(d, pattern = "_clinical_clean\\.csv$", full.names = TRUE)
  fm <- list.files(d, pattern = "^cohort_manifest_.*\\.csv$", full.names = TRUE)
  if (!length(fe) || !length(fc)) { cat("  !! 数据文件不全，跳过\n"); next }

  expr <- readRDS(fe[1]); clin <- fread(fc[1])
  old  <- if (length(fm)) fread(fm[1]) else data.table()
  sup  <- if (!is.null(SUPP[[acc]])) SUPP[[acc]] else list()

  # ---------- 从数据现算（不引用旧 manifest 的表面值）
  smp <- expr[seq_len(min(3000, nrow(expr))), seq_len(min(30, ncol(expr))), drop = FALSE]
  e_med <- median(smp, na.rm = TRUE); e_max <- max(smp, na.rm = TRUE)
  is_log <- e_max < 25
  id_col <- intersect(c("patient_id","gsm","sample_id"), colnames(clin))[1]
  stopifnot(identical(colnames(expr), as.character(clin[[id_col]])))

  has_surv <- "OS_time" %in% colnames(clin) && any(!is.na(clin$OS_time))
  n_ev <- if (has_surv) sum(clin$OS_event, na.rm = TRUE) else NA_integer_
  mfu  <- if (has_surv) round(median(clin$OS_time, na.rm = TRUE), 1) else NA_real_
  n_st <- if ("stage" %in% colnames(clin)) sum(!is.na(clin$stage)) else NA_integer_

  r <- data.table(
    cohort            = as.character(sup$cohort %||% pick(old, "cohort", acc)),
    platform          = as.character(sup$platform %||% pick(old, "platform")),
    n_raw             = as.integer(sup$n_raw %||% pick(old, "n_raw")),
    n_after_tumor     = as.integer(sup$n_after_tumor %||% pick(old, "n_after_tumor")),
    n_after_histology = as.integer(sup$n_after_histology %||% pick(old, "n_after_histology")),
    n_final           = ncol(expr),                       # ★ 现算
    n_events          = as.integer(n_ev),                 # ★ 现算
    event_rate        = if (has_surv) sprintf("%.1f%%", 100*n_ev/ncol(expr)) else NA_character_,
    n_genes           = nrow(expr),                       # ★ 现算
    expr_scale        = ifelse(is_log, "log2", "linear"), # ★ 现算
    expr_median       = round(e_med, 3),                  # ★ 现算
    expr_max          = round(e_max, 1),                  # ★ 现算
    need_log          = !is_log,                          # ★ 现算
    has_survival      = has_surv,                         # ★ 现算
    time_key          = as.character(sup$time_key %||% pick(old, "time_key")),
    time_src_unit     = as.character(sup$time_src_unit %||% pick(old, "time_src_unit")),
    time_mult         = as.numeric(sup$time_mult %||% pick(old, "time_mult")),
    time_unit         = ifelse(has_surv, "days", NA_character_),
    median_followup_days = mfu,                           # ★ 现算
    stage_type        = as.character(sup$stage_type %||% pick(old, "stage_type")),
    n_with_stage      = as.integer(n_st),                 # ★ 现算
    histology_filter  = as.character(sup$histology_filter %||%
                          pick(old, "histology_filter", "无（全部为 LUAD）")),
    biotype_filter    = as.character(sup$biotype_filter %||% pick(old, "biotype", "无")),
    assay_used        = as.character(sup$assay_used %||% pick(old, "assay")),
    sample_dedup      = as.character(sup$sample_dedup %||%
                          pick(old, "dedup_rule", "读入即按 ID 去重")),
    gene_collapse     = as.character(sup$gene_collapse %||%
                          pick(old, "gene_collapse", "同名取均表达最大")),
    anno              = as.character(sup$anno %||% pick(old, "anno")),
    source_md5        = as.character(pick(old, "source_md5",
                          pick(old, "source_md5_expr"))),
    built             = as.character(pick(old, "built")),
    harmonized        = format(Sys.time(), "%Y-%m-%d %H:%M")
  )
  setcolorder(r, SCHEMA)

  # ---------- 与旧 manifest 的差异（现算值 vs 旧记录）
  cat("  现算 vs 旧 manifest：\n")
  for (f in c("n_final","n_events","n_genes","need_log","has_survival","n_with_stage")) {
    o <- pick(old, f, pick(old, "n_samples"))
    n <- r[[f]][1]
    flag <- if (is.na(o)) "旧值缺失" else if (identical(as.character(o), as.character(n))) "一致" else "!! 不一致"
    cat(sprintf("    %-16s 旧=%-14s 新=%-14s %s\n", f,
                ifelse(is.na(o),"-",as.character(o)), as.character(n), flag))
  }

  # ---------- 备份并写出
  if (length(fm)) {
    bak <- paste0(fm[1], ".bak_", STAMP)
    if (!file.exists(bak)) file.copy(fm[1], bak)
    cat("  已备份旧 manifest -> ", basename(bak), "\n", sep="")
  }
  out_f <- file.path(d, paste0("cohort_manifest_", sub("_LUAD$", "", acc), ".csv"))
  # 保持与既有文件名一致，避免 load_cohort 匹配到两个
  if (length(fm)) out_f <- fm[1]
  write.csv(r, out_f, row.names = FALSE)
  cat("  已写出 -> ", basename(out_f), "\n", sep="")
  rows[[acc]] <- r
}

# ============================================================ registry
hr("合并 registry")
reg <- rbindlist(rows, fill = TRUE)
setcolorder(reg, SCHEMA)
f <- file.path(DER, "cohort_registry.csv")
write.csv(reg, f, row.names = FALSE)
kv("写出", f)

cat("\n  核心字段一览：\n\n")
show <- c("cohort","n_final","n_events","event_rate","n_genes",
          "expr_scale","need_log","time_src_unit","has_survival","n_with_stage")
print(reg[, ..show], row.names = FALSE)

hr("汇总")
sv <- reg[has_survival == TRUE]
kv("有生存的队列数", nrow(sv))
kv("合计样本", sum(sv$n_final))
kv("合计事件", sum(sv$n_events))
kv("事件率区间", sprintf("%.1f%% – %.1f%%",
    100*min(sv$n_events/sv$n_final), 100*max(sv$n_events/sv$n_final)))
kv("需 log 变换的队列", paste(reg$cohort[reg$need_log], collapse=", "))
kv("时间源单位", paste(unique(na.omit(sv$time_src_unit)), collapse=" / "))

na_cnt <- colSums(is.na(reg))
if (any(na_cnt > 0)) {
  cat("\n  仍有缺失的字段（缺失队列数）：\n")
  for (n in names(na_cnt)[na_cnt > 0]) cat(sprintf("    %-22s %d\n", n, na_cnt[n]))
  cat("  -> CPTAC 的时间相关字段缺失是预期的（无生存数据）\n")
}

hr("验证：协议现在应读声明而非推断")
cat("  重跑自检确认：\n")
cat("    Rscript ~/thesis_LUAD_CAF/01_bases/A_bulk_OS/scripts/scoring_protocol.R\n")
cat("  第 1 节 TCGA_LUAD 的『log 变换』应仍为『是』，但来源已从推断改为 manifest 声明。\n")

cat("\n输出: ", OUT, "\n", sep="")
