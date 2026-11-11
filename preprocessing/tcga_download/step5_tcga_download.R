# =============================================================================
# step5_tcga_download.R
# 下载 TCGA-LUAD RNA-seq 数据及临床信息
# =============================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

DATA_DIR <- path.expand("~/caf_spatial_project/data/TCGA_LUAD")
dir.create(DATA_DIR, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# STEP 1: 下载 RNA-seq 数据（STAR counts）
# =============================================================================
cat("\n[STEP 1] 查询 TCGA-LUAD RNA-seq 数据...\n")

query_rna <- GDCquery(
  project           = "TCGA-LUAD",
  data.category     = "Transcriptome Profiling",
  data.type         = "Gene Expression Quantification",
  workflow.type     = "STAR - Counts",
  sample.type       = "Primary Tumor"
)

cat(sprintf("  找到样本数: %d\n", nrow(getResults(query_rna))))

cat("\n[STEP 2] 下载数据...\n")
GDCdownload(
  query     = query_rna,
  directory = DATA_DIR,
  method    = "api"
)

cat("\n[STEP 3] 整理数据为 SummarizedExperiment...\n")
se <- GDCprepare(
  query     = query_rna,
  directory = DATA_DIR
)

cat(sprintf("  样本数: %d\n", ncol(se)))
cat(sprintf("  基因数: %d\n", nrow(se)))

saveRDS(se, file.path(DATA_DIR, "TCGA_LUAD_SE.rds"))
cat("\n[STEP 4] 下载临床数据...\n")

clinical <- GDCquery_clinic(
  project = "TCGA-LUAD",
  type    = "clinical"
)

cat(sprintf("  临床样本数: %d\n", nrow(clinical)))
write.csv(clinical,
          file.path(DATA_DIR, "TCGA_LUAD_clinical.csv"),
          row.names=FALSE)

cat("\n✅ TCGA-LUAD 数据下载完成！\n")
cat(sprintf("  → %s/TCGA_LUAD_SE.rds\n",      DATA_DIR))
cat(sprintf("  → %s/TCGA_LUAD_clinical.csv\n", DATA_DIR))
