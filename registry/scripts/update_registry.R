#!/usr/bin/env Rscript
# ==============================================================================
# update_registry.R —— 更新 signature 收录表与排除表
#
# 字段设计说明：
#   evidence_level  权重/基因列表的证据强度
#     1 = 从作者提供的数据反解（R2=1，最强）
#     2 = 包内文件 + md5 锁定
#     3 = 原文正文/表格，自洽性检查通过
#     4 = 原文正文，但存在文本矛盾（已标注）
#     5 = 无权重需求（A 类等权 ssGSEA）
#   caf_definition  该 signature 如何定义 CAF —— 用于结果分组解读
#   overlap_cohorts 原文用过的、与底座 A 重叠的队列 —— 结果解读时需带上
#
# 只写 03_signature_registry/ 之下。
# ==============================================================================

d <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")

reg <- data.frame(stringsAsFactors = FALSE,
  sig_id = sprintf("SIG-%03d", 1:10),
  name = c(
    "Ren 2023 CAF 4-gene model",
    "MCP-counter Fibroblasts",
    "ESTIMATE StromalSignature",
    "Zhang & Cheng 2023 CAFRS",
    "Navab 2011 CAF prognostic signature",
    "Si 2024 CAF 8-gene model",
    "Wang(LDDD) 2026 CAF 8-gene model",
    "Wang 2023 CAF 5-gene model",
    "Ruhle 2023 CAF MARKER",
    "Ruhle 2023 CAF (EPIC)"),
  class = c("B","B","A","B","A","B","B","B","A","A"),
  n_genes = c(4L, 8L, 141L, 7L, 11L, 8L, 8L, 5L, 27L, 20L),
  caf_definition = c(
    "TCGA 挖掘（单细胞 CAF cluster 相关基因 -> Cox/LASSO）",
    "解卷积标志基因（转录组特异性排序）",
    "泛癌基质基因集（不区分细胞类型）",
    "TCGA 挖掘（CAF 相关基因 -> Cox/LASSO）",
    "原代 CAF vs 配对正常成纤维细胞培养对比",
    "TCGA 挖掘（CAF 相关基因 -> Cox/LASSO）",
    "TCGA 挖掘（CAF 相关基因 -> 机器学习）",
    "WGCNA 模块基因 -> Cox/LASSO",
    "综述总结的公认 CAF 标志基因",
    "解卷积参比谱衍生的 CAF 特异基因"),
  source = c(
    "Front Immunol 2023;14:1201573",
    "MCPcounter GitHub master Signatures/genes.txt",
    "estimate 1.0.13 extdata/SI_geneset.gmt",
    "Aging (Albany NY) 2023;15(11):4986-5006",
    "PNAS 2011;108(17):7160-7165",
    "Sci Rep 2024;14:23243",
    "Lett Drug Des Discov 2026 (PII S1570180826000631)",
    "Cancer Immunol Immunother 2023;72:2423-2442",
    "OncoImmunology 2023;12(1):2274130 Table S1",
    "OncoImmunology 2023;12(1):2274130 Table S1"),
  doi = c(
    "10.3389/fimmu.2023.1201573", "-", "-",
    "10.18632/aging.204774", "10.1073/pnas.1014506108",
    "10.1038/s41598-024-74336-1", "待核实",
    "10.1007/s00262-023-03428-0",
    "10.1080/2162402X.2023.2274130", "10.1080/2162402X.2023.2274130"),
  evidence_level = c(1L, 2L, 2L, 4L, 3L, 4L, 3L, 3L, 5L, 5L),
  provenance = c(
    "权重由作者 GitHub 提供的 riskscore 反解，R2=1.0000000000，残差 1.27e-14",
    "md5 944af31812c3f2aa61526074604e1c77，下载于 2026-08-30",
    "md5 9628130acbd541d8a33f21b4ab1e864f",
    "权重取自正文 Results 公式",
    "基因列表取自正文 Table 1（13 探针对应 11 基因）",
    "权重取自正文公式",
    "权重取自正文公式",
    "权重取自正文 Results 公式（beta 系数）",
    "基因列表取自补充材料 Table S1（嵌入图，已人工转录）",
    "基因列表取自补充材料 Table S1（嵌入图，已人工转录）"),
  overlap_cohorts = c(
    "TCGA", "-", "-", "TCGA, GSE41271", "-",
    "TCGA, GSE31210, GSE37745, GSE50081, GSE68465, GSE3141",
    "待核实", "TCGA, GSE41271, GSE37745, GSE42127",
    "TCGA", "TCGA"),
  note = c(
    "原文 Methods 与 Results 的 CLEC3B 符号相反；反解证明 Results 版(-0.123)正确。风险基因 EXO1/CCNB1 均为增殖基因",
    "等权，权重全 1。包内不含基因表，运行时从 GitHub 拉取，故须锁定 md5",
    "★ 不得用作实验二主分析的纯度协变量（SAP 7.2）。同时作为被测对象与潜在混杂读数，两个角色不可混用",
    "正文写 SERINE2，补充图 1 写 SERPINE2，以后者为准；正文方向描述与三个负系数矛盾（写反）",
    "唯一有生物学锚点者：15 例配对原代培养对比。全平台 100% 覆盖（在 HG-U133 上开发）",
    "正文公式 CLEC3B=+0.077，但图 4E 标为 Protective 且条形在负轴；照正文原样使用。与 SIG-001 共享 CLEC3B/EXO1/CCNB1",
    "期刊为 Letters in Drug Design & Discovery（Bentham），学科定位为药物设计，非肿瘤/生信领域。LATS2 被列为风险基因(+0.336)，与其公认抑癌基因身份相反",
    "原文使用跨队列 ComBat 校正 + 固定切点 2.01527；本研究改为队列内 Z 标准化（SAP 4.1）。训练集报 AUC 0.64(2y)/0.67(4y)",
    "教科书式 CAF 标志，非统计筛选产物。与 SIG-010 同出一表但仅 1 个基因重叠",
    "该基因集使我们此前因'需自定拆分规则'而排除的 EPIC 变为可纳入 —— 此处的拆分由原作者完成并发表")
)

exc <- data.frame(stringsAsFactors = FALSE,
  name = c("Cords 2024","EPIC CAF (包内)","xCell Fibroblasts","Kong 2024","Chen 2025"),
  source = c(
    "Cancer Cell 2024;42(3):396-412",
    "EPIC 1.1.7 (TRef$sigGenes)",
    "xCell 1.1.0 (xCell.data$signatures)",
    "Mol Genet Genomics 2024;299:70",
    "Int J Med Sci 2025;22(12):2956-2972"),
  criterion = c("E1","I1","I1","E2","E1"),
  reason = c(
    "成像质谱流式蛋白分型，无可提取的转录组基因列表",
    "sigGenes(98) 为七种细胞类型共用；按 refProfiles 拆出 CAF 专属(20)需自定规则，非已发表列表。注：Ruhle 2023 已发表其拆分结果，见 SIG-010",
    "九套特征集(ENCODE/FANTOM/HPCA 各 3)，九套交集为 0，基因数 17-136；任何合并规则均为本研究自定",
    "仅 2 个基因(COL5A2, COL6A2)，基因数 < 3",
    "28 基因模型，完整名单未公开。正文仅可视化 top5/top10(Fig 7A-B，图中仅可读出 9 个)；补充材料含 4 图 2 表，无基因列表；数据可用性声明称代码需向通讯作者索取。已核查全文正文、全部图表与补充材料。已发邮件索取，待回复")
)

write.csv(reg, file.path(d, "registry_notes.csv"), row.names = FALSE)
write.csv(exc, file.path(d, "excluded_notes.csv"), row.names = FALSE)

hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")

hr("已收录")
print(reg[, c("sig_id","name","class","n_genes","evidence_level")], row.names = FALSE)

hr("按类别")
cat("  A 类（等权 ssGSEA，用零分布 A/B）: ", sum(reg$class=="A"), " 个，基因数 ",
    paste(reg$n_genes[reg$class=="A"], collapse=" / "), "\n", sep="")
cat("  B 类（原文权重加权和，用 5.3 复合零分布）: ", sum(reg$class=="B"), " 个，基因数 ",
    paste(reg$n_genes[reg$class=="B"], collapse=" / "), "\n", sep="")

hr("按 CAF 定义方式（供结果分组解读）")
grp <- c("原代细胞对比"="SIG-005", "教科书标志"="SIG-009",
         "解卷积衍生"="SIG-002, SIG-010", "泛基质"="SIG-003",
         "统计筛选（TCGA/WGCNA 挖掘）"="SIG-001, SIG-004, SIG-006, SIG-007, SIG-008")
for (n in names(grp)) cat(sprintf("  %-28s %s\n", n, grp[n]))

hr("证据等级分布")
lv <- c("1 反解(R2=1)","2 包内文件+md5","3 正文，自洽","4 正文，有矛盾","5 无权重需求(A类)")
for (i in 1:5) {
  ids <- reg$sig_id[reg$evidence_level==i]
  cat(sprintf("  %-18s %s\n", lv[i], ifelse(length(ids), paste(ids, collapse=", "), "—")))
}

hr("已排除")
print(exc[, c("name","criterion")], row.names = FALSE)

hr("汇总")
cat("  候选总数 ", nrow(reg)+nrow(exc), "  |  纳入 ", nrow(reg), "  |  排除 ", nrow(exc), "\n", sep="")
cat("  SAP 3.3 目标区间 12-20，当前 ", nrow(reg), "，", 
    ifelse(nrow(reg)>=12, "已达下限", paste0("距下限还差 ", 12-nrow(reg), " 个")), "\n", sep="")
cat("\n  写出:\n    registry_notes.csv\n    excluded_notes.csv\n")
