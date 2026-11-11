# 阶段 3 事实清单（20261025，sci2_stage3_facts.R 现算）

来源文件：
- stage3_readouts_long_20261024.csv（md5 349081cc890715851e6d09a2e6fba1d3）
- stage3_cox_20261024.csv（md5 d802746107a03ff80e6a77b947a2e9e4）
- stage3_consistency_20261024.csv（md5 41ad99bd5426bfb3b84057714ed53fe5）
- stage3_readout_notes_20261024.csv（md5 f1bf07d8fc738a6a894f8d4a09633939）

## G1 相关结构的跨队列一致性（共有读数，15 对队列）
- cophenetic 相关中位 0.76（0.57 至 0.86）；相关矩阵上三角 Spearman 中位 0.90（0.86 至 0.95）。最低的 3 对：TCGA_LUAD–GSE68465 0.57；GSE68465–GSE30219 0.61；GSE31210–GSE30219 0.61。
- 不可写：『签名分为 X 类』（未设类数 k）。叶序相邻不等于同类。

## G2 合理性检查
- SIG-002（MCP-counter Fibroblasts 基因）与 MCP_Fibroblasts 的 Spearman：0.98 至 1.00。
- TCGA：MCP_Fibroblasts 与 EPIC_CAFs 0.92；xCell_Fibroblasts 与 EPIC_CAFs 0.37。

## G3 两个『成纤维细胞』读数
- MCP_Fibroblasts 与 xCell_Fibroblasts 的 Spearman：0.19 至 0.58（GSE30219 0.42，GSE31210 0.19，GSE50081 0.48，GSE68465 0.41，GSE72094 0.58，TCGA_LUAD 0.49）。
- 与增殖指数：MCP_Fibroblasts 0.03 至 0.37；xCell_Fibroblasts -0.54 至 -0.32（为负的队列 6/6）。xCell_Fibroblasts 与 xCell_StromaScore 0.74 至 0.92。
- 单变量 Cox z：MCP_Fibroblasts -0.39 至 2.46（>0 的队列 5/6，|z|≥1.96 的 2）；xCell_Fibroblasts -2.33 至 -0.46（<0 的队列 6/6，|z|≥1.96 的 1）。
- 可写：两种成纤维细胞读数中度相关，但与增殖的相关方向不同，预后关联方向相反。
- 不可写：『成纤维细胞具有保护 / 促进作用』；『xCell 的负向关联由增殖导致』（未做校正检验）。

## G4 基质含量类签名（SIG-002/003/009/010）
- 与 MCP_Fibroblasts 的 Spearman 0.64 至 1.00；与增殖指数 -0.32 至 0.29；TCGA 中与 ABSOLUTE 纯度 -0.63 至 -0.45。
- 6 队列中 |z|≥1.96 的队列数：SIG-002 1，SIG-003 0，SIG-009 1，SIG-010 0。
- SIG-003 与 estimate 包 StromalSignature 共有 136/141 基因（阶段 3 日志第 0 节）。

## G5 有信号签名与增殖指数（阶段 2 F3 的阶段 3 复核）
- SIG-001/004/006/007/008/011/012 与增殖指数的 Spearman（6 队列）：-0.66 至 0.92。

## G6 预后横向比较：增殖指数作为参照读数
- 增殖指数单变量 z：2.54 至 3.94，|z|≥1.96 的队列 6/6。
- GEO 5 队列中，上述 7 个签名的可评估组合 34 个，|z| 超过同队列增殖指数的 13 个。
- TCGA 中 7 个签名 |z| 超过增殖指数（3.94）的 5 个：SIG-001 4.99，SIG-004 6.75，SIG-006 4.50，SIG-007 8.04，SIG-011 7.27。TCGA 是 SIG-001/004/006/011/012 的开发队列。
- 不可写：『签名不优于 / 优于增殖指数』（SAP 未预设 z 的比较检验，仅描述）。
