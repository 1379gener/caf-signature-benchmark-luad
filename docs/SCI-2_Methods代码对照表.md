# SCI-2 Methods–代码对照表（v5，2026-10-25；⚠️ 已清零）

> 依据：`/tmp/sci2_methods_evidence_20261025_1348.txt`（取证 1）与 `/tmp/sci2_methods_evidence2_20261025_1442.txt`（取证 2）。行号均指取证中所列脚本版本（md5 见取证；null_engine.R 7ef1e0b5…，03_baseA_tcga.R 2c477c58…，SAP_v1.2.txt af0b8664…）。
> 状态：✅ 有代码证据且与拟写 Methods 一致；⚠️ 证据不足 / 代码与 SAP 不一致 / 计划了但未执行。
> 放置：`~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/text/SCI-2_Methods代码对照表.md`。投稿前 ⚠️ 必须清零。

## 1. 队列与数据底座

| # | Methods 拟写（要点） | 证据 | 状态 |
|---|---|---|---|
| 1.1 | TCGA-LUAD 来源为 GDC 的 SummarizedExperiment；表达用 tpm_unstrand，log2(TPM+1) | 03_baseA_tcga.R L11、L57；08_harmonize L36（assay_used = tpm_unstrand）；scoring_protocol L69–73（need_log 时 log2(x+1)） | ✅ |
| 1.2 | TCGA 只保留 gene_type 为 protein_coding、且 gene_name 非空的行，再按 symbol 合并 | 03_baseA_tcga.R L29、L105–110（过滤）、L118–125（合并）；manifest n_genes 19,938（取证 2 C） | ✅ |
| 1.3 | 同一病人多个 aliquot 取文库最大者；去重后断言无重复 | 03_baseA_tcga.R L25、L84、L92 | ✅ |
| 1.4 | 同名基因取平均表达最高的一行 | TCGA L28、L118–122；GEO 05_rev3 L342–344；CPTAC 06 L92–95 | ✅ |
| 1.5 | OS：死亡者用 days_to_death，其余用 days_to_last_follow_up | 03_baseA_tcga.R L142–144 | ✅ |
| 1.6 | 纳入全分期 | 03_baseA_tcga.R L30 | ✅ |
| 1.7 | TCGA：GDC SE 共 540 个原发肿瘤 aliquot，来自 517 名病人（12 人有多个 aliquot）；每人保留文库最大的一个 → 517 人；排除 OS 时间缺失 70 人（存活 66、死亡 4）和 OS ≤ 0 天 4 人 → 443 人、182 事件 | 03_baseA_tcga.R L59–92、L163、L188–194；build_20260827_2333.log L22/65/100–103/118；取证 3 B2 从 SE 独立复算一致，同一病人各 aliquot 的生存字段无不一致。**注意**：manifest 把 517 记在 n_after_tumor 列，实际是去重后人数，写 Methods 和流程图时不要照搬这个列名 | ✅ |
| 1.8 | GEO：characteristics 按键名解析；时间统一换算为天（月 ×30.4375，年 ×365.25） | 05_rev3 L11–12、L44/55/70/81/92、L268 | ✅ |
| 1.9 | 组织学精确匹配：GSE50081 为 "adenocarcinoma"，GSE30219 为 "ADC"；GSE68465 按 disease_state 取 "Lung Adenocarcinoma" | 05_rev3 L79、L90、L41 | ✅ |
| 1.10 | GSE31210：246 样本中 226 例为肿瘤；遵从原作者 exclude 标志排除 22 例 → 204 例、30 事件 | 05_rev3 L58–59、L297–306；geo_build_20260905_0009.log L19、L33；verify_ch04 L68–69 | ✅ |
| 1.11 | GEO 有效样本：通过筛选、未被 exclude、OS>0、生存状态非缺失 | 05_rev3 L306 | ✅ |
| 1.12 | 探针到基因：GSE68465（GPL96）用 hgu133a.db；GSE31210 / GSE50081 / GSE30219（GPL570）用 hgu133plus2.db；GSE72094（GPL15048）用 GPL SOFT 注释表；CPTAC 源表已是基因层面。去掉无符号探针，每个探针只映射一次 | 05_rev3 L180–189、L199；manifest anno 列（取证 2 C） | ✅ |
| 1.13 | log 尺度判定：最大值 ≥ 25 视为线性尺度，做 log2(x+1) | 05_rev3 L332–335、L373–374；scoring_protocol L69–73 | ✅ |
| 1.14 | CPTAC 没有可用生存数据，不进入任何生存分析 | SAP L206–212；scoring_protocol L75–77 | ✅ |
| 1.15 | 读入每个队列时断言表达矩阵列与临床表行一一对应、无重复 | scoring_protocol L65–66 | ✅ |
| 1.16 | GSE72094 n = 398、113 事件（与 SCI-1 一致） | 05_rev3 L73、L312–317 | ✅ |
| 1.17 | 各 GEO 队列样本流转：GSE68465 462 → 组织类型筛选 443 → 有可用 OS 442（236 事件）；GSE72094 442 → 442 → 398（113）；GSE31210 246 → 226 → exclude 标志 −22 → 204（30）；GSE50081 181 → 严格匹配 adenocarcinoma 127 → 127（51）；GSE30219 307 → 严格匹配 ADC 85 → 83（43）。「有可用 OS」= 时间非缺失且 > 0、生存状态非缺失 | geo_build_20260829_1315.log（L7–33、69–99、105–141、147–184）；geo_build_20260905_0009.log L19/33；05_rev3 L306 | ✅ |

## 2. 签名注册与锁定

| # | Methods 拟写 | 证据 | 状态 |
|---|---|---|---|
| 2.1 | 12 个签名，在任何零分布运行前锁定；锁定表 md5 4ad3218931d95dfb6ee67a7ba07b9a75（2026-09-03） | lock_registry.R L3–21；锁定表名与 md5 见取证 L518；run_benchmark L64 拒绝在未锁定时运行 | ✅ |
| 2.2 | 分类规则：A 类 ≥10 个基因，等权 ssGSEA；B 类 3–9 个基因，必须有原文权重 | lock_registry.R L11–12 | ✅ |
| 2.3 | 任意两个签名重叠 >90% 时只保留其一 | lock_registry.R L13、L132–141 | ✅ |
| 2.4 | 基因符号统一到现行 HGNC：只在别名唯一映射到一个 Entrez ID 时改写 | normalize_symbols.R L12、L53–73；symbol_mapping_20260903.csv | ✅ |
| 2.5 | B 类权重：能拿到作者数据的，用最小二乘从风险分反解（R² = 1）；其余取原文数值 | recover_weights.R L109–113；recover_ren2023_coef.R L44–45；锁定表 evidence_level 字段 | ✅ |
| 2.6 | SIG-012 原文权重敏感性分析（事先登记）：原文权重逐基因 Z 加权和，§5.3 加权复合零分布 1000 次；主分析不变。结果按队列角色报告：强于随机 = TCGA（原文训练队列）+ GSE30219/GSE31210/GSE72094（原文验证队列）；原文未用过的 GSE68465、GSE50081 两种打分均为不特异 | sci2_sens_sig012_origweights.R（md5 14abe9dd…）→ sens_SIG012_origweights_20261025.csv；闸门：等权 z 与阶段 1 一致；角色取自 train_cohort_annotation.csv 与锁定表 overlap_cohorts | ✅（DEVIATIONS 2026-10-25 D） |
| 2.7 | 解卷积来源的基因集：MCP-counter（GitHub master）、ESTIMATE（estimate 1.0.13 SI_geneset.gmt）、EPIC 衍生 | update_registry_v2 L52–53、L111–112；export_deconv_signatures.R | ✅ |

## 3. 打分与 Cox

| # | Methods 拟写 | 证据 | 状态 |
|---|---|---|---|
| 3.1 | A 类：GSVA 1.48.3 的 ssGSEA，ssgsea.norm = FALSE，在 log2 表达矩阵上打分 | null_engine L218–225；scoring_protocol L28–29、L154–163；版本见取证 F | ✅ |
| 3.2 | B 类：各基因先在队列内 Z 标准化，再按原文权重加权求和；缺失基因不重新分配权重 | null_engine L229–235 | ✅ |
| 3.3 | 打分基因 = 签名基因（精确符号）∩ 候选池；别名映射只用于判断覆盖率 | run_benchmark L158；DEVIATIONS 2026-10-16 第 2 条 | ✅（已登记偏差） |
| 3.4 | 覆盖率（含别名回退）< 70% 的组合判为不可评估 | scoring_protocol L33、L132–145；run_benchmark L141–148 | ✅ |
| 3.5 | 单变量 Cox：分数在分析样本内标准化（per SD），报告 Wald z；有效样本为 OS>0 且结局非缺失 | null_engine L252–259 | ✅ |
| 3.6 | 事件 < 25 的队列降为敏感性分析队列（事先规定）；本研究无队列触发 | run_benchmark L97–98、L147；verify_ch04 L70–71 | ✅ |

## 4. 零分布（阶段 1）

| # | Methods 拟写 | 证据 | 状态 |
|---|---|---|---|
| 4.1 | 候选池：在 ≥50% 样本中表达值 > 0 的基因 | null_engine L28、L46–47。**事实：芯片队列 log2 值全部 > 0，此条只在 TCGA 起作用（19,938 → 17,362）**，Methods 须如实写明 | ✅（写法须改） |
| 4.2 | 在 log2 尺度上按平均表达和方差各分 10 个等频箱，组成联合箱 | null_engine L27、L49–61 | ✅ |
| 4.3 | 零分布 A：为每个签名基因从同一联合箱不放回抽一个替代基因（不取签名自身基因）；该箱不够时依次退到平均表达 ±1、±2 箱（方差箱不变）；重复 1000 次。实测 71 个组合均未发生退让 | null_engine L77–117；verify_nullengine R3/E1（相邻箱 0、全池兜底 0）。代码另有 SAP 未写的全池兜底分支（L111），从未触发 | ✅（DEVIATIONS 2026-10-25 (1)） |
| 4.4 | 零分布 B：随机取签名中一个基因的联合箱抽种子基因，取与其 Spearman 相关最高的 k−1 个基因组成模块（种子与模块均不取签名自身基因）；模块平均表达与平均方差须在真实签名 ±20% 内，否则重抽，至多 50 次；第 51 次直接采用。71,000 个成员中 142 个（0.20%）为第 51 次采用，其中 118 个超出 ±35%；按最坏情形计算，剔除这些成员不可能改变任何组合 p_emp_B < 0.05 的判定 | null_engine L131–167；verify_nullengine E2/E3、R7；sci2_nullB_relax_bound.R → nullB_relax_bound_20261025.csv（两种口径均 0/71） | ✅（DEVIATIONS 2026-10-25 (2)） |
| 4.5 | 零分布 B 质量验证：零分布模块内平均相关的中位数 ≥ 真实签名的 70%。实现用前 200 个模块；用全部 1000 个复核，71 个组合判定无一改变 | null_engine L174–179；verify_nullengine R5/R6/F1 | ✅（DEVIATIONS 2026-10-25 (3)） |
| 4.6 | 共表达用行标准化的秩矩阵计算（等价于 Spearman） | null_engine L66–70、L122–126 | ✅ |
| 4.7 | 零分布 B 只对 A 类签名报告；B 类签名只报零分布 A（即 §5.3 第 1 条加权复合零分布，verify_nullengine G1/G2 验证）；B 类的零分布 B 结果放补充材料并注明 SAP 未规定、权重按模块位置分配；A、B 类分表 | DEVIATIONS 2026-10-25 G；SAP L301–305、L405–409、L491–498 | ✅（Limitations 须写 SAP 前后不一致） |
| 4.8 | 经验 p = (1 + #{\|z_null\| ≥ \|z_real\|}) / (1 + N)，双侧 | null_engine L263–266 | ✅ |
| 4.9 | seed 由（队列, 签名, 零分布类型）确定，基数 20260830 | null_engine L33–39 | ✅ |
| 4.10 | 分类：\|z\| < 1.96 → 无信号；否则 p_emp_A ≥ 0.05 → 不特异；否则 → 强于随机 | diagnose_results L93–94 | ✅ |
| 4.11 | pass_B_main 只计入共表达质量验证通过的组合；pass_B_sens 不过滤 | diagnose_results L147–149 | ✅ |
| 4.12 | GSE30219 × SIG-001/002 的零分布次数从 50 补足到 1000 | rerun_gse30219_n1000.sh；DEVIATIONS 2026-10-16 | ✅ |
| 4.13 | SAP §12：随机抽 20 个零分布成员（抽样 seed 20261025），从表达矩阵独立复算 Cox z，与引擎输出一致（最大差 0，容差 1e-6） | verify_nullengine H1；verify_nullengine_members_20261025.csv | ✅ |

## 5. 队列角色与过拟合（事后探索性）

| # | Methods 拟写 | 证据 | 状态 |
|---|---|---|---|
| 5.1 | 训练队列取自各原文 Methods 的人工标注；「原文用过」取自锁定前登记的 overlap_cohorts；SIG-007 两项锁定时为「待核实」，2026-10-25 据原文补正（训练 TCGA，验证 GSE30219） | analyze_cohort_role_20261025.R（L49、L90）；train_cohort_annotation_20261025.csv；DEVIATIONS 2026-10-25 SIG-007 | ✅ |
| 5.2 | 过拟合比 = \|开发队列 z\| / \|外部队列 z 的中位数\|；结果 3.17（SIG-007）/ 3.06 / 2.67 / 1.85 / 1.35 / 1.26；通过率 开发 6/6、原文验证 3/13、完全独立 6/52 | analyze_cohort_role_20261025.R；results_by_cohort_role_20261025.csv、overfit_gradient_20261025.csv（旧 20261016 版作废，SIG-007 被误计为完全独立） | ✅ |
| 5.3 | 该拆分为 SAP 之外的事后分析 | analyze_cohort_role L9–16；SAP §14 修订 5 | ✅ |
| 5.4 | SIG-007 出处：Li W, Gui T, Hou Y, Lett Drug Des Discov, doi 10.1016/j.lddd.2026.100333；TCGA 训练、GSE30219 验证；8 个系数与 weights 文件一致。锁定表名称「Wang」、note「Bentham」有误，锁定表不改 | 原文 §2.1、§3.4；DEVIATIONS 2026-10-25 SIG-007 | ✅（出版年【待核实】：页眉 2025，DOI 2026） |

## 6. 混杂分解（阶段 2）

| # | Methods 拟写 | 证据 | 状态 |
|---|---|---|---|
| 6.1 | M1–M4 的协变量组成 | stage2 L52–54、L91、L93–98 | ✅ |
| 6.2 | 分数、纯度、增殖、年龄在各模型自己的分析样本内标准化；性别作因子；分期为二分类 | stage2 L48–51、L77；stage_coding L16–18 | ✅ |
| 6.3 | 纯度：Aran 2015 Supplementary Data 1 的 ABSOLUTE，LUAD 原发灶，按病人去重 | purity_tcga L27–44 | ✅ |
| 6.4 | 增殖指数：131 个 meta-PCNA 基因 log2 表达的逐样本中位数；精确匹配加唯一别名回退 | prolif_index L41–45 | ✅ |
| 6.5 | TCGA 分期：主分析用 GDC；敏感性分析用 paper_Tumor.stage 补缺失，GDC 已有的值不改 | tcga_sens L34–38 | ✅ |
| 6.6 | 衰减幅度用同一样本上重新拟合的 M1 作分母；只在 M1 的 \|z\| ≥ 1.96 时报告 | stage2 L86–124、L122 | ✅ |
| 6.7 | 全样本 M1 与阶段 1 的 z_real 一致（闸门） | stage2 L109–111 | ✅ |
| 6.8 | GEO 队列不计算 M2/M4（没有与表达无关的纯度来源） | stage2 L84；SAP §7.2 | ✅ |
| 6.9 | 共线性诊断：VIF = 1/(1−r²)，并报告实测 SE 比 | diagnostics L54–64 | ✅ |

## 7. 冗余与相关（阶段 3）及事后探索

| # | Methods 拟写 | 证据 | 状态 |
|---|---|---|---|
| 7.1 | 外部读数：MCP-counter 10 类（按 HUGO 符号）；xCell 4 项（TCGA 设 rnaseq = TRUE，芯片设 FALSE）；ESTIMATE ImmuneSignature 用同一 ssGSEA 路径；增殖指数 | redundancy L75–83 | ✅ |
| 7.2 | EPIC 只在 TCGA 计算，输入为 2^x − 1（还原到 TPM） | redundancy L86–87 | ✅ |
| 7.3 | 两两 Spearman（成对完整观测）；剔除无变异的读数 | redundancy L104–106 | ✅ |
| 7.4 | 层次聚类：距离 1 − \|ρ\|，ward.D2 | redundancy L107 | ✅ |
| 7.5 | 跨队列一致性：在两队列共有的读数上重新聚类，算 cophenetic 相关和相关矩阵上三角的 Spearman；EPIC 和纯度不参与 | redundancy L116–122 | ✅ |
| 7.6 | 每个读数在每个队列做单变量 Cox（per SD） | redundancy L130–137（取证只到 L137） | ✅ |
| 7.7 | 事后探索：两种成纤维细胞读数在单变量和加增殖指数后的 Cox | explore_fibro L20–40；DEVIATIONS 2026-10-25 | ✅ |

## 8. 软件版本

R 4.3.1；GSVA 1.48.3；survival 3.5.7；data.table 1.18.2.1；org.Hs.eg.db 3.17.0（别名映射）；AnnotationDbi 1.62.2；MCPcounter 1.2.0；xCell 1.1.0；EPIC 1.1.7；estimate 1.0.13；dendextend 1.17.1；readxl 1.4.3；GEOquery 2.68.0；ggplot2 4.0.2；cowplot 1.1.1 ✅

---

## ⚠️ 清单与关闭条件

| 编号 | 问题 | 怎么关闭 |
|---|---|---|
| A | ~~TCGA protein_coding 过滤~~ | 已关闭（取证 2） |
| B | ~~样本流转~~ | 已关闭（取证 3 + geo_build_20260829_1315.log）。GEO 最后一步只记录「可用于生存分析」的总数，未分别记录时间缺失、时间 ≤ 0、状态缺失各几例；Methods 合并写为「无可用 OS 信息」，不拆分 |
| C | ~~探针注释~~ | 已关闭（取证 2） |
| D | ~~SIG-012 原文权重敏感性分析~~ | 已完成（DEVIATIONS 2026-10-25 D）。SIG-005 原文是否区分上调/下调基因【待核实】，只影响 Limitations 的措辞，不补做分析 |
| E1 | ~~零分布 A 全池兜底~~ | 已关闭：从未触发（DEVIATIONS (1)） |
| E2 | ~~零分布 B 放宽后未检查 ±35%~~ | 已关闭：判定翻转上界两种口径均 0/71（DEVIATIONS (2)） |
| F | ~~QC 只用前 200 个模块~~ | 已关闭：1000 个复核判定不变（DEVIATIONS (3)） |
| G | ~~B 类签名的零分布 B~~ | 已决定按 §5.3 处理（DEVIATIONS 2026-10-25 G）；涉及 pass_B 的图须重排 |
| H | ~~20 个零分布成员独立复算~~ | 已关闭：20/20 一致 |
| I | ~~SIG-007 的训练队列~~ | 已关闭：TCGA 训练、GSE30219 验证；角色重算 6/6、3/13、6/52 |
| — | §5.3 第 2 条（复合指数）：注册表中没有复合指数（取证 3 G0），本条不适用，Methods 不写 | 已关闭 |
| 新 | SAP §5.3 末条「三类结果分别成表，不混排」 | 现有 Fig 与汇总表把 A、B 类放在同一张里（按类别分区）。成稿时 A 类、B 类分表；是否需要拆图，看图的分区是否足够清楚 |

### G：取证 2 之后的判断

**1. 零分布 A 用于 B 类 = §5.3 第 1 条。** §5.3 第 1 条原文（SAP L491–493）：「随机抽取 k 个基因，逐个与原模型对应基因做联合箱匹配；各随机基因在队列内 Z 标准化后套用原模型的同一组权重加权求和」。代码：sample_null_A 按签名基因顺序逐个联合箱匹配（L102–113），run_benchmark L193 把 `sg$weights[gp]` 按同一顺序配给替代基因，score_sets_weighted 逐基因 Z 后加权求和。第 j 个替代基因拿到的就是第 j 个真实基因的权重，和 §5.3 一致。sample_null_weighted（L188–210）与 sample_null_A 的抽样代码只差在箱不够时的退让方式；箱够用时两者在同一 seed 下给出完全相同的集合。**已验证（verify_nullengine G1/G2 均为 0）**：7 个 B 类签名的全部组合中，两个函数在同一 seed 下逐成员给出相同基因集和权重。这部分不是偏离，Methods 照 §5.3 写（DEVIATIONS 2026-10-25 (4)）。

**2. 零分布 B 用于 B 类：SAP 前后不一致，代码选了其中一种读法。**
- SAP 表格（L301–305）和 §5.3 只给 B 类规定了加权复合零分布，没有零分布 B。
- 但 §4.4（L405–409）写的是「对每个 signature，报告其通过零分布 A 与零分布 B 的队列数……A 类与 B 类分开汇总」，字面上包括 B 类。
- 代码对 B 类也跑了零分布 B。模块由种子加相关最高的 k−1 个基因组成，按相关性排序；权重按位置分配，所以权重和基因之间没有对应关系，第 1 个权重总是给种子基因。
- **影响有多大（注册表已确认，取证 3 G0）**：B 类 7 个为 SIG-001/002/004/006/007/008/011，A 类 5 个为 SIG-003/005/009/010/012。pass_B_main 不为 0 的 SIG-004（3）、SIG-006（2）、SIG-007（4）、SIG-008（2）、SIG-011（1），以及 pass_B_sens 的 SIG-001（5），**全部是 B 类**。A 类 5 个签名的 pass_B 全部为 0。所以 pass_B 的全部非零结果都落在 SAP 有歧义的这一块。

**两种处理（你选；都不改动任何已有数字，都写进 DEVIATIONS）：**
- **(1) 按 §4.4 保留。** B 类的 pass_B 照报，Methods 写明：B 类零分布 B 的权重按模块内位置分配，与基因没有对应关系；SAP 在这一点上前后不一致。
- **(2) 按表格和 §5.3 处理。** 主结果中 B 类只报零分布 A（= §5.3），pass_B 只对 A 类报告（5 个签名全部为 0）；B 类的零分布 B 结果放补充材料，注明「SAP 未为 B 类规定此项」。

两种读法下结论的写法不同：(1) 下是「部分 B 类签名在共表达匹配对照下仍强于随机」；(2) 下是「主分析中没有签名通过零分布 B」。
