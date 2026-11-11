
## 2026-10-16 GSE30219 × SIG-001/SIG-002 零分布次数补足
- 发现：verify_ch3.R / 盘点显示这两个组合 n_null = 50（2026-09-03 23:50 --test 运行产物），SAP 规定 1000
- 原因：run_benchmark.R 对已存在的 RDS 跳过，全量运行未覆盖
- 处置：旧 RDS 移至 null_benchmark_superseded_20261016/；以 --force 重跑，seed 规则不变
- 影响：只影响这两个组合的 p_emp_A / p_emp_B 及其下游计数；真实 z 不受影响

## 2026-10-16 覆盖率口径与实际打分基因不一致（实现层面，结果产出后发现）
- 发现：verify_ch3.R 第 5 节打分覆盖审计（scoring_coverage_audit_20261016.csv）
- 现象：gene_coverage() 按 SAP §4.1 把别名回退命中的基因计入覆盖率；run_benchmark.R 打分时用
  gp = intersect(签名基因原名, 候选池)，不使用别名映射，因此回退命中的基因未进入打分
- 涉及 5/71 个可评估组合、8 个基因（均为现行 HGNC 符号 → 平台旧名）：
  GSE72094 × SIG-003：CCN4→WISP1, ADGRA2→GPR124, PLPPR4→LPPR4, NME8→TXNDC3, TENM4→ODZ4（实际 135/141）
  GSE72094 × SIG-006：TMT1A→METTL7A（7/8）；TCGA_LUAD × SIG-006：TMT1A→METTL7A（7/8）
  GSE72094 × SIG-007：CEMIP→KIAA1199（7/8）；GSE72094 × SIG-009：C5AR2→GPR77（26/27）
- 影响：上述组合实际覆盖率 0.875–0.963，均高于 0.70，可评估性不变；全部已报数字本来就按实际打分
  基因计算（verify_ch3 逐位复现），数值不变
- 处置：不重跑（结果产出后不改打分方式）；正文与补充表的覆盖率一律报实际纳入打分的基因数（n_in_pool），
  Methods 写明别名映射只用于判断可评估性、不用于打分；B 类模型在缺失基因时按剩余基因及其原权重计算

## 2026-10-17 GSE30219 分期推导「待修」标注经核实为误判，关闭（阶段 2 任何模型运行之前）
- 来源：SAP v1.1/v1.2 §2.4 及 scoring_protocol.R 待办「85 例中 81 例判为 I 期，与原始 TNM 分布不符，须在 M4 运行前修正」
- 核查：sci2_gse30219_stage_verify.R 从 raw/GSE30219_series_matrix.txt.gz（md5 7d7f2a0056397cde7262a54f83f9f8fe）
  按键名读取 pT/pN/pM，以显式规则独立映射，不调用整备脚本的 tnm2stage()
- 结果：腺癌 85 例中 81 例为 I 期，该数字准确；底座 83 例 TNM 字符串 83/83 一致、分期 83/83 一致（I 79 / II 3 / III 1）
- 误判原因：将腺癌子集与全队列（含 SCC、LCNE 等）的 TNM 分布对比。腺癌 T3 2 例、T4 0 例、N2/N3 0 例；
  非腺癌 T3 29、T4 21、N2/N3 40
- 处置：不修改任何数据；该待办关闭
- 局限：pT 未区分 T2a/T2b，T2N0 10 例统一归 I 期；按 AJCC 第 7/8 版，其中 T2b 者应为 IIA，本数据无法区分

## 2026-10-17 SAP §14 追加：M4 中 stage 的编码方式（阶段 2 任何模型运行之前决定）
- 原文：SAP v1.2 §7 规定 M4 = M1 + 纯度 + 增殖 + stage + age + sex，未规定 stage 的编码
- 问题：GSE30219 腺癌 83 例中非 I 期仅 4 例（II 3 / III 1）；GSE31210 与 GSE50081 无 III/IV 期。
  多级因子在这些队列中估计不稳或不收敛
- 决定：所有队列统一编码为二分类 stage_ge_II（0 = I 期，1 = II/III/IV 期）；A/B/C 亚分并入所属分期组；
  非罗马数字分期组的取值记为缺失。M4 按完整病例分析，各队列 stage 缺失比例照 SAP §7 报告
- 代价：不区分 II 期与 III/IV 期；在 TCGA 等分期范围较宽的队列中，残余混杂可能多于多级编码
- 实现：scripts/sci2_m4_stage_coding.R → data_derived/m4_stage_binary_<日期>.csv；M4 只读该文件
- 性质：方案细化，发生于结果产出之前

## 2026-10-17 SAP §14 追加：TCGA 分期缺失的处理（阶段 2 任何模型运行之前决定）
- 现象：TCGA_LUAD 443 例中 64 例（14.4%）分期缺失；核查确认源文件 TCGA_LUAD_SE.rds 的 GDC 字段
  ajcc_pathologic_stage 对这 64 人本即为 NA，底座整备无漏取
- 第二来源：SE colData 的 paper_Tumor.stage（TCGA 肺腺癌标志性论文的临床注释）。两源均有值的 158 人中，
  分期组一致 153 人，二分类（I vs ≥II）一致 153 人（96.8%，Cohen kappa 0.935）；5 例不一致全部落在 I/II 交界处
  （GDC I → 论文 II 3 例，GDC II → 论文 I 2 例）
- 决定：主分析 M4 只用 GDC 分期，完整病例分析，如实报告缺失 14.4%；
  事前登记一个敏感性分析：GDC 缺失者用 paper_Tumor.stage 补充（可补 29 人：I 9 / II 10 / III 7 / IV 3；
  仍缺 35 人，7.9%），GDC 有值者一律不改。两版结果并列报告，结论以主分析为准
- 实现：scripts/sci2_m4_stage_tcga_sens.R → data_derived/m4_stage_binary_TCGA_sens_<日期>.csv
- 性质：方案细化，发生于结果产出之前；不修改底座

## 2026-10-17 SAP §14 追加：TCGA 纯度来源落实、分期缺失的结局相关性、衰减幅度的同样本计算（均在 M2/M4 运行之前）
### (1) 纯度来源（按 SAP §7.2 原文执行，无偏离）
- Aran D, Sirota M, Butte AJ. Nat Commun 2015;6:8971 Supplementary Data 1 的 ABSOLUTE 列；
  文件 raw/purity/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx（md5 c459e6a965789b96860fc77bd346c681）
- 原表缺失写作字符串 NaN，读入时按缺失处理；LUAD 原发灶有 ABSOLUTE 者 354 人，同一病人无多管取值冲突
- 对底座 443 人覆盖 301 人（67.9%）；有/无纯度两组事件率 0.415 / 0.401
- 与 SE 的 paper_Purity.ABSOLUTE.calls（逗号小数，转换后）在共同 191 人中 Spearman 1.000，189 人完全相同
- 实现：scripts/sci2_purity_tcga.R → data_derived/m2_purity_TCGA_20261017.csv
### (2) TCGA 分期缺失与结局强相关（数据特征，经核实非编码错误）
- GDC 分期缺失 64 人中 59 人死亡（92.2%），有分期的 379 人中 123 人死亡（32.5%）；
  底座 OS_event / OS_time 与 SE 的 vital_status / days_to_* 逐人一致（64/64、379/379）
- 后果：M4 主分析（完整病例）剔除一个事件率极高的亚组，存在选择偏倚；敏感性版本补回其中 29 人（27 事件）。
  两版并列报告，讨论中须明示此缺失机制
### (3) 衰减幅度在同一完整病例样本上计算
- SAP §7 仅规定 M4 采用完整病例分析，未规定 M1 的样本。为避免把样本变化混入衰减幅度：
  Attenuation 一律以「在 M4 完整病例样本上重新拟合的 M1」为分母；M2 同理（M1 在 M2 样本上重新拟合）。
  全体样本的 M1 另行单列报告
- 性质：方案细化，发生于结果产出之前

## 2026-10-17 SAP §14 追加：meta-PCNA 重叠数据更正、增殖指数的计算细则（均在 M3 运行之前）
### (1) SAP v1.2 §7.2 第 5 条「实测」数字有误，更正
- 原文：SIG-001 重叠 2/4（50%）、SIG-006 重叠 3/8（37.5%），据此两者标注为结构性重叠
- 复算（sci2_prolif_diag.R，按符号与 Entrez ID 双重核对 confounders/metaPCNA_table.csv）：
  SIG-001 1/4（25%）、SIG-006 1/8（12.5%），重叠基因均仅为 CCNB1；EXO1（Entrez 9156）与 CDC25C（Entrez 995）
  不在 meta-PCNA 131 基因中。原数字疑为按生物学常识计入「增殖相关基因」，而非按 meta-PCNA 名单
- 后果：按 §7.2 第 5 条规则（>30%），无签名构成结构性重叠；SIG-001 与 SIG-006 的 M3 结果计入主汇总
- 连带更正：SAP 称 15 个基因改名，metaPCNA_table.csv 实为 14 个；「SIG-001 一半基因属 meta-PCNA」的表述不成立，不得使用
### (2) 增殖指数计算细则
- 定义按 SAP §7.2：131 基因 log2 表达（load_cohort 协议）的逐样本中位数
- 基因匹配：先按现行 HGNC 符号精确匹配，未命中者用 resolve_aliases() 唯一别名回退。与签名打分（只用精确匹配）口径不同，
  理由：增殖指数是混杂因子、非被测对象，覆盖越全，校正越准；Methods 中写明
- 只对 6 个有生存数据的队列计算；CPTAC 不进入 M3/M4，其矩阵中 11 个 meta-PCNA 基因含 134 个 NA，不做处理
- 实现：scripts/sci2_prolif_index.R → data_derived/m3_prolif_metaPCNA_<日期>.csv
- 性质：数据更正与方案细化，发生于 M3 运行之前

## 2026-10-17 SAP §14 追加：衰减幅度的报告条件（M1–M4 运行之前决定）
- 问题：Attenuation = [log HR_M1 − log HR_Mk] / log HR_M1，M1 效应接近零时分母趋零，数值不稳且无解释意义
- 决定：仅当同一分析样本上的 M1 满足 |z| ≥ 1.96 时报告衰减幅度并纳入汇总；其余组合只报告各模型的 HR 与 95% CI，
  衰减幅度标为「无可衰减信号」。原始计算值保留在输出文件中（attenuation 列），汇总只用 attenuation_reported 列
- 其他实现细节：分数与连续协变量（纯度、增殖、年龄）在各模型自身的分析样本内 Z 标准化；sex 作因子；
  stage 为二分类 stage_ge_II；GEO 队列 M2/M4 按 SAP §7.2 不计算；全样本 M1 必须与阶段 1 的 z_real 逐个一致（闸门）
- 实现：scripts/sci2_stage2_confounding.R
- 性质：方案细化，发生于结果产出之前

## 2026-10-17 SAP §14 追加：衰减幅度的报告条件（M1–M4 运行之前决定）
- 问题：Attenuation = [log HR_M1 − log HR_Mk] / log HR_M1，M1 效应接近零时分母趋零，数值不稳且无解释意义
- 决定：仅当同一分析样本上的 M1 满足 |z| ≥ 1.96 时报告衰减幅度并纳入汇总；其余组合只报告各模型的 HR 与 95% CI，
  衰减幅度标为「无可衰减信号」。原始计算值保留在输出文件中（attenuation 列），汇总只用 attenuation_reported 列
- 其他实现细节：分数与连续协变量（纯度、增殖、年龄）在各模型自身的分析样本内 Z 标准化；sex 作因子；
  stage 为二分类 stage_ge_II；GEO 队列 M2/M4 按 SAP §7.2 不计算；全样本 M1 必须与阶段 1 的 z_real 逐个一致（闸门）
- 实现：scripts/sci2_stage2_confounding.R
- 性质：方案细化，发生于结果产出之前

## 2026-10-18 SAP §14 追加：阶段 3（§8）的读数与一致性指标（阶段 3 运行之前决定）
- §8 第 2 条「通用基质 / 免疫解卷积读数」落实为：MCP-counter 全部 10 类；xCell 的 Fibroblasts / ImmuneScore /
  StromaScore / MicroenvironmentScore；ESTIMATE ImmuneSignature（estimate 包 SI_geneset，与签名同一 ssGSEA 路径）；
  meta-PCNA 增殖指数。以上 6 队列均计算
- EPIC（CAFs / Endothelial）仅在 TCGA 计算：EPIC 要求线性 TPM 尺度，芯片数据不满足其设计前提；TCGA 由 log2(TPM+1) 还原。
  ABSOLUTE 纯度仅 TCGA。两者不进入跨队列一致性比较
- §8 第 3 条「聚类结构在各队列间的一致性」落实为：两两队列在共有读数上分别以 1−|rho|、ward.D2 重新聚类，
  主指标为 cophenetic 相关（dendextend::cor_cophenetic），辅指标为相关矩阵上三角的 Spearman；报中位与范围。不设类数 k
- §8 第 4 条：全部读数在各队列做单变量 Cox（per SD），与阶段 1/2 同一有效样本定义；外部读数作参照，单独标注
- 实现：scripts/sci2_stage3_redundancy.R
- 性质：方案细化，发生于阶段 3 结果产出之前

## 2026-10-25 事后探索性分析：成纤维细胞读数校正增殖后的预后关联（运行之前登记）
- 动机：阶段 3 发现 MCP_Fibroblasts 与 xCell_Fibroblasts 预后方向相反（z>0 5/6；z<0 6/6），
  且 xCell_Fibroblasts 与增殖指数 6/6 负相关（−0.54 至 −0.32），MCP_Fibroblasts 与增殖近乎不相关（0.03–0.37）
- 性质：SAP 之外的事后探索性分析，发生于阶段 3 结果产出之后；不构成预注册检验，论文中须如此标注
- 设计（事前固定，不再调整）：仅限 MCP_Fibroblasts 与 xCell_Fibroblasts 两个读数；6 个生存队列各拟合
  Surv ~ scale(读数) 与 Surv ~ scale(读数) + scale(meta-PCNA 增殖指数)；有效样本定义同阶段 1–3；
  报告两个模型的 HR、95% CI、z，以及各模型方向一致的队列数与 |z|≥1.96 的队列数；不做多重比较校正，不计算衰减幅度
- 无论结果如何均如实报告；不据此修改阶段 1–3 的任何结论
- 实现：scripts/sci2_stage3_explore_fibro.R → data_derived/stage3_explore_fibro_prolif_<日期>.csv

## 2026-10-25 零分布引擎复核（verify_nullengine；均为实现层面，阶段 1 结果产出之后发现）

依据：scripts/sci2_verify_nullengine.R（md5 62d0a5f8b3175b9024432866a01b1061），日志 logs/verify_nullengine_20261025_1456.log，
输出 data_derived/verify_nullengine_combo_20261025.csv、verify_nullengine_members_20261025.csv。
按 RDS 中记录的 seed 重新生成全部 71 个组合的零分布；重现性闸门 H1、R1–R7 共 8 条全部通过
（含 SAP §12 要求的 20 个零分布成员从表达矩阵独立复算 Cox z，最大差 0）。以下计数针对的即为阶段 1 实际使用的零分布。

### (1) 零分布 A 的第三层退让（SAP §5.1 未规定）
null_engine.R L111：同箱与平均表达 ±1、±2 箱均无可用基因时，从全池（除签名基因）任取，且不计入 n_fallback_A。
实测：71 个组合中相邻箱退让 0 次、全池兜底 0 次。该分支从未触发，对结果无影响。Methods 只写 SAP 的两层规则，并注明实测退让为 0。

### (2) 零分布 B「放宽至 ±35%」未被检查（与 SAP §5.2 第 3 条字面不符）
null_engine.R L151/L167：连续 50 次未落在 ±20% 内时，第 51 次把容差设为 ±35%，但无论该次模块是否落在 ±35% 内都直接采用。
实测：71,000 个零分布 B 成员中 142 个（0.20%）为第 51 次直接采用，其中 118 个（0.17%）的平均表达或平均方差偏离真实签名超过 ±35%；
分布于 20 个组合，集中在 TCGA × SIG-012（64 个，其中 57 个超出 ±35%）；最大偏离：平均表达 2.33 倍（TCGA × SIG-004）、平均方差 3.43 倍（GSE30219 × SIG-004）。
处置（看到 p 值影响之前决定）：不重抽、不改动任何已有 p_emp_B。用 scripts/sci2_nullB_relax_bound.R 计算判定翻转的上界：
不知道受影响成员各自的 z，取所有可能情形，若剔除这些成员后 p_emp_B 的整个可能区间仍与原值位于 0.05 同一侧，则判定不可能改变。
分别按「超出 ±35%」和「第 51 次采用的全部成员」两种口径计算。
- 若两种口径下均无组合可能翻转：Methods 如实报告上述实现与计数，结论不变。
- 若有组合可能翻转：列出这些组合；另行登记成员级敏感性分析（剔除超出 ±35% 的成员后重算 p_emp_B）后再运行，主结果仍为原值。

### (3) 共表达质量验证只用前 200 个模块（SAP §5.2 写「各模块」）
null_engine.R L175。实测：用全部 1000 个模块重算，71 个组合的 QC 判定无一改变。Methods 写明实现为前 200 个，并报告用 1000 个复核判定一致。

### (4) B 类签名的零分布 A 与 SAP §5.3 第 1 条等价
§5.3 专用函数 sample_null_weighted 从未被调用；B 类调用的是 sample_null_A，并把权重按签名基因顺序配给对应替代基因。
实测：7 个 B 类签名的全部组合中，sample_null_weighted 与 sample_null_A 在同一 seed 下给出的基因集逐成员相同，权重相同。
因此 B 类的零分布 A 即 §5.3 第 1 条规定的加权复合零分布，不构成偏离。B 类的零分布 B 如何报告另行决定（见对照表 G）。
- 结果（2026-10-25，sci2_nullB_relax_bound.R → data_derived/nullB_relax_bound_20261025.csv）：两种口径下可能翻转 p_emp_B < 0.05 判定的组合均为 0/71。
  离阈值最近者 GSE50081 × SIG-004（p = 0.044，剔除后区间上界 0.0440）。按 (2) 的事先规定，不做成员级敏感性分析，Methods 如实报告实现与计数。

## 2026-10-25 B 类签名零分布 B 的报告方式（对照表 G；不改动任何已有数字）

事实：SAP v1.2 表格（L301–305）与 §5.3 只为 B 类签名规定了加权复合零分布；§4.4（L405–409）「对每个 signature 报告通过零分布 A 与零分布 B 的队列数」字面上包括 B 类，两处不一致。
代码对 B 类也运行了零分布 B，权重按模块内位置分配（第 1 个权重恒给种子基因），权重与基因无对应关系。
阶段 1 中 pass_B 的全部非零结果（SIG-004/006/007/008/011 的 pass_B_main，SIG-001 的 pass_B_sens）均来自 B 类签名；A 类 5 个签名的 pass_B 全部为 0。
决定：按 SAP 表格与 §5.3 执行。
- 主结果：零分布 B（pass_B_main / pass_B_sens）只对 A 类签名报告；B 类签名只报零分布 A（经 verify_nullengine G1/G2 验证即 §5.3 第 1 条的加权复合零分布）。
- B 类签名的零分布 B 结果移至补充材料，注明「SAP 未为 B 类规定此项；权重按模块位置分配」，不作解释。
- 按 §5.3 末条，A 类、B 类结果分表报告。涉及 pass_B 的图须按此重排（数字不变）。
- Limitations 中写明 SAP 这一处前后不一致，以及上述处理在阶段 1 结果产出之后决定。

## 2026-10-25 SIG-012 原文权重敏感性分析（对照表 D；运行之前登记）

事实：注册时（update_registry_v2.R L83）写明「原文权重另存 SIG-012_weights_original.csv 作敏感性分析」，此前未执行。
文件 genes/SIG-012_weights_original.csv（md5 c91c10cce0c0cd73f7d5808292f79a9b，11 个基因，与 SIG-012_genes.txt 一致）。
主分析不变：SIG-012 按 ≥10 基因规则归 A 类，等权 ssGSEA。
敏感性分析做法（scripts/sci2_sens_sig012_origweights.R）：
- 打分基因 = SIG-012_genes.txt ∩ 候选池（精确符号，同主分析）；各基因队列内 Z 后按原文权重加权求和，缺失基因不重新分配权重
- 单变量 Cox（per SD），Wald z
- 零分布：§5.3 第 1 条加权复合零分布（逐基因联合箱匹配 + 同一组权重），1000 次；seed = make_seed(队列, "SIG-012", "W")
- 不做零分布 B（与上条 G 的处理一致）
- 报告 z_w、HR_w、p_emp_W、零分布 |z| 95% 分位、与主分析分类的一致性、两种分数的 Spearman 相关
- 仅作敏感性分析报告；不论结果如何，主分析与阶段 1 汇总数字不变
- 结果（2026-10-25，data_derived/sens_SIG012_origweights_20261025.csv；闸门：等权 z 与阶段 1 一致）：
  原文权重分数与等权 ssGSEA 分数在 6 个队列均负相关（Spearman −0.52 至 −0.06），z 方向全部相反（等权 −4.17 至 −2.28；原文权重 2.13 至 7.19）。
  原因：11 个原文权重中 7 个为负；等权 ssGSEA 不区分权重符号。两种打分读出的方向相反，只比较 |z|。
  分类一致 4/6；GSE31210、GSE30219 由「不特异」变为「强于随机」（p_emp_W 0.002、0.010）；按原文权重「强于随机」4/6，主分析 2/6。
  主分析与阶段 1 汇总数字不变。Limitations 须写明：A 类等权规则用于正负权重混合的原文模型时，所测不是原文模型本身。
  解读前须核对 SIG-012 的训练/验证队列（原文权重在训练队列上拟合）。

## 2026-10-25 SIG-007 出处更正与队列角色重算（I 项；锁定表不改）

依据：原文全文 Li W, Gui T, Hou Y. Development of a RiskScore model related to cancer-associated fibroblasts for assessing patient prognosis
and drug sensitivity in lung adenocarcinoma. Lett Drug Des Discov. doi:10.1016/j.lddd.2026.100333（PII S1570180826000631；页眉印「22 (2025) 100333」，出版年【待核实】）。
- §2.1：TCGA-LUAD（500 肿瘤 + 59 对照）为训练集；GSE30219（83 例）为验证集；另在 CPTAC 验证（Figure S1；底座 A 的 CPTAC 无生存数据，不影响角色）。
- §3.4 公式与 genes/SIG-007_weights.csv 8 个系数逐一相符（2026-10-25 核对）。
锁定表（md5 4ad3218931d95dfb6ee67a7ba07b9a75）中三处记录有误，锁定表不改，以本条为准：
  (1) 名称「Wang(LDDD)」：第一作者为 Wan Li；(2) note「Bentham」：期刊主页为 KeAi Publishing，平台 ScienceDirect；
  (3) doi「待核实」、overlap_cohorts「待核实」：分别为上述 DOI 与「TCGA, GSE30219」。
影响：analyze_cohort_role.R 把「待核实」按「原文未用过」处理，SIG-007 的 6 个组合全部计入「完全独立队列」，其中 TCGA（训练）与 GSE30219（原文验证）各通过 1 次。
处置：新建 scripts/analyze_cohort_role_20261025.R（复制原脚本，仅改 SIG-007 的 train_cohort 为 TCGA、原文用过的队列为 TCGA 与 GSE30219），
输出新日期文件，旧文件保留。阶段 1 的 z、p_emp 与分类不变；按角色的通过率与过拟合梯度以新文件为准。
- 结果（2026-10-25，results_by_cohort_role_20261025.csv、overfit_gradient_20261025.csv；train_cohort_annotation_20261025.csv）：
  通过率 开发 6/6、原文验证 3/13、完全独立 6/52（旧 5/5、2/12、8/54）；仅预后模型类：6/6、3/11、6/24。
  过拟合比 3.17（SIG-007）/ 3.06 / 2.67 / 1.85 / 1.35 / 1.26。20261016 版两个文件作废，不删除。

## 2026-10-25 SAP §3.4 补执行：SCANI myCAF（运行之前登记；发生于阶段 1–3 结果产出之后）

SAP §3.4：「SCANI 的 myCAF ssGSEA 组分作为 A 类的一个基因集纳入，适用完全相同的零分布；在结果表中不置于首位，不作特殊标注。」
此前未执行，且无排除记录（见本日「SAP §6、§3.3、§3.4」条 (3)）。作者决定补执行（2026-10-25）。
基因集：~/caf_spatial_project/results/01_scRNA/CAF_subtype_signatures.rds 的 $myCAF（md5 17d3458f1a9f55c0e90316aed443cd3e；
  即 SCI-1 step6e_final_model.R L14、L28 所用）。不做符号改写（与锁定表签名一致：打分精确匹配，别名只用于覆盖率）。
做法：scripts/sci2_scani_supp.R，与 run_benchmark.R 的 A 类相同（覆盖率 ≥ 70%；ssGSEA ssgsea.norm = FALSE；单变量 Cox per SD；
  零分布 A、B 各 1000 次，seed = make_seed(队列, "SCANI_myCAF", 类型)；p_emp 双侧；共表达 QC）。结果另存，已锁定的 71 个组合不改动。
队列角色：SCI-1 的开发队列为 TCGA-LUAD（II–IV 期 166 例），外部验证为 GSE72094；其余 4 个 GEO 队列为完全独立队列。
报告：不论结果如何均报告；纳入后 A 类每队列 BH 家族由 5 变为 6（SAP §6），主结果按含 SCANI 的家族计算，不含 SCANI 的结果并列；
  论文中注明该签名系补执行、晚于其他结果，并披露作者与 SCANI（SCI-1）的关系。

## 2026-10-25 SAP §6、§3.3、§3.4 未执行或未遵守（Methods–代码对照第二轮发现；阶段 1 结果产出之后）

依据：SAP_v1.2.txt §3、§6 原文（2026-10-25 取证）；diagnose_results.R、run_benchmark.R 中无 p.adjust。
(1) §6 多重比较：SAP 规定「跨 signature（同一队列、同一类别内）做 BH FDR，q < 0.05」；跨队列与零分布 A/B 之间不校正。
    实现中所有通过判定均为未校正 p_emp < 0.05，此前汇总（31/25/15、pass_A、按角色通过率）均为未校正结果。
    处置（看结果之前决定）：scripts/sci2_sap6_bh.R 按 §6 补算 BH；按 SAP，BH 校正后的结果为主结果，未校正结果并列报告。
    零分布 B 的 BH 只对 A 类（见 2026-10-25 G）。不改动任何 z 与 p_emp。
(2) §3.3「两类的结果不得合并计算通过队列数的总体比例」：此前汇总与稿件 v0 将 A、B 类合并（71 个组合、按角色 6/6 等）。
    处置：所有比例按类别分开报告；稿件 v1 据此改写。
(3) §3.4 规定 SCANI 的 myCAF ssGSEA 组分作为 A 类基因集纳入；锁定表（12 个）与 excluded_notes.csv（6 个）中均无 SCANI，
    无任何排除记录。处置：作者决定补执行（见本日「SAP §3.4 补执行：SCANI myCAF」条）。
(4) 候选签名的来源：registry 中无检索策略记录（数据库、检索式、日期）。18 个候选、纳入 12、排除 6 的标准见 SAP §3.1–3.2。
    论文中写为非系统性文献检索，并列入 Limitations。

## 2026-10-25 SAP 的性质与论文措辞（作者决定）
SAP v1.2 封面写「预注册」，定稿日期与签署栏空白；未在任何公共平台登记，也未以邮件等方式取得外部时间戳（导师授权作者自行决定）。
论文中不得使用「pre-registered」「time-stamped」；统一写为「analysis plan written before the benchmark was run; not publicly registered」，
可核实的依据为锁定表 MD5（4ad3218931d95dfb6ee67a7ba07b9a75，2026-09-03）与本偏差记录。SAP 全文作为补充文件提交。

## 2026-10-25 SCANI 并入 A 类 BH；事后检查：方向、HR 版过拟合比、通过率 CI（运行之前登记）
SCANI 补跑完成（scani_supp_summary_20261025.csv；6 队列均可评估）。scripts/sci2_classA_scani_and_checks.R，只读已有结果：
(1) SAP §3.4 + §6：SCANI_myCAF 并入 A 类，每队列 BH 家族 = 6，q_A、q_B 重算；主结果 = 含 SCANI，不含 SCANI 的结果（sap6_bh_combo_20261025.csv）并列。
    SCANI 角色：TCGA = 开发（SCI-1 开发样本为 II–IV 期子集）；GSE72094 = SCI-1 外部验证；其余 4 个 = 完全独立。
    SCANI 零分布 B「第 51 次直接采用」成员按 relaxed 口径做最坏情形上界。
(2) 事后、描述性（SAP 未规定）：
    (a) 方向：B 类与 SCANI 预期 z > 0；其余 A 类原文未给等权分数方向，记为未定义、不判定。
    (b) 过拟合比的 HR 版（|log HR| per SD），与 z 版并列；用于回应 TCGA 检验力/平台混杂，不论结果均报告。
    (c) 「类别 × 队列角色」通过数的精确二项 95% CI。

## 2026-10-25 SIG-009/SIG-010 出处作者更正（锁定表不改；结果产出之后发现）
- 发现：锁定表 name 记为「Ruhle 2023」，但所记出处 OncoImmunology 2023;12(1):2274130 Table S1、DOI 10.1080/2162402X.2023.2274130 均正确，该文第一作者为 Koeck S（作者列表中无 Rühle）
- 处置：锁定表不改；稿件 Table 2、参考文献与补充表 S2 显示名改为「Koeck 2023」；补充表 S1 第 28 条
- 影响：基因列表与全部结果不变

## 2026-11-09 审稿准备：事后描述性汇总（结果后，运行前登记）
均为对已有结果文件的汇总，不重跑零分布、不改变任何主结果或判定：
1. B 类 8 个通过组合在零分布 B（Table S3，SAP 未规定）下的结果汇总。
2. 各队列 × 类别零分布 A 的 p95(|z|) 与事件数；据此近似估算最小可检测 HR per SD（se ≈ 1/√D）。
3. 各队列零分布 A 成员中 |z_null| ≥ 1.96 的比例（按类别，中位数与范围）。
4. 敏感性：对全部 77 个组合的 p_emp 做一次全局 BH；报告通过数是否变化。
5. 各签名基因与 131 个 meta-PCNA 基因的重叠数。
6. TCGA 中各签名与增殖指数、MCP-counter / EPIC 成纤维细胞读数的 Spearman ρ（核对 Discussion 解释段）。
7. Table S6 中 xCell / MCP-counter 单变量 |z| ≥ 1.96 的队列数；增殖指数 6 队列 z 值的出处。

## 2026-11-09 补充表新增 Table S8；Table S6 补增殖指数（结果后，运行前登记）
- 新增 Table S8：由 data_derived/review_summaries_20261109/ 的 02_nullA_summary_by_cohort.csv 与 04_global_bh.csv 生成（随机成员 |z|≥1.96 比例、最小可检测 HR、三种 BH 族定义的通过数与逐组合 q 值）。
- Table S6 增加增殖指数（Prolif_metaPCNA）在 6 个队列的单变量 Cox 结果，来源 stage3_cox_20261024.csv。
- 不改变 S1–S7 其他内容和任何主结果；补充表 xlsx 以新日期输出，旧版移入 superseded。

## 2026-11-09 补充表新增 Table S8；Table S6 补增殖指数（结果后，运行前登记）
- 新增 Table S8：由 data_derived/review_summaries_20261109/ 的 02_nullA_summary_by_cohort.csv 与 04_global_bh.csv 生成（随机成员 |z|≥1.96 比例、最小可检测 HR、三种 BH 族定义的通过数与逐组合 q 值）。
- Table S6 增加增殖指数（Prolif_metaPCNA）在 6 个队列的单变量 Cox 结果，来源 stage3_cox_20261024.csv。
- 不改变 S1–S7 其他内容和任何主结果；补充表 xlsx 以新日期输出，旧版移入 superseded。

## 2026-11-09 补充表新增 Table S8；Table S6 补增殖指数（结果后，运行前登记）
- 新增 Table S8：由 data_derived/review_summaries_20261109/ 的 02_nullA_summary_by_cohort.csv 与 04_global_bh.csv 生成（随机成员 |z|≥1.96 比例、最小可检测 HR、三种 BH 族定义的通过数与逐组合 q 值）。
- Table S6 增加增殖指数（Prolif_metaPCNA）在 6 个队列的单变量 Cox 结果，来源 stage3_cox_20261024.csv。
- 不改变 S1–S7 其他内容和任何主结果；补充表 xlsx 以新日期输出，旧版移入 superseded。

## 2026-11-10 签名来源途径重建；Table S2 增 found_via 列（结果后，运行前登记）
- 依据：作者核对建注册表时的对话记录（无数据库、检索式、检索日期记录；10-25 已定为「非系统性检索」）。
- 18 个候选的来源：起点文献（SCI-1 已引用）4 个——Ren 2023、Si 2024、Cords 2024、Chen 2025；起点基质读数 4 个——ESTIMATE、MCP-counter、EPIC、xCell；核对过程中追加 7 个——Zhang Y 2023、Navab 2011、Li 2025(SIG-007)、Wang S 2023、Zhang H 2023、Huang 2023、Peng 2023；未记录 3 个——SIG-009、SIG-010（Koeck 2023）、Kong 2024。
- 起点清单中的 Liu 2025 与 Bagaev TME 分型未作为候选评估，原因未记录；在 S2 中列出、状态为 not assessed。
- 稿件 2.3 节改为如实描述；Bagaev 2021 加入参考文献。不改变任何结果。

## 2026-11-10 比例风险假定检验（结果后，运行前登记）
- 对主分析的 77 个签名–队列组合（含 SCANI），用与主分析相同的标准化分数和单变量 Cox 模型（OS），以 survival::cox.zph（KM 变换）检验比例风险假定。
- 报告：p < 0.05 的组合数与清单；对违背者另报按随访时间分段（≤3 年 / >3 年）的 HR 作描述。
- 描述性敏感性分析，不改变任何判定；结果写入 Table S8（新增 S8c）与 Limitations。

## 2026-11-10 SAP 中未执行的两项；File S1 编制（结果后记录）
- SAP §4.1「跨队列共同基因池（11,311 个）作敏感性分析」与 §9「循环特征选择佐证小节」均未执行，原因未记录。
- File S1 = 英文说明页（性质、版本经过、中期数字、勘误、未按计划执行的条款）+ 英文译文 + 中文原文原样（Paper2_SAP_v1.2.docx md5 9eee91d86a25da1facb213e8704983e7，未改动）。
- S1 第 11 条补入第 79 行的连带更正（SAP 称 15 个基因改名，实为 14 个）。

## 2026-11-11 · SIG-007 publication year in Table S2 (bibliographic correction; no analysis)
- Finding: the locked registry and Table S2 give the SIG-007 year as 2026 (taken from the DOI/PII). The publisher's page says Lett Drug Des Discov Vol. 22, Issue 12, December 2025, article 100333; manuscript ref 9 already reads 2025, 22, 100333.
- Action: the locked registry is not edited. In Table S2, the SIG-007 `name` becomes "Li (LDDD) 2025 …" (keeping the original registry label in brackets) and `source` becomes "Lett Drug Des Discov 2025;22(12):100333". The DOI is unchanged. Add Table S1 No. 34.
- Effect on results: none.
- Script: scripts/sci2_supp_S2year_20261111.R; output data_derived/supp_tables_20261111/. Input 20261110c is not overwritten.
