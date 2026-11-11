#!/usr/bin/env Rscript
# ==============================================================================
# analyze_cohort_role.R —— 按「开发队列 vs 独立队列」拆分结果
#
# 动机：统计筛选类 signature 的通过率 37%，但 15 次通过中有 6 次发生在 TCGA，
#   而这 6 个 signature 全部在 TCGA 上开发。若不拆分，会把「过拟合到开发队列」
#   误读为「该 signature 特异」。
#
# ★ 本分析的性质
#   SAP §10 路径 (b) 规定「检验基因数、类别、亚型来源、推导方式与通过率的关联」，
#   未明文包含「开发队列 vs 独立队列」这一维度。
#   因此本分析属【SAP 之外的事后探索性分析】，须在 §14 记录并在论文中如此标注。
#
#   但需注意：各 signature 用过哪些队列，是在【清单锁定前】就登记在
#   registry 的 overlap_cohorts 字段中的（锁定文件 md5 可证）。
#   即：数据是事前记录的，分析是事后想到的。这两点都要如实说明。
#
# 只读，不修改任何已有结果。
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

CH    <- path.expand("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark")
REG   <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")
STAMP <- format(Sys.Date(), "%Y%m%d")
OUT   <- file.path("/tmp", format(Sys.time(), "cohort_role_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 74), "\n", t, "\n", strrep("=", 74), "\n", sep = "")
sh <- function(t) cat("\n  --- ", t, "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-32s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("按队列角色拆分结果\n时间: ", format(Sys.time()), "\n", sep = "")

# ============================================================ 0. 训练队列人工标注
hr("0. 各 signature 的训练队列（人工标注）")
cat("  来源：各原文 Methods 的明文表述。overlap_cohorts 只记录『用过』，\n")
cat("  不区分训练与验证，故训练队列须单独标注。\n\n")

TRAIN <- data.table(
  sig_id = sprintf("SIG-%03d", 1:12),
  train_cohort = c(
    "TCGA",        # SIG-001 Ren 2023：TCGA 训练
    NA,            # SIG-002 MCP-counter：非预后模型，无训练队列
    NA,            # SIG-003 ESTIMATE：非预后模型
    "TCGA",        # SIG-004 Zhang Y：TCGA 内部拆分训练/测试
    NA,            # SIG-005 Navab 2011：原代细胞培养对比，非队列训练
    "TCGA",        # SIG-006 Si 2024：TCGA 训练
    "待核实",      # SIG-007 LDDD 2026：全文获取受限，训练队列未确认
    "GSE41271",    # SIG-008 Wang S 2023：★ GSE41271 训练，TCGA 是验证集
    NA,            # SIG-009 Rühle CAF MARKER：综述汇总，非预后模型
    NA,            # SIG-010 Rühle EPIC CAF：解卷积衍生，非预后模型
    "TCGA",        # SIG-011 Zhang H 2023：TCGA 训练
    "TCGA"))       # SIG-012 Huang 2023：TCGA 训练
print(TRAIN, row.names = FALSE)
cat("\n  ★ SIG-008 的训练队列是 GSE41271（不在底座 A 中），TCGA 对它是【验证集】。\n")
cat("     这一点很重要：若 SIG-008 在 TCGA 上同样表现突出，其性质与其余五个不同。\n")
cat("  ★ SIG-007 的训练队列待核实，本分析中按『未知』处理，单独列出。\n")

# ============================================================ 1. 载入
hr("1. 载入结果与登记表")
fe <- sort(list.files(file.path(CH, "data_derived"),
                      pattern = "^results_extended_.*\\.csv$", full.names = TRUE),
           decreasing = TRUE)
if (!length(fe)) stop("未找到 results_extended_*.csv，请先跑 diagnose_results.R")
E <- fread(fe[1]); kv("扩展汇总表", basename(fe[1]))

fr <- sort(list.files(REG, pattern = "^signature_registry_locked_.*\\.csv$",
                      full.names = TRUE), decreasing = TRUE)
L <- fread(fr[1]); kv("锁定表", basename(fr[1]))
kv("锁定表 md5", unname(tools::md5sum(fr[1])))

if (!"overlap_cohorts" %in% names(L)) stop("锁定表缺 overlap_cohorts 字段")
E <- E[evaluable == TRUE]
kv("可评估组合", nrow(E))

# ============================================================ 2. 判定队列角色
hr("2. 判定每个组合的队列角色")

# 队列名归一：底座 A 的命名 vs 原文表述
norm_co <- function(x) {
  x <- toupper(gsub("[^A-Za-z0-9]", "", x))
  x <- sub("^TCGALUAD$", "TCGA", x)
  x
}
used_map <- setNames(lapply(L$overlap_cohorts, function(s) {
  if (is.na(s) || s %in% c("-", "待核实")) return(character(0))
  norm_co(strsplit(s, ",")[[1]])
}), L$sig_id)

E <- merge(E, TRAIN, by = "sig_id", all.x = TRUE)
E[, co_norm := norm_co(cohort)]
E[, used_in_original := mapply(function(s, c0) c0 %in% used_map[[s]], sig_id, co_norm)]
E[, is_train := !is.na(train_cohort) & train_cohort != "待核实" &
                norm_co(train_cohort) == co_norm]
E[, role := fifelse(is_train, "开发队列",
             fifelse(used_in_original, "原文验证队列", "完全独立队列"))]
# 训练队列未知者单独标记
E[!is.na(train_cohort) & train_cohort == "待核实" & used_in_original,
  role := "原文用过（训练/验证未分）"]

sh("角色分布")
print(E[, .N, by = role][order(-N)], row.names = FALSE)

sh("逐 signature 的队列角色")
print(dcast(E, sig_id + class ~ role, fun.aggregate = length, value.var = "cohort"),
      row.names = FALSE)

# ============================================================ 3. 按角色的通过率
hr("3. 通过率：开发队列 vs 独立队列")
agg <- E[, .(n = .N,
             pass_A = sum(p_emp_A < 0.05, na.rm = TRUE),
             rate_A = sprintf("%.0f%%", 100 * mean(p_emp_A < 0.05, na.rm = TRUE)),
             median_abs_z = round(median(abs(z_real), na.rm = TRUE), 2),
             median_zA_p95 = round(median(zA_p95, na.rm = TRUE), 2)),
         by = role][order(-pass_A)]
print(agg, row.names = FALSE)

sh("仅限有预后模型的 signature（排除非模型类：SIG-002/003/005/009/010）")
MODEL <- setdiff(sprintf("SIG-%03d", 1:12), c("SIG-002","SIG-003","SIG-005","SIG-009","SIG-010"))
M <- E[sig_id %in% MODEL]
kv("纳入 signature", MODEL)
agg2 <- M[, .(n = .N,
              pass_A = sum(p_emp_A < 0.05, na.rm = TRUE),
              rate_A = sprintf("%.0f%%", 100 * mean(p_emp_A < 0.05, na.rm = TRUE)),
              median_abs_z = round(median(abs(z_real), na.rm = TRUE), 2)),
          by = role][order(-pass_A)]
print(agg2, row.names = FALSE)
cat("\n  非模型类（解卷积标志、泛基质、教科书标志、原代对比）本就不为预后设计，\n")
cat("  纳入角色对比会稀释效应，故单列。\n")

# ============================================================ 4. 过拟合梯度
hr("4. 过拟合梯度：开发队列上的 z 相对外部队列")
G <- M[!is.na(train_cohort) & train_cohort != "待核实",
       .(z_train = z_real[role == "开发队列"][1],
         z_ext_median = median(z_real[role != "开发队列"], na.rm = TRUE),
         z_ext_max = z_real[role != "开发队列"][which.max(abs(z_real[role != "开发队列"]))],
         n_ext = sum(role != "开发队列")), by = .(sig_id, class, train_cohort)]
G[, ratio := round(abs(z_train) / abs(z_ext_median), 2)]
G[, `:=`(z_train = round(z_train, 2), z_ext_median = round(z_ext_median, 2),
         z_ext_max = round(z_ext_max, 2))]
print(G[order(-ratio)], row.names = FALSE)
cat("\n  ratio = |开发队列 z| / |外部队列 z 的中位|\n")
cat("  ratio 明显 > 1 表示该 signature 在其开发队列上的表现远超外部 ——\n")
cat("  这是过拟合到开发队列的直接证据，非该 signature 的特异性。\n")

# ============================================================ 5. 关键对照
hr("5. 关键对照：SIG-008")
S8 <- E[sig_id == "SIG-008", .(cohort, role, z_real = round(z_real, 2),
                               p_emp_A = round(p_emp_A, 4))]
print(S8[order(role, cohort)], row.names = FALSE)
cat("\n  SIG-008 的训练队列是 GSE41271，不在底座 A 中。\n")
cat("  因此它在本研究的六个队列上【全部是外部验证】。\n")
cat("  若它的通过数与那些在 TCGA 上训练的 signature 相当，\n")
cat("  说明后者的高通过率不能全部归因于过拟合。反之则支持过拟合解释。\n")

# ============================================================ 6. 落盘
hr("6. 输出")
f1 <- file.path(CH, "data_derived", sprintf("results_by_cohort_role_%s.csv", STAMP))
fwrite(E[, .(cohort, sig_id, class, role, used_in_original, is_train,
             train_cohort, z_real, zA_p95, p_emp_A, p_emp_B, coexpr_qc_pass)], f1)
kv("逐组合角色标注", f1)
f2 <- file.path(CH, "data_derived", sprintf("overfit_gradient_%s.csv", STAMP))
fwrite(G, f2); kv("过拟合梯度", f2)

hr("7. 须记入 SAP §14 的内容")
cat("  # 4 | ", STAMP, " | 新增事后探索性分析：按『开发队列 vs 独立队列』拆分通过率\n", sep = "")
cat("  #      理由：统计筛选类的 15 次通过中有 6 次发生在其开发队列（TCGA），\n")
cat("  #            不拆分会把过拟合误读为特异性。\n")
cat("  #      性质：SAP §10 路径 (b) 未明文包含该维度，属 SAP 之外的事后分析。\n")
cat("  #            但各 signature 用过哪些队列，已于清单锁定前登记在\n")
cat("  #            overlap_cohorts 字段中（锁定表 md5 可证）。\n")
cat("  #      结果产出前/后：★ 之后\n")
cat("\n  论文中须写明：该拆分为事后探索性分析，不构成预注册的检验。\n")

cat("\n输出: ", OUT, "\n", sep = "")
