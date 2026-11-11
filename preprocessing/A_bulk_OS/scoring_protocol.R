#!/usr/bin/env Rscript
# ==============================================================================
# scoring_protocol.R —— 底座 A 统一打分协议（可 source 的模块 + 自检）
#
# 目的：把打分口径写成【代码】而不是文档。下游一律通过本模块读数据、打分，
#       不自行 readRDS、不自行调 gsva —— 口径由代码强制，不靠人记。
#
# 依据：《口径取证备忘_20260829.md》§4。取证已证实当年四个队列用了四套
#       塌缩规则与基因宇宙，GSE72094 的 SCANI 分数与 TCGA 口径的 Spearman
#       相关仅 0.53。统一协议是 Paper 2 的前提。
#
# 用法：
#   source("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/scripts/scoring_protocol.R")
#   co <- load_cohort("GSE72094")
#   cv <- gene_coverage(my_genes, co$expr)
#   sc <- score_ssgsea(co$expr, list(sig = my_genes))
#
# 直接运行则执行自检：Rscript scoring_protocol.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(GSVA)
})

# ============================== 协议常量（改这里，不改下面） ==================
PROTOCOL <- list(
  version        = "1.0 (2026-08-29)",
  gsva_version   = "1.48.3",          # 锁定；旧式接口
  ssgsea_norm    = FALSE,             # 归一化对 Z 标准化后的指数无影响，但仍锁定
  gene_collapse  = "同名取均表达最大（已在整备阶段完成）",
  gene_universe  = "各队列自有基因池（该队列检出的全部基因）",
  log_scale      = "统一 log2(x+1)；manifest 中 need_log=FALSE 的队列不再变换",
  min_coverage   = 0.70,              # signature 基因覆盖率阈值
  dedup_then_z   = TRUE,              # 先去重再 Z 标准化（整备阶段已去重）
  drop_constant  = FALSE,           # 不预先剔除常数行，交由 GSVA 处理
  alias_fallback = TRUE              # 不预先剔除常数行，交由 GSVA 处理并记录
)

DER <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived")

COHORTS <- c("TCGA_LUAD", "GSE68465", "GSE72094", "GSE31210",
             "GSE50081", "GSE30219", "CPTAC_LUAD")

# ============================== 数据读取 ======================================
#' 读入一个队列，返回已按协议处理好的表达矩阵与临床表
#' @return list(acc, expr, clin, man, id_col, n, n_events, has_survival)
load_cohort <- function(acc, root = DER, verbose = TRUE) {
  d <- file.path(root, acc)
  if (!dir.exists(d)) stop("队列目录不存在: ", d)

  fe <- list.files(d, pattern = "_expr_symbol.*\\.rds$", full.names = TRUE)
  fc <- list.files(d, pattern = "_clinical_clean\\.csv$", full.names = TRUE)
  fm <- list.files(d, pattern = "^cohort_manifest_.*\\.csv$", full.names = TRUE)
  if (!length(fe) || !length(fc)) stop("队列文件不全: ", acc)

  expr <- readRDS(fe[1])
  clin <- fread(fc[1])
  man  <- if (length(fm)) fread(fm[1]) else data.table()

  # ID 列：各队列命名不同（patient_id / gsm / sample_id）
  id_col <- intersect(c("patient_id", "gsm", "sample_id"), colnames(clin))[1]
  if (is.na(id_col)) stop("未能识别 ID 列: ", acc)

  # 断言：表达矩阵列与临床表行一一对应，顺序一致
  stopifnot(identical(colnames(expr), as.character(clin[[id_col]])))
  stopifnot(!anyDuplicated(colnames(expr)))

  # 对数变换：以 manifest 为准；manifest 缺失时按取值域推断
  need_log <- if ("need_log" %in% colnames(man)) as.logical(man$need_log[1]) else {
    mx <- max(expr[seq_len(min(2000, nrow(expr))), , drop = FALSE], na.rm = TRUE)
    mx >= 25
  }
  if (isTRUE(need_log)) expr <- log2(expr + 1)

  has_surv <- if ("has_survival" %in% colnames(man)) as.logical(man$has_survival[1]) else
                all(c("OS_time", "OS_event") %in% colnames(clin)) &&
                any(!is.na(clin$OS_time))

  if (verbose) {
    cat(sprintf("  %-12s %5d 样本 x %6d 基因  | log 变换 %-3s | 生存 %-3s | 事件 %s\n",
        acc, ncol(expr), nrow(expr), ifelse(isTRUE(need_log), "是", "否"),
        ifelse(isTRUE(has_surv), "有", "无"),
        ifelse(isTRUE(has_surv), sum(clin$OS_event, na.rm = TRUE), "-")))
  }

  list(acc = acc, expr = expr, clin = clin, man = man, id_col = id_col,
       n = ncol(expr), n_genes = nrow(expr),
       n_events = if (isTRUE(has_surv)) sum(clin$OS_event, na.rm = TRUE) else NA_integer_,
       has_survival = isTRUE(has_surv), logged = isTRUE(need_log))
}

# ============================== 覆盖率 ========================================

# ============================== 别名回退（SAP §4.1 补充条款）==================
#' signature 基因在矩阵中缺失时，尝试用其别名取值
#' 保守规则：恰有一个别名存在于矩阵中才回退；零个或多个均视为真实缺失。
#' @return list(present, rescued, missing, map)
#'   present 可用于取值的行名（原符号 + 回退到的别名）
#'   map     data.frame(gene, used) 记录全部回退映射
resolve_aliases <- function(genes, rownames_expr) {
  genes <- unique(as.character(genes))
  direct <- intersect(genes, rownames_expr)
  miss   <- setdiff(genes, rownames_expr)
  if (!length(miss) || !requireNamespace("org.Hs.eg.db", quietly = TRUE))
    return(list(present = direct, rescued = character(0), missing = miss,
                map = data.frame(gene = character(0), used = character(0),
                                 stringsAsFactors = FALSE)))
  suppressPackageStartupMessages({
    requireNamespace("AnnotationDbi", quietly = TRUE) })
  rescued <- character(0); used <- character(0); still <- character(0)
  for (g in miss) {
    eid <- tryCatch(suppressMessages(AnnotationDbi::select(
             org.Hs.eg.db::org.Hs.eg.db, keys = g, columns = "ENTREZID",
             keytype = "SYMBOL")$ENTREZID[1]), error = function(e) NA)
    al <- if (is.na(eid)) character(0) else tryCatch(suppressMessages(
            AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = eid,
              columns = "ALIAS", keytype = "ENTREZID")$ALIAS),
            error = function(e) character(0))
    hit <- setdiff(intersect(unique(al), rownames_expr), c(direct, used))
    if (length(hit) == 1L) {                 # 恰有一个才回退
      rescued <- c(rescued, g); used <- c(used, hit)
    } else {
      still <- c(still, g)
    }
  }
  list(present = c(direct, used), rescued = rescued, missing = still,
       map = data.frame(gene = rescued, used = used, stringsAsFactors = FALSE))
}

#' signature 基因在该队列的覆盖率（含别名回退，SAP §4.1）
#' 低于 PROTOCOL$min_coverage 的组合标记为不可评估，不给经验 p 值
gene_coverage <- function(genes, expr, use_alias = TRUE) {
  genes <- unique(as.character(genes))
  if (use_alias) {
    r <- resolve_aliases(genes, rownames(expr))
    present <- r$present; missing <- r$missing
    rescued <- r$rescued; map <- r$map
  } else {
    present <- intersect(genes, rownames(expr)); missing <- setdiff(genes, present)
    rescued <- character(0)
    map <- data.frame(gene = character(0), used = character(0), stringsAsFactors = FALSE)
  }
  frac <- length(present) / length(genes)
  list(n_total = length(genes), n_present = length(present),
       frac = frac, evaluable = frac >= PROTOCOL$min_coverage,
       present = present, missing = missing,
       n_rescued = length(rescued), alias_map = map)
}

# ============================== 打分 ==========================================
#' ssGSEA 打分，GSVA 1.48.3 旧式接口，参数锁定
#' @param expr 已按协议处理的表达矩阵（log2 尺度）
#' @param gene_sets 具名 list
score_ssgsea <- function(expr, gene_sets, verbose = FALSE) {
  v <- as.character(packageVersion("GSVA"))
  if (package_version(v) >= "1.50")
    stop("GSVA 版本为 ", v, "，本协议锁定 1.48.x 旧式接口。",
         "混用两套接口会产生不可比的分数。")
  gs <- lapply(gene_sets, function(g) resolve_aliases(g, rownames(expr))$present)
  empty <- names(gs)[vapply(gs, length, 0L) == 0L]
  if (length(empty)) stop("以下基因集在该队列中无任何基因: ", paste(empty, collapse = ", "))
  GSVA::gsva(expr, gs, method = "ssgsea",
             ssgsea.norm = PROTOCOL$ssgsea_norm, verbose = verbose)
}

#' 队列内 Z 标准化
zscore <- function(x) as.numeric(scale(as.numeric(x)))

#' 合成 SCANI（示例：复合指数的标准写法）
compose_scani <- function(expr, mycaf_genes) {
  stopifnot(all(c("TIMP1", "PDPN") %in% rownames(expr)))
  s <- score_ssgsea(expr, list(myCAF = mycaf_genes))["myCAF", ]
  zscore(s) + zscore(expr["TIMP1", ]) - zscore(expr["PDPN", ])
}

# ============================== 汇总 manifest =================================
#' 把七个队列的 manifest 合并成一张总表
build_registry <- function(root = DER, cohorts = COHORTS) {
  rows <- list()
  for (a in cohorts) {
    fm <- list.files(file.path(root, a), pattern = "^cohort_manifest_.*\\.csv$",
                     full.names = TRUE)
    if (!length(fm)) next
    m <- fread(fm[1])
    rows[[a]] <- m
  }
  rbindlist(rows, fill = TRUE)
}

# ============================== 自检 ==========================================
if (sys.nframe() == 0L) {

  OUT <- file.path("/tmp", format(Sys.time(), "protocol_selftest_%Y%m%d_%H%M.txt"))
  sink(OUT, split = TRUE)
  on.exit(sink(), add = TRUE)

  hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep = "")

  cat("统一打分协议 · 自检\n时间: ", format(Sys.time()), "\n输出: ", OUT, "\n", sep = "")

  hr("0. 协议常量")
  for (n in names(PROTOCOL)) cat(sprintf("  %-16s %s\n", n, PROTOCOL[[n]]))
  cat(sprintf("\n  当前环境 GSVA: %s  %s\n", as.character(packageVersion("GSVA")),
      ifelse(as.character(packageVersion("GSVA")) == PROTOCOL$gsva_version,
             "✔ 与协议一致", "!! 与协议不符")))

  hr("1. 载入七个队列")
  CO <- list()
  for (a in COHORTS) {
    CO[[a]] <- tryCatch(load_cohort(a), error = function(e) {
      cat(sprintf("  %-12s !! %s\n", a, conditionMessage(e))); NULL })
  }
  CO <- CO[!vapply(CO, is.null, TRUE)]
  cat(sprintf("\n  成功载入 %d / %d 个队列\n", length(CO), length(COHORTS)))

  surv <- CO[vapply(CO, function(x) x$has_survival, TRUE)]
  cat(sprintf("  其中有生存数据 %d 个，合计 %d 例、%d 事件\n", length(surv),
      sum(vapply(surv, function(x) x$n, 0L)),
      sum(vapply(surv, function(x) x$n_events, 0L))))

  hr("2. 关键基因齐备性")
  keyg <- c("TIMP1", "PDPN", "PCNA", "MKI67", "ACTA2", "COL1A1", "FAP", "POSTN", "DCN")
  cat(sprintf("  %-10s %s\n", "基因", paste(sprintf("%-11s", names(CO)), collapse = "")))
  for (g in keyg)
    cat(sprintf("  %-10s %s\n", g,
        paste(sprintf("%-11s", vapply(CO, function(x)
          ifelse(g %in% rownames(x$expr), "有", "!! 缺"), "")), collapse = "")))

  hr("3. myCAF signature 的逐队列覆盖率（阈值 70%）")
  gsf <- path.expand("~/caf_spatial_project/results/01_scRNA/CAF_subtype_signatures.rds")
  if (file.exists(gsf)) {
    sigs <- readRDS(gsf)
    for (sn in names(sigs)) {
      cat(sprintf("\n  [%s]  %d 基因\n", sn, length(sigs[[sn]])))
      for (a in names(CO)) {
        cv <- gene_coverage(sigs[[sn]], CO[[a]]$expr)
        cat(sprintf("    %-12s %2d/%2d  %5.1f%%  %s%s\n", a, cv$n_present, cv$n_total,
            100 * cv$frac, ifelse(cv$evaluable, "可评估", "!! 不可评估"),
            ifelse(length(cv$missing) && length(cv$missing) <= 6,
                   paste0("  缺: ", paste(cv$missing, collapse = ",")), "")))
      }
    }
  } else cat("  !! 未找到基因集文件\n")

  hr("4. 冒烟测试：在最小队列上跑一次完整打分")
  small <- names(CO)[which.min(vapply(CO, function(x) x$n, 0L))]
  cat("  选用队列: ", small, "\n", sep = "")
  if (file.exists(gsf)) {
    t0 <- Sys.time()
    sc <- compose_scani(CO[[small]]$expr, sigs$myCAF)
    cat(sprintf("  compose_scani 完成，耗时 %.1f 秒\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    cat(sprintf("  分数：n=%d  范围 [%.3f, %.3f]  均值 %.4f  标准差 %.4f\n",
        length(sc), min(sc), max(sc), mean(sc), sd(sc)))
    cat("  （均值应约为 0，标准差约为 sqrt(3)≈1.73 量级 —— 三个 Z 分数的线性组合）\n")
  }

  hr("5. 合并 manifest")
  reg <- build_registry()
  if (nrow(reg)) {
    keep <- intersect(c("cohort","platform","n_final","n_samples","n_events","n_genes",
                        "expr_scale","need_log","time_src_unit","has_survival",
                        "n_with_stage","median_followup_days"), colnames(reg))
    print(reg[, ..keep], row.names = FALSE)
    f <- file.path(DER, "cohort_registry.csv")
    write.csv(reg, f, row.names = FALSE)
    cat("\n  已写出总表: ", f, "\n", sep = "")
  } else cat("  !! 无 manifest 可合并\n")

  hr("6. 协议自检结论")
  cat("  本模块可被下游 source。下游脚本不得自行 readRDS 表达矩阵、\n")
  cat("  不得自行调用 gsva —— 一律通过 load_cohort() / score_ssgsea()。\n")
  cat("\n  待办（不阻塞实验一）：\n")
  cat("    - GSE30219 的分期推导有误（85 例中 81 例判为 I 期），只影响 M4 模型\n")
  cat("    - CPTAC 无生存数据，has_survival=FALSE，自动排除于生存分析之外\n")
  cat("    - 【待定】signature 收录是否分 A/B 两类（见备忘 §7 第 1 条）\n")
  cat("    - 【待定】GSE31210 是否遵从 exclude 标志（22 例）\n")

  cat("\n输出: ", OUT, "\n", sep = "")
}
