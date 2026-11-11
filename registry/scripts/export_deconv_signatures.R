#!/usr/bin/env Rscript
# ==============================================================================
# export_deconv_signatures.R —— 从解卷积包导出基质相关基因集
#
# 目标：把 MCP-counter / EPIC / xCell / ESTIMATE 中与成纤维细胞、基质相关的
#   基因集导出为 SIG-00X_genes.txt，纳入 Paper 2 的被测清单。
#
# 这四个是本领域公认的通用读数，审稿人会期待看到它们作为对照。
#
# 两段式：
#   A. 勘察 —— 打印每个包里可用的签名名称、基因数、示例基因
#   B. 导出 —— 仅在提供 --write 时落盘
#
# 用法：
#   Rscript export_deconv_signatures.R            # 只勘察
#   Rscript export_deconv_signatures.R --write    # 勘察 + 导出
# ==============================================================================

WRITE  <- "--write" %in% commandArgs(trailingOnly = TRUE)
OUTDIR <- path.expand("~/thesis_LUAD_CAF/03_signature_registry/genes")
LOGDIR <- path.expand("~/thesis_LUAD_CAF/03_signature_registry")
OUT <- file.path("/tmp", format(Sys.time(), "deconv_export_%Y%m%d_%H%M.txt"))
sink(OUT, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep = "")
sh <- function(t) cat("\n  --- ", t, "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-30s %s\n", paste0(k, ":"), paste(v, collapse = ", ")))

cat("解卷积包基因集导出\n时间: ", format(Sys.time()),
    "\n模式: ", ifelse(WRITE, "勘察 + 导出", "仅勘察（加 --write 导出）"), "\n", sep = "")

COLLECT <- list()   # name -> list(genes, source, note)

add <- function(id, name, genes, source, note = "") {
  genes <- unique(toupper(trimws(as.character(genes))))
  genes <- genes[nzchar(genes) & !is.na(genes) & genes != "NA"]
  COLLECT[[id]] <<- list(name = name, genes = genes, source = source, note = note)
  cat(sprintf("    [%s] %-34s %4d 基因   %s\n", id, name, length(genes),
              ifelse(length(genes) >= 10, "A 类", "B 类（等权，权重取全 1）")))
  cat(sprintf("           示例: %s\n", paste(head(genes, 10), collapse = ", ")))
}

# ============================================================ MCP-counter
hr("1. MCP-counter")
if (requireNamespace("MCPcounter", quietly = TRUE)) {
  kv("版本", as.character(packageVersion("MCPcounter")))
  e <- asNamespace("MCPcounter")
  objs <- ls(e)
  kv("包内对象", objs)
  # 标志基因表通常名为 MCPcounter.genes / genes
  tb <- NULL
  for (o in objs) {
    x <- tryCatch(get(o, envir = e), error = function(err) NULL)
    if (is.data.frame(x) && nrow(x) > 100) { tb <- x; kv("使用对象", o); break }
  }
  if (is.null(tb)) {
    d <- tryCatch(data(package = "MCPcounter")$results[, "Item"], error = function(e) character(0))
    kv("包内数据集", d)
    for (o in d) { x <- tryCatch({ eval(parse(text = paste0("MCPcounter::", o))) },
                                 error = function(e) NULL)
      if (is.data.frame(x) && nrow(x) > 100) { tb <- x; kv("使用数据集", o); break } }
  }
  if (!is.null(tb)) {
    kv("标志基因表维度", paste(dim(tb), collapse = " x "))
    kv("列名", colnames(tb))
    cc <- grep("Cell.?population|population|cell_type", colnames(tb), ignore.case = TRUE, value = TRUE)
    gc <- grep("HUGO|symbol|gene", colnames(tb), ignore.case = TRUE, value = TRUE)
    if (length(cc) && length(gc)) {
      cat("\n  细胞类型分布：\n")
      tt <- sort(table(tb[[cc[1]]]), decreasing = TRUE)
      for (i in seq_along(tt)) cat(sprintf("    %-30s %4d\n", names(tt)[i], tt[i]))
      fb <- grep("fibroblast", names(tt), ignore.case = TRUE, value = TRUE)
      sh(paste0("成纤维细胞相关：", paste(fb, collapse = " | ")))
      for (f in fb) add(paste0("MCP_", gsub("[^A-Za-z]", "", f)),
                        paste("MCP-counter", f), tb[[gc[1]]][tb[[cc[1]]] == f],
                        "MCPcounter 包内标志基因表")
    } else cat("  !! 未识别细胞类型列或基因列，列名见上\n")
  } else cat("  !! 未找到标志基因表\n")
} else cat("  未安装\n")

# ============================================================ EPIC
hr("2. EPIC")
if (requireNamespace("EPIC", quietly = TRUE)) {
  kv("版本", as.character(packageVersion("EPIC")))
  for (ref in c("TRef", "BRef")) {
    x <- tryCatch(get(ref, envir = asNamespace("EPIC")), error = function(e) NULL)
    if (is.null(x)) next
    sh(paste0("参比谱 ", ref))
    kv("  组成", names(x))
    if (!is.null(x$refProfiles)) kv("  细胞类型", colnames(x$refProfiles))
    if (!is.null(x$sigGenes)) {
      kv("  sigGenes 总数", length(x$sigGenes))
      if (ref == "TRef") add("EPIC_TRef_sig", "EPIC TRef signature genes（含 CAF）",
                             x$sigGenes, "EPIC::TRef$sigGenes",
                             "TRef 的细胞类型含 CAFs；sigGenes 为全部细胞类型共用的标志基因")
    }
  }
  cat("\n  注：EPIC 的 sigGenes 是全部细胞类型共用的标志基因集，非 CAF 专属。\n")
  cat("      若需 CAF 专属，须从 refProfiles 中按 CAFs 列的相对高表达自行定义 ——\n")
  cat("      那属于本研究自建 signature，不符合 SAP §3.1 第 1 条（须已发表）。\n")
  cat("      故此处按原样纳入 sigGenes，并在收录表中注明其非 CAF 专属。\n")
} else cat("  未安装\n")

# ============================================================ xCell
hr("3. xCell")
if (requireNamespace("xCell", quietly = TRUE)) {
  kv("版本", as.character(packageVersion("xCell")))
  sig <- tryCatch(get("xCell.data", envir = asNamespace("xCell")), error = function(e) NULL)
  if (is.null(sig)) {
    d <- tryCatch(data(package = "xCell")$results[, "Item"], error = function(e) character(0))
    kv("包内数据集", d)
  }
  if (!is.null(sig)) {
    kv("xCell.data 组成", names(sig))
    gs <- sig$signatures
    if (!is.null(gs)) {
      nm <- if (inherits(gs, "GeneSetCollection")) {
        if (requireNamespace("GSEABase", quietly = TRUE))
          vapply(gs, GSEABase::setName, "") else names(gs)
      } else names(gs)
      kv("签名总数", length(nm))
      fb <- grep("^Fibroblast", nm, ignore.case = TRUE, value = TRUE)
      sh(paste0("Fibroblast 相关签名（", length(fb), " 套）"))
      for (f in fb) cat("    ", f, "\n", sep = "")
      cat("\n  ★ xCell 对每个细胞类型提供 3 套来自不同数据源的签名。\n")
      cat("    SAP 未规定取哪一套。本脚本采取的规则：取三套的【并集】，\n")
      cat("    作为单一 signature 纳入，并在收录表中注明。\n")
      cat("    理由：任取其一皆为任意选择；并集是唯一无需额外判断的确定规则。\n")
      if (length(fb)) {
        gl <- unlist(lapply(fb, function(f) {
          i <- which(nm == f)
          if (inherits(gs, "GeneSetCollection")) GSEABase::geneIds(gs[[i]]) else gs[[i]] }))
        add("xCell_Fibroblasts", "xCell Fibroblasts（3 套并集）", gl,
            "xCell::xCell.data$signatures",
            paste0("并集自: ", paste(fb, collapse = " | ")))
      }
    } else cat("  !! xCell.data 中无 signatures\n")
  } else cat("  !! 未找到 xCell.data\n")
} else cat("  未安装\n")

# ============================================================ ESTIMATE
hr("4. ESTIMATE")
if (requireNamespace("estimate", quietly = TRUE)) {
  kv("版本", as.character(packageVersion("estimate")))
  sc <- tryCatch(get("SI_geneset", envir = asNamespace("estimate")), error = function(e) NULL)
  if (is.null(sc)) {
    f <- system.file("extdata", "SI_geneset.gmt", package = "estimate")
    kv("extdata 路径", f)
    if (nzchar(f) && file.exists(f)) {
      L <- readLines(f, warn = FALSE)
      for (l in L) {
        p <- strsplit(l, "\t")[[1]]
        nm2 <- p[1]; g <- p[-(1:2)]
        cat(sprintf("\n    集合 %s : %d 基因\n", nm2, length(g)))
        if (grepl("stromal", nm2, ignore.case = TRUE))
          add("ESTIMATE_stromal", "ESTIMATE StromalSignature", g,
              "estimate 包 SI_geneset.gmt",
              "★ 不得用作实验二主分析的纯度协变量（SAP §7.2）")
      }
    } else cat("  !! 未找到 SI_geneset\n")
  } else {
    kv("SI_geneset 维度", paste(dim(sc), collapse = " x "))
    print(utils::head(sc, 3))
    cat("  -> 结构非预期，请把上面几行贴回\n")
  }
} else cat("  未安装\n")

# ============================================================ 汇总与导出
hr("5. 汇总")
if (!length(COLLECT)) {
  cat("  未采集到任何基因集。请把上面各节的对象名 / 列名贴回，我据实修改。\n")
} else {
  cat(sprintf("  共采集 %d 个基因集：\n\n", length(COLLECT)))
  cat(sprintf("  %-22s %-36s %6s  %s\n", "ID", "名称", "基因数", "类别"))
  for (id in names(COLLECT)) {
    x <- COLLECT[[id]]
    cat(sprintf("  %-22s %-36s %6d  %s\n", id, substr(x$name, 1, 36), length(x$genes),
                ifelse(length(x$genes) >= 10, "A", "B")))
  }

  sh("两两重叠（提示 SAP §3.2 的 >90% 重叠排除规则）")
  ids <- names(COLLECT)
  if (length(ids) > 1) for (i in seq_along(ids)) for (j in seq_along(ids)) {
    if (j <= i) next
    a <- COLLECT[[ids[i]]]$genes; b <- COLLECT[[ids[j]]]$genes
    ov <- length(intersect(a, b))
    cat(sprintf("    %-20s vs %-20s 交集 %3d  (%.0f%% / %.0f%%)\n",
        ids[i], ids[j], ov, 100*ov/length(a), 100*ov/length(b)))
  }

  if (WRITE) {
    dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
    exist <- list.files(OUTDIR, pattern = "^SIG-[0-9]+_genes\\.txt$")
    nxt <- if (length(exist)) max(as.integer(sub("^SIG-([0-9]+).*", "\\1", exist))) + 1L else 1L
    sh("导出")
    meta <- list()
    for (id in names(COLLECT)) {
      x <- COLLECT[[id]]
      sid <- sprintf("SIG-%03d", nxt)
      writeLines(x$genes, file.path(OUTDIR, paste0(sid, "_genes.txt")))
      if (length(x$genes) < 10)
        write.csv(data.frame(gene = x$genes, weight = 1),
                  file.path(OUTDIR, paste0(sid, "_weights.csv")), row.names = FALSE, quote = FALSE)
      cat(sprintf("    %s  <- %-34s %4d 基因\n", sid, x$name, length(x$genes)))
      meta[[sid]] <- data.frame(sig_id = sid, name = x$name, n_genes = length(x$genes),
        class = ifelse(length(x$genes) >= 10, "A", "B"),
        source = x$source, note = x$note, exported = format(Sys.time(), "%Y-%m-%d %H:%M"),
        stringsAsFactors = FALSE)
      nxt <- nxt + 1L
    }
    m <- do.call(rbind, meta)
    f <- file.path(LOGDIR, "deconv_signatures_meta.csv")
    write.csv(m, f, row.names = FALSE)
    cat("\n  元数据 -> ", f, "\n", sep = "")
    cat("  请据此填写收录表对应行。\n")
  } else {
    cat("\n  （未加 --write，本次只勘察不导出）\n")
  }
}

hr("结束")
cat("输出: ", OUT, "\n", sep = "")
