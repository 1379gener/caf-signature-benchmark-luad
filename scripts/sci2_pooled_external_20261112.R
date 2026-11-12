# sci2_pooled_external_20261112.R —— 每个签名在外部队列上的合并检验（Stouffer，权重 sqrt(死亡数)）；零分布 A 成员按序号跨队列配对合并
# 登记：DEVIATIONS.md 2026-11-12「外部队列合并检验」（运行之前）。只读主分析 rds；输出 pooled_external_20261112/ 与 supp_tables_20261112b/
source("~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/scripts/sci2_common_20261111.R")
suppressMessages(library(openxlsx))
OD <- "data_derived/pooled_external_20261112"; if (dir.exists(OD)) stop("已存在，不覆盖：", OD); dir.create(OD)
sink("logs/pooled_external_20261112.log", split = TRUE); cat("run:", format(Sys.time()), "\n")

## ---- 主分析 rds：z_real、z_null_A、事件数
f1 <- list.files("data_derived/null_benchmark", "\\.rds$", full.names = TRUE)
f2 <- list.files("data_derived/null_benchmark_scani_supp", "\\.rds$", full.names = TRUE)
R <- list(); for (f in c(f1, f2)) { x <- readRDS(f); if (!isTRUE(x$evaluable)) next
  ev <- if (!is.null(x$ev_used)) x$ev_used else x$n_events
  R[[paste(x$cohort, x$sig_id)]] <- list(cohort = x$cohort, sig_id = x$sig_id, class = if (is.null(x$class)) "A" else x$class,
    z = x$z_real, zn = x$z_null_A, ev = ev, p = x$p_emp_A) }
stopifnot(length(R) == 77, all(vapply(R, function(r) length(r$zn) == 1000, TRUE)))
D <- rbindlist(lapply(R, function(r) data.table(cohort = r$cohort, sig_id = r$sig_id, class = r$class, z = r$z, ev = r$ev, p = r$p)))

## ---- 队列角色（取自 2026-10-25 的队列角色分析输出）
rf <- list.files("data_derived", "role", recursive = TRUE, full.names = TRUE); rf <- rf[grepl("\\.csv$", rf)]
ok <- vapply(rf, function(f) all(c("cohort", "sig_id", "role") %in% names(fread(f, nrows = 1))), TRUE); rf <- rf[ok]
rf <- sort(rf, decreasing = TRUE); cat("角色文件候选：\n"); print(rf); cat("使用：", rf[1], "\n")
RL <- fread(rf[1])[, .(cohort = gsub("-", "_", as.character(cohort)), sig_id, role)]
RL[cohort == "TCGA", cohort := "TCGA_LUAD"]
sc <- data.table(cohort = names(NEXP), sig_id = "SCANI_myCAF", role = fifelse(names(NEXP) == "TCGA_LUAD", "development", fifelse(names(NEXP) == "GSE72094", "validation", "independent")))
RL <- rbind(RL[sig_id != "SCANI_myCAF"], sc); stopifnot(!anyDuplicated(RL[, .(cohort, sig_id)]))
RL[, role3 := fifelse(grepl("开发|develop", role, ignore.case = TRUE), "development",
              fifelse(grepl("验证|valid", role, ignore.case = TRUE), "validation", "independent"))]
print(RL[, .N, by = .(role, role3)])
D <- merge(D, RL[, .(cohort, sig_id, role3)], by = c("cohort", "sig_id"), all.x = TRUE)
if (anyNA(D$role3)) { cat("!! 无角色的组合：\n"); print(D[is.na(role3)]) }
stopifnot(!anyNA(D$role3))
print(dcast(D[, .N, by = .(class, role3)], class ~ role3, value.var = "N"))
stopifnot(D[class == "B", .N, by = role3][order(role3), N] == c(5L, 28L, 8L), D[class == "A", .N, by = role3][order(role3), N] == c(2L, 28L, 6L))

## ---- 合并检验
pool <- function(keys) {
  w <- sqrt(vapply(keys, function(k) R[[k]]$ev, 0)); z <- vapply(keys, function(k) R[[k]]$z, 0)
  ZN <- sapply(keys, function(k) R[[k]]$zn); okm <- stats::complete.cases(ZN)
  Z <- sum(w * z) / sqrt(sum(w^2)); Zn <- as.numeric(ZN[okm, , drop = FALSE] %*% w) / sqrt(sum(w^2))
  list(Z = Z, p = (1 + sum(abs(Zn) >= abs(Z))) / (1 + length(Zn)), p95 = unname(quantile(abs(Zn), .95)), n_null = length(Zn),
       n_coh = length(keys), deaths = sum(w^2)) }
out <- list()
for (s in unique(D$sig_id)) for (set in c("external", "independent")) {
  d <- D[sig_id == s & (if (set == "external") role3 != "development" else role3 == "independent")]
  if (nrow(d) < 2) next
  r <- pool(paste(d$cohort, d$sig_id))
  out[[paste(s, set)]] <- data.table(sig_id = s, class = d$class[1], cohort_set = set, cohorts = paste(sort(d$cohort), collapse = ", "),
    n_cohorts = r$n_coh, deaths = r$deaths, z_pooled = r$Z, abs_z_null_p95 = r$p95, p_emp_pooled = r$p, n_null = r$n_null,
    n_single_cohort_pass_unadj = sum(d$p < 0.05)) }
P <- rbindlist(out)
P[, q_pooled := p.adjust(p_emp_pooled, "BH"), by = .(class, cohort_set)]
setorder(P, cohort_set, class, sig_id); print(P[, .(cohort_set, class, sig_id, n_cohorts, deaths, z = round(z_pooled, 2), p95 = round(abs_z_null_p95, 2), p = round(p_emp_pooled, 3), q = round(q_pooled, 3))])
fwrite(P, file.path(OD, "pooled_external_20261112.csv")); fwrite(D, file.path(OD, "pooled_inputs_roles_20261112.csv"))
for (set in c("external", "independent")) for (cl in c("A", "B")) { x <- P[cohort_set == set & class == cl]
  cat(sprintf("%s, class %s: %d signatures; |z_pooled| >= 1.96: %d; p_emp < 0.05: %d; q < 0.05: %d (%s)\n", set, cl, nrow(x),
    sum(abs(x$z_pooled) >= 1.96), sum(x$p_emp_pooled < 0.05), sum(x$q_pooled < 0.05), paste(x[q_pooled < 0.05, sig_id], collapse = ", "))) }

## ---- 补充表 20261112b：S8g + S1 No. 38
inF <- "data_derived/supp_tables_20261112/SCI-2_Supplementary_Tables_S1-S9_EN.xlsx"
outD <- "data_derived/supp_tables_20261112b"; outF <- file.path(outD, basename(inF)); stopifnot(file.exists(inF), !dir.exists(outD)); dir.create(outD)
wb <- loadWorkbook(inF); r4 <- function(x) round(x, 4)
addWorksheet(wb, "S8g_pooled_external")
writeData(wb, "S8g_pooled_external", P[, .(sig_id, class, cohort_set, cohorts, n_cohorts, deaths, z_pooled = r4(z_pooled),
  abs_z_null_A_p95 = r4(abs_z_null_p95), p_emp_pooled = r4(p_emp_pooled), q_pooled_within_class = r4(q_pooled), n_null,
  n_single_cohort_pass_unadjusted = n_single_cohort_pass_unadj)], headerStyle = createStyle(textDecoration = "bold"))
s1 <- read.xlsx(inF, "S1_deviation_log", check.names = FALSE); stopifnot(nrow(s1) == 37)
r <- s1[1, ]; r[1, ] <- NA; r[1, 1] <- 38L; r[1, 2] <- "2026-11-12"; r[1, 3] <- "Pooled test in external cohorts (not in the plan)"
r[1, 4] <- "For each signature, z was combined across cohorts other than its development cohort (and, separately, across fully independent cohorts) by weighted Stouffer (weights: square root of deaths); null-A members were combined across cohorts by index."
r[1, 5] <- "After"; r[1, 6] <- "Post hoc analysis registered before it was run; reported in Table S8g; main analysis unchanged."
xe <- P[cohort_set == "external"]; xi <- P[cohort_set == "independent"]
r[1, 7] <- sprintf("External cohorts: %d of %d signatures exceeded the pooled null after false-discovery-rate control within class; fully independent cohorts: %d of %d (Table S8g).",
  sum(xe$q_pooled < 0.05), nrow(xe), sum(xi$q_pooled < 0.05), nrow(xi))
writeData(wb, "S1_deviation_log", r, startRow = nrow(s1) + 2, colNames = FALSE); saveWorkbook(wb, outF)
for (n in getSheetNames(inF)) { a <- read.xlsx(inF, n, check.names = FALSE); b <- read.xlsx(outF, n, check.names = FALSE)
  cat(sprintf("%-28s %s\n", n, if (identical(a, b)) "identical" else sprintf("CHANGED (rows %d -> %d)", nrow(a), nrow(b)))) }
b <- read.xlsx(outF, "S8g_pooled_external"); cat(sprintf("S8g_pooled_external NEW %d rows; CJK: %s\n", nrow(b), any(grepl("[一-鿿]", unlist(b)))))
cat("S1 No. 38 effect:", r[1, 7], "\n"); cat("md5", tools::md5sum(outF), "\n"); cat("done:", format(Sys.time()), "\n"); sink()
