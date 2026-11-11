# SIG-007 year correction in Table S2 + Table S1 No. 34. Reads 20261110c, writes 20261111. No analysis.
suppressMessages(library(openxlsx))
inF  <- "data_derived/supp_tables_20261110c/SCI-2_Supplementary_Tables_S1-S8_EN.xlsx"
outD <- "data_derived/supp_tables_20261111"; outF <- file.path(outD, basename(inF))
stopifnot(file.exists(inF), !file.exists(outF)); dir.create(outD, showWarnings = FALSE)
wb <- loadWorkbook(inF)

## S2
s2 <- read.xlsx(inF, sheet = "S2_registry", check.names = FALSE)
i  <- which(s2$id == "SIG-007"); stopifnot(length(i) == 1)
stopifnot(grepl("2026", s2$name[i]), grepl("2026", s2$source[i]))
new_name <- "Li (LDDD) 2025 CAF 8-gene model [registry label 'Wang(LDDD) 2026'; corrected, see Table S1 Nos. 21 and 34]"
new_src  <- "Lett Drug Des Discov 2025;22(12):100333 (PII S1570180826000631)"
writeData(wb, "S2_registry", new_name, startCol = match("name",   names(s2)), startRow = i + 1, colNames = FALSE)
writeData(wb, "S2_registry", new_src,  startCol = match("source", names(s2)), startRow = i + 1, colNames = FALSE)

## S1 No. 34
s1 <- read.xlsx(inF, sheet = "S1_deviation_log", check.names = FALSE)
stopifnot(nrow(s1) == 33, max(as.integer(s1[[1]])) == 33)
row34 <- s1[1, ]; row34[1, ] <- NA
row34[1, 1] <- 34L
row34[1, 2] <- "2026-11-11"
row34[1, 3] <- "SIG-007 registry record"
row34[1, 4] <- "Registry and Table S2 gave the publication year as 2026, taken from the DOI and PII. The published version is Vol. 22, Issue 12, December 2025, article 100333 (publisher's page)."
row34[1, 5] <- "After"
row34[1, 6] <- "Locked registry not edited; Table S2 name and source corrected to 2025; DOI unchanged."
row34[1, 7] <- "None."
writeData(wb, "S1_deviation_log", row34, startRow = nrow(s1) + 2, colNames = FALSE)
saveWorkbook(wb, outF)

## checks: only S1 and S2 differ
for (n in getSheetNames(inF)) {
  a <- read.xlsx(inF, sheet = n, check.names = FALSE); b <- read.xlsx(outF, sheet = n, check.names = FALSE)
  cat(sprintf("%-28s %s\n", n, if (identical(a, b)) "identical" else sprintf("CHANGED (rows %d -> %d)", nrow(a), nrow(b))))
}
b2 <- read.xlsx(outF, sheet = "S2_registry", check.names = FALSE); print(t(b2[b2$id == "SIG-007", c("name","source","doi")]))
b1 <- read.xlsx(outF, sheet = "S1_deviation_log", check.names = FALSE); print(t(tail(b1, 1)))
cat("md5", tools::md5sum(outF), "\n")
