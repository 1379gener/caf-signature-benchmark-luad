#!/usr/bin/env Rscript
# 安装实验三所需的解卷积包
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

try_install <- function(name, expr) {
  cat("\n=== ", name, " ===\n", sep = "")
  if (requireNamespace(name, quietly = TRUE)) {
    cat("  已安装 ", as.character(packageVersion(name)), "\n", sep = ""); return(invisible())
  }
  r <- tryCatch({ eval(expr); TRUE },
                error = function(e) { cat("  !! 失败: ", conditionMessage(e), "\n", sep = ""); FALSE })
  if (r && requireNamespace(name, quietly = TRUE))
    cat("  安装成功 ", as.character(packageVersion(name)), "\n", sep = "")
}

try_install("MCPcounter",
  quote(remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source", upgrade = "never")))
try_install("EPIC",
  quote(remotes::install_github("GfellerLab/EPIC", upgrade = "never")))
try_install("xCell",
  quote(remotes::install_github("dviraran/xCell", upgrade = "never")))
try_install("estimate",
  quote(install.packages("estimate", repos = "http://r-forge.r-project.org")))

cat("\n=== 汇总 ===\n")
for (p in c("MCPcounter", "EPIC", "xCell", "estimate")) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  cat(sprintf("  %-12s %s\n", p, ifelse(is.na(v), "未安装", v)))
}
