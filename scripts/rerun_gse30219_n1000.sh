#!/usr/bin/env bash
# ==============================================================================
# rerun_gse30219_n1000.sh —— 补跑 n_null=50 的两个组合（GSE30219 × SIG-001 / SIG-002）
#   原因：2026-09-03 23:50 的 --test 运行产物（n_null=50）一直留在 null_benchmark/，
#         run_benchmark.R 跳过已有 RDS，故从未按 SAP 的 1000 次重跑。
#   做法：先登记偏差 → 旧 RDS 移到 superseded/（不删）→ --force 只重跑这两个组合
#   seed 不变（make_seed 由 队列|签名|类型 决定），只是重抽次数从 50 到 1000
# 用法：bash ~/thesis_LUAD_CAF/02_chapters/ch04_benchmark/rerun_gse30219_n1000.sh
# ==============================================================================
set -euo pipefail
CH=~/thesis_LUAD_CAF/02_chapters/ch04_benchmark
NB=$CH/data_derived/null_benchmark
SUP=$CH/data_derived/null_benchmark_superseded_$(date +%Y%m%d)
DEV=$CH/DEVIATIONS.md
mkdir -p "$SUP"

for s in SIG-001 SIG-002; do
  f=$NB/GSE30219__$s.rds
  [ -f "$f" ] || { echo "缺 $f"; exit 1; }
  n=$(Rscript -e "cat(readRDS('$f')\$n_null)")
  [ "$n" = "50" ] || { echo "$f 的 n_null=$n，不是 50，停止"; exit 1; }
done

cat >> "$DEV" <<EOF

## $(date +%Y-%m-%d) GSE30219 × SIG-001/SIG-002 零分布次数补足
- 发现：verify_ch3.R / 盘点显示这两个组合 n_null = 50（2026-09-03 23:50 --test 运行产物），SAP 规定 1000
- 原因：run_benchmark.R 对已存在的 RDS 跳过，全量运行未覆盖
- 处置：旧 RDS 移至 $(basename "$SUP")/；以 --force 重跑，seed 规则不变
- 影响：只影响这两个组合的 p_emp_A / p_emp_B 及其下游计数；真实 z 不受影响
EOF
echo "已登记：$DEV"

mv "$NB/GSE30219__SIG-001.rds" "$NB/GSE30219__SIG-002.rds" "$SUP/"
md5sum "$SUP"/*.rds > "$SUP/MD5.txt"

cd "$CH/scripts"
Rscript run_benchmark.R --cohort=GSE30219 --sig=SIG-001,SIG-002 --force

Rscript -e '
for (s in c("SIG-001","SIG-002")) {
  o <- readRDS(sprintf("'"$SUP"'/GSE30219__%s.rds", s)); n <- readRDS(sprintf("'"$NB"'/GSE30219__%s.rds", s))
  cat(sprintf("%s  n_null %d -> %d | z_real %.4f -> %.4f | p_emp_A %.4f -> %.4f | p_emp_B %.4f -> %.4f | QC %s -> %s\n",
      s, o$n_null, n$n_null, o$z_real, n$z_real, o$p_emp_A, n$p_emp_A, o$p_emp_B, n$p_emp_B,
      o$coexpr_qc_pass, n$coexpr_qc_pass)) }'
