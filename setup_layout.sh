#!/usr/bin/env bash
# Recreate the original server paths as symlinks so the scripts run unchanged.
# Scripts refer to ~/thesis_LUAD_CAF/... and ~/caf_spatial_project/...; run this where those paths do not exist.
# Extract the release asset into inputs/derived/ first (see README).
set -euo pipefail
R=$(cd "$(dirname "$0")" && pwd); T=~/thesis_LUAD_CAF; P=~/caf_spatial_project
[ -e "$T" ] && { echo "STOP: $T already exists"; exit 1; }
[ -e "$P" ] && { echo "STOP: $P already exists"; exit 1; }
ls "$R"/inputs/derived/*/*_expr_symbol*.rds >/dev/null 2>&1 || echo "NOTE: expression matrices not found in inputs/derived (extract the release asset)"

A=$T/01_bases/A_bulk_OS; C=$T/02_chapters/ch04_benchmark
mkdir -p "$T"/00_admin "$T"/05_manuscripts/SCI-1_IJMS/tables "$A" "$C"/figures "$C"/text "$R"/inputs/raw/purity
for f in "$R"/plan/*; do ln -s "$f" "$T"/00_admin/; done
ln -s "$R"/registry "$T"/03_signature_registry
ln -s "$R"/scani_genes/S1_Table_myCAF_signature_50genes.csv "$T"/05_manuscripts/SCI-1_IJMS/tables/
ln -s "$R"/preprocessing/A_bulk_OS "$A"/scripts
ln -s "$R"/inputs/derived "$A"/derived
ln -s "$R"/inputs/raw "$A"/raw
ln -s "$R"/scripts "$C"/scripts; ln -s "$R"/results "$C"/data_derived; ln -s "$R"/logs "$C"/logs
ln -s "$R"/docs/DEVIATIONS.md "$C"/DEVIATIONS.md
for f in verify_ch04.R verify_ch3.R rerun_gse30219_n1000.sh; do ln -s "$R"/scripts/$f "$C"/$f; done
mkdir -p "$P"/results/01_scRNA "$P"/data/TCGA_LUAD "$P"/scripts
ln -s "$R"/scani_genes/CAF_subtype_signatures.rds "$P"/results/01_scRNA/
ln -s "$R"/preprocessing/tcga_download/step5_tcga_download.R "$P"/scripts/
echo "Layout created under $T and $P."
echo "TCGA_LUAD_SE.rds (needed by sci2_purity_tcga.R and sci2_m4_stage_tcga_sens.R) is produced by $P/scripts/step5_tcga_download.R."
echo "Place the Aran et al. purity table in $R/inputs/raw/purity/."
