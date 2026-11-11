#!/usr/bin/env bash
# Recreate the original server paths as symlinks so the scripts run unchanged.
# Scripts refer to ~/thesis_LUAD_CAF/...; run this on a machine where that path does not yet exist.
set -euo pipefail
R=$(cd "$(dirname "$0")" && pwd); T=~/thesis_LUAD_CAF
[ -e "$T" ] && { echo "STOP: $T already exists"; exit 1; }
C=$T/02_chapters/ch04_benchmark; mkdir -p "$T"/00_admin "$T"/05_manuscripts/SCI-1_IJMS/tables "$C" "$T"/01_bases
for f in "$R"/plan/*; do ln -s "$f" "$T"/00_admin/; done
ln -s "$R"/registry "$T"/03_signature_registry
ln -s "$R"/scani_genes/S1_Table_myCAF_signature_50genes.csv "$T"/05_manuscripts/SCI-1_IJMS/tables/
ln -s "$R"/scripts "$C"/scripts; ln -s "$R"/results "$C"/data_derived; ln -s "$R"/logs "$C"/logs
mkdir -p "$C"/figures "$C"/text; ln -s "$R"/docs/DEVIATIONS.md "$C"/DEVIATIONS.md
for f in verify_ch04.R verify_ch3.R rerun_gse30219_n1000.sh; do ln -s "$R"/scripts/$f "$C"/$f; done
echo "Layout created under $T. Place processed inputs under $T/01_bases/ and ~/caf_spatial_project/data/ (see README)."
