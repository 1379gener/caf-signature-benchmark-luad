# CAF signature benchmark in lung adenocarcinoma

Code and results for: *Prognostic specificity of cancer-associated fibroblast signatures in lung adenocarcinoma: a benchmark against expression-matched random gene sets* (Wang Y et al., submitted to *Genes*).

## Contents
| Folder | Content |
|---|---|
| `plan/` | Internal statistical analysis plan v1.2 (Chinese original, `.docx` and `.txt`) and File S1 (note, English translation, original) |
| `registry/` | Locked signature registry (`signature_registry_locked_20260903.csv`, MD5 `4ad3218931d95dfb6ee67a7ba07b9a75`), gene lists and published weights, symbol mapping, meta-PCNA genes, registry scripts. Source PDFs are not included (copyright). |
| `scani_genes/` | SCANI myCAF gene list, from https://github.com/1379gener/scani-luad (tag v9.1) |
| `scripts/` | Benchmark, null distributions, confounding models, figures, supplementary tables, verification |
| `results/` | All result files. Current supplementary tables: `supp_tables_20261111/`. Earlier dated files are kept because the deviation log refers to them. |
| `figures/` | Figures as submitted |
| `logs/` | Run logs |
| `docs/` | Working deviation log `DEVIATIONS.md` (mostly Chinese; the English version is Table S1) and the Methods-to-code correspondence table |

## Rerunning
Scripts use the paths of the original analysis server (`~/thesis_LUAD_CAF/...`) and are kept exactly as run. `setup_layout.sh` recreates those paths as symlinks to this repository. Then add the input data below.

Input data (not included):
- GEO: GSE68465, GSE72094, GSE31210, GSE50081, GSE30219
- TCGA-LUAD expression and clinical data (Genomic Data Commons; TCGA-CDR)
- CPTAC-LUAD (LinkedOmics)

Processed inputs are expected under `~/thesis_LUAD_CAF/01_bases/` and `~/caf_spatial_project/data/`.

## Software
R 4.3.1; GSVA 1.48.3, survival 3.5-7, data.table 1.18.2.1, org.Hs.eg.db 3.17.0, MCPcounter 1.2.0, xCell 1.1.0, EPIC 1.1.7, estimate 1.0.13, dendextend 1.17.1, ggplot2 4.0.2.

## Integrity
`MD5SUMS.txt` lists MD5 checksums of all files.
