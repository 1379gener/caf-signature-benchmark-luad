# CAF signature benchmark in lung adenocarcinoma

Code and results for: *Prognostic specificity of cancer-associated fibroblast signatures in lung adenocarcinoma: a benchmark against expression-matched random gene sets* (Wang Y et al., submitted to *Genes*).

## Contents
| Folder | Content |
|---|---|
| `preprocessing/A_bulk_OS/` | Cohort preparation: TCGA-LUAD (`03_baseA_tcga.R`), GEO download (`07_download_geo.R`) and preparation (`05_baseA_geo_rev3.R`; `rev2` kept for the record), CPTAC-LUAD (`06_baseA_cptac_rev2.R`), manifests (`08_harmonize_manifests.R`), the GSE31210 exclusion fix and `scoring_protocol.R`, which every analysis script sources (cohort loading, assertions, scoring rules) |
| `preprocessing/tcga_download/` | Download of TCGA-LUAD from the Genomic Data Commons (TCGAbiolinks; GDC Data Release 45.0, 4 December 2025) |
| `inputs/derived/` | Per-cohort clinical tables (`*_clinical_clean.csv`), manifests, build logs and `cohort_registry.csv`. The processed expression matrices are in the release asset (below). |
| `plan/` | Internal statistical analysis plan v1.2 (Chinese original, `.docx` and `.txt`) and File S1 (note, English translation, original) |
| `registry/` | Locked signature registry (`signature_registry_locked_20260903.csv`, MD5 `4ad3218931d95dfb6ee67a7ba07b9a75`), gene lists and published weights, symbol mapping, meta-PCNA genes, registry scripts. Source PDFs are not included (copyright). |
| `scani_genes/` | SCANI myCAF gene list (`S1_Table_myCAF_signature_50genes.csv`, from https://github.com/1379gener/scani-luad, tag v9.1) and the R object read by the scripts (`CAF_subtype_signatures.rds`) |
| `scripts/` | Benchmark (`run_benchmark.R`, `null_engine.R`), confounding models, readouts, figures, supplementary tables, verification |
| `results/` | All result files. Current supplementary tables: `supp_tables_20261111/`. Earlier dated files are kept because the deviation log refers to them. |
| `figures/` | Figures as submitted |
| `logs/` | Run logs |
| `docs/` | Working deviation log `DEVIATIONS.md` (mostly Chinese; the English version is Table S1) and the Methods-to-code correspondence table |
| `environment/` | `sessionInfo.txt` (R 4.3.1, Bioconductor 3.17, package versions, RNG kind) |

## Data sources
- GEO series matrices (downloaded with `preprocessing/A_bulk_OS/07_download_geo.R`): GSE68465 (GPL96), GSE72094 (GPL15048, with the GPL15048 annotation table), GSE31210, GSE50081, GSE30219 (GPL570).
- TCGA-LUAD: GDC STAR-Counts `tpm_unstrand`, Data Release 45.0. Overall survival from the GDC fields `vital_status`, `days_to_death` and `days_to_last_follow_up`.
- CPTAC-LUAD: LinkedOmics RNA-seq (`.cct`), used for correlations only.
- Tumour purity: ABSOLUTE estimates from Aran et al., Nat Commun 2015 (Supplementary Data 1, `41467_2015_BFncomms9971_MOESM1236_ESM.xlsx`), placed in `inputs/raw/purity/`.

## Rerunning
1. Download the release asset `SCI-2_processed_expression_20261111.tar` (check it against the `.md5` file) and extract it into `inputs/derived/` (`tar -xf SCI-2_processed_expression_20261111.tar -C inputs/derived`). Alternatively rebuild the matrices from the sources above with the scripts in `preprocessing/`.
2. Run `setup_layout.sh`. The scripts are kept exactly as run and use the paths of the original analysis server (`~/thesis_LUAD_CAF/...`, `~/caf_spatial_project/...`); the script recreates those paths as symlinks to this repository.
3. Run `scripts/run_benchmark.R`, then the scripts named in `docs/` and Table S1 in date order. `scripts/sci2_verify_nullengine.R` and `verify_ch04.R` recompute every reported z and empirical p value.

Null distributions are reproducible from fixed seeds: seed = 20260830 + (Σ position × UTF-8 code of "cohort|signature|null type") mod 100000, with R's default RNG (Mersenne-Twister, Inversion, Rejection).

## Integrity
`MD5SUMS.txt` lists MD5 checksums of all files in the repository.
