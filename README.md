# Pan-Pediatric Cancer lncRNA Pipeline

Computational pipeline for bulk RNA-seq differential expression, functional enrichment, gene regulatory network analysis, and weighted gene co-expression network analysis (WGCNA) across 34 pediatric cancer types.

**PI:** Khan Lab, National Cancer Institute (NCI/NIH)  
**Analyst:** Abid Al Reza  
**HPC Environment:** NIH Biowulf (SLURM)  
**R version:** 4.4.2 | **Bioconductor:** 3.20

---

## Table of Contents

1. [Cohort Overview](#cohort-overview)
2. [Project Structure](#project-structure)
3. [Required Input Files](#required-input-files)
4. [Setup](#setup)
5. [Configuration](#configuration)
6. [Pipeline Execution](#pipeline-execution)
7. [SLURM Resource Allocation](#slurm-resource-allocation)
8. [Analysis Modules](#analysis-modules)
9. [Output Structure](#output-structure)
10. [Key Parameters](#key-parameters)
11. [Package Versions](#package-versions)
12. [Troubleshooting](#troubleshooting)

---

## Cohort Overview

| Attribute | Value |
|---|---|
| Total samples | 2,881 |
| Tumor types analyzed | 34 (NB.MYCN.Unk excluded; NL used as reference) |
| Normal (NL) samples | 147 |
| Library types | PolyA (1,358) · RiboZero (1,523) |
| Study centers | StJude (1,565) · NCI (880) · BeatCC (204) · BCCA (195) · CHOP (25) · MSKCC (12) |
| Annotation | GENCODE v42 · Ensembl v108 |
| Alignment/Quantification | RSEM |

**Tumor types (34 active):**

| Code | Tumor | n | Code | Tumor | n |
|---|---|---|---|---|---|
| GLI | Glioma (HGG/LGG/GNOS merged) | 360 | MBL.G3 | Medulloblastoma Group 3 | 41 |
| NB.MYCN.NA | Neuroblastoma MYCN non-amp | 243 | CPC | Choroid Plexus Carcinoma | 36 |
| WT | Wilms Tumor | 211 | MEL | Melanoma | 36 |
| OS | Osteosarcoma | 210 | WTB | Wilms Tumor Bilateral | 31 |
| EPMT | Ependymoma | 202 | ACC | Adrenocortical Carcinoma | 29 |
| EWS | Ewing Sarcoma | 178 | MBL.U | Medulloblastoma Unclassified | 29 |
| RMS.FN | RMS Fusion-Negative | 167 | THPA | Thyroid Papillary Carcinoma | 29 |
| RMS.FP | RMS Fusion-Positive | 101 | ASPS | Alveolar Soft Part Sarcoma | 21 |
| NB.MYCN.A | Neuroblastoma MYCN-amp | 99 | HGNET | High-Grade Neuroepithelial Tumor | 20 |
| RT | Rhabdoid Tumor | 74 | SPZM | Spitzoid Melanoma | 20 |
| MBL.G4 | Medulloblastoma Group 4 | 67 | CCSK | Clear Cell Sarcoma of Kidney | 19 |
| DSRCT | Desmoplastic Small Round Cell | 59 | MBL.WNT | Medulloblastoma WNT | 19 |
| MBL.SHH | Medulloblastoma SHH | 58 | PAST | Papillary Astrocytoma | 17 |
| HBL | Hepatoblastoma | 52 | GB | Glioblastoma | 14 |
| SS | Synovial Sarcoma | 50 | GNG | Ganglioglioma | 11 |
| ACPG | Adamantinomatous Craniopharyngioma | 44 | RMS.SS | RMS Sclerosing/Spindle | 7 |
| ATRT | Atypical Teratoid/Rhabdoid Tumor | 44 | RBL | Retinoblastoma | 42 |

---

## Project Structure

```
pan_pediatric_pipeline/
│
├── scripts/                        # Shared reusable function modules
│   ├── 00_packages.R               # Centralised package loader with pinned versions
│   ├── 00_config.R                 # All paths, filenames, and analysis parameters
│   ├── data_filtering_functions.R  # Sample subsetting, count filtering, annotation
│   ├── deseq_functions.R           # DESeq2 wrappers, PCA, DGE annotation, volcano
│   ├── functional_analysis_functions.R  # GO, KEGG, GSEA, Reactome, NCI enrichment
│   └── TF_enrichment_function.R    # Hypergeometric TF enrichment (unified)
│
├── analysis/                       # Numbered analysis drivers — run in order
│   ├── 1_library_type_batch_correction.R   # ComBat-seq library-type correction
│   ├── 2_analyze_samples_2025_02_21.R      # Per-tumour DESeq2 (one-vs-all)
│   ├── 3_TPM_calculation_RSEM_f.R          # TPM calculation + cross-validation
│   ├── 4_GO_analysis.R                     # GO / DO enrichment per tumour
│   ├── 5_TFs_enrichments_heatmap.R         # TF enrichment heatmap across tumours
│   ├── 6_GRN_TF_enrichment.R               # GENIE3 GRN TF enrichment per tumour
│   └── 7_WGCNA.R                           # Full WGCNA co-expression analysis
│
└── slurm/                          # HPC job submission scripts
    ├── run_sbatch.sh               # Master launcher (submit this to run everything)
    ├── analyze_samples_2025_02_21.sh   # DGE worker (called by run_sbatch.sh)
    ├── grn_tf_enrichment.sh            # GRN worker (called by run_sbatch.sh)
    └── wgcna.sh                        # WGCNA worker (called by run_sbatch.sh)
```

---

## Required Input Files

All input files must be present in `raw_data/` under `PROJECT_ROOT` before running. The project root is defined in `scripts/00_config.R`.

| File | Description |
|---|---|
| `sample_sheet_final_2024_04_12.csv` | Master sample metadata (2,881 rows × 30 columns) |
| `exp_count.txt` | Raw RSEM gene-level count matrix (genes × samples) |
| `exp_gene_length.txt` | RSEM effective gene lengths (matched to count matrix) |
| `exp_tpm.txt` | Pre-calculated RSEM TPM matrix (used for cross-validation) |
| `gencode_v42_annotation.genes.csv` | GENCODE v42 gene annotation |
| `raw_counts_data__post_bc.tsv` | Batch-corrected counts (output of step 1; input to steps 2, 7) |
| `NCI_GeneSet_v36.gmt` | NCI pathway gene sets in GMT format |
| `GRN_genie3_threshold_0.02.txt` | GENIE3 GRN edge table (threshold 0.02; for step 6) |
| `GRN_genie3_threshold_0.03.txt` | GENIE3 GRN edge table (threshold 0.03; for step 7) |
| `DatabaseExtract_v_1.01.xlsx` | TF database from humantfs.ccbr.utoronto.ca |
| `RSEM_no_filter_DGE_master_list.xlsx` | Master DGE list across all tumours (sheet: `vs_all_tumor`) |

> The file `raw_counts_data__post_bc.tsv` is generated by `1_library_type_batch_correction.R` and must exist before running any downstream steps.

---

## Setup

### 1. Clone / transfer the pipeline

```bash
# Transfer the pipeline folder to your project directory on Biowulf
scp -r pan_pediatric_pipeline/ \
  user@biowulf.nih.gov:/gpfs/gsfs10/users/khanlab2/lncRNA_Project___2024_02_29/
```

### 2. Edit the project root path

Open `scripts/00_config.R` and update the single `PROJECT_ROOT` variable to match your HPC path:

```r
PROJECT_ROOT <- "/gpfs/gsfs10/users/khanlab2/lncRNA_Project___2024_02_29/abid_WIP_2024_02_29/2024_04_12_RSEM"
```

This is the **only** path you should ever need to change. All other paths are derived from it.

### 3. Create the logs directory

```bash
mkdir -p /path/to/pan_pediatric_pipeline/logs
```

### 4. Install R packages (first run only)

Packages are installed automatically the first time `00_packages.R` is sourced. To pre-install manually on an interactive node:

```bash
sinteractive --mem=16g --cpus-per-task=4
module load R/4.4.2
Rscript -e "source('scripts/00_packages.R')"
```

If any Bioconductor package install fails due to write-protected system paths, install to a user library:

```r
# Inside R
.libPaths(c("~/R/4.4", .libPaths()))
BiocManager::install("PackageName", lib = "~/R/4.4")
```

---

## Configuration

All parameters live in `scripts/00_config.R`. The table below covers the most commonly adjusted settings.

### Paths

| Variable | Default | Purpose |
|---|---|---|
| `PROJECT_ROOT` | `/gpfs/.../2024_04_12_RSEM` | Single root for all paths |
| `PATHS$output` | `0__output_directory/` | DESeq2 and general outputs |
| `PATHS$wgcna_output` | `WGCNA/WGCNA_output/` | All WGCNA outputs |
| `PATHS$go_output` | `GoPlot_MDS/functionalAnalysis/` | GO enrichment plots |

### DESeq2 / DGE parameters

| Parameter | Default | Description |
|---|---|---|
| `DGE_PARAMS$padj_threshold` | `0.01` | BH-adjusted p-value cut-off |
| `DGE_PARAMS$log2fold_threshold` | `1` | \|log2FC\| cut-off (2-fold) |
| `DGE_PARAMS$lfc_shrink_type` | `"ashr"` | LFC shrinkage estimator |
| `DGE_PARAMS$min_count` | `10` | edgeR filterByExpr min counts |
| `DGE_PARAMS$min_prop` | `0.75` | Fraction of samples requiring min counts |

### WGCNA parameters

| Parameter | Default | Description |
|---|---|---|
| `WGCNA_PARAMS$power` | `12` | Soft-threshold power — **review `03_network_topology_analysis.pdf` and adjust before full run** |
| `WGCNA_PARAMS$network_type` | `"signed"` | Signed network (preserves sign of correlation) |
| `WGCNA_PARAMS$min_module_size` | `30` | Minimum genes per module |
| `WGCNA_PARAMS$merge_threshold` | `0.2` | ME dissimilarity cutoff for merging modules |
| `WGCNA_PARAMS$t_value_cutoff` | `5` | Min mixed-model t-value for module–tumour association |

### Enrichment parameters

| Parameter | Default | Description |
|---|---|---|
| `ENRICH_PARAMS$pvalue_cutoff` | `0.05` | ORA p-value threshold |
| `ENRICH_PARAMS$qvalue_cutoff` | `0.05` | ORA q-value threshold |
| `ENRICH_PARAMS$go_simplify_cutoff` | `0.7` | Wang semantic similarity for GO simplify() |
| `TF_PARAMS$cutoff` | `5` | Minimum TF-target / DEG overlap for hypergeometric test |
| `GRN_PARAMS$fdr_cutoff` | `0.05` | Bonferroni FDR cutoff for master TF enrichment table |

---

## Pipeline Execution

### Full automated run (recommended)

From the `slurm/` directory, a single command submits all stages with correct SLURM dependencies:

```bash
cd /path/to/pan_pediatric_pipeline
bash slurm/run_sbatch.sh
```

This will print a job ID summary:

```
[STAGE 2] Submitting DGE jobs (34 tumours, batch size 5)...
  Job 1 submitted (ID 12345678) | tumours: ACC,ACPG,ASPS,ATRT,CCSK
  ...
  Job 7 submitted (ID 12345684) | tumours: SS,THPA,WT,WTB
[STAGE 3] Submitting GRN TF enrichment (after DGE)...
  GRN job submitted (ID 12345685)
[STAGE 4] Submitting WGCNA (independent)...
  WGCNA job submitted (ID 12345686)
```

### Stage 1 prerequisite (run manually once)

Batch correction must be run before stages 2–4 as it produces `raw_counts_data__post_bc.tsv`:

```bash
sbatch --mem=256g --cpus-per-task=8 --gres=lscratch:50 --time=12:00:00 \
  --wrap="cd analysis && Rscript 1_library_type_batch_correction.R"
```

Or interactively:

```bash
sinteractive --mem=64g --cpus-per-task=8
module load R/4.4.2
cd analysis
Rscript 1_library_type_batch_correction.R
```

### Running individual stages manually

```bash
module load R/4.4.2
cd analysis

# Run any single script
Rscript 2_analyze_samples_2025_02_21.R "ACC,ATRT,EWS"   # DGE for 3 tumours
Rscript 3_TPM_calculation_RSEM_f.R                       # TPM (no arguments)
Rscript 4_GO_analysis.R
Rscript 5_TFs_enrichments_heatmap.R
Rscript 6_GRN_TF_enrichment.R
Rscript 7_WGCNA.R
```

### SLURM dependency diagram

```
Stage 1 (manual)
    └── 1_library_type_batch_correction.R
            │
            ▼
Stage 2 (7 parallel jobs)
    ├── DGE: ACC–CCSK
    ├── DGE: CPC–GB
    ├── DGE: GLI–MBL.G3
    ├── DGE: MBL.G4–MEL
    ├── DGE: NB.MYCN.A–RBL
    ├── DGE: RMS.FN–SPZM
    └── DGE: SS–WTB
            │  (afterok: all 7)
            ▼
Stage 3     6_GRN_TF_enrichment.R
            (4_GO_analysis.R and 5_TFs_enrichments_heatmap.R run separately)

Stage 4     7_WGCNA.R  ◄── independent, runs in parallel with stages 2–3
```

### SLURM resource allocation

Resource values in `slurm/run_sbatch.sh` were set as follows. **These are estimates, not empirically measured** — validate and tune them after your first successful run.

| Job | `--mem` | `--cpus-per-task` | `--time` | Basis |
|---|---|---|---|---|
| Stage 1 — batch correction | 256g | 8 | 12:00:00 | Copied from your original `run___Sbatch.sh` |
| Stage 2 — DGE (×7) | 256g | 8 | 12:00:00 | Copied from your original `run___Sbatch.sh` |
| Stage 3 — GRN enrichment | 16g | 2 | 02:00:00 | Pure R loop, no large matrices, no parallelism |
| Stage 4 — WGCNA | 256g | 16 | 24:00:00 | TOM scales O(genes²); doubled CPUs for WGCNA threads + DESeq2 BiocParallel |

#### How to measure actual usage after a run

```bash
# Replace <JOBID> with the job ID printed by run_sbatch.sh
sacct -j <JOBID> --format=JobID,MaxRSS,Elapsed,CPUTime
```

| Field | What it tells you |
|---|---|
| `MaxRSS` | Peak memory actually used — if well under `--mem`, reduce the allocation |
| `Elapsed` | Actual wall time — if consistently under limit, reduce `--time` |
| `CPUTime` | CPU × time product — if low relative to allocation, you may be over-requesting CPUs |

#### Tuning guidelines

- **DGE jobs (Stage 2):** Memory scales with genes retained after filtering and samples in the contrast. If `MaxRSS` is consistently below 64 GB, reduce `DGE_MEM` to `"64g"` in `run_sbatch.sh`.
- **WGCNA (Stage 4):** TOM construction is the bottleneck — memory is proportional to genes². With ~20,000 genes after filtering, 256 GB is the safe floor. Wall time of 24h is conservative; typical runs finish in 6–12h depending on cluster load.
- **CPUs vs. BiocParallel cores:** `VIS_PARAMS$n_cores` in `scripts/00_config.R` controls how many cores BiocParallel and WGCNA threads use. This value must always be ≤ `--cpus-per-task` in the SLURM script. If you reduce CPUs in the SLURM script, reduce `n_cores` in config to match, or jobs will be throttled by the scheduler.

#### Updating resource values

All resource variables are defined at the top of `slurm/run_sbatch.sh` in one place:

```bash
DGE_MEM="256g"; DGE_CPUS=8;  DGE_SCRATCH="lscratch:50";  DGE_TIME="12:00:00"
GRN_MEM="16g";  GRN_CPUS=2;  GRN_TIME="02:00:00"
WGC_MEM="256g"; WGC_CPUS=16; WGC_SCRATCH="lscratch:100"; WGC_TIME="24:00:00"
```

Edit those six lines only — do not edit the `sbatch` calls below them.


---

## Analysis Modules

### `1_library_type_batch_correction.R`

Corrects for library-type (PolyA vs RiboZero) batch effects using `sva::ComBat_seq`, which operates in count space preserving integer values. Generates pre- and post-correction VST matrices and 3-D PCA plots coloured by tumor type and shaped by library type or study centre.

**Key outputs:**  `raw_counts_data__post_bc.tsv`, `0_assay_vst.csv`, `0_post_bc_assay_vst.csv`, `1_pre_bc_PCA.html`, `1_post_bc_PCA.html`

---

### `2_analyze_samples_2025_02_21.R`

Runs one-versus-all-other-tumours DESeq2 for each tumor type. Accepts a comma-separated list of tumor codes as a command-line argument, enabling parallel SLURM execution across 7 jobs. Includes batch effect correction (StudyCenters) in the design formula when the model matrix is full-rank.

**Per-tumour outputs:**

| File | Description |
|---|---|
| `1_samples.csv` | Sample metadata subset |
| `2_assay_vst.csv` | VST-normalised expression |
| `3_PCA_varience.png` | Scree plot |
| `4_PC1_vs_PC2_using_vst.png` | Static 2-D PCA |
| `4_1_PCA_plotly.html` | Interactive 3-D PCA |
| `5_MA_plot.png` | MA plot (ashr-shrunk LFC) |
| `6_DGElist_full_with_lfc_shrinkage.csv` | All genes annotated |
| `7_DGElist_filtered_with_lfc_shrinkage.csv` | Significant DEGs only |
| `8_volcano_plot.html` | Interactive volcano |
| `9_volcano_plot.png` | Static volcano |
| `10_object_dimension.txt` | Dimension sanity check |

---

### `3_TPM_calculation_RSEM_f.R`

Calculates TPM two ways — manual formula and `DGEobj.utils::convertCounts()` — and cross-validates both against the pre-calculated RSEM TPM file. Exports TPM and raw count tables for a user-defined list of biologically queried genes.

---

### `4_GO_analysis.R`

Iterates over all 34 tumour types. For each, runs GO ORA (all ontologies), GO simplification (Wang semantic similarity), GO GSEA, and Disease Ontology (DO) enrichment on up- and down-regulated gene sets separately.

---

### `5_TFs_enrichments_heatmap.R`

Reads the master TF enrichment result file (produced by step 6) and renders a `-log(FDR)` ComplexHeatmap across all tumour types for TFs with overlap ≥ 25.

---

### `6_GRN_TF_enrichment.R`

For each tumour type, tests whether up-regulated DEGs (log2FC ≥ 2, padj < 0.05) are over-represented among the target genes of each TF in the GENIE3 network (threshold 0.02) using a hypergeometric test with Bonferroni correction. Writes per-tumour results and a master annotated table filtered to FDR < 0.05.

> Note: unlike the DEG-level TF analysis in step 5, this script does not require the TF to itself be a DEG (`require_tf_in_list = FALSE`).

---

### `7_WGCNA.R`

Full WGCNA pipeline in 19 numbered sections:

| Section | Description |
|---|---|
| 1 | Data loading and PAR_Y gene removal |
| 2 | QC and outlier detection on raw counts |
| 3 | VST normalisation (DESeq2 intercept-only model) |
| 4 | Soft-threshold power selection (scale-free topology) |
| 5 | TOM construction and dynamic tree cut module detection |
| 6 | Module merging by eigengene dissimilarity |
| 7–8 | Module membership calculation and intra-modular connectivity |
| 9 | Cytoscape edge/node export (pre- and post-merge) |
| 10 | Module–tumour Pearson correlation heatmap |
| 11 | Module–tumour mixed-model (lmerTest; StudyCenters as random effect) |
| 12 | Filtered module identification (t-value + BH adj-p cutoffs) |
| 13 | Module–module correlation (corrplot + MDS + dendrogram) |
| 14 | Module membership violin plot |
| 15–16 | Hub gene identification (top 1, top 5, top 25%) |
| 17 | lncRNA / non-lncRNA proportion per module |
| 18 | Per-module GO, KEGG, NCI, and TF enrichment |
| 19 | TF enrichment heatmaps for modules and hub genes |

> **Important:** After step 4 completes, inspect `03_network_topology_analysis.pdf` and confirm `WGCNA_PARAMS$power = 12` is appropriate for your data before a full run. Adjust in `scripts/00_config.R` if needed.

---

## Output Structure

```
0__output_directory/
├── 0_assay_vst.csv                     # Pre-BC VST matrix
├── 0_post_bc_assay_vst.csv             # Post-BC VST matrix
├── 1_pre_bc_PCA.html / _studycenters.html
├── 1_post_bc_PCA.html / _studycenters.html
└── DESeq2_analysis_2025_02_21/
    └── specific_v_all_other_tumors/
        ├── ACC/
        │   ├── 1_samples.csv ... 10_object_dimension.txt
        │   └── Grn-enrichments.ACC.*.txt
        ├── ATRT/  ...
        └── Genie3_GRN_enrichment_DEGs-up_enrichmentfdr-0.05.txt

GoPlot_MDS/functionalAnalysis/
└── <TUMOUR>/
    ├── up_genes/   1_DGE_up_genes.csv ... 15_Do_enrichment_dot_plot.png
    └── down_genes/ ...

WGCNA/
├── WGCNA_workspace_power_12.RData
└── WGCNA_output/
    ├── 01_sampleClustering_hclust_raw.pdf
    ├── 03_network_topology_analysis.pdf         ← review before full run
    ├── 05–07_module_detection_and_merging.pdf
    ├── 12_moduleTrait_correlation.txt
    ├── 14_module_significance_BH_mixed_model.txt
    ├── 14_2_filtered_module_t-values_mixed_model.png
    ├── 20_module_TF_heatmap.pdf / .png
    ├── cytoscape_inputs/oldMEs/ & newMEs/
    └── filtered_mod_memberships/
        └── <MODULE_COLOR>/
            ├── 1_DGE_list.csv
            ├── 2_GO_terms.csv
            ├── 3_KEGG_pathways.csv
            ├── 5_TF_DGE_list.csv
            └── hub_genes/ ...
```

---

## Key Parameters

The three parameters most likely to need adjustment for a new dataset:

1. **`PROJECT_ROOT`** in `scripts/00_config.R` — change once to relocate the whole pipeline.

2. **`WGCNA_PARAMS$power`** in `scripts/00_config.R` — run WGCNA step 4 first (soft-threshold plot), inspect `03_network_topology_analysis.pdf`, then set the power where `signed R² ≥ 0.80` before proceeding with the full run.

3. **`DGE_PARAMS$padj_threshold`** and **`DGE_PARAMS$log2fold_threshold`** — currently `0.01` and `1` respectively. The GRN enrichment script uses its own `GRN_PARAMS$fc_threshold = 2` and `GRN_PARAMS$alpha = 0.05` for DEG definition, which are intentionally more lenient to capture more TF targets.

---

## Package Versions

All packages are loaded and version-checked via `scripts/00_packages.R`. Versions were observed directly in SLURM output logs where available, and cross-referenced against the Bioconductor 3.20 / R 4.4.2 release manifest.

### CRAN

| Package | Version |
|---|---|
| tidyverse | 2.0.0 |
| dplyr | 1.1.4 |
| ggplot2 | 3.5.1 |
| readr | 2.1.5 |
| tibble | 3.2.1 |
| tidyr | 1.3.1 |
| purrr | 1.0.4 |
| forcats | 1.0.0 |
| lubridate | 1.9.4 |
| stringr | 1.5.1 |
| WGCNA | 1.73 |
| lmerTest | 3.1.3 |
| lme4 | 1.1.35 |
| corrplot | 0.92 |
| ComplexHeatmap | 2.22.0 |
| circlize | 0.4.16 |
| pheatmap | 1.0.12 |
| plotly | 4.10.4 |
| readxl | 1.4.3 |
| reshape2 | 1.4.4 |
| ashr | 2.2.63 |

### Bioconductor (release 3.20)

| Package | Version |
|---|---|
| clusterProfiler | 4.14.4 |
| DOSE | 4.0.0 |
| enrichplot | 1.26.6 |
| pcaExplorer | 3.0.0 |
| DESeq2 | 1.46.0 |
| edgeR | 4.4.0 |
| limma | 3.62.0 |
| sva | 3.54.0 |
| BiocParallel | 1.40.0 |
| biomaRt | 2.62.0 |
| org.Hs.eg.db | 3.20.0 |
| ReactomePA | 1.50.0 |

---

## Troubleshooting

**`source("../scripts/00_packages.R")` fails with "No such file"**  
Scripts must be run from the `analysis/` directory. Either `cd analysis` first, or call `Rscript analysis/2_analyze_samples_2025_02_21.R` from the project root only via the SLURM wrapper scripts which handle the `cd` automatically.

**biomaRt warning: "package version same as or greater than current"**  
This is a benign warning from BiocManager when `biomaRt` is already up to date. The package loads correctly.

**WGCNA TOM construction is very slow or runs out of memory**  
TOM construction scales as O(genes²). With ~20,000 genes and 500 samples, 256 GB RAM is required. If memory is insufficient, reduce the gene set further by tightening the `filterByExpr` parameters in `data_filtering_functions.R`.

**`lmerTest::lmer` convergence warnings in `7_WGCNA.R`**  
Mixed models for modules with very few samples in a tumour group may not converge. The warning is recorded but the pipeline continues — those entries will appear as `NA` in the significance table.

**SLURM job fails in Stage 3 with "DGE file not found"**  
The GRN enrichment script expects DGE files at `5_DGElist_full.csv` per tumour directory. Confirm all 34 Stage 2 DGE jobs completed successfully with `sacct -j <jobid> --format=JobID,State,ExitCode` before investigating.

**`getTRN_enrichment` returns an empty data frame for all tumours**  
Check that the TRN file columns are named `tfs` and `targets` after loading. The scripts rename columns `[1:2]` on load but column order in the source file must match.

---

*Pipeline maintained by Abid Al Reza · NCI/NIH · Axle Informatics contract*  
*Last updated: 2025-02-21*
