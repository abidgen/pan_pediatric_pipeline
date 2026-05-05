## =============================================================================
## FILE: 00_config.R
## PURPOSE: Single source of truth for all project paths and analysis
## =============================================================================

# ── Project root ───────────────────────────────────────────────────────────────
# Adjust this one line to move the whole project to a different location.
PROJECT_ROOT <- "/gpfs/gsfs10/users/khanlab2/lncRNA_Project___2024_02_29/abid_WIP_2024_02_29/2024_04_12_RSEM"

# ── Sub-directories ────────────────────────────────────────────────────────────
PATHS <- list(
  raw_data   = file.path(PROJECT_ROOT, "raw_data"),
  scripts    = file.path(PROJECT_ROOT, "scripts"),
  output     = file.path(PROJECT_ROOT, "0__output_directory"),
  dge_lists  = file.path(PROJECT_ROOT, "DGE_lists"),
  go_output  = file.path(PROJECT_ROOT, "GoPlot_MDS", "functionalAnalysis"),
  adhoc      = file.path(PROJECT_ROOT, "other_adhoc_analysis")
)

# ── Raw data files ─────────────────────────────────────────────────────────────
FILES <- list(
  sample_sheet       = file.path(PATHS$raw_data, "sample_sheet_final_2024_04_12.csv"),
  raw_counts         = file.path(PATHS$raw_data, "exp_count.txt"),
  raw_lengths        = file.path(PATHS$raw_data, "exp_gene_length.txt"),
  raw_tpm            = file.path(PATHS$raw_data, "exp_tpm.txt"),
  annotation_ensembl = file.path(PATHS$raw_data, "gencode_v42_annotation.genes.csv"),
  nci_geneset        = file.path(PATHS$raw_data, "NCI_GeneSet_v36.gmt"),
  post_bc_counts     = file.path(PATHS$raw_data, "raw_counts_data__post_bc.tsv")
)

# ── BioMart parameters ─────────────────────────────────────────────────────────
BIOMART <- list(
  ensembl_version = 108,    # Ensembl version for biomaRt queries
  dataset         = "hsapiens_gene_ensembl"
)

# ── DESeq2 / DGE thresholds ────────────────────────────────────────────────────
DGE_PARAMS <- list(
  padj_threshold      = 0.01,   # Adjusted p-value cut-off
  log2fold_threshold  = 1,      # |log2FC| cut-off (= 2-fold change)
  lfc_shrink_type     = "ashr", # LFC shrinkage method
  min_count           = 10,     # edgeR filterByExpr: min counts
  min_prop            = 0.75    # edgeR filterByExpr: min sample proportion
)

# ── GO / pathway enrichment thresholds ────────────────────────────────────────
ENRICH_PARAMS <- list(
  pvalue_cutoff   = 0.05,
  qvalue_cutoff   = 0.05,
  padj_method     = "fdr",
  show_categories = 10,
  go_simplify_cutoff = 0.7,
  do_pvalue_cutoff = 0.01,
  do_qvalue_cutoff = 0.20
)

# ── TF enrichment parameters ──────────────────────────────────────────────────
TF_PARAMS <- list(
  cutoff = 5,     # Minimum TF-target / DEG overlap
  FC     = 2,     # log2 fold-change threshold
  alpha  = 0.05   # p-value threshold
)

# ── Batch correction parameters ───────────────────────────────────────────────
BATCH_PARAMS <- list(
  batch_col = "Enrichement_Step_2",   # Column used as batch variable
  group_col = "Recoding_172"          # Biological group column
)

# ── PCA / visualisation ───────────────────────────────────────────────────────
VIS_PARAMS <- list(
  marker_size    = 5,
  marker_opacity = 0.8,
  n_cores        = 8    # BiocParallel cores (must match SLURM --cpus-per-task)
)

# ── Sample exclusions ─────────────────────────────────────────────────────────
EXCLUDED_SAMPLES <- c("NB.MYCN.Unk")

# ── Tumor subtype groups (used for sibling-exclusion logic) ──────────────────
GROUPED_SUBTYPES <- c("RMS", "NB", "MBL")


# ── WGCNA parameters ──────────────────────────────────────────────────────────
WGCNA_PARAMS <- list(
  network_type    = "signed",   # Network directionality
  power           = 12,         # Soft-threshold power (set after sft analysis)
  min_module_size = 30,         # Minimum genes per module (cutreeDynamic)
  deep_split      = 4,          # Tree cut sensitivity (0–4)
  merge_threshold = 0.2,        # MEDiss threshold for merging close modules
  t_value_cutoff  = 5,          # Min t-value to call a module-tumour association
  mm_cutoff       = 0.5         # Module membership value threshold
)

# ── WGCNA output paths ────────────────────────────────────────────────────────
PATHS$wgcna_root   <- file.path(PROJECT_ROOT, "WGCNA")
PATHS$wgcna_output <- file.path(PROJECT_ROOT, "WGCNA", "WGCNA_output")

# ── GRN / TF enrichment parameters ───────────────────────────────────────────
GRN_FILES <- list(
  grn_network  = file.path(PATHS$raw_data, "GRN_genie3_threshold_0.02.txt"),
  grn_wgcna    = file.path(PATHS$raw_data, "GRN_genie3_threshold_0.03.txt"),
  tf_db        = file.path(PATHS$raw_data, "DatabaseExtract_v_1.01.xlsx"),
  dge_master   = file.path(PATHS$dge_lists, "RSEM_no_filter_DGE_master_list.xlsx")
)

GRN_PARAMS <- list(
  fdr_cutoff      = 0.05,  # FDR threshold for master results table
  overlap_heatmap = 25,    # Min overlap to include TF in heatmap
  fc_threshold    = 2,     # log2FC cutoff for DEG definition in GRN script
  alpha           = 0.05   # p-value cutoff for DEG definition in GRN script
)

message("00_config.R loaded successfully.")

# ── Tumor-type color palette ───────────────────────────────────────────────────
# Named color vector matching all 35 Recoding_172 levels (34 tumours + NL).
# Used for PCA plots and WGCNA eigengene heatmaps. Sourced from RColorBrewer
# Set1/Set3 palettes to maximise visual separation across 35 categories.
MY_COLORS <- c(
  ACC="#E41A1C", ACPG="#377EB8", ASPS="#4DAF4A", ATRT="#984EA3",
  CCSK="#FF7F00", CPC="#A65628", DSRCT="#F781BF", EPMT="#999999",
  EWS="#66C2A5", GB="#FC8D62",  GLI="#8DA0CB", GNG="#E78AC3",
  HBL="#A6D854", HGNET="#FFD92F", MBL.G3="#E5C494", MBL.G4="#B3B3B3",
  MBL.SHH="#1B9E77", MBL.U="#D95F02", MBL.WNT="#7570B3", MEL="#E7298A",
  NB.MYCN.A="#66A61E", NB.MYCN.NA="#E6AB02", NB.MYCN.Unk="#A6761D",
  OS="#666666", PAST="#8DD3C7", RBL="#FFFFB3", RMS.FN="#BEBADA",
  RMS.FP="#FB8072", RMS.SS="#80B1D3", RT="#FDB462", SPZM="#B3DE69",
  SS="#FCCDE5", THPA="#D9D9D9", WT="#BC80BD", WTB="#CCEBC5",
  NL="#000000"
)
