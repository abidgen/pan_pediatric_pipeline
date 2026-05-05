## =============================================================================
## FILE: TFs_enrichments_heatmap.R
## PURPOSE: Visualise TF enrichment results across all tumor types as a
##          ComplexHeatmap where cells = -log(FDR).
##
## INPUT:  Genie3_GRN_enrichment_DEGs-up_enrichmentfdr-0.05.txt
##           – Output of getTRN_enrichment() merged across all comparisons.
##           Expected columns: TF (or HGNC.symbol), set, fdr, overlap.
##
## DEPENDENCIES: 00_packages.R (ComplexHeatmap 2.22.0, circlize 0.4.16,
##                               reshape2 1.4.4, plyr 1.8.9)
##               00_config.R
## =============================================================================

rm(list = ls())

source("../scripts/00_packages.R")
source("../scripts/00_config.R")


# ── Parameters ─────────────────────────────────────────────────────────────────
INPUT_FILE   <- file.path(
  PATHS$output, "DESeq2_analysis", "specific_v_all_other_tumors",
  "Genie3_GRN_enrichment_DEGs-up_enrichmentfdr-0.05.txt"
)
MIN_OVERLAP  <- 25L     # Only keep TF rows with overlap >= this threshold
MIN_FREQ     <- 1L      # Minimum number of tumor types a TF must appear in
EPSILON      <- 1e-100  # Replace exact-zero FDR to avoid log(0)


# ── Load and filter TF enrichment results ─────────────────────────────────────
TFs_res  <- utils::read.delim(INPUT_FILE)
TFs_res  <- TFs_res[TFs_res$overlap >= MIN_OVERLAP, ]

# Optional: further filter to TFs enriched in at least MIN_FREQ tumor types
TFs_freq <- as.data.frame(table(TFs_res$TF))
colnames(TFs_freq) <- c("TF", "Freq")
# Uncomment the line below to apply frequency filter:
# TFs_res <- TFs_res[TFs_res$TF %in% TFs_freq$TF[TFs_freq$Freq >= MIN_FREQ], ]


# ── Reshape to TF × tumor-type FDR matrix ────────────────────────────────────
# acast produces a wide matrix: rows = TF (HGNC.symbol), cols = set (tumor type)
heatmap_mat <- reshape2::acast(TFs_res, HGNC.symbol ~ set, value.var = "fdr")

# Replace NA (TF not enriched in that tumor) with 1 (no enrichment)
heatmap_mat[is.na(heatmap_mat)] <- 1

# Clamp exact zeros to avoid -Inf after log transformation
heatmap_mat[heatmap_mat == 0] <- EPSILON

# Transform to -log(FDR): higher = more significant
heatmap_mat <- -log(heatmap_mat)


# ── Colour scale ───────────────────────────────────────────────────────────────
# White → red → dark red maps increasing significance
col_fun <- circlize::colorRamp2(
  breaks = c(0, 50, 100),
  colors = c("white", "red", "darkred")
)


# ── Draw heatmap ───────────────────────────────────────────────────────────────
ht <- ComplexHeatmap::Heatmap(
  matrix               = heatmap_mat,
  col                  = col_fun,
  name                 = "-log(FDR)",
  row_title            = "Transcription Factors",
  column_title         = "Tumor Types",
  row_names_gp         = grid::gpar(fontsize = 8),
  column_names_gp      = grid::gpar(fontsize = 8),
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  show_row_dend        = TRUE,
  show_column_dend     = TRUE,
  heatmap_legend_param = list(title = "-log(FDR)")
)

ComplexHeatmap::draw(ht)
