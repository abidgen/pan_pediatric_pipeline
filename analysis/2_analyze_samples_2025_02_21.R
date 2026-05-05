## =============================================================================
## FILE: analyze_samples_2025_02_21.R
## PURPOSE: Main analysis driver — one-vs-all DESeq2 differential expression
##          for each tumor type against all other tumors.
##
## USAGE (via SLURM):
##   Rscript analyze_samples_2025_02_21.R "ACC,ACPG,ATRT"
##   The first argument is a comma-separated list of tumor-type labels.
##
## OUTPUTS (per tumor type, written to output_dir_master/specific_v_all_other_tumors/<TUMOR>):
##   1_samples.csv                        – Sample metadata subset
##   2_assay_vst.csv                      – VST-normalised expression
##   3_PCA_varience.png                   – Scree plot
##   4_PC1_vs_PC2_using_vst.png           – Static 2-D PCA
##   4_1_PCA_plotly.html                  – Interactive 3-D PCA
##   5_MA_plot.png                        – MA plot (LFC shrinkage)
##   6_DGElist_full_with_lfc_shrinkage.csv
##   7_DGElist_filtered_with_lfc_shrinkage.csv
##   8_volcano_plot.html                  – Interactive volcano
##   9_volcano_plot.png                   – Static volcano
##   10_object_dimension.txt              – Dimension sanity check
##
## DEPENDENCIES: 00_packages.R, 00_config.R,
##               data_filtering_functions.R, deseq_functions.R
## =============================================================================

rm(list = ls())

# ── Load shared modules ────────────────────────────────────────────────────────
source("../scripts/00_packages.R")
source("../scripts/00_config.R")
source("../scripts/data_filtering_functions.R")
source("../scripts/deseq_functions.R")
setwd(PROJECT_ROOT)          # Defined in 00_config.R


# ── Parameters ─────────────────────────────────────────────────────────────────
# Reference / comparison contrasts; extend the list to run more comparisons.
n_v_a_set <- c("all_other_tumors")   # Active contrast directions
n_v_a <- list(
  Normal           = c("Normal",          "all_other_tumors"),
  all_other_tumors = c("all_other_tumors", "Normal")
)


# ── Create master output directory ────────────────────────────────────────────
output_dir_master <- file.path(PATHS$output, "DESeq2_analysis_2025_02_21")
if (!dir.exists(output_dir_master)) {
  dir.create(output_dir_master, recursive = TRUE)
}
message("Output root: ", output_dir_master)


# ── Load annotation ────────────────────────────────────────────────────────────
annotation_ensembl  <- utils::read.csv(FILES$annotation_ensembl)
annotation_enhanced <- create_annotation_enhanced()  # biomaRt query
annotation_enhanced <- dplyr::left_join(
  annotation_ensembl,
  annotation_enhanced,
  by = c("gene_id" = "ensembl_gene_id_version")
)
rm(annotation_ensembl)


# ── Load and clean sample sheet ───────────────────────────────────────────────
samples_master <- utils::read.delim(
  FILES$sample_sheet,
  header     = TRUE,
  sep        = ",",
  check.names = FALSE,
  row.names  = 1L
)
# Drop ambiguous NB.MYCN.Unk samples (confirmed exclusion in 00_config.R)
samples_master <- samples_master[
  !samples_master$Recoding_172 %in% EXCLUDED_SAMPLES, ]


# ── Load batch-corrected raw counts ───────────────────────────────────────────
raw_counts <- utils::read.delim(
  FILES$post_bc_counts,
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE,
  row.names   = 1L
)
# Pre-filter to samples in master sheet
# DGE_PARAMS$samples_col is not set in config; default to "Recoding_172"
samples_col    <- "Recoding_172"
raw_counts_data <- filter_raw_counts(raw_counts, samples_master, samples_col)


# ── Parse tumor-type argument from SLURM command line ─────────────────────────
args        <- base::commandArgs(trailingOnly = TRUE)
samples_set <- base::strsplit(args[1L], ",")[[1L]]
message("Processing tumor types: ", paste(samples_set, collapse = ", "))


# ── Main analysis loop ────────────────────────────────────────────────────────
for (contrast_dir in n_v_a_set) {

  reference_sample <- n_v_a[[contrast_dir]][1L]
  exclude          <- n_v_a[[contrast_dir]][2L]

  for (specific_sample in samples_set) {

    message("\n>>> Tumor: ", specific_sample,
            "  |  Reference: ", reference_sample)

    # ── Subset samples ──────────────────────────────────────────────────────
    samples <- specific_v_reference(
      data            = samples_master,
      samples_col     = "Recoding_172",
      specific_sample = specific_sample,
      exclude         = exclude
    )

    # ── Subset and filter raw counts ────────────────────────────────────────
    raw_counts_filtered <- filter_raw_counts(
      raw_counts  = raw_counts_data,
      samples     = samples,
      samples_col = "Recoding_172"
    )

    # ── Create per-tumor output directory ───────────────────────────────────
    output_dir <- file.path(
      output_dir_master,
      paste0("specific_v_", reference_sample),
      specific_sample
    )
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # ── Write sample metadata ────────────────────────────────────────────────
    utils::write.csv(as.data.frame(samples),
                     file.path(output_dir, "1_samples.csv"))

    # ── DESeq2 ──────────────────────────────────────────────────────────────
    dds <- create_dds(raw_counts_filtered, samples, reference_sample)

    # ── VST normalisation ────────────────────────────────────────────────────
    vst_mat <- DESeq2::assay(DESeq2::varianceStabilizingTransformation(dds))
    utils::write.csv(
      as.data.frame(vst_mat) %>% dplyr::mutate_if(is.numeric, round_2),
      file.path(output_dir, "2_assay_vst.csv")
    )

    # ── PCA ──────────────────────────────────────────────────────────────────
    PCA <- create_PCA_data(vst_mat, samples, "specific_v_all")

    # Scree plot
    png(file.path(output_dir, "3_PCA_varience.png"),
        units = "in", width = 15, height = 4, res = 300)
    print(pcaExplorer::pcascree(
      PCA$pcaobj,
      type  = "pev",
      title = "Proportion of explained variance"
    ))
    dev.off()

    # Static PC1 vs PC2
    png(file.path(output_dir, "4_PC1_vs_PC2_using_vst.png"),
        units = "in", width = 5, height = 3.5, res = 300)
    print(create_PCA_plot(PCA$pca_df, PCA$pve, "specific_v_all", "red", "blue"))
    dev.off()

    # Interactive 3-D PCA
    htmlwidgets::saveWidget(
      create_PCA_plotly(PCA, specific_sample, "specific_v_all"),
      file.path(output_dir, "4_1_PCA_plotly.html"),
      selfcontained = TRUE
    )

    # ── DESeq2 results + LFC shrinkage (ashr) ────────────────────────────────
    results_raw <- DESeq2::results(
      dds,
      contrast = c("specific_v_all", specific_sample, reference_sample),
      alpha    = 0.05,
      parallel = TRUE
    )
    results_shrunk <- DESeq2::lfcShrink(
      dds,
      contrast = c("specific_v_all", specific_sample, reference_sample),
      type     = DGE_PARAMS$lfc_shrink_type,
      alpha    = 0.05,
      parallel = TRUE,
      res      = results_raw
    )

    # MA plot (uses shrunk estimates)
    png(file.path(output_dir, "5_MA_plot.png"),
        units = "in", width = 10, height = 10, res = 300)
    DESeq2::plotMA(results_shrunk, ylim = c(-2, 2))
    dev.off()

    # ── Annotate and write DGE list ──────────────────────────────────────────
    DGE <- create_annotated_DGElist(
      results             = as.data.frame(results_shrunk),
      annotation_ensembl  = annotation_enhanced,
      padj_threshold      = DGE_PARAMS$padj_threshold,
      log2fold_threshold  = DGE_PARAMS$log2fold_threshold
    )
    utils::write.csv(DGE$full,
                     file.path(output_dir, "6_DGElist_full_with_lfc_shrinkage.csv"))
    utils::write.csv(DGE$filtered,
                     file.path(output_dir, "7_DGElist_filtered_with_lfc_shrinkage.csv"))

    # ── Volcano plots ────────────────────────────────────────────────────────
    volcano_plot <- create_volcano_plot(DGE$full)

    htmltools::save_html(
      plotly::ggplotly(volcano_plot),
      file = file.path(output_dir, "8_volcano_plot.html")
    )
    png(file.path(output_dir, "9_volcano_plot.png"),
        units = "in", width = 4, height = 3.5, res = 300)
    print(volcano_plot)
    dev.off()

    # ── Object dimension log ─────────────────────────────────────────────────
    files_list <- list(
      samples             = samples,
      raw_counts          = raw_counts,
      raw_counts_filtered = raw_counts_filtered,
      vst                 = as.data.frame(vst_mat),
      DGE_full            = DGE$full,
      DGE_anno_only       = DGE$anno_only,
      DGE_filtered        = DGE$filtered
    )
    sink(file.path(output_dir, "10_object_dimension.txt"))
    print(lapply(files_list, find_dim))
    sink()

    message("<<< Done: ", specific_sample)
  }
}
