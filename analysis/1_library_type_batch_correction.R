## =============================================================================
## FILE: 1_library_type_batch_correction.R
## PURPOSE: Two-stage batch correction for Pan-Pediatric Cancer RNA-seq data.
##   Stage 1 – Remove library-type (PolyA vs RiboZero) batch effect via
##             sva::ComBat_seq (count-space correction).
##   Stage 2 – Downstream DESeq2 / VST / PCA on pre- and post-corrected data.
##
## DEPENDENCIES: 00_packages.R (limma 3.62.0, edgeR 4.4.0, sva 3.54.0,
##                               DESeq2 1.46.0, tidyverse 2.0.0,
##                               plotly 4.10.4, htmlwidgets 1.6.4)
##               00_config.R, data_filtering_functions.R, deseq_functions.R
## =============================================================================

rm(list = ls())

# ── Load shared modules ────────────────────────────────────────────────────────
source("../scripts/00_packages.R")
source("../scripts/00_config.R")
source("../scripts/data_filtering_functions.R")
source("../scripts/deseq_functions.R")

setwd(PROJECT_ROOT)


# ── Helper: 3-D PCA coloured by tumour type, shaped by library type ───────────
#' @param PCA_object  Output of create_PCA_data().
#' @param label       String appended to the plot title.
.plot_pca_enrichment <- function(PCA_object, label) {
  total_var <- 100 * sum(
    summary(PCA_object$pcaobj)$importance["Proportion of Variance", 1:3]
  )
  title_text <- sprintf("%s  —  Total Explained Variance = %.2f%%",
                        label, total_var)

  plotly::plot_ly(
    PCA_object$pca_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color  = ~Recoding_172,
    colors = MY_COLORS,
    marker = list(size    = VIS_PARAMS$marker_size,
                  opacity = VIS_PARAMS$marker_opacity)
  ) %>%
    plotly::add_markers(
      symbol  = ~Enrichement_Step_2,
      symbols = c("circle", "cross")
    ) %>%
    plotly::layout(
      legend = list(title = list(text = "Tumor Types")),
      title  = title_text,
      scene  = list(
        bgcolor = "white",
        xaxis   = list(title = paste0("PC1 (", PCA_object$pve[1], "%)")),
        yaxis   = list(title = paste0("PC2 (", PCA_object$pve[2], "%)")),
        zaxis   = list(title = paste0("PC3 (", PCA_object$pve[3], "%)"))
      )
    )
}

# Helper: 3-D PCA shaped by study centre instead of library type
.plot_pca_studycentre <- function(PCA_object, label) {
  total_var <- 100 * sum(
    summary(PCA_object$pcaobj)$importance["Proportion of Variance", 1:3]
  )
  plotly::plot_ly(
    PCA_object$pca_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color  = ~Recoding_172,
    colors = MY_COLORS,
    marker = list(size    = VIS_PARAMS$marker_size,
                  opacity = VIS_PARAMS$marker_opacity)
  ) %>%
    plotly::add_markers(
      symbol  = ~StudyCenters,
      symbols = c("square-open", "circle-open", "circle",
                  "diamond", "diamond-open", "cross")
    ) %>%
    plotly::layout(
      legend = list(title = list(text = "Tumor Types")),
      title  = sprintf("%s  —  %.2f%% variance", label, total_var),
      scene  = list(
        bgcolor = "white",
        xaxis   = list(title = paste0("PC1 (", PCA_object$pve[1], "%)")),
        yaxis   = list(title = paste0("PC2 (", PCA_object$pve[2], "%)")),
        zaxis   = list(title = paste0("PC3 (", PCA_object$pve[3], "%)"))
      )
    )
}


# ── Create output directory ────────────────────────────────────────────────────
output_dir <- file.path(PATHS$output)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ── Load sample sheet ─────────────────────────────────────────────────────────
samples <- utils::read.delim(
  FILES$sample_sheet,
  header      = TRUE,
  sep         = ",",
  check.names = FALSE,
  row.names   = 1L
)


# ── Load raw counts ───────────────────────────────────────────────────────────
raw_counts_data <- utils::read.table(
  FILES$raw_counts,
  header    = TRUE,
  row.names = 1L,
  sep       = "\t",
  encoding  = "ISO-8859-1"
)

# Standardise column names to match sample sheet (e.g. X123.456 → 123-456)
sample_names <- gsub("[.]", "-", colnames(raw_counts_data))
sample_names <- gsub("^X", "",  sample_names)
colnames(raw_counts_data) <- sample_names

# Guard: all count columns should be present in the sample sheet
stopifnot(sum(!(colnames(raw_counts_data) %in% samples$sample_names)) == 0L)

# Re-order count columns to match sample sheet row order
raw_counts_data <- raw_counts_data[, samples$sample_names]
stopifnot(all(colnames(raw_counts_data) == samples$sample_names))

saveRDS(raw_counts_data, file.path(PATHS$raw_data, "raw_counts_data.rds"))


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 1 – PRE-BATCH-CORRECTION
# ═══════════════════════════════════════════════════════════════════════════════

# ── Low-expression filtering (edgeR filterByExpr) ────────────────────────────
pre_y <- edgeR::DGEList(counts = raw_counts_data, samples = samples)
keep  <- edgeR::filterByExpr(
  pre_y,
  group     = pre_y$samples$Recoding_172,
  min.count = DGE_PARAMS$min_count,
  min.prop  = DGE_PARAMS$min_prop
)
pre_y <- pre_y[keep, ]
message("Pre-BC: ", nrow(pre_y$counts), " genes retained after filterByExpr.")

# ── DESeq2 object (pre-batch-correction) ────────────────────────────────────
dds_pre <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(pre_y$counts),
  colData   = pre_y$samples,
  design    = ~ Recoding_172
)
dds_pre <- DESeq2::estimateSizeFactors(dds_pre)
dds_pre$Recoding_172 <- stats::relevel(dds_pre$Recoding_172, ref = "NL")
dds_pre <- DESeq2::DESeq(dds_pre, parallel = TRUE)
saveRDS(dds_pre, file.path(PATHS$raw_data, "pre_bc_vs_normal_dds.rds"))

# ── VST (pre-batch-correction) ────────────────────────────────────────────────
pre_vst_file <- file.path(output_dir, "0_assay_vst.csv")
pre_vst_mat  <- DESeq2::varianceStabilizingTransformation(dds_pre)

utils::write.csv(
  as.data.frame(DESeq2::assay(pre_vst_mat)) %>%
    dplyr::mutate_if(is.numeric, round_2),
  pre_vst_file
)
# Re-read for exact numeric reproducibility (avoids floating-point drift)
pre_vst <- utils::read.delim(pre_vst_file,
                              header = TRUE, sep = ",",
                              check.names = FALSE, row.names = 1L)

# ── PCA (pre-batch-correction) ────────────────────────────────────────────────
pre_bc_PCA <- create_PCA_data(pre_vst, samples,
                               "Recoding_172", "Enrichement_Step_2", "StudyCenters")

htmlwidgets::saveWidget(
  .plot_pca_enrichment(pre_bc_PCA, "Pre-BC — all samples"),
  file.path(output_dir, "1_pre_bc_PCA.html"),
  selfcontained = TRUE
)
htmlwidgets::saveWidget(
  .plot_pca_studycentre(pre_bc_PCA, "Pre-BC — study centres"),
  file.path(output_dir, "1_pre_bc_PCA_studycenters.html"),
  selfcontained = TRUE
)


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 2 – LIBRARY-TYPE BATCH CORRECTION (ComBat_seq)
# ═══════════════════════════════════════════════════════════════════════════════

# Cast factor columns before ComBat_seq
pre_y$samples$Enrichement_Step_2 <- as.factor(pre_y$samples$Enrichement_Step_2)
pre_y$samples$Recoding_172       <- as.factor(pre_y$samples$Recoding_172)

# ComBat_seq operates in count space and preserves integer counts
raw_counts_post_bc <- sva::ComBat_seq(
  counts = as.matrix(pre_y$counts),
  batch  = pre_y$samples[[BATCH_PARAMS$batch_col]],  # Library type
  group  = pre_y$samples[[BATCH_PARAMS$group_col]]   # Tumour type
)
raw_counts_post_bc <- as.data.frame(raw_counts_post_bc)
stopifnot(all(colnames(raw_counts_post_bc) == rownames(samples)))

# Write batch-corrected counts for downstream scripts
readr::write_tsv(
  raw_counts_post_bc %>% tibble::rownames_to_column(var = "gene_ID"),
  FILES$post_bc_counts
)


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 3 – POST-BATCH-CORRECTION DESEQ2 / VST / PCA
# ═══════════════════════════════════════════════════════════════════════════════

# ── Low-expression filtering (post-BC, relaxed: min.count = 1) ──────────────
post_y <- edgeR::DGEList(counts = raw_counts_post_bc, samples = samples)
keep   <- edgeR::filterByExpr(
  post_y,
  group     = post_y$samples$Recoding_172,
  min.count = 1L,
  min.prop  = DGE_PARAMS$min_prop
)
post_y <- post_y[keep, ]
message("Post-BC: ", nrow(post_y$counts), " genes retained after filterByExpr.")
gc()

# ── DESeq2 object (post-batch-correction) ────────────────────────────────────
dds_post <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(post_y$counts),
  colData   = post_y$samples,
  design    = ~ Recoding_172
)
dds_post <- DESeq2::estimateSizeFactors(dds_post)
dds_post$Recoding_172 <- stats::relevel(dds_post$Recoding_172, ref = "NL")
gc()
dds_post <- DESeq2::DESeq(dds_post, parallel = TRUE)
saveRDS(dds_post, file.path(PATHS$raw_data, "post_bc_vs_normal_dds.rds"))

# ── VST (post-batch-correction) ───────────────────────────────────────────────
gc()
post_vst_file <- file.path(output_dir, "0_post_bc_assay_vst.csv")
post_vst_mat  <- DESeq2::varianceStabilizingTransformation(dds_post)

utils::write.csv(
  as.data.frame(DESeq2::assay(post_vst_mat)) %>%
    dplyr::mutate_if(is.numeric, round_2),
  post_vst_file
)
post_bc_vst <- utils::read.delim(post_vst_file,
                                  header = TRUE, sep = ",",
                                  check.names = FALSE, row.names = 1L)

# ── PCA (post-batch-correction) ───────────────────────────────────────────────
gc()
post_bc_PCA <- create_PCA_data(post_bc_vst, samples,
                                "Recoding_172", "Enrichement_Step_2", "StudyCenters")

htmlwidgets::saveWidget(
  .plot_pca_enrichment(post_bc_PCA, "Post-BC — all samples"),
  file.path(output_dir, "1_post_bc_PCA.html"),
  selfcontained = TRUE
)
htmlwidgets::saveWidget(
  .plot_pca_studycentre(post_bc_PCA, "Post-BC — study centres"),
  file.path(output_dir, "1_post_bc_PCA_studycenters.html"),
  selfcontained = TRUE
)

message("Batch correction complete. Outputs written to: ", output_dir)
