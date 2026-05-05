## =============================================================================
## FILE: deseq_functions.R
## PURPOSE: DESeq2 dataset construction, PCA helpers, DGE list annotation,
##          and volcano plot generation.
##
## DEPENDENCIES: 00_packages.R (DESeq2 1.46.0, pcaExplorer 3.0.0,
##                               limma 3.62.0, plotly 4.10.4,
##                               tidyverse 2.0.0, htmltools 0.5.8)
##               00_config.R   (DGE_PARAMS, VIS_PARAMS)
## =============================================================================


# ── 1. DESeq2 dataset builder ──────────────────────────────────────────────────

#' Build and run a DESeq2 analysis object.
#'
#' Automatically tests whether the batch covariate (StudyCenters) can be
#' included in the design formula (full-rank check), falls back to the
#' simpler design if not.
#'
#' @param raw_counts_filtered genes × samples count matrix; columns must
#'   match rownames of `samples`.
#' @param samples             Sample-sheet data.frame.
#' @param reference_sample    Level in `specific_v_all` to use as baseline.
#' @return Fitted DESeqDataSet object ready for results().
create_dds <- function(raw_counts_filtered, samples, reference_sample) {

  # Guard: sample IDs must align between counts and metadata
  if (!identical(colnames(raw_counts_filtered), rownames(samples))) {
    stop("Column names of raw_counts_filtered must match rownames of samples.")
  }

  if (!reference_sample %in% samples$specific_v_all) {
    stop("reference_sample '", reference_sample,
         "' not found in samples$specific_v_all.")
  }

  # Choose design formula based on rank of model matrix
  mm <- stats::model.matrix(~ StudyCenters + specific_v_all, data = samples)
  if (limma::is.fullrank(mm)) {
    message("Full-rank model matrix: including StudyCenters batch covariate.")
    design_formula <- ~ StudyCenters + specific_v_all
  } else {
    message("Rank-deficient model matrix: dropping StudyCenters from design.")
    design_formula <- ~ specific_v_all
  }

  # Construct DESeqDataSet (counts must be integer)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(raw_counts_filtered),
    colData   = samples,
    design    = design_formula
  )

  dds <- DESeq2::estimateSizeFactors(dds)

  # Set reference level for contrast computation
  dds$specific_v_all <- stats::relevel(dds$specific_v_all,
                                       ref = reference_sample)

  # Run full DESeq2 pipeline with multi-core support
  DESeq2::DESeq(dds, parallel = TRUE)
}


# ── 2. PCA data preparation ────────────────────────────────────────────────────

#' Run PCA on VST data and attach sample metadata.
#'
#' @param vst        genes × samples VST matrix (numeric).
#' @param samples    Sample-sheet data.frame (rows ordered identically to VST cols).
#' @param PCA_col    Column name from `samples` to attach to the PCA data frame.
#' @return Named list:
#'   \itemize{
#'     \item pcaobj  – prcomp result
#'     \item pca_df  – data.frame of PC coordinates + `PCA_col`
#'     \item pve     – numeric vector of % variance explained per PC
#'   }
create_PCA_data <- function(vst, samples, PCA_col, ...) {
  # Extra column names passed via ... are attached to the PCA data frame.
  # Allows batch_correction to attach Enrichement_Step_2 and StudyCenters.

  # prcomp expects samples-as-rows, genes-as-columns
  pcaobj <- stats::prcomp(t(vst))

  pca_df <- as.data.frame(pcaobj$x) %>%
    dplyr::mutate(!!PCA_col := samples[[PCA_col]])

  pve <- round(pcaobj$sdev^2 / sum(pcaobj$sdev^2) * 100, 2L)


  # Attach any additional metadata columns requested via ...
  extra_cols <- as.character(c(...))
  for (col in extra_cols) {
    if (col %in% colnames(samples)) {
      pca_df[[col]] <- samples[[col]]
    } else {
      warning("create_PCA_data: column '", col, "' not found in samples.")
    }
  }

  list(pcaobj = pcaobj, pca_df = pca_df, pve = pve)
}


# ── 3. Static 2-D PCA plot (ggplot2) ─────────────────────────────────────────

#' Generate a 2-D PCA scatter plot with group ellipses.
#'
#' @param pca_df   PCA data.frame from create_PCA_data()$pca_df.
#' @param pve      Percent variance explained vector.
#' @param PCA_col  Grouping column name (must exist in pca_df).
#' @param c1,c2    Hex/named colours for the two comparison groups.
#' @return ggplot object.
create_PCA_plot <- function(pca_df, pve, PCA_col, c1, c2) {
  ggplot2::ggplot(pca_df,
                  ggplot2::aes(PC1, PC2,
                               colour = .data[[PCA_col]])) +
    ggplot2::geom_point(size = 2L, alpha = 0.2) +
    ggplot2::geom_polygon(stat = "ellipse",
                          ggplot2::aes(fill = .data[[PCA_col]]),
                          alpha = 0.3) +
    ggplot2::scale_color_manual(values = c(c1, c2)) +
    ggplot2::scale_fill_manual(values  = c(c1, c2)) +
    ggplot2::theme_bw() +
    ggplot2::xlab(paste0("PC1 (", pve[1], "%)")) +
    ggplot2::ylab(paste0("PC2 (", pve[2], "%)")) +
    ggplot2::labs(colour = "Sample Type", fill = "Sample Type")
}


# ── 4. Interactive 3-D PCA plot (plotly) ──────────────────────────────────────

#' Generate an interactive 3-D PCA scatter plot.
#'
#' @param PCA_object     Output list from create_PCA_data().
#' @param sub_sample_type String label used in the plot title.
#' @param PCA_col        Column in pca_df used for colour grouping.
#' @return plotly object (can be saved with htmlwidgets::saveWidget).
create_PCA_plotly <- function(PCA_object, sub_sample_type, PCA_col) {

  # Sum variance explained by the first three PCs for the title
  total_var <- 100 * sum(
    summary(PCA_object$pcaobj)$importance["Proportion of Variance", 1:3]
  )
  title_text <- sprintf("%s  —  Total Explained Variance = %.2f%%",
                        sub_sample_type, total_var)

  # plotly uses formula NSE; build the color formula at runtime from PCA_col
  color_formula <- stats::as.formula(paste0("~", PCA_col))

  plotly::plot_ly(
    PCA_object$pca_df,
    x      = ~PC1, y = ~PC2, z = ~PC3,
    color  = color_formula,
    colors = c("red", "blue"),
    marker = list(size    = VIS_PARAMS$marker_size,
                  opacity = VIS_PARAMS$marker_opacity)
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


# ── 5. Annotated DGE list builder ─────────────────────────────────────────────

#' Annotate DESeq2 results and split into useful subsets.
#'
#' Steps:
#'   1. Drops rows with any NA (incomplete cases).
#'   2. Adds threshold flags and novel-gene indicator.
#'   3. Left-joins with Ensembl annotation.
#'   4. Returns five subsets for downstream use.
#'
#' @param results            as.data.frame(DESeq2::lfcShrink(...)) output.
#' @param annotation_ensembl Gene annotation data.frame; must have a `gene_id` column.
#' @param padj_threshold     Adjusted p-value cut-off (from DGE_PARAMS).
#' @param log2fold_threshold |log2FC| cut-off (from DGE_PARAMS).
#' @return Named list with elements:
#'   full, anno_only, filtered, filtered_novel, filtered_anno
create_annotated_DGElist <- function(results,
                                     annotation_ensembl,
                                     padj_threshold    = DGE_PARAMS$padj_threshold,
                                     log2fold_threshold = DGE_PARAMS$log2fold_threshold) {

  annotated <- results[stats::complete.cases(results), ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    dplyr::mutate(
      # TRUE if gene passes both statistical and effect-size thresholds
      threshold_pass = padj < padj_threshold &
                       abs(log2FoldChange) > log2fold_threshold,
      # Novel genes follow "MST…" naming convention
      Novel_gene     = grepl("^MST", gene_id)
    ) %>%
    dplyr::left_join(annotation_ensembl, by = "gene_id") %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))

  list(
    full           = annotated,
    anno_only      = dplyr::filter(annotated, !Novel_gene),
    filtered       = dplyr::filter(annotated,  threshold_pass),
    filtered_novel = dplyr::filter(annotated,  threshold_pass & Novel_gene),
    filtered_anno  = dplyr::filter(annotated,  threshold_pass & !Novel_gene)
  )
}


# ── 6. Volcano plot ────────────────────────────────────────────────────────────

#' Build a volcano plot from annotated DGE results.
#'
#' Points are coloured by the pre-computed `threshold_pass` flag.
#' Dashed reference lines mark ±2 log2FC and FDR = 0.05.
#'
#' @param DGE_list data.frame from create_annotated_DGElist()$full.
#' @return ggplot object.
create_volcano_plot <- function(DGE_list) {
  ggplot2::ggplot(DGE_list,
                  ggplot2::aes(log2FoldChange, -log10(padj),
                               name = gene_name)) +
    ggplot2::geom_point(ggplot2::aes(colour = threshold_pass),
                        size  = 1L,
                        alpha = 0.3) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black",
                                           "TRUE"  = "red")) +
    ggplot2::geom_vline(xintercept = c(-2, 2),
                        colour     = "green",
                        linetype   = 3L) +
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        colour     = "blue",
                        linetype   = 3L) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}
