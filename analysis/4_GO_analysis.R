## =============================================================================
## FILE: GO_analysis.R
## PURPOSE: Per-tumor GO enrichment (ORA + GSEA) and Disease Ontology (DO)
##          analysis.  Iterates over every tumor type defined in the sample
##          sheet, processing up-regulated and down-regulated gene sets
##          separately.
##
## OUTPUTS (per tumor, per direction, under PATHS$go_output):
##   1_DGE_<direction>.csv
##   2_enriched_GO_terms.csv          – raw enrichGO
##   3_enriched_GO_terms_simplified.csv
##   4/5_go_enrichment_{bar,dot}_plot.png
##   6/7_go_enrichment_simplified_{bar,dot}_plot.png
##   8_GSEA_Go.csv
##   9_GSEA_Go.png
##   10_GSEA_Go_simplified.csv
##   11_GSEA_Go_simplified.png
##   12_enriched_DO_terms.csv
##   14/15_Do_enrichment_{bar,dot}_plot.png
##
## DEPENDENCIES: 00_packages.R (clusterProfiler 4.14.4, DOSE 4.0.0,
##                               enrichplot 1.26.6, org.Hs.eg.db 3.20.0,
##                               BiocParallel 1.40.0, tidyverse 2.0.0)
##               00_config.R, data_filtering_functions.R,
##               functional_analysis_functions.R
## =============================================================================

rm(list = ls())

source("../scripts/00_packages.R")
source("../scripts/00_config.R")
source("../scripts/data_filtering_functions.R")
source("../scripts/functional_analysis_functions.R")

setwd(PROJECT_ROOT)


# ── GO simplification parameters ──────────────────────────────────────────────
# clusterProfiler::simplify() removes semantically redundant GO terms.
SIMPLIFY_ARGS <- list(
  cutoff     = ENRICH_PARAMS$go_simplify_cutoff,
  by         = "p.adjust",
  select_fun = min,
  measure    = "Wang",
  semData    = NULL
)


# ── Load sample sheet ─────────────────────────────────────────────────────────
samples_master <- utils::read.delim(
  FILES$sample_sheet,
  header      = TRUE,
  sep         = ",",
  check.names = FALSE,
  row.names   = 1L
)

# Build the set of tumor types to process (exclude Normal)
samples_set <- unique(samples_master$Recoding_172)
samples_set <- samples_set[samples_set != "NL"]


# ── Main loop: one iteration per tumor type ────────────────────────────────────
for (sample in samples_set) {

  message("\n>>> GO analysis: ", sample)

  # ── Create output directory ────────────────────────────────────────────────
  output_dir_sample <- file.path(PATHS$go_output, sample)
  if (!dir.exists(output_dir_sample)) {
    dir.create(output_dir_sample, recursive = TRUE)
  }

  # ── Load DGE results for this tumor type ──────────────────────────────────
  dge_dir <- file.path(PATHS$output, "DESeq2_analysis",
                       "specific_v_all_other_tumors", sample)

  DGE_filtered_anno <- utils::read.delim(
    file.path(dge_dir, "6_DGElist_filtered_with_lfc_shrinkage.csv"),
    header = TRUE, sep = ",", check.names = FALSE, row.names = 1L
  )
  DGE_full <- utils::read.delim(
    file.path(dge_dir, "5_DGElist_full_with_lfc_shrinkage.csv"),
    header = TRUE, sep = ",", check.names = FALSE, row.names = 1L
  )

  # Background universe: all ENTREZ IDs from the full DGE table
  ent_uni            <- as.character(DGE_full$entrezgene_id)
  ent_uni_gene_names <- DGE_full$gene_name

  # Split into up / down gene sets
  DGE_func <- create_func_DGE(DGE_filtered_anno)

  # ── Inner loop: up-regulated and down-regulated sets ──────────────────────
  for (direction in c("up_genes", "down_genes")) {

    output_dir_0 <- file.path(output_dir_sample, direction)
    if (!dir.exists(output_dir_0)) dir.create(output_dir_0, recursive = TRUE)

    # Write gene list
    utils::write.csv(DGE_func[[direction]],
                     file.path(output_dir_0, paste0("1_DGE_", direction, ".csv")))

    ent_gene <- as.character(DGE_func[[direction]]$entrezgene_id)

    # ── GO over-representation analysis ──────────────────────────────────
    eGO     <- try(create_eGO(ent_gene, ent_uni))
    try(utils::write.csv(eGO, file.path(output_dir_0, "2_enriched_GO_terms.csv")))

    eGoSimp <- try(do.call(clusterProfiler::simplify,
                            c(list(x = eGO), SIMPLIFY_ARGS)))
    try(utils::write.csv(eGoSimp,
                          file.path(output_dir_0, "3_enriched_GO_terms_simplified.csv")))

    # Bar and dot plots — raw GO
    .save_png <- function(path, w, h, expr) {
      grDevices::png(path, units = "in", width = w, height = h, res = 300)
      try(print(expr))
      grDevices::dev.off()
    }

    .save_png(file.path(output_dir_0, "4_go_enrichment_bar_plot.png"), 8, 5,
              enrichplot::barplot(eGO, showCategory = ENRICH_PARAMS$show_categories))
    .save_png(file.path(output_dir_0, "5_go_enrichment_dot_plot.png"), 6, 7,
              enrichplot::dotplot(eGO, showCategory = ENRICH_PARAMS$show_categories))

    # Bar and dot plots — simplified GO
    .save_png(file.path(output_dir_0, "6_go_enrichment_simplified_bar_plot.png"), 8, 5,
              enrichplot::barplot(eGoSimp, showCategory = ENRICH_PARAMS$show_categories))
    .save_png(file.path(output_dir_0, "7_go_enrichment_simplified_dot_plot.png"), 6, 7,
              enrichplot::dotplot(eGoSimp, showCategory = ENRICH_PARAMS$show_categories))

    # ── GSEA (GO) ────────────────────────────────────────────────────────
    GSEA_GO <- try(create_GSEA_GO(DGE_func[[direction]]))
    try(utils::write.csv(GSEA_GO, file.path(output_dir_0, "8_GSEA_Go.csv")))

    .save_png(file.path(output_dir_0, "9_GSEA_Go.png"), 10, 10,
              enrichplot::dotplot(GSEA_GO,
                                  showCategory = ENRICH_PARAMS$show_categories,
                                  split        = ".sign") +
                ggplot2::facet_grid(. ~ .sign))

    GSEA_GoSimp <- try(do.call(clusterProfiler::simplify,
                                c(list(x = GSEA_GO), SIMPLIFY_ARGS)))
    try(utils::write.csv(GSEA_GoSimp,
                          file.path(output_dir_0, "10_GSEA_Go_simplified.csv")))

    .save_png(file.path(output_dir_0, "11_GSEA_Go_simplified.png"), 10, 10,
              enrichplot::dotplot(GSEA_GoSimp,
                                  showCategory = ENRICH_PARAMS$show_categories,
                                  split        = ".sign") +
                ggplot2::facet_grid(. ~ .sign))

    # ── Disease Ontology (DO) enrichment ──────────────────────────────────
    eDO <- try(DOSE::enrichDGN(
      ent_gene,
      pvalueCutoff  = ENRICH_PARAMS$do_pvalue_cutoff,
      pAdjustMethod = "BH",
      universe      = ent_uni,
      minGSSize     = 10L,
      maxGSSize     = 500L,
      qvalueCutoff  = ENRICH_PARAMS$do_qvalue_cutoff,
      readable      = FALSE
    ))
    try(utils::write.csv(eDO, file.path(output_dir_0, "12_enriched_DO_terms.csv")))

    .save_png(file.path(output_dir_0, "14_Do_enrichment_bar_plot.png"), 8, 5,
              enrichplot::barplot(eDO, showCategory = ENRICH_PARAMS$show_categories))
    .save_png(file.path(output_dir_0, "15_Do_enrichment_dot_plot.png"), 6, 7,
              enrichplot::dotplot(eDO, showCategory = ENRICH_PARAMS$show_categories))

    message("  Done: ", sample, " / ", direction)
  }
}

message("\nGO_analysis.R completed.")
