## =============================================================================
## FILE: analysis/6_GRN_TF_enrichment.R
## PURPOSE: Per-tumour TF enrichment using the GENIE3 gene regulatory network.
##          Tests whether DEGs (up-regulated by FC and p-value) are
##          over-represented among each TF's target genes.
##
## INPUTS  (defined in 00_config.R):
##   GRN_FILES$grn_network   – GENIE3 GRN edge table (tab-separated)
##   GRN_FILES$tf_db         – TF database (DatabaseExtract_v_1.01.xlsx)
##   PATHS$output/.../       – Per-tumour DGE full result CSVs
##
## OUTPUTS (written relative to each tumour's DGE output directory):
##   Grn-enrichments.<tumour>.overlap_cutoff.<n>.alpha.<a>.DEGs_fc.<fc>.<date>.txt
##   Genie3_GRN_enrichment_DEGs-up_enrichmentfdr-<fdr>.txt  (master table)
##
## REDUNDANCY REMOVED:
##   getTRN_enrichment() previously duplicated here; now sourced from
##   scripts/TF_enrichment_function.R with require_tf_in_list = FALSE.
##
## DEPENDENCIES: 00_packages.R (stringr 1.5.1, readxl 1.4.3, plyr 1.8.9,
##                               dplyr 1.1.4)
##               00_config.R, TF_enrichment_function.R
## =============================================================================

rm(list = ls())

source("../scripts/00_packages.R")
source("../scripts/00_config.R")
source("../scripts/TF_enrichment_function.R")

setwd(PROJECT_ROOT)


# ── Locate per-tumour DGE directories ─────────────────────────────────────────
dge_root  <- file.path(PATHS$output, "DESeq2_analysis", "specific_v_all_other_tumors")
dge_dirs  <- list.dirs(dge_root, full.names = TRUE, recursive = FALSE)
# Exclude the root itself and keep only immediate tumour subdirectories
dge_dirs  <- dge_dirs[file.info(dge_dirs)$isdir]
stopifnot(length(dge_dirs) > 0L)
message("Found ", length(dge_dirs), " tumour DGE directories.")


# ── Load and reshape GRN ──────────────────────────────────────────────────────
# GRN has one row per TF-target pair; aggregate targets per TF to a string list
TRN <- utils::read.delim(GRN_FILES$grn_network, sep = "\t")
colnames(TRN)[1:2] <- c("tfs", "targets")
TRN <- stats::aggregate(targets ~ tfs, TRN, toString)
message("GRN loaded: ", nrow(TRN), " TFs.")


# ── Per-tumour enrichment loop ────────────────────────────────────────────────
master_res <- NULL

for (dge_dir in dge_dirs) {

  tumour_name <- basename(dge_dir)
  full_dge_path <- file.path(dge_dir, "5_DGElist_full.csv")   # full (unfiltered) DGE

  if (!file.exists(full_dge_path)) {
    warning("Skipping ", tumour_name, " — DGE file not found: ", full_dge_path)
    next
  }

  all_genes <- utils::read.delim(full_dge_path, sep = ",")
  universe  <- nrow(all_genes)

  # Up-regulated DEGs only: FC >= threshold AND padj < alpha
  DEGs <- all_genes[
    !is.na(all_genes$log2FoldChange) &
    !is.na(all_genes$padj) &
    all_genes$log2FoldChange >= GRN_PARAMS$fc_threshold &
    all_genes$padj            <  GRN_PARAMS$alpha,
    "ensembl_gene_id"
  ]

  message("  ", tumour_name, ": ", length(DEGs), " up-regulated DEGs")

  # Run TF enrichment (TF need not be a DEG itself: require_tf_in_list = FALSE)
  res <- getTRN_enrichment(
    TRN                = TRN,
    GeneList           = DEGs,
    universe           = universe,
    filename           = tumour_name,
    require_tf_in_list = FALSE
  )

  # Write per-tumour result
  out_file <- file.path(
    dge_dir,
    paste("Grn-enrichments", tumour_name,
          "overlap_cutoff", TF_PARAMS$cutoff,
          "alpha", GRN_PARAMS$alpha,
          "DEGs_fc", GRN_PARAMS$fc_threshold,
          "date_created", Sys.Date(),
          "txt", sep = ".")
  )
  utils::write.table(res, out_file, row.names = FALSE, quote = FALSE, sep = "\t")

  master_res <- rbind(master_res, res)
}

stopifnot(!is.null(master_res))


# ── Annotate master results with TF metadata ──────────────────────────────────
TFs_desc <- as.data.frame(readxl::read_excel(GRN_FILES$tf_db))
TFs_desc <- TFs_desc[, c("Ensembl ID", "HGNC symbol", "EntrezGene Description")]
colnames(TFs_desc)[1L] <- "TF"
TFs_desc <- TFs_desc[TFs_desc$TF %in% master_res$TF, ]

# plyr::join preserves row order (unlike merge)
master_res <- plyr::join(master_res, TFs_desc, by = "TF")
master_res <- master_res[master_res$fdr < GRN_PARAMS$fdr_cutoff, ]


# ── Write master output ───────────────────────────────────────────────────────
master_out <- file.path(
  dge_root,
  paste0("Genie3_GRN_enrichment_DEGs-up_enrichmentfdr-",
         GRN_PARAMS$fdr_cutoff, ".txt")
)
utils::write.table(master_res, master_out, row.names = FALSE, quote = FALSE)
message("Master TF enrichment written to: ", master_out)
