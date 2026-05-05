## =============================================================================
## FILE: TPM_calculation_RSEM_f.R
## PURPOSE: Calculate TPM from raw counts + effective gene lengths, then
##          cross-validate against pre-calculated RSEM TPM output.
##
## TWO METHODS:
##   (A) Manual formula  — tpm_gen() — matches RSEM's own calculation.
##   (B) DGEobj.utils::convertCounts() — independent library cross-check.
##
## OUTPUTS (in raw_data/):
##   tpm_data.tsv                             – TPM via manual formula
##   queried_genes_tpm_using_counts_and_lengths.csv
##   queried_genes_tpm_pre_calculated.csv
##   queried_genes_raw_count.csv
##
## DEPENDENCIES: 00_packages.R (edgeR 4.4.0, DGEobj.utils 1.0.6,
##                               tidyverse 2.0.0)
##               00_config.R, data_filtering_functions.R
## =============================================================================

rm(list = ls())

source("../scripts/00_packages.R")
source("../scripts/00_config.R")
source("../scripts/data_filtering_functions.R")

setwd(PROJECT_ROOT)


# ── 1. Load annotation ─────────────────────────────────────────────────────────
annotation_ensembl  <- utils::read.csv(FILES$annotation_ensembl)
annotation_enhanced <- create_annotation_enhanced()
annotation_enhanced <- dplyr::left_join(
  annotation_ensembl,
  annotation_enhanced,
  by = c("gene_id" = "ensembl_gene_id_version")
)
rownames(annotation_enhanced) <- annotation_enhanced$gene_id
rm(annotation_ensembl)


# ── 2. Load and align sample sheet ────────────────────────────────────────────
samples <- utils::read.delim(
  "raw_data/sample_sheet_final_2024_09_17.csv",
  header      = TRUE,
  sep         = ",",
  check.names = FALSE,
  row.names   = 1L
)
colnames(samples)[3L] <- "Enrichement_Step_2"
samples$sample_names  <- rownames(samples)

# Collapse GLI glioma subtypes
samples$Recoding_172[samples$Recoding_172 %in% c("LGG", "HGG", "GNOS")] <- "GLI"

# Sort for downstream consistency
samples <- samples[
  with(samples, order(Recoding_172, Enrichement_Step_2, sample_names)), ]


# ── 3. Load raw counts and gene lengths ───────────────────────────────────────
raw_counts_data <- utils::read.table(
  FILES$raw_counts,
  header    = TRUE,
  row.names = 1L,
  sep       = "\t"
)
raw_lengths_data <- utils::read.table(
  FILES$raw_lengths,
  header    = TRUE,
  row.names = 1L,
  sep       = "\t"
)

# Confirm row and column order match between counts and lengths
stopifnot(
  all(colnames(raw_counts_data) == colnames(raw_lengths_data)),
  all(rownames(raw_counts_data) == rownames(raw_lengths_data))
)

# Standardise column names to match sample sheet
sample_names <- gsub("[.]", "-", colnames(raw_counts_data))
sample_names <- gsub("^X",  "",  sample_names)
colnames(raw_counts_data)  <- sample_names
colnames(raw_lengths_data) <- sample_names
stopifnot(sum(!(colnames(raw_counts_data) %in% samples$sample_names)) == 0L)


# ── 4. TPM calculation (Method A — manual formula) ───────────────────────────

#' Standard TPM formula: normalise by gene length, then by library size.
#'
#' @param counts  genes × samples raw count matrix.
#' @param lengths genes × samples effective-length matrix (matching dimensions).
#' @return genes × samples TPM matrix.
tpm_gen <- function(counts, lengths) {
  rate_per_kb <- counts / lengths             # RPK: reads per kilobase
  t(t(rate_per_kb) * 1e6 / colSums(rate_per_kb))  # Scale to per-million
}

tpm_data <- as.data.frame(tpm_gen(raw_counts_data, raw_lengths_data)) %>%
  dplyr::mutate_if(is.numeric, round_2)

readr::write_tsv(
  tpm_data %>% tibble::rownames_to_column(var = "gene_ID"),
  file.path(PATHS$raw_data, "tpm_data.tsv")
)
message("TPM (method A) written to tpm_data.tsv")


# ── 5. TPM cross-validation (Method B — DGEobj.utils) ────────────────────────
tpm_data_b <- DGEobj.utils::convertCounts(
  as.matrix(raw_counts_data),
  unit      = "TPM",
  geneLength = as.matrix(raw_lengths_data),
  log       = FALSE,
  normalize = "none",
  prior.count = NULL
)
tpm_data_b <- as.data.frame(tpm_data_b) %>% dplyr::mutate_if(is.numeric, round_2)
message("TPM methods A and B match: ",
        all.equal(tpm_data, tpm_data_b, check.attributes = FALSE))


# ── 6. Cross-validate against pre-calculated RSEM TPM ────────────────────────
precalc_tpm <- utils::read.table(
  FILES$raw_tpm,
  header    = TRUE,
  row.names = 1L,
  sep       = "\t"
)
colnames(precalc_tpm) <- sample_names

stopifnot(
  all(colnames(precalc_tpm) == colnames(tpm_data)),
  all(rownames(precalc_tpm) == rownames(tpm_data))
)
message("Pre-calculated TPM matches computed TPM: ",
        isTRUE(all.equal(precalc_tpm, tpm_data, check.attributes = FALSE)))


# ── 7. Export queried-gene tables ─────────────────────────────────────────────
# Genes of biological interest for this project
QUERY_GENES <- c("MYCN", "MYC", "MYOD1", "PHOX2B", "KDSR", "FGFR4", "INAFM2")

queried_gene_ids <- annotation_enhanced[
  annotation_enhanced$gene_name %in% QUERY_GENES, "gene_id"]

.extract_queried <- function(tpm_matrix, gene_ids, samples_df, label) {
  mat <- as.data.frame(t(tpm_matrix[gene_ids, ]))
  colnames(mat) <- annotation_enhanced[colnames(mat), "gene_name"]
  mat %>%
    tibble::rownames_to_column("sample_names") %>%
    dplyr::left_join(
      samples_df[, c("sample_names", "Enrichement_Step_2", "Recoding_172")],
      by = "sample_names"
    ) %>%
    utils::write.csv(file.path(PATHS$raw_data, label))
  message("Written: ", label)
}

.extract_queried(tpm_data,     queried_gene_ids, samples,
                 "queried_genes_tpm_using_counts_and_lengths.csv")
.extract_queried(precalc_tpm,  queried_gene_ids, samples,
                 "queried_genes_tpm_pre_calculated.csv")

# Raw counts for the same queried genes
raw_q <- as.data.frame(t(raw_counts_data[queried_gene_ids, ]))
colnames(raw_q) <- annotation_enhanced[colnames(raw_q), "gene_name"]
raw_q %>%
  tibble::rownames_to_column("sample_names") %>%
  dplyr::left_join(
    samples[, c("sample_names", "Enrichement_Step_2", "Recoding_172")],
    by = "sample_names"
  ) %>%
  utils::write.csv(file.path(PATHS$raw_data, "queried_genes_raw_count.csv"))
message("Written: queried_genes_raw_count.csv")
