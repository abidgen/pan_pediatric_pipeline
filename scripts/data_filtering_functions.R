## =============================================================================
## FILE: data_filtering_functions.R
## PURPOSE: Sample-sheet manipulation, raw-count filtering, and utility
##          functions shared across all analysis scripts.
##
## DEPENDENCIES: 00_packages.R (tidyverse, biomaRt, AnnotationHub)
##               00_config.R   (BIOMART, DGE_PARAMS)
## =============================================================================


# ── 1. Ensembl annotation ──────────────────────────────────────────────────────

#' Retrieve enhanced Ensembl gene annotation via biomaRt.
#'
#' Fetches five key attributes for every human gene from Ensembl,
#' then deduplicates on ensembl_gene_id_version so each row is unique.
#'
#' @param version  Integer. Ensembl archive version (default from BIOMART config).
#' @return data.frame with columns:
#'   ensembl_gene_id_version, ensembl_gene_id, entrezgene_id,
#'   description, external_gene_name
create_annotation_enhanced <- function(version = BIOMART$ensembl_version) {

  ensembl <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = BIOMART$dataset,
    version = version
  )

  annotation_df <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id_version",
      "ensembl_gene_id",
      "entrezgene_id",
      "description",
      "external_gene_name"
    ),
    uniqueRows = TRUE,
    values     = "1",       # Dummy filter; returns full table
    mart       = ensembl
  )

  # Remove duplicates introduced by multi-mapping transcripts
  annotation_df[!duplicated(annotation_df$ensembl_gene_id_version), ]
}


# ── 2. Enrichment-step selector ───────────────────────────────────────────────

#' Filter a data frame to a single library-type enrichment category.
#'
#' @param data            data.frame to filter (usually the sample sheet).
#' @param enrichment_col  Column name containing enrichment labels.
#' @param enrichment_type String to match (e.g. "polya", "ribozero").
#' @return Filtered data.frame.
select_enrichment <- function(data, enrichment_col, enrichment_type) {
  data[data[[enrichment_col]] == enrichment_type, ]
}


# ── 3. Category-recoding helpers ───────────────────────────────────────────────
# These three functions map a raw tumor-type string to a three-level label:
#   <specific>  |  Normal  |  all_other_tumors

#' Prefix match: groups any subtype starting with `specific` as that category.
specific_cat_v_all_col <- function(x, specific) {
  if (startsWith(x, specific)) specific
  else if (x == "NL")          "Normal"
  else                         "all_other_tumors"
}

#' Exact match: only the exact `specific` label maps to itself.
specific_v_all_col <- function(x, specific) {
  if (x == specific) x
  else if (x == "NL") "Normal"
  else                "all_other_tumors"
}

#' Exact match variant that folds "NL" into all_other_tumors.
specific_v_all_col_with_NL <- function(x, specific) {
  if (x == specific) x else "all_other_tumors"
}


# ── 4. Sample subsetting ───────────────────────────────────────────────────────

#' Build a "specific tumor vs all other tumors" sample subset.
#'
#' For grouped tumor families (RMS, NB, MBL) the function removes
#' sibling subtypes so they do not contaminate the "all other" group.
#'
#' @param data            Sample-sheet data.frame.
#' @param samples_col     Column containing tumor-type labels.
#' @param specific_sample Tumor type to call "specific".
#' @param exclude         Category to drop after recoding
#'                        (either "Normal" or "all_other_tumors").
#' @return data.frame with a new `specific_v_all` factor column.
specific_v_reference <- function(data, samples_col, specific_sample, exclude) {

  # Remove sibling subtypes for grouped families (e.g. RMS.1, RMS.2)
  family_prefix <- unlist(strsplit(specific_sample, "[.]"))[1]
  if (family_prefix %in% GROUPED_SUBTYPES) {
    siblings <- unique(data[[samples_col]])
    siblings <- siblings[startsWith(siblings, family_prefix) &
                         siblings != specific_sample]
    data <- data[!(data[[samples_col]] %in% siblings), ]
  }

  # Recode tumor labels into three-level contrast variable
  data$specific_v_all <- vapply(
    data[[samples_col]],
    FUN       = function(x) specific_v_all_col(x, specific = specific_sample),
    FUN.VALUE = character(1L)
  )

  # Drop unwanted category (e.g., exclude Normal when doing tumor-vs-tumor)
  data <- data[data$specific_v_all != exclude, ]

  # Set factor levels: reference group first, specific sample last
  lvls <- unique(data$specific_v_all)
  data$specific_v_all <- factor(
    data$specific_v_all,
    levels = c(lvls[lvls != specific_sample], specific_sample)
  )

  data
}


#' Variant of specific_v_reference that retains Normal samples.
#'
#' @inheritParams specific_v_reference
#' @return data.frame with `specific_v_all` factor (Normal kept as its own level).
specific_v_reference_with_NL <- function(data, samples_col, specific_sample) {

  family_prefix <- unlist(strsplit(specific_sample, "[.]"))[1]
  if (family_prefix %in% GROUPED_SUBTYPES) {
    siblings <- unique(data[[samples_col]])
    siblings <- siblings[startsWith(siblings, family_prefix) &
                         siblings != specific_sample]
    data <- data[!(data[[samples_col]] %in% siblings), ]
  }

  data$specific_v_all <- vapply(
    data[[samples_col]],
    FUN       = function(x) specific_v_all_col_with_NL(x, specific = specific_sample),
    FUN.VALUE = character(1L)
  )

  lvls <- unique(data$specific_v_all)
  data$specific_v_all <- factor(
    data$specific_v_all,
    levels = c(lvls[lvls != specific_sample], specific_sample)
  )

  data
}


#' Prefix-based category comparison (e.g., all RMS subtypes vs rest).
#'
#' @inheritParams specific_v_reference
specific_cat_v__reference <- function(data, samples_col, specific_sample, exclude) {

  data$specific_v_all <- vapply(
    data[[samples_col]],
    FUN       = function(x) specific_cat_v_all_col(x, specific = specific_sample),
    FUN.VALUE = character(1L)
  )

  data <- data[data$specific_v_all != exclude, ]
  lvls <- unique(data$specific_v_all)
  data$specific_v_all <- factor(
    data$specific_v_all,
    levels = c(lvls[lvls != specific_sample], specific_sample)
  )

  data
}


# ── 5. Row-number index table ─────────────────────────────────────────────────

#' Build a per-group row index table for vectorised subsetting.
#'
#' Returns a data.frame with first row, last row, and sample count for
#' each unique tumor group.  Used inside filter_raw_counts().
#'
#' @param data        data.frame sorted by samples_col.
#' @param samples_col Column containing group labels.
#' @return data.frame(group, f_row_number, l_row_number, sample_count).
row_num_table <- function(data, samples_col) {
  data$f_row_number <- seq_len(nrow(data))
  idx <- data[!duplicated(data[[samples_col]]), c(samples_col, "f_row_number")]
  rownames(idx) <- NULL

  idx <- idx %>%
    dplyr::mutate(
      l_row_number = dplyr::lead(f_row_number, 1L) - 1L,
      sample_count = NA_integer_
    )

  # Fix the last group's end boundary
  idx$l_row_number[nrow(idx)] <- nrow(data)
  idx$sample_count <- idx$l_row_number - idx$f_row_number + 1L
  idx
}


# ── 6. Raw count filtering ────────────────────────────────────────────────────

#' Retain only genes expressed above threshold in ≥75 % of samples
#' in at least one tumor group.
#'
#' Threshold: ≥15 counts in ≥75 % of per-group samples.
#' This is stricter than DESeq2's default filterByExpr and removes
#' very lowly expressed genes before normalisation.
#'
#' @param raw_counts   genes × samples count matrix / data.frame.
#' @param samples      Sample-sheet data.frame (rownames = sample IDs).
#' @param samples_col  Column containing tumor-type group labels.
#' @return Filtered count matrix containing only samples in `samples`.
filter_raw_counts <- function(raw_counts, samples, samples_col) {

  samples_row_counts <- row_num_table(samples, samples_col)

  # Restrict to the requested sample set
  raw_counts_sub <- raw_counts[, rownames(samples), drop = FALSE]

  # Per-group logical pass/fail matrix (one column per group)
  group_pass <- vapply(
    seq_len(nrow(samples_row_counts)),
    FUN = function(i) {
      col_idx <- samples_row_counts$f_row_number[i]:samples_row_counts$l_row_number[i]
      rowSums(raw_counts_sub[, col_idx] >= 15L) >=
        0.75 * samples_row_counts$sample_count[i]
    },
    FUN.VALUE = logical(nrow(raw_counts_sub))
  )

  # Keep gene if it passes the threshold in ANY group
  keep <- apply(group_pass, 1L, any)
  raw_counts_sub[keep, , drop = FALSE]
}


# ── 7. Utility functions ───────────────────────────────────────────────────────

#' Return dimensions of an object as a named list.
find_dim <- function(x) list(dim = dim(x))

#' Round a numeric value to 2 decimal places.
round_2 <- function(x) round(x, 2L)
