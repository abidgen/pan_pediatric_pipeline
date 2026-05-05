## =============================================================================
## FILE: scripts/TF_enrichment_function.R
## PURPOSE: Hypergeometric TF regulatory enrichment against a GRN.
##          Shared by:
##            analysis/6_GRN_TF_enrichment.R   (require_tf_in_list = FALSE)
##            analysis/7_WGCNA.R               (require_tf_in_list = FALSE)
##            analysis/4_GO_analysis.R (DEG-level, require_tf_in_list = TRUE)
##
## DEPENDENCIES: 00_packages.R (stringr 1.5.1)
##               00_config.R   (TF_PARAMS)
## =============================================================================

# NOTE: 00_packages.R and 00_config.R must be sourced by the calling script
# before sourcing this file. Do not source them here to avoid double-loading.


#' Hypergeometric TF enrichment against a transcriptional regulatory network.
#'
#' For each TF row in `TRN`:
#'   1. Optionally checks that the TF is itself in the query gene list.
#'   2. Computes the overlap between the gene list and that TF's targets.
#'   3. One-sided hypergeometric test  P(X >= overlap).
#'   4. Bonferroni FDR correction across all tested TFs.
#'
#' The `require_tf_in_list` flag reconciles two original analysis scripts:
#'   - DEG-level analysis required the TF to be a DEG itself     (TRUE).
#'   - GRN / WGCNA analysis tests all TFs in the network         (FALSE).
#'
#' @param TRN               data.frame with columns `tfs` and `targets`
#'                          (comma-separated target IDs per TF row).
#' @param GeneList          Character vector of gene IDs to test.
#' @param universe          Integer total background gene count.
#' @param filename          Label written to the `set` column (e.g. tumour name).
#' @param require_tf_in_list Logical. If TRUE (default) skip TFs not in GeneList.
#'                          Set FALSE for GRN / WGCNA enrichment.
#'
#' @return data.frame(TF, set, overlap, p, gene.size, tf.size, universe, fdr)
#'   sorted by ascending p-value.
getTRN_enrichment <- function(TRN,
                              GeneList,
                              universe,
                              filename,
                              require_tf_in_list = TRUE) {

  n_gene <- length(GeneList)

  res <- data.frame(
    TF = character(), set = character(), overlap = integer(),
    p = double(), gene.size = integer(), tf.size = integer(),
    universe = integer(), stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(TRN))) {

    # Parse comma-separated targets; sanitise whitespace and deduplicate
    targets <- unique(stringr::str_trim(strsplit(TRN$targets[i], ",")[[1L]]))
    n_tf    <- length(targets)

    if (n_tf   < TF_PARAMS$cutoff) next
    if (n_gene < TF_PARAMS$cutoff) next
    if (require_tf_in_list && !(TRN$tfs[i] %in% GeneList)) next

    overlap <- length(intersect(GeneList, targets))
    if (overlap < TF_PARAMS$cutoff) next

    p_val <- stats::phyper(
      q = overlap - 1L, m = n_gene, n = universe - n_gene,
      k = n_tf, lower.tail = FALSE
    )

    res[nrow(res) + 1L, ] <- list(
      TF = TRN$tfs[i], set = filename, overlap = overlap,
      p = p_val, gene.size = n_gene, tf.size = n_tf, universe = universe
    )
  }

  if (nrow(res) == 0L) return(res)

  res$p   <- as.double(res$p)
  res     <- res[order(res$p), ]
  res$fdr <- stats::p.adjust(res$p, method = "bonferroni")
  res
}
