## =============================================================================
## FILE: functional_analysis_functions.R
## PURPOSE: GO, KEGG, GSEA, Reactome, and NCI gene-set enrichment wrappers.
##
## DEPENDENCIES: 00_packages.R (clusterProfiler 4.14.4, DOSE 4.0.0,
##                               enrichplot 1.26.6, org.Hs.eg.db 3.20.0,
##                               ReactomePA 1.50.0, tidyverse 2.0.0)
##               00_config.R   (ENRICH_PARAMS, FILES)
## =============================================================================


# ── 0. NCI gene-set (loaded once at source time) ───────────────────────────────
# Loaded when this module is sourced. If the GMT file is absent the variable is
# set to NULL and a warning is issued; scripts using create_NCI_enrichment() or
# create_NCI_GSEA() will receive NULL results and should handle them with try().
if (file.exists(FILES$nci_geneset)) {
  NCI_GeneSet <- clusterProfiler::read.gmt(FILES$nci_geneset)
} else {
  NCI_GeneSet <- NULL
  warning("NCI_GeneSet GMT file not found at: ", FILES$nci_geneset,
          "\n  create_NCI_enrichment() and create_NCI_GSEA() will return NULL.")
}


# ── 1. Shared helper: ranked gene list ────────────────────────────────────────

#' Prepare a named numeric vector of log2FC values suitable for GSEA.
#'
#' Removes NA values, deduplicates on `id_col`, then calls tibble::deframe()
#' to produce a named vector sorted by the caller (not sorted here so the
#' caller can apply pre-GSEA sorting if needed).
#'
#' @param df     data.frame with at least an ID column and a fold-change column.
#' @param id_col Bare column name for gene identifiers.
#' @param fc_col Bare column name for log2 fold change (default = log2FoldChange).
#' @return Named numeric vector (names = gene IDs, values = log2FC).
prep_ranked_genes <- function(df, id_col, fc_col = log2FoldChange) {
  df %>%
    tidyr::drop_na({{ id_col }}, {{ fc_col }}) %>%
    dplyr::distinct({{ id_col }}, .keep_all = TRUE) %>%
    dplyr::select({{ id_col }}, {{ fc_col }}) %>%
    tibble::deframe()
}


# ── 2. DGE directional split ──────────────────────────────────────────────────

#' Split filtered, annotated DGE results into up / down / all.
#'
#' @param DGE_filtered_anno_only data.frame; typically DGE$filtered_anno
#'   from create_annotated_DGElist().
#' @return Named list: all_genes, up_genes, down_genes.
create_func_DGE <- function(DGE_filtered_anno_only) {
  list(
    all_genes  = DGE_filtered_anno_only,
    up_genes   = dplyr::filter(DGE_filtered_anno_only, log2FoldChange > 0),
    down_genes = dplyr::filter(DGE_filtered_anno_only, log2FoldChange < 0)
  )
}


# ── 3. GO over-representation analysis ────────────────────────────────────────

#' Run GO enrichment (all ontologies: BP, CC, MF) against a background set.
#'
#' @param ent_gene Character vector of ENTREZ IDs for the gene set of interest.
#' @param ent_uni  Character vector of ENTREZ IDs for the background universe.
#' @return enrichResult from clusterProfiler::enrichGO().
create_eGO <- function(ent_gene, ent_uni) {
  clusterProfiler::enrichGO(
    gene          = ent_gene,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pvalueCutoff  = ENRICH_PARAMS$pvalue_cutoff,
    qvalueCutoff  = ENRICH_PARAMS$qvalue_cutoff,
    pAdjustMethod = ENRICH_PARAMS$padj_method,
    universe      = ent_uni,
    readable      = TRUE   # Convert ENTREZ IDs to gene symbols in output
  )
}


# ── 4. KEGG over-representation analysis ──────────────────────────────────────

#' Run KEGG pathway enrichment for Homo sapiens.
#'
#' @param ent_gene Character vector of ENTREZ IDs.
#' @return enrichResult with readable gene symbols.
create_eKEGG <- function(ent_gene) {
  eKEGG <- clusterProfiler::enrichKEGG(
    gene          = ent_gene,
    organism      = "hsa",
    pvalueCutoff  = ENRICH_PARAMS$pvalue_cutoff,
    qvalueCutoff  = ENRICH_PARAMS$qvalue_cutoff
  )
  # Convert internal IDs to readable gene symbols when result is non-null
  if (!is.null(eKEGG)) {
    clusterProfiler::setReadable(eKEGG,
                                 OrgDb   = org.Hs.eg.db,
                                 keyType = "ENTREZID")
  } else {
    eKEGG
  }
}


# ── 5. GO gene set enrichment analysis (GSEA) ─────────────────────────────────

#' Run GSEA against GO terms using gene SYMBOL identifiers.
#'
#' @param DGE_filtered_anno_only data.frame with a `gene_name` column and
#'   a `log2FoldChange` column.
#' @return gseaResult from clusterProfiler::gseGO().
create_GSEA_GO <- function(DGE_filtered_anno_only) {
  gene_list <- prep_ranked_genes(DGE_filtered_anno_only, gene_name)
  gene_list <- sort(gene_list, decreasing = TRUE)  # GSEA requires sorted input

  clusterProfiler::gseGO(
    geneList = gene_list,
    ont      = "ALL",
    keyType  = "SYMBOL",
    OrgDb    = org.Hs.eg.db,
    eps      = 1e-10   # Minimum p-value precision
  )
}


# ── 6. KEGG GSEA ──────────────────────────────────────────────────────────────

#' Run GSEA against KEGG pathways using ENTREZ IDs.
#'
#' @param DGE_filtered_anno_only data.frame with an `entrezgene_id` column.
#' @return gseaResult from clusterProfiler::gseKEGG().
create_GSEA_KEGG <- function(DGE_filtered_anno_only) {
  gene_list <- prep_ranked_genes(DGE_filtered_anno_only, entrezgene_id)
  gene_list <- sort(gene_list, decreasing = TRUE)

  clusterProfiler::gseKEGG(
    geneList     = gene_list,
    organism     = "hsa",
    pvalueCutoff = ENRICH_PARAMS$pvalue_cutoff,
    verbose      = FALSE
  )
}


# ── 7. NCI gene-set over-representation ───────────────────────────────────────

#' Run enrichment against the NCI gene-set collection (GMT format).
#'
#' @param ent_gene_names     Character vector of gene SYMBOLS of interest.
#' @param ent_uni_gene_names Character vector of background gene SYMBOLS.
#' @return enrichResult from clusterProfiler::enricher().
create_NCI_enrichment <- function(ent_gene_names, ent_uni_gene_names) {
  if (is.null(NCI_GeneSet)) {
    warning("create_NCI_enrichment: NCI_GeneSet is NULL — skipping.")
    return(NULL)
  }
  clusterProfiler::enricher(
    gene          = ent_gene_names,
    universe      = ent_uni_gene_names,
    TERM2GENE     = NCI_GeneSet,
    pvalueCutoff  = ENRICH_PARAMS$pvalue_cutoff,
    qvalueCutoff  = ENRICH_PARAMS$qvalue_cutoff,
    pAdjustMethod = "BH",   # Benjamini-Hochberg
    minGSSize     = 10L,
    maxGSSize     = 500L
  )
}


# ── 8. NCI GSEA ───────────────────────────────────────────────────────────────

#' Run GSEA against the NCI gene-set collection.
#'
#' @param DGE_filtered_anno_only data.frame with `gene_name` and `log2FoldChange`.
#' @return gseaResult from clusterProfiler::GSEA().
create_NCI_GSEA <- function(DGE_filtered_anno_only) {
  if (is.null(NCI_GeneSet)) {
    warning("create_NCI_GSEA: NCI_GeneSet is NULL — skipping.")
    return(NULL)
  }
  gene_list <- prep_ranked_genes(DGE_filtered_anno_only, gene_name)
  gene_list <- sort(gene_list, decreasing = TRUE)

  clusterProfiler::GSEA(
    geneList  = gene_list,
    TERM2GENE = NCI_GeneSet
  )
}


# ── 9. Reactome over-representation ───────────────────────────────────────────

#' Run Reactome pathway enrichment using ENTREZ IDs.
#'
#' @param ent_gene ENTREZ IDs for the gene set of interest.
#' @param ent_uni  ENTREZ IDs for the background universe.
#' @return enrichResult from ReactomePA::enrichPathway().
create_eRO <- function(ent_gene, ent_uni) {
  ReactomePA::enrichPathway(
    gene          = ent_gene,
    organism      = "human",
    universe      = ent_uni,
    pvalueCutoff  = ENRICH_PARAMS$pvalue_cutoff,
    qvalueCutoff  = ENRICH_PARAMS$qvalue_cutoff,
    pAdjustMethod = ENRICH_PARAMS$padj_method,
    minGSSize     = 10L,
    maxGSSize     = 500L,
    readable      = TRUE
  )
}


# ── 10. Reactome GSEA ─────────────────────────────────────────────────────────

#' Run GSEA against Reactome pathways.
#'
#' @param DGE_filtered_anno_only data.frame with `entrezgene_id` and `log2FoldChange`.
#' @return gseaResult from ReactomePA::gsePathway().
create_GSEA_RO <- function(DGE_filtered_anno_only) {
  gene_list <- prep_ranked_genes(DGE_filtered_anno_only, entrezgene_id)
  gene_list <- sort(gene_list, decreasing = TRUE)

  ReactomePA::gsePathway(
    geneList     = gene_list,
    organism     = "human",
    pvalueCutoff = ENRICH_PARAMS$pvalue_cutoff,
    eps          = 1e-10,
    verbose      = FALSE
  )
}
