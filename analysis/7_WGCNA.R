## =============================================================================
## FILE: analysis/7_WGCNA.R
## PURPOSE: Weighted Gene Co-expression Network Analysis (WGCNA) for the
##          Pan-Paediatric Cancer cohort.
##
## PIPELINE STAGES (numbered to match output file prefixes):
##   01  QC on raw counts         — sample clustering + outlier removal
##   02  Normalisation (VST)
##   03  Soft-threshold selection — scale-free topology fit
##   04  TOM construction         — topological overlap matrix
##   05  Module detection         — dynamic tree cut
##   06  Module merging           — eigengene dissimilarity
##   07–09  Module eigengenes + membership
##   10  Intra-modular connectivity
##   11–13  Module–tumour correlation (Pearson + mixed model)
##   14  Filtered module heatmaps
##   15–18  Module–module correlation, MDS, dendrogram
##   19  Module membership violin plot
##   20  TF enrichment heatmap
##   21–23  Hub gene identification
##   24  lncRNA / non-lncRNA proportion
##   99  Final summary table
##
## REDUNDANCY REMOVED vs original:
##   - specific_v_reference__() removed; uses shared specific_v_reference()
##   - All package library() calls removed; loaded via 00_packages.R
##   - All View() calls removed (not SLURM-safe)
##   - create_NCI_enrichemnt typo fixed → create_NCI_enrichment
##   - Hardcoded paths replaced with PATHS / GRN_FILES / WGCNA_PARAMS
##
## DEPENDENCIES: 00_packages.R (WGCNA 1.73, CorLevelPlot 0.99.2,
##                               lmerTest 3.1.3, lme4 1.1.35,
##                               corrplot 0.92, ComplexHeatmap 2.22.0,
##                               circlize 0.4.16, pheatmap 1.0.12,
##                               RColorBrewer 1.1.3, reshape2 1.4.4,
##                               plyr 1.8.9, tidyverse 2.0.0,
##                               DESeq2 1.46.0, BiocParallel 1.40.0)
##               00_config.R (WGCNA_PARAMS, GRN_FILES, PATHS)
##               data_filtering_functions.R, functional_analysis_functions.R,
##               TF_enrichment_function.R
## =============================================================================

rm(list = ls())

source("../scripts/00_packages.R")
source("../scripts/00_config.R")
source("../scripts/data_filtering_functions.R")
source("../scripts/functional_analysis_functions.R")
source("../scripts/TF_enrichment_function.R")

setwd(PROJECT_ROOT)

# Maximize WGCNA thread usage (respects SLURM --cpus-per-task allocation)
WGCNA::enableWGCNAThreads(nThreads = VIS_PARAMS$n_cores)
options(stringsAsFactors = FALSE)

# ── Create master output directory ────────────────────────────────────────────
if (!dir.exists(PATHS$wgcna_output)) {
  dir.create(PATHS$wgcna_output, recursive = TRUE)
}


# =============================================================================
# SECTION 1: Load data
# =============================================================================

# ── Annotation ────────────────────────────────────────────────────────────────
annotation_ensembl  <- utils::read.csv(FILES$annotation_ensembl)
annotation_enhanced <- create_annotation_enhanced()
annotation_enhanced <- dplyr::left_join(
  annotation_ensembl, annotation_enhanced,
  by = c("gene_id" = "ensembl_gene_id_version")
)
rownames(annotation_enhanced) <- annotation_enhanced$gene_id
rm(annotation_ensembl)

# ── Sample sheet (tumour samples only; exclude Normal) ────────────────────────
samples_with_NL <- utils::read.delim(
  FILES$sample_sheet, header = TRUE, sep = ",",
  check.names = FALSE, row.names = 1L
)
samples <- samples_with_NL[samples_with_NL$Recoding_172 != "NL", ]
message("Samples loaded: ", nrow(samples))

# ── Batch-corrected raw counts ────────────────────────────────────────────────
raw_counts <- utils::read.delim(
  FILES$post_bc_counts, header = TRUE, sep = "\t",
  check.names = FALSE, row.names = 1L
)

# Filter to tumour samples and remove lowly expressed genes
raw_counts_data <- filter_raw_counts(raw_counts, samples, "Recoding_172")

# Remove PAR_Y pseudo-autosomal genes (can cause spurious correlations)
raw_counts_filtered <- raw_counts_data %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  dplyr::filter(!stringr::str_detect(gene_id, "PAR_Y")) %>%
  tibble::column_to_rownames(var = "gene_id")

# Guard: sample order must be consistent
stopifnot(all(samples$sample_names == colnames(raw_counts_filtered)))
message("Counts filtered: ", nrow(raw_counts_filtered), " genes x ",
        ncol(raw_counts_filtered), " samples")


# =============================================================================
# SECTION 2: QC — outlier detection on raw counts
# =============================================================================

raw_counts_t <- t(raw_counts_filtered)

gsg <- WGCNA::goodSamplesGenes(raw_counts_t)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)   > 0L) message("Removing genes: ",
    paste(names(raw_counts_t)[!gsg$goodGenes], collapse = ", "))
  if (sum(!gsg$goodSamples) > 0L) message("Removing samples: ",
    paste(rownames(raw_counts_t)[!gsg$goodSamples], collapse = ", "))
  raw_counts_t <- raw_counts_t[gsg$goodSamples, gsg$goodGenes]
} else {
  message("QC (raw): all samples and genes passed.")
}

# Hierarchical clustering dendrogram to visually confirm outlier removal
htree_raw <- stats::hclust(stats::dist(raw_counts_t), method = "average")
grDevices::pdf(file.path(PATHS$wgcna_output,
                         "01_sampleClustering_hclust_raw.pdf"),
               width = 400, height = 10)
graphics::par(cex = 1.3, mar = c(0, 4, 2, 0))
plot(htree_raw,
     main = "Sample clustering — raw counts (outlier detection)",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
grDevices::dev.off()


# =============================================================================
# SECTION 3: Normalisation (DESeq2 VST)
# =============================================================================

dds_wgcna <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(raw_counts_filtered),
  colData   = samples,
  design    = ~ 1   # Intercept only; WGCNA uses expression patterns not contrasts
)
dds_wgcna <- DESeq2::estimateSizeFactors(dds_wgcna)
dds_wgcna <- DESeq2::DESeq(dds_wgcna, parallel = TRUE)

norm_counts_t <- DESeq2::assay(
  DESeq2::varianceStabilizingTransformation(dds_wgcna)
) %>% t()

# QC on normalised counts
gsg2 <- WGCNA::goodSamplesGenes(norm_counts_t)
if (!gsg2$allOK) {
  norm_counts_t <- norm_counts_t[gsg2$goodSamples, gsg2$goodGenes]
  message("Post-VST QC removed outliers.")
} else {
  message("QC (VST): all samples and genes passed.")
}

htree_norm <- stats::hclust(stats::dist(norm_counts_t), method = "average")
grDevices::pdf(file.path(PATHS$wgcna_output,
                         "02_sampleClustering_hclust_normalised.pdf"),
               width = 400, height = 10)
graphics::par(cex = 1.3, mar = c(0, 4, 2, 0))
plot(htree_norm,
     main = "Sample clustering — VST normalised (outlier detection)",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
grDevices::dev.off()


# =============================================================================
# SECTION 4: Soft-threshold selection
# =============================================================================

powers      <- c(1:20, seq(22, 30, by = 2))
networkType <- WGCNA_PARAMS$network_type

sft <- WGCNA::pickSoftThreshold(
  norm_counts_t,
  powerVector = powers,
  verbose     = 5,
  networkType = networkType
)
message("Estimated soft-threshold power: ", sft$powerEstimate)

grDevices::pdf(file.path(PATHS$wgcna_output, "03_network_topology_analysis.pdf"),
               width = 10, height = 15)
graphics::par(mfrow = c(2, 1))
cex1 <- 0.9

# Scale-free topology fit (R^2)
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
graphics::abline(h = 0.80, col = "red")

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")
grDevices::dev.off()


# =============================================================================
# SECTION 5: TOM construction and module detection
# =============================================================================

power <- WGCNA_PARAMS$power   # Set in 00_config.R after reviewing sft plot

adjacency_matrix <- WGCNA::adjacency(norm_counts_t, power = power,
                                     type = networkType)
TOM     <- WGCNA::TOMsimilarity(adjacency_matrix, TOMType = networkType)
dissTOM <- 1 - TOM

# Gene dendrogram
geneTree <- stats::hclust(stats::as.dist(dissTOM), method = "average")

grDevices::pdf(file.path(PATHS$wgcna_output,
                         "04_gene_clustering_TOM_dissimilarity.pdf"),
               width = 12, height = 9)
plot(geneTree, xlab = "", sub = "",
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
grDevices::dev.off()

# Dynamic tree cut to define initial modules
dynamicMods <- WGCNA::cutreeDynamic(
  dendro          = geneTree,
  distM           = dissTOM,
  deepSplit       = WGCNA_PARAMS$deep_split,
  pamRespectsDendro = FALSE,
  minClusterSize  = WGCNA_PARAMS$min_module_size
)
dynamicColors <- WGCNA::labels2colors(dynamicMods)
message("Initial modules: ", length(table(dynamicColors)))

grDevices::pdf(file.path(PATHS$wgcna_output,
                         "05_gene_dendrogram_module_colors.pdf"),
               width = 8, height = 6)
WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                           dendroLabels = FALSE, hang = 0.03,
                           addGuide = TRUE, guideHang = 0.05,
                           main = "Gene dendrogram and module colors")
grDevices::dev.off()


# =============================================================================
# SECTION 6: Module merging
# =============================================================================

MEList <- WGCNA::moduleEigengenes(norm_counts_t, colors = dynamicColors)
MEs    <- MEList$eigengenes
MEDiss <- 1 - stats::cor(MEs)
METree <- stats::hclust(stats::as.dist(MEDiss), method = "average")

grDevices::pdf(file.path(PATHS$wgcna_output,
                         "06_clustering_module_eigengenes.pdf"),
               width = 12, height = 12)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
graphics::abline(h = WGCNA_PARAMS$merge_threshold, col = "red")
grDevices::dev.off()

merge        <- WGCNA::mergeCloseModules(norm_counts_t, dynamicColors,
                                         cutHeight = WGCNA_PARAMS$merge_threshold,
                                         verbose = 3)
mergedColors <- merge$colors
mergedMEs    <- merge$newMEs
message("Merged modules: ", length(table(mergedColors)))

grDevices::pdf(file.path(PATHS$wgcna_output, "07_merged_module_tree.pdf"),
               width = 12, height = 9)
WGCNA::plotDendroAndColors(
  geneTree,
  cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
grDevices::dev.off()

utils::write.table(merge$oldMEs,
                   file.path(PATHS$wgcna_output, "08_oldMEs.txt"))
utils::write.table(merge$newMEs,
                   file.path(PATHS$wgcna_output, "09_newMEs_merged.txt"))


# =============================================================================
# SECTION 7: Module membership calculation
# =============================================================================

unmergedMEs <- merge$oldMEs

module_membership <- data.frame(
  ID                     = colnames(norm_counts_t),
  module                 = dynamicColors,
  module_name            = paste0("ME", dynamicColors),
  color                  = dynamicColors,
  module_membership_value = NA_real_,
  stringsAsFactors = FALSE
)
rownames(module_membership) <- colnames(norm_counts_t)

for (i in seq_len(nrow(module_membership))) {
  gene_id  <- module_membership$ID[i]
  me_name  <- module_membership$module_name[i]
  module_membership$module_membership_value[i] <-
    stats::cor(unmergedMEs[, me_name], norm_counts_t[, gene_id])
}

utils::write.table(
  module_membership,
  file.path(PATHS$wgcna_output,
            paste0("10_module_membership_power_", power,
                   "_networktype_", networkType, ".txt")),
  sep = "\t"
)

# Write per-module gene lists (for downstream enrichment tools)
old_mod_dir <- file.path(PATHS$wgcna_output, "mod_memberships", "old_modules")
if (!dir.exists(old_mod_dir)) dir.create(old_mod_dir, recursive = TRUE)

for (i in seq_along(names(table(module_membership$module)))) {
  mod <- names(table(module_membership$module))[i]
  genes <- module_membership[module_membership$module == mod, "ID"]
  utils::write.table(
    genes,
    file.path(old_mod_dir, paste("sb", i, mod, "genes.txt", sep = "_")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

utils::write.table(
  module_membership$ID,
  file.path(PATHS$wgcna_output, "module_membership_universe.txt"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)


# =============================================================================
# SECTION 8: Intra-modular connectivity
# =============================================================================

intra_connectivity <- WGCNA::intramodularConnectivity.fromExpr(
  norm_counts_t, dynamicColors,
  networkType = networkType, power = power,
  getWholeNetworkConnectivity = TRUE
)
rownames(intra_connectivity) <- colnames(norm_counts_t)
intra_connectivity$module    <- dynamicColors

utils::write.table(
  intra_connectivity,
  file.path(PATHS$wgcna_output,
            paste0("10_intra_modular_connectivity_power_", power,
                   "_networktype_", networkType, ".txt")),
  sep = "\t"
)


# =============================================================================
# SECTION 9: Cytoscape export
# =============================================================================

.export_cytoscape <- function(ME_list, colors, label, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  for (i in seq_along(ME_list)) {
    module_color <- substring(names(ME_list)[i], 3L)
    if (module_color == "grey") next   # Grey = unassigned genes; skip
    inModule <- is.finite(match(colors, module_color))
    modGenes <- colnames(norm_counts_t)[inModule]
    modTOM   <- TOM[inModule, inModule]
    dimnames(modTOM) <- list(modGenes, modGenes)
    WGCNA::exportNetworkToCytoscape(
      modTOM,
      edgeFile = file.path(out_dir, paste0(label, "_edges-", module_color, ".txt")),
      nodeFile = file.path(out_dir, paste0(label, "_nodes-", module_color, ".txt")),
      weighted = TRUE, threshold = -1,
      nodeNames = modGenes, nodeAttr = colors[inModule]
    )
  }
}

.export_cytoscape(merge$oldMEs, dynamicColors, "origin",
                  file.path(PATHS$wgcna_output, "cytoscape_inputs", "oldMEs"))
.export_cytoscape(merge$newMEs, mergedColors,  "merged",
                  file.path(PATHS$wgcna_output, "cytoscape_inputs", "newMEs"))


# =============================================================================
# SECTION 10: Module–tumour Pearson correlation
# =============================================================================

# Build binary tumour-type indicator matrix (samples × tumour types)
tumor_types         <- names(table(samples$Recoding_172))
tumor_binary        <- as.data.frame(
  vapply(tumor_types,
         function(tt) as.integer(samples$Recoding_172 == tt),
         integer(nrow(samples)))
)
rownames(tumor_binary) <- rownames(samples)

# Recalculate ordered MEs for correlation
MEs0 <- WGCNA::moduleEigengenes(norm_counts_t, dynamicColors)$eigengenes
MEs  <- WGCNA::orderMEs(MEs0)
stopifnot(all(rownames(MEs) == rownames(tumor_binary)))

moduleTraitCor    <- stats::cor(MEs, tumor_binary, use = "p")
moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(norm_counts_t))

utils::write.table(moduleTraitCor,
                   file.path(PATHS$wgcna_output, "12_moduleTrait_correlation.txt"))
utils::write.table(moduleTraitPvalue,
                   file.path(PATHS$wgcna_output, "12_moduleTrait_pValue.txt"))

# Mask non-significant entries for heatmap
moduleTraitCor_masked <- moduleTraitCor
moduleTraitCor_masked[moduleTraitPvalue > 0.05] <- NA

textMatrix <- paste0(round_2(moduleTraitCor), "\n(",
                     signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

WGCNA_COLORS <- c("#0000CD", "#2133D4", "#4367DC", "#659AE3", "#87CEEB",
                  "white", "#FFC0CB", "#F29098", "#E66065", "#D93032", "#CD0000")

grDevices::pdf(file.path(PATHS$wgcna_output,
                         "12_Module-tumour_relationships_Pearson.pdf"),
               width = 30, height = 40)
graphics::par(mar = c(15, 12, 3, 3))
WGCNA::labeledHeatmap(
  Matrix       = moduleTraitCor_masked,
  xLabels      = colnames(tumor_binary),
  yLabels      = colnames(MEs),
  ySymbols     = colnames(MEs),
  colorLabels  = FALSE,
  colors       = WGCNA_COLORS,
  textMatrix   = textMatrix,
  setStdMargins = FALSE,
  cex.text     = 0.9,
  zlim         = c(-1, 1),
  main         = "Module–tumour relationships (Pearson)"
)
grDevices::dev.off()


# =============================================================================
# SECTION 11: Module–tumour mixed-model correlation (StudyCenters as RE)
# =============================================================================

module_significance_data <- as.data.frame(
  matrix(NA_real_, ncol = length(tumor_types), nrow = ncol(MEs),
         dimnames = list(colnames(MEs), tumor_types))
)
module_correlation_data  <- module_significance_data

for (tumor_type in tumor_types) {

  # Subset samples: specific tumour vs all others (excluding Normal)
  samples_specific <- specific_v_reference(
    data            = samples,
    samples_col     = "Recoding_172",
    specific_sample = tumor_type,
    exclude         = "Normal"
  )

  tumor_factor  <- samples_specific$specific_v_all
  batch_effect  <- samples_specific$StudyCenters
  sample_MEs    <- MEs[rownames(samples_specific), ]

  # Fit linear mixed model per module; extract t-value and BH-adjusted p
  temp <- t(sapply(seq_len(ncol(sample_MEs)), function(x) {
    fit    <- lmerTest::lmer(sample_MEs[, x] ~ tumor_factor + (1 | batch_effect))
    coefs  <- summary(fit)$coefficients
    p_vals <- coefs[, 5L]   # p-value
    t_vals <- coefs[, 4L]   # t-value (used as effect size proxy)
    c(p_vals, t_vals)
  }))
  rownames(temp) <- colnames(sample_MEs)

  # Second row = specific-tumour coefficient (first = intercept)
  p_col  <- as.data.frame(temp)[2L]
  p_adj  <- stats::p.adjust(p_col[, 1L], method = "BH")
  t_col  <- as.data.frame(temp)[4L]

  module_significance_data[rownames(module_significance_data), tumor_type] <- p_adj
  module_correlation_data[rownames(module_correlation_data),   tumor_type] <- t_col[rownames(module_correlation_data), ]
}

utils::write.table(module_significance_data,
                   file.path(PATHS$wgcna_output,
                             "14_module_significance_BH_mixed_model.txt"))
utils::write.table(module_correlation_data,
                   file.path(PATHS$wgcna_output,
                             "14_module_t-values_mixed_model.txt"))

# Mask entries with adj-p > 0.05 for plotting
module_correlation_masked <- module_correlation_data
module_correlation_masked[module_significance_data > 0.05] <- NA

textMatrix_mm <- matrix(
  paste0(round(data.matrix(module_correlation_data), 3L), "\n(",
         signif(data.matrix(module_significance_data), 1L), ")"),
  nrow(module_correlation_data), ncol(module_correlation_data)
)

grDevices::pdf(file.path(PATHS$wgcna_output,
                         "14_module_significance_t-values_mixed_model.pdf"),
               width = 30, height = 40)
graphics::par(mar = c(15, 12, 3, 3))
WGCNA::labeledHeatmap(
  Matrix       = module_correlation_masked,
  xLabels      = colnames(module_correlation_masked),
  yLabels      = rownames(module_correlation_masked),
  ySymbols     = rownames(module_correlation_masked),
  colorLabels  = FALSE, colors = WGCNA_COLORS,
  textMatrix   = textMatrix_mm,
  setStdMargins = FALSE, cex.text = 0.9, zlim = c(-100, 100),
  main = "Module–tumour t-values (mixed model); adj-p > 0.05 = grey"
)
grDevices::dev.off()


# =============================================================================
# SECTION 12: Filter significant high-correlation modules
# =============================================================================

# Keep modules where at least one tumour exceeds the t-value cutoff
t_cutoff    <- WGCNA_PARAMS$t_value_cutoff
cor_filtered <- module_correlation_masked
cor_filtered[cor_filtered < t_cutoff] <- NA

filtered_module_list <- rownames(
  module_correlation_data[rowSums(!is.na(cor_filtered)) > 0, ]
)
message("Filtered modules (t > ", t_cutoff, "): ", length(filtered_module_list))

filtered_corr_NA  <- module_correlation_data[filtered_module_list, ]
filtered_corr_NA[filtered_corr_NA < t_cutoff] <- NA
filtered_sig      <- module_significance_data[filtered_module_list, ]
filtered_sig[filtered_sig <= 1e-100] <- 1e-100

filtered_module_membership <- module_membership[
  module_membership$module_name %in% filtered_module_list, ]
utils::write.table(filtered_module_membership,
                   file.path(PATHS$wgcna_output,
                             "19_0_filtered_module_membership.txt"))

textMatrix_filt <- matrix(
  paste0(round(data.matrix(filtered_corr_NA), 3L), "\n(",
         signif(data.matrix(filtered_sig), 1L), ")"),
  nrow(filtered_corr_NA), ncol(filtered_corr_NA)
)

grDevices::png(file.path(PATHS$wgcna_output,
                         "14_2_filtered_module_t-values_mixed_model.png"),
               units = "in", width = 30, height = 25, res = 300)
graphics::par(mar = c(7, 18, 5, 5))
WGCNA::labeledHeatmap(
  Matrix       = filtered_corr_NA,
  xLabels      = colnames(filtered_corr_NA),
  yLabels      = rownames(filtered_corr_NA),
  ySymbols     = rownames(filtered_corr_NA),
  font.lab.x = 2, font.lab.y = 2, cex.lab = 1.5,
  colorLabels  = FALSE, colors = WGCNA_COLORS,
  textMatrix   = textMatrix_filt,
  setStdMargins = FALSE, cex.text = 0.9, zlim = c(-100, 100),
  main = paste0("Filtered modules  |  t < ", t_cutoff, " or adj-p > 0.05 = grey")
)
grDevices::dev.off()


# =============================================================================
# SECTION 13: Module–module correlation
# =============================================================================

# Only significant modules: keep modules where at least one tumour has adj-p < 0.05
sig_module_names <- rownames(filtered_sig)[
  rowSums(filtered_sig < 0.05, na.rm = TRUE) > 0
]
ME_sig <- MEs[, sig_module_names, drop = FALSE]

# p-value matrix for corrplot
.cor_pmat <- function(mat) {
  mat  <- as.matrix(mat)
  n    <- ncol(mat)
  pmat <- matrix(NA_real_, n, n, dimnames = list(colnames(mat), colnames(mat)))
  diag(pmat) <- 0
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      tmp      <- stats::cor.test(mat[, i], mat[, j])
      pmat[i, j] <- pmat[j, i] <- tmp$p.value
    }
  }
  pmat
}
p_mat <- .cor_pmat(ME_sig)

col_fun <- grDevices::colorRampPalette(
  c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

grDevices::pdf(file.path(PATHS$wgcna_output,
                         "15_module-module_correlation_colour.pdf"),
               width = 30, height = 40)
corrplot::corrplot(
  stats::cor(ME_sig), method = "color", col = col_fun(200),
  type = "upper", order = "hclust", addCoef.col = "black",
  tl.col = "black", tl.srt = 45,
  p.mat = p_mat, sig.level = 0.01, insig = "blank", diag = FALSE
)
grDevices::dev.off()

# Euclidean distance-based dendrogram and MDS
distance <- stats::dist(t(MEs), method = "euclidean")
cluster  <- stats::hclust(distance, method = "average")

grDevices::pdf(file.path(PATHS$wgcna_output, "17_module_dendrogram.pdf"),
               width = 10, height = 20)
plot(cluster)
grDevices::dev.off()

MDS <- stats::cmdscale(distance, k = 2L)
grDevices::pdf(file.path(PATHS$wgcna_output, "18_module_MDS.pdf"),
               width = 10, height = 10)
plot(MDS, col = gsub("ME", "", rownames(MDS)))
graphics::text(MDS, labels = rownames(MDS), cex = 0.7, pos = 3L)
grDevices::dev.off()


# =============================================================================
# SECTION 14: Module membership violin plot
# =============================================================================

violin_data <- filtered_module_membership

p_violin <- ggplot2::ggplot(
  violin_data,
  ggplot2::aes(x = module_name, y = module_membership_value,
               fill = module_name)
) +
  ggplot2::geom_violin(alpha = 0.5) +
  ggplot2::geom_boxplot(width = 0.1, fill = "white") +
  ggplot2::geom_hline(yintercept = WGCNA_PARAMS$mm_cutoff, col = "red") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "none",
    text = ggplot2::element_text(size = 20),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
  )

x_labels <- ggplot2::ggplot_build(p_violin)$layout$panel_params[[1]]$x$get_labels()

grDevices::pdf(file.path(PATHS$wgcna_output,
                         "19_module_membership_violin_plot.pdf"),
               width = 40, height = 10)
print(p_violin +
        ggplot2::scale_fill_manual(values = gsub("ME", "", x_labels)))
grDevices::dev.off()


# =============================================================================
# SECTION 15: Hub gene identification
# =============================================================================

# Top hub gene per module
module_hub_genes <- data.frame(
  WGCNA::chooseTopHubInEachModule(
    norm_counts_t, dynamicColors,
    omitColors = "grey", power = power, type = networkType
  )
)
colnames(module_hub_genes) <- "gene_id"
module_hub_genes$module    <- rownames(module_hub_genes)
module_hub_genes           <- merge(module_hub_genes, annotation_enhanced,
                                    by = "gene_id", all.x = TRUE)

utils::write.table(
  module_hub_genes,
  file.path(PATHS$wgcna_output,
            paste0("21_module_hub_genes_power_", power,
                   "_networktype_", networkType, ".txt")),
  sep = "\t"
)

# All hub genes ranked by intra-modular adjacency (custom topHubs function)
topHubs <- function(datExpr, colorh, omitColors = "grey",
                    power = 2, type = "signed") {
  modules <- names(table(colorh))
  if (!is.na(omitColors)[1L]) modules <- modules[!modules %in% omitColors]

  connectivity_table <- data.frame(
    gene = character(), connectivity_rowSums_adj = double(),
    module = character()
  )
  for (m in modules) {
    adj    <- WGCNA::adjacency(datExpr[, colorh == m], power = power, type = type)
    sorted <- sort(rowSums(adj), decreasing = TRUE)
    connectivity_table <- rbind(
      connectivity_table,
      data.frame(
        gene                   = names(sorted),
        connectivity_rowSums_adj = sorted,
        module                 = m
      )
    )
  }
  connectivity_table
}

connectivity_table <- topHubs(
  norm_counts_t, dynamicColors,
  omitColors = "grey", power = power, type = networkType
)
colnames(connectivity_table)[1L] <- "gene_id"
connectivity_table <- merge(connectivity_table, annotation_enhanced,
                            by = "gene_id", all.x = TRUE) %>%
  dplyr::group_by(module) %>%
  dplyr::arrange(dplyr::desc(connectivity_rowSums_adj), .by_group = TRUE)

utils::write.table(
  connectivity_table,
  file.path(PATHS$wgcna_output,
            paste0("22_all_hub_genes_power_", power,
                   "_networktype_", networkType, ".txt")),
  sep = "\t"
)

filtered_mod_dir <- file.path(PATHS$wgcna_output, "filtered_mod_memberships")
if (!dir.exists(filtered_mod_dir)) dir.create(filtered_mod_dir, recursive = TRUE)

# Top 5 and top 25% hub genes for filtered modules
module_top5 <- connectivity_table %>%
  dplyr::group_by(module) %>%
  dplyr::top_n(5L, wt = connectivity_rowSums_adj)

module_top25pct <- connectivity_table %>%
  dplyr::group_by(module) %>%
  dplyr::filter(connectivity_rowSums_adj >
                  stats::quantile(connectivity_rowSums_adj, 0.75))

filt_mods_bare <- sub("ME", "", filtered_module_list)

filtered_top5    <- dplyr::filter(module_top5,    module %in% filt_mods_bare) %>%
  dplyr::mutate(module_name = paste0("ME", module))
filtered_top25   <- dplyr::filter(module_top25pct, module %in% filt_mods_bare) %>%
  dplyr::mutate(module_name = paste0("ME", module))

utils::write.table(filtered_top5,   file.path(PATHS$wgcna_output,
  paste0("22_filtered_top5_hub_genes_power_",    power, "_networktype_", networkType, ".txt")), sep = "\t")
utils::write.table(filtered_top25,  file.path(PATHS$wgcna_output,
  paste0("23_filtered_top25pct_hub_genes_power_", power, "_networktype_", networkType, ".txt")), sep = "\t")


# =============================================================================
# SECTION 16: lncRNA proportion per module
# =============================================================================

.count_gene_types <- function(tbl, label_col, power_str, network_str) {
  tbl2 <- tbl
  tbl2$gene_type[tbl2$gene_type != "lncRNA"] <- "non-lncRNA"
  counts <- tbl2 %>%
    dplyr::count(module, gene_type) %>%
    dplyr::mutate(per = 100 * n / sum(n))
  utils::write.table(counts,
    file.path(PATHS$wgcna_output,
              paste0(label_col, "_genetype_count_power_", power_str,
                     "_networktype_", network_str, ".txt")),
    sep = "\t")
  counts
}

genetype_all  <- .count_gene_types(
  dplyr::filter(connectivity_table, module %in% filt_mods_bare),
  "24_module_specific", power, networkType)
genetype_hub  <- .count_gene_types(
  filtered_top25, "24_module_hub", power, networkType)

.barplot_genetype <- function(data, title_str, file_path) {
  grDevices::png(file_path, units = "in", width = 15, height = 5, res = 300)
  print(
    ggplot2::ggplot(data.frame(data),
                    ggplot2::aes(fill = gene_type, y = per, x = module)) +
      ggplot2::geom_bar(position = "fill", stat = "identity") +
      ggplot2::labs(title = title_str) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  )
  grDevices::dev.off()
}

.barplot_genetype(genetype_all, "Module gene-type proportion",
  file.path(PATHS$wgcna_output,
            paste0("24_module_genetype_pct_power_", power,
                   "_networktype_", networkType, ".png")))
.barplot_genetype(genetype_hub, "Hub gene-type proportion",
  file.path(PATHS$wgcna_output,
            paste0("24_hub_genetype_pct_power_", power,
                   "_networktype_", networkType, ".png")))


# =============================================================================
# SECTION 17: Per-module enrichment (GO, KEGG, NCI, TF)
# =============================================================================

# Load DGE master list and universe
DGE_master <- readxl::read_excel(GRN_FILES$dge_master, sheet = "vs_all_tumor")
module_universe_raw <- utils::read.delim(
  file.path(PATHS$wgcna_output, "module_membership_universe.txt"), header = FALSE)
module_universe <- annotation_enhanced[
  annotation_enhanced$gene_id %in% module_universe_raw[, 1L],
  c("gene_id", "gene_type", "gene_name", "ensembl_gene_id", "entrezgene_id")
]
ent_universe <- module_universe$entrezgene_id

# Load GRN for TF enrichment (threshold 0.03 for WGCNA)
TRN_wgcna <- utils::read.delim(GRN_FILES$grn_wgcna, sep = "\t")
colnames(TRN_wgcna)[1:2] <- c("tfs", "targets")
TRN_wgcna <- stats::aggregate(targets ~ tfs, TRN_wgcna, toString)
tf_universe <- length(unique(module_universe$ensembl_gene_id))

DGE_common_cols <- c("gene_id", "gene_name", "ensembl_gene_id",
                     "entrezgene_id", "gene_type", "Novel_gene")
modules_to_run  <- names(table(filtered_module_membership$module))

master_TF_res      <- NULL
hub_TF_res         <- NULL

.save_enrichment_plots <- function(obj, dir, prefix, title_suffix, showCat) {
  if (is.null(obj) || length(obj@result$ID) == 0L) return()
  grDevices::png(file.path(dir, paste0(prefix, "_bar.png")),
                 units = "in", width = 8, height = 5, res = 300)
  try(print(enrichplot::barplot(obj, showCategory = showCat) +
              ggplot2::ggtitle(title_suffix)))
  grDevices::dev.off()
  grDevices::png(file.path(dir, paste0(prefix, "_dot.png")),
                 units = "in", width = 6, height = 7, res = 300)
  try(print(enrichplot::dotplot(obj, showCategory = showCat) +
              ggplot2::ggtitle(title_suffix)))
  grDevices::dev.off()
}

for (mod_color in modules_to_run) {

  mod_dir <- file.path(filtered_mod_dir, mod_color)
  if (!dir.exists(mod_dir)) dir.create(mod_dir, recursive = TRUE)

  # Relevant tumour types for this module (t > cutoff)
  relevant_tumors <- colnames(filtered_corr_NA)[
    !is.na(filtered_corr_NA[paste0("ME", mod_color), ]) &
    filtered_corr_NA[paste0("ME", mod_color), ] > t_cutoff
  ]

  # Module gene IDs
  mod_gene_ids <- filtered_module_membership[
    filtered_module_membership$module == mod_color, "ID"]
  utils::write.table(mod_gene_ids,
    file.path(mod_dir, paste("sb", mod_color, "genes.txt", sep = "_")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Module-specific DGE rows
  mod_DGE <- DGE_master[DGE_master$gene_id %in% mod_gene_ids,
                         c(DGE_common_cols, relevant_tumors)]
  utils::write.csv(as.data.frame(mod_DGE), file.path(mod_dir, "1_DGE_list.csv"))

  # Universe subset for this module
  mod_members <- module_universe[module_universe$gene_id %in% mod_gene_ids, ]

  # ── GO enrichment ────────────────────────────────────────────────────────
  mod_eGO <- try(create_eGO(mod_members$entrezgene_id, ent_universe))
  if (!inherits(mod_eGO, "try-error") && !is.null(mod_eGO)) {
    mod_eGO <- mod_eGO %>% dplyr::filter(Count >= 5L)
    if (length(mod_eGO@result$ID) > 0L) {
      utils::write.csv(mod_eGO, file.path(mod_dir, "2_GO_terms.csv"))
      utils::write.csv(mod_eGO@result[, "p.adjust", drop = FALSE],
                       file.path(mod_dir, "2_GO_terms_for_revigo.csv"))
      .save_enrichment_plots(mod_eGO, mod_dir, "2_GO",
                             paste("GO —", mod_color), ENRICH_PARAMS$show_categories)
    }
  }

  # ── KEGG enrichment ───────────────────────────────────────────────────────
  mod_eKEGG <- try(create_eKEGG(mod_members$entrezgene_id))
  if (!inherits(mod_eKEGG, "try-error") && !is.null(mod_eKEGG)) {
    mod_eKEGG <- mod_eKEGG %>% dplyr::filter(Count >= 5L)
    if (length(mod_eKEGG@result$ID) > 0L) {
      utils::write.csv(mod_eKEGG, file.path(mod_dir, "3_KEGG_pathways.csv"))
      .save_enrichment_plots(mod_eKEGG, mod_dir, "3_KEGG",
                             paste("KEGG —", mod_color), ENRICH_PARAMS$show_categories)
    }
  }

  # ── NCI gene-set enrichment ───────────────────────────────────────────────
  mod_NCI <- try(create_NCI_enrichment(mod_members$gene_name,
                                       module_universe$gene_name))
  if (!inherits(mod_NCI, "try-error") && !is.null(mod_NCI)) {
    mod_NCI <- mod_NCI %>% dplyr::filter(Count >= 5L)
    if (length(mod_NCI@result$ID) > 0L) {
      utils::write.csv(mod_NCI, file.path(mod_dir, "4_NCI_terms.csv"))
      .save_enrichment_plots(mod_NCI, mod_dir, "4_NCI",
                             paste("NCI —", mod_color), ENRICH_PARAMS$show_categories)
    }
  }

  # ── TF enrichment (all modules; require_tf_in_list = FALSE) ───────────────
  tf_res <- getTRN_enrichment(TRN_wgcna, mod_members$ensembl_gene_id,
                              tf_universe, mod_color,
                              require_tf_in_list = FALSE)
  if (nrow(tf_res) > 0L) {
    tf_DGE <- DGE_master[DGE_master$ensembl_gene_id %in% tf_res$TF,
                          c(DGE_common_cols, relevant_tumors)]
    utils::write.csv(as.data.frame(tf_DGE), file.path(mod_dir, "5_TF_DGE_list.csv"))
  }
  master_TF_res <- rbind(master_TF_res, tf_res)

  # ── Hub gene sub-analysis ─────────────────────────────────────────────────
  hub_dir <- file.path(mod_dir, "hub_genes")
  if (!dir.exists(hub_dir)) dir.create(hub_dir)

  hub_ids <- filtered_top25[filtered_top25$module == mod_color, "gene_id"]
  utils::write.table(hub_ids,
    file.path(hub_dir, paste("sb", mod_color, "hub_genes.txt", sep = "_")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  hub_members <- module_universe[module_universe$gene_id %in% hub_ids$gene_id, ]

  hub_eGO <- try(create_eGO(hub_members$entrezgene_id, ent_universe))
  if (!inherits(hub_eGO, "try-error") && !is.null(hub_eGO)) {
    hub_eGO <- hub_eGO %>% dplyr::filter(Count >= 5L)
    if (length(hub_eGO@result$ID) > 0L) {
      utils::write.csv(hub_eGO, file.path(hub_dir, "2_GO_terms.csv"))
      .save_enrichment_plots(hub_eGO, hub_dir, "2_GO",
                             paste("Hub GO —", mod_color), ENRICH_PARAMS$show_categories)
    }
  }

  hub_tf_res <- getTRN_enrichment(TRN_wgcna, hub_members$ensembl_gene_id,
                                  tf_universe, mod_color,
                                  require_tf_in_list = FALSE)
  hub_TF_res <- rbind(hub_TF_res, hub_tf_res)

  message("  Module ", mod_color, " done.")
}


# =============================================================================
# SECTION 18: TF enrichment heatmaps (modules and hub genes)
# =============================================================================

.write_tf_master_and_heatmap <- function(tf_res_raw, tf_db_path,
                                         out_prefix, fdr_cut, overlap_cut) {
  if (is.null(tf_res_raw) || nrow(tf_res_raw) == 0L) return(invisible(NULL))

  TF_db <- as.data.frame(readxl::read_excel(tf_db_path))
  TF_db <- TF_db[TF_db[["Is TF?"]] == "Yes",
                 c("Ensembl ID", "HGNC symbol", "EntrezGene Description")]
  colnames(TF_db)[1:2] <- c("TF", "HGNC.symbol")

  tf_annot <- plyr::join(
    tf_res_raw[as.double(tf_res_raw$p) |> is.finite(), ],
    TF_db[TF_db$TF %in% tf_res_raw$TF, ],
    by = "TF"
  )
  tf_annot$overlap   <- as.numeric(tf_annot$overlap)
  tf_annot$tf.size   <- as.numeric(tf_annot$tf.size)
  tf_annot$universe  <- as.numeric(tf_annot$universe)
  tf_filt <- tf_annot[!is.na(tf_annot$fdr) & tf_annot$fdr < fdr_cut, ]

  utils::write.table(tf_filt, paste0(out_prefix, "_TF_enrichment.txt"),
                     row.names = FALSE, quote = FALSE)

  # Heatmap of -log(FDR) for TFs with overlap >= overlap_cut
  heat_data <- tf_filt[tf_filt$overlap >= overlap_cut, ]
  if (nrow(heat_data) == 0L) return(invisible(NULL))

  mat <- reshape2::acast(heat_data, HGNC.symbol ~ set, value.var = "fdr")
  mat[is.na(mat)] <- 1
  mat[mat == 0]   <- 1e-100
  mat <- -log(mat)

  col_fun <- circlize::colorRamp2(c(0, 50, 100),
                                   c("white", "red", "darkred"))

  for (fmt in c("pdf", "png")) {
    if (fmt == "pdf") {
      grDevices::pdf(paste0(out_prefix, "_TF_heatmap.pdf"), width = 10, height = 15)
    } else {
      grDevices::png(paste0(out_prefix, "_TF_heatmap.png"),
                     units = "in", width = 10, height = 15, res = 300)
    }
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        mat, col = col_fun,
        row_names_gp    = grid::gpar(fontsize = 18),
        column_names_gp = grid::gpar(fontsize = 18)
      )
    )
    grDevices::dev.off()
  }
}

.write_tf_master_and_heatmap(
  master_TF_res, GRN_FILES$tf_db,
  file.path(PATHS$wgcna_output, "20_module"),
  GRN_PARAMS$fdr_cutoff, GRN_PARAMS$overlap_heatmap
)
.write_tf_master_and_heatmap(
  hub_TF_res, GRN_FILES$tf_db,
  file.path(PATHS$wgcna_output, "20_module_hub_genes"),
  GRN_PARAMS$fdr_cutoff, GRN_PARAMS$overlap_heatmap
)


# =============================================================================
# SECTION 19: Save workspace
# =============================================================================

save.image(file = file.path(PATHS$wgcna_root,
                            paste0("WGCNA_workspace_power_", power, ".RData")))
message("WGCNA analysis complete. Workspace saved.")
