## =============================================================================
## FILE: 00_packages.R
## PURPOSE: Central package loader with pinned versions.
##          Source this file ONCE from any downstream script instead of
##          scattering library() calls across multiple files.
##
## ENVIRONMENT (confirmed from SLURM logs):
##   R            4.4.2  (2024-10-31)
##   Bioconductor 3.20
##   BiocManager  1.30.25
##
## VERSIONS EXTRACTED FROM SLURM LOGS (slurm-4870927{2-9}.out):
##   Directly observed:  tidyverse 2.0.0, dplyr 1.1.4, ggplot2 3.5.1,
##                       readr 2.1.5, tibble 3.2.1, tidyr 1.3.1,
##                       purrr 1.0.4, forcats 1.0.0, lubridate 1.9.4,
##                       stringr 1.5.1, BiocManager 1.30.25,
##                       clusterProfiler 4.14.4, DOSE 4.0.0,
##                       enrichplot 1.26.6, pcaExplorer 3.0.0
##   Inferred from Bioc 3.20 / R 4.4.2 environment:
##                       DESeq2 1.46.0, edgeR 4.4.0, limma 3.62.0,
##                       sva 3.54.0, BiocParallel 1.40.0,
##                       biomaRt 2.62.0, AnnotationHub 3.14.0,
##                       org.Hs.eg.db 3.20.0, ReactomePA 1.50.0,
##                       ComplexHeatmap 2.22.0
## =============================================================================

# ── Bootstrap BiocManager ──────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Lock Bioconductor release to the version used during development
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)

# ── Version registry ───────────────────────────────────────────────────────────
# Centralise all pinned versions here; downstream code never hard-codes them.
PKG_VERSIONS <- list(
  # CRAN
  tidyverse    = "2.0.0",
  dplyr        = "1.1.4",
  ggplot2      = "3.5.1",
  readr        = "2.1.5",
  tibble       = "3.2.1",
  tidyr        = "1.3.1",
  purrr        = "1.0.4",
  forcats      = "1.0.0",
  lubridate    = "1.9.4",
  stringr      = "1.5.1",
  plotly       = "4.10.4",
  htmltools    = "0.5.8",
  htmlwidgets  = "1.6.4",
  pheatmap     = "1.0.12",
  readxl       = "1.4.3",
  reshape2     = "1.4.4",
  plyr         = "1.8.9",
  circlize     = "0.4.16",
  gprofiler2   = "0.2.3",
  ashr         = "2.2.63",
  # WGCNA + mixed-model stack (CRAN)
  WGCNA        = "1.73",
  CorLevelPlot = "0.99.2",
  gridExtra    = "2.3",
  MASS         = "7.3.60",
  lmerTest     = "3.1.3",
  lme4         = "1.1.35",
  RColorBrewer = "1.1.3",
  dichromat    = "2.0.0.1",
  corrplot     = "0.92",
  # Bioconductor (Bioc 3.20)
  DESeq2           = "1.46.0",
  edgeR            = "4.4.0",
  limma            = "3.62.0",
  sva              = "3.54.0",
  BiocParallel     = "1.40.0",
  biomaRt          = "2.62.0",
  AnnotationHub    = "3.14.0",
  pcaExplorer      = "3.0.0",
  clusterProfiler  = "4.14.4",
  DOSE             = "4.0.0",
  enrichplot       = "1.26.6",
  org.Hs.eg.db     = "3.20.0",
  ReactomePA       = "1.50.0",
  ComplexHeatmap   = "2.22.0",
  DGEobj.utils     = "1.0.6"
)


# ── Install helper ─────────────────────────────────────────────────────────────
# Installs a package only if it is absent; does NOT downgrade an existing
# installation so that the HPC system library is respected.
.install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing) > 0L) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    installer(missing)
  }
}

# ── CRAN packages ──────────────────────────────────────────────────────────────
cran_pkgs <- c(
  "tidyverse", "plotly", "htmltools", "htmlwidgets",
  "pheatmap", "readxl", "reshape2", "plyr",
  "circlize", "gprofiler2", "ashr", "stringr",
  "WGCNA", "CorLevelPlot", "gridExtra", "MASS",
  "lmerTest", "lme4", "RColorBrewer", "dichromat", "corrplot"
)
.install_if_missing(cran_pkgs, install.packages)

# ── Bioconductor packages ──────────────────────────────────────────────────────
bioc_pkgs <- c(
  "DESeq2", "edgeR", "limma", "sva", "BiocParallel",
  "biomaRt", "AnnotationHub", "pcaExplorer",
  "clusterProfiler", "DOSE", "enrichplot",
  "org.Hs.eg.db", "ReactomePA", "ComplexHeatmap",
  "DGEobj.utils"
)
.install_if_missing(bioc_pkgs, BiocManager::install)

# ── Version assertion helper ───────────────────────────────────────────────────
# Called after loading to warn if the runtime version differs from the pin.
assert_pkg_version <- function(pkg, expected) {
  actual <- as.character(packageVersion(pkg))
  if (actual != expected) {
    warning(sprintf(
      "[version mismatch] %s: expected %s, got %s", pkg, expected, actual
    ))
  }
}

# ── Load all packages ──────────────────────────────────────────────────────────
all_pkgs <- c(cran_pkgs, bioc_pkgs)
invisible(lapply(all_pkgs, library, character.only = TRUE))

# ── Assert versions post-load ──────────────────────────────────────────────────
# Only checks packages whose versions were directly confirmed in SLURM logs.
for (pkg in c("dplyr", "ggplot2", "readr", "tibble", "tidyr",
              "purrr", "forcats", "lubridate", "stringr",
              "clusterProfiler", "DOSE", "enrichplot", "pcaExplorer")) {
  assert_pkg_version(pkg, PKG_VERSIONS[[pkg]])
}

# ── Parallelisation ────────────────────────────────────────────────────────────
# Register 8 cores for BiocParallel (matches SLURM job allocation).
BiocParallel::register(BiocParallel::MulticoreParam(8))

message("00_packages.R loaded successfully.")
