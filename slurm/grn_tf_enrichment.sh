#!/bin/bash
## =============================================================================
## FILE: slurm/grn_tf_enrichment.sh
## PURPOSE: SLURM worker for GENIE3 GRN TF enrichment (6_GRN_TF_enrichment.R).
##          Must run AFTER all DGE jobs (2_analyze_samples_2025_02_21.R) finish.
##
## RESOURCES: Single-threaded R loop over tumour directories; 16 GB is enough.
##
## USAGE (called from run_sbatch.sh — do not submit directly):
##   sbatch --mem=16g --cpus-per-task=2 --time=02:00:00 \
##          slurm/grn_tf_enrichment.sh
## =============================================================================
#SBATCH --job-name=grn_tf
#SBATCH --output=../logs/slurm-grn-%j.out
#SBATCH --error=../logs/slurm-grn-%j.err

set -euo pipefail

module load R/4.4.2

printf '\n[%s] GRN TF enrichment job %s started\n' \
  "$(date '+%Y-%m-%d %H:%M:%S')" "$SLURM_JOB_ID"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../analysis"

Rscript 6_GRN_TF_enrichment.R

printf '\n[%s] GRN TF enrichment job %s completed.\n' \
  "$(date '+%Y-%m-%d %H:%M:%S')" "$SLURM_JOB_ID"
