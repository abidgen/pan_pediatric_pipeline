#!/bin/bash
## =============================================================================
## FILE: slurm/analyze_samples_2025_02_21.sh
## PURPOSE: SLURM worker — runs analyze_samples_2025_02_21.R for one batch
##          of comma-separated tumor types.
##
## CALLED BY: run_sbatch.sh (never submit this script directly).
##
## ARGUMENT:
##   $1  Comma-separated list of tumor-type labels, e.g. "ACC,ACPG,ATRT"
##
## RESOURCE DEFAULTS (overridden by run_sbatch.sh via sbatch flags):
##   --mem=256g  --cpus-per-task=8  --gres=lscratch:50  --time=12:00:00
##
## ENVIRONMENT (confirmed from SLURM logs):
##   R 4.4.2  |  Bioconductor 3.20  |  8 cores
## =============================================================================
#SBATCH --job-name=dge_analysis
#SBATCH --output=../logs/slurm-%j.out
#SBATCH --error=../logs/slurm-%j.err

set -euo pipefail   # Exit on error, unset variable, or pipe failure

# ── Load pinned R module ───────────────────────────────────────────────────────
module load R/4.4.2

# ── Validate argument ──────────────────────────────────────────────────────────
TUMOR_LIST="${1:-}"
if [[ -z "$TUMOR_LIST" ]]; then
  echo "ERROR: No tumor-type list supplied." >&2
  echo "Usage: sbatch $0 \"ACC,ACPG,ATRT\"" >&2
  exit 1
fi

printf '\n[%s] Job %s started  |  tumors: %s\n' \
  "$(date '+%Y-%m-%d %H:%M:%S')" "$SLURM_JOB_ID" "$TUMOR_LIST"

# ── Run analysis (cd so that relative source() paths resolve correctly) ────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../analysis"

Rscript 2_analyze_samples_2025_02_21.R "$TUMOR_LIST"

printf '\n[%s] Job %s completed successfully.\n' \
  "$(date '+%Y-%m-%d %H:%M:%S')" "$SLURM_JOB_ID"
