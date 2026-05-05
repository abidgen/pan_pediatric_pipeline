#!/bin/bash
## =============================================================================
## FILE: slurm/wgcna.sh
## PURPOSE: SLURM worker for the WGCNA analysis (7_WGCNA.R).
##          Single-job, no batching required — WGCNA is a sequential workflow.
##
## RESOURCES: 256 GB RAM for TOM construction on ~20k genes x ~500 samples.
##            16 CPUs: 8 for WGCNA threads + 8 for DESeq2 BiocParallel.
##
## USAGE (called from run_sbatch.sh — do not submit directly):
##   sbatch --mem=256g --cpus-per-task=16 --gres=lscratch:100 \
##          --time=24:00:00 slurm/wgcna.sh
## =============================================================================
#SBATCH --job-name=wgcna
#SBATCH --output=../logs/slurm-wgcna-%j.out
#SBATCH --error=../logs/slurm-wgcna-%j.err

set -euo pipefail

module load R/4.4.2

printf '\n[%s] WGCNA job %s started\n' \
  "$(date '+%Y-%m-%d %H:%M:%S')" "$SLURM_JOB_ID"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../analysis"

Rscript 7_WGCNA.R

printf '\n[%s] WGCNA job %s completed.\n' \
  "$(date '+%Y-%m-%d %H:%M:%S')" "$SLURM_JOB_ID"
