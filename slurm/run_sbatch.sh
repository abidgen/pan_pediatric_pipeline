#!/bin/bash
## =============================================================================
## FILE: slurm/run_sbatch.sh
## PURPOSE: Master launcher — submits the full pipeline in dependency order.
##
## PIPELINE ORDER:
##   Stage 1 (manual prerequisite, run once):
##     1_library_type_batch_correction.R   — interactive or separate sbatch
##
##   Stage 2 — DGE (7 parallel jobs, 5 tumours each):
##     2_analyze_samples_2025_02_21.R      — via analyze_samples_2025_02_21.sh
##
##   Stage 3 — runs after ALL Stage 2 jobs succeed:
##     3_TPM_calculation_RSEM_f.R  (standalone, no DGE dependency)
##     4_GO_analysis.R
##     5_TFs_enrichments_heatmap.R
##     6_GRN_TF_enrichment.R       — via grn_tf_enrichment.sh
##
##   Stage 4 — independent long-running job:
##     7_WGCNA.R                   — via wgcna.sh
##
## USAGE:
##   bash slurm/run_sbatch.sh
##
## SLURM dependencies use --dependency=afterok:<jobid> to enforce ordering.
## =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/../logs"
mkdir -p "$LOG_DIR"

# ── Resources ─────────────────────────────────────────────────────────────────
DGE_MEM="256g"; DGE_CPUS=8;  DGE_SCRATCH="lscratch:50";  DGE_TIME="12:00:00"
GRN_MEM="16g";  GRN_CPUS=2;  GRN_TIME="02:00:00"
WGC_MEM="256g"; WGC_CPUS=16; WGC_SCRATCH="lscratch:100"; WGC_TIME="24:00:00"

# ── Tumor-type master list (34 types) ─────────────────────────────────────────
SAMPLES=(
    "ACC"       "ACPG"      "ASPS"      "ATRT"      "CCSK"
    "CPC"       "DSRCT"     "EPMT"      "EWS"       "GB"
    "GLI"       "GNG"       "HBL"       "HGNET"     "MBL.G3"
    "MBL.G4"    "MBL.SHH"   "MBL.U"     "MBL.WNT"   "MEL"
    "NB.MYCN.A" "NB.MYCN.NA" "OS"       "PAST"      "RBL"
    "RMS.FN"    "RMS.FP"    "RMS.SS"    "RT"        "SPZM"
    "SS"        "THPA"      "WT"        "WTB"
)
TOTAL=${#SAMPLES[@]}    # 34
BATCH_SIZE=5

# ── Helper: extract job ID from sbatch stdout ──────────────────────────────────
get_jobid() { echo "$1" | grep -oP '\d+$'; }

# =============================================================================
# STAGE 2: DGE — 7 parallel jobs
# =============================================================================
printf '\n[STAGE 2] Submitting DGE jobs (%d tumours, batch size %d)...\n' \
  "$TOTAL" "$BATCH_SIZE"

DGE_JOB_IDS=()
JOB_COUNT=0

for (( i = 0; i < TOTAL; i += BATCH_SIZE )); do
    BATCH=("${SAMPLES[@]:$i:$BATCH_SIZE}")
    LIST_STRING=$(IFS=,; echo "${BATCH[*]}")
    JOB_COUNT=$(( JOB_COUNT + 1 ))

    OUT=$(sbatch \
        --mem="$DGE_MEM" --cpus-per-task="$DGE_CPUS" \
        --gres="$DGE_SCRATCH" --time="$DGE_TIME" \
        "$SCRIPT_DIR/analyze_samples_2025_02_21.sh" "$LIST_STRING")

    JID=$(get_jobid "$OUT")
    DGE_JOB_IDS+=("$JID")
    printf '  Job %d submitted (ID %s) | tumours: %s\n' \
      "$JOB_COUNT" "$JID" "$LIST_STRING"
done

# Build afterok dependency string: afterok:id1:id2:...
DGE_DEP="afterok:$(IFS=:; echo "${DGE_JOB_IDS[*]}")"
printf '[STAGE 2] %d DGE jobs submitted. Downstream dependency: %s\n' \
  "${#DGE_JOB_IDS[@]}" "$DGE_DEP"

# =============================================================================
# STAGE 3: GRN TF enrichment — depends on all Stage 2 jobs
# =============================================================================
printf '\n[STAGE 3] Submitting GRN TF enrichment (after DGE)...\n'

GRN_OUT=$(sbatch \
    --dependency="$DGE_DEP" \
    --mem="$GRN_MEM" --cpus-per-task="$GRN_CPUS" --time="$GRN_TIME" \
    "$SCRIPT_DIR/grn_tf_enrichment.sh")
GRN_JID=$(get_jobid "$GRN_OUT")
printf '  GRN job submitted (ID %s)\n' "$GRN_JID"

# =============================================================================
# STAGE 4: WGCNA — independent, long-running (no DGE dependency)
# =============================================================================
printf '\n[STAGE 4] Submitting WGCNA (independent)...\n'

WGC_OUT=$(sbatch \
    --mem="$WGC_MEM" --cpus-per-task="$WGC_CPUS" \
    --gres="$WGC_SCRATCH" --time="$WGC_TIME" \
    "$SCRIPT_DIR/wgcna.sh")
WGC_JID=$(get_jobid "$WGC_OUT")
printf '  WGCNA job submitted (ID %s)\n' "$WGC_JID"

# =============================================================================
# Summary
# =============================================================================
printf '\n=== Pipeline submitted ===\n'
printf '  Stage 2 DGE  : %d jobs (%s)\n' \
  "${#DGE_JOB_IDS[@]}" "$(IFS=,; echo "${DGE_JOB_IDS[*]}")"
printf '  Stage 3 GRN  : job %s (after DGE)\n' "$GRN_JID"
printf '  Stage 4 WGCNA: job %s (independent)\n' "$WGC_JID"
