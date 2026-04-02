#!/bin/bash
set -e

# Load configuration file
CONFIG_FILE="$1"
source "$CONFIG_FILE"
shift

# Required variables check
: "${SPLIT_PDBS_SCRIPT:?not set}"
: "${SUPERIMPOSE_SCRIPT:?not set}"
: "${FPOCKET_BIN:?not set}"
: "${MAIN_ANALYSIS_DIR:?not set}"
# fallback
: "${N_JOBS:=1}"

# Load conda
# adjust this line depending on your conda installation
source ~/anaconda3/etc/profile.d/conda.sh

# Optional argument (residue range for secondary PRS)
RANGE_AB="$1"

# Activate main environment
conda activate ensembleflex

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

BASE_DIR=$(pwd)
OUT_DIR="${BASE_DIR}/out"
mkdir -p "$OUT_DIR"

# Detect whether input contains multiple models
MODEL_COUNT=$(grep -h "^MODEL" *.pdb 2>/dev/null | wc -l)

if [ "$MODEL_COUNT" -le 1 ]; then
  SINGLE_MODEL=true
else
  SINGLE_MODEL=false
fi

# =============================
# STEP1: Split PDB (if multi-model)
# =============================
if [ "$SINGLE_MODEL" = false ]; then
  mkdir -p "$OUT_DIR/split_pdbs"
  Rscript "$SPLIT_PDBS_SCRIPT" \
    -i "$BASE_DIR/" \
    -o "$OUT_DIR/split_pdbs/"
fi

# =============================
# STEP2: Structural superimposition
# =============================
mkdir -p "$OUT_DIR/superimposed"

if [ "$SINGLE_MODEL" = false ]; then
  Rscript "$SUPERIMPOSE_SCRIPT" \
    -i "$OUT_DIR/split_pdbs/" \
    -o "$OUT_DIR/"
else
  cp *.pdb "$OUT_DIR/superimposed/"
fi

# =============================
# STEP3: Conformer analysis
# =============================
if [ "$SINGLE_MODEL" = false ]; then
  Rscript "$MAIN_ANALYSIS_DIR/1_run_conformer.R" \
    "$OUT_DIR/superimposed/" \
    "$OUT_DIR/result_1_conformer/" \
    "$OUT_DIR/"
fi

# =============================
# STEP4: ANM analysis
# =============================
if [ "$SINGLE_MODEL" = false ]; then
  python "$MAIN_ANALYSIS_DIR/2_run_ANM.py" \
    "$OUT_DIR/superimposed/" \
    "$OUT_DIR/result_2_ANM"
fi

# =============================
# STEP5: PRS analysis
# =============================
python "$MAIN_ANALYSIS_DIR/3_run_PRS.py" \
  "$OUT_DIR/superimposed" \
  "$OUT_DIR/result_3_PRS" \
  --n_jobs "$N_JOBS"

# =============================
# STEP5.5: Pocket detection (fpocket)
# =============================
echo "[STEP 5.5] fpocket"

conda deactivate
conda activate freesasa

export FPOCKET_BIN

bash "$MAIN_ANALYSIS_DIR/4_bash_fpocket_multi.sh" \
  "$OUT_DIR/superimposed" \
  "$OUT_DIR/result_4_fpockets"

# =============================
# STEP6: SASA calculation (FreeSASA)
# =============================
TMP="$OUT_DIR/tmp.tsv"
FINAL="$OUT_DIR/result_5_freesasa.tsv"

echo -e "conformer\tchain\tresid\tresname\tASA" > "$FINAL"
> "$TMP"

for pdb in "$OUT_DIR/superimposed/"*.pdb; do
    base=$(basename "$pdb" .pdb)

    fpocket_dir="$OUT_DIR/result_4_fpockets/${base}_out"
    if [ ! -d "$fpocket_dir" ]; then
        continue
    fi

    python "$MAIN_ANALYSIS_DIR/5a_run_freesasa.py" "$pdb" \
      --fpocket \
      --fpocket_out "$fpocket_dir" >> "$TMP"

done

cat "$TMP" >> "$FINAL"
rm "$TMP"

# =============================
# STEP7: Integrate PRS and SASA
# =============================
python "$MAIN_ANALYSIS_DIR/6_run_prs_sasa.py" \
  "$OUT_DIR/result_3_PRS/PRS_per_conformer.tsv" \
  "$FINAL" \
  --outdir "$OUT_DIR"

# =============================
# STEP8: Secondary PRS analysis
# =============================
conda deactivate
conda activate ensembleflex

if [ -n "$RANGE_AB" ]; then
  python "$MAIN_ANALYSIS_DIR/7_run_secondary_prs.py" \
    "$OUT_DIR/superimposed/" \
    "$OUT_DIR/result_3_PRS/PRS_per_conformer.tsv" \
    "$OUT_DIR/" \
    -r "$RANGE_AB"
else
  python "$MAIN_ANALYSIS_DIR/7_run_secondary_prs.py" \
    "$OUT_DIR/superimposed/" \
    "$OUT_DIR/result_3_PRS/PRS_per_conformer.tsv" \
    "$OUT_DIR/"
fi

echo "DONE"