#!/bin/bash
# Usage: bash 4_bash_fpocket_multi.sh <input_dir> <output_dir>

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Require fpocket binary path
: "${FPOCKET_BIN:?not set}"

mkdir -p "$OUTPUT_DIR"

# Temporary directory for fpocket execution
TMP_DIR="${OUTPUT_DIR}/tmp_fpocket"
mkdir -p "$TMP_DIR"

for pdb in "$INPUT_DIR"/*.pdb; do
    base=$(basename "$pdb")

    cp "$pdb" "$TMP_DIR/"

    "$FPOCKET_BIN" -f "$TMP_DIR/$base"

    mv "$TMP_DIR/${base%.pdb}_out" "$OUTPUT_DIR/"

    rm "$TMP_DIR/$base"
done

# Remove temporary directory
rmdir "$TMP_DIR"