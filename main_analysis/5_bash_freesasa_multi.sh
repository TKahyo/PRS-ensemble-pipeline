#!/bin/bash
# Usage: bash 5_bash_freesasa_multi.sh <pdb_dir> <mol2_dir> <out_dir>

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <pdb_dir> <mol2_dir> <out_dir>"
    exit 1
fi

PDB_DIR="$1"
MOL2_DIR="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

# Temporary file for aggregated results
tmp_all="$OUT_DIR/tmp_freesasa.tsv"
> "$tmp_all"

# Loop over all PDB files
for pdb_file in "$PDB_DIR"/*.pdb; do
    base=$(basename "$pdb_file" .pdb)
    mol2_file="$MOL2_DIR/$base/pocket0.mol2"

    # Skip if MOL2 file is missing
    if [ ! -f "$mol2_file" ]; then
        echo "Warning: MOL2 file not found for $base, skipping"
        continue
    fi

    echo "Processing $base ..."

    # Append FreeSASA results
    python3 ./main_analysis/5a_run_freesasa.py \
        "$pdb_file" "$mol2_file" >> "$tmp_all"
done

# Write final TSV
final_tsv="$OUT_DIR/result_freesasa.tsv"
echo -e "conformer\tchain\tresid\tresname\tASA" > "$final_tsv"
cat "$tmp_all" >> "$final_tsv"

# Cleanup
rm "$tmp_all"

echo "All done. Final TSV: $final_tsv"