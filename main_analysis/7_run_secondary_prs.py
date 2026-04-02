#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
from glob import glob
from Bio.PDB import PDBParser, DSSP
import matplotlib.pyplot as plt
import seaborn as sns


def classify_ss(ss_code: str) -> str:
    """Map DSSP secondary structure code into Helix / Sheet / Other."""
    if ss_code in ("H", "G", "I"):
        return "Helix"
    elif ss_code in ("E", "B"):
        return "Sheet"
    else:
        return "Other"


def aggregate_for_one_conformer(conformer, pdb_file, prs_df, parser, col_resnums_set=None):
    """
    Process one conformer:
      - assign DSSP secondary structure
      - aggregate PRS by secondary structure
      - optionally restrict residues by range
    """
    if conformer not in prs_df.index:
        print(f"Warning: {conformer} not in PRS TSV. Skipping.")
        return None, []

    structure = parser.get_structure(conformer, pdb_file)
    model = structure[0]

    # DSSP (requires mkdssp)
    dssp = DSSP(model, pdb_file, dssp="mkdssp")

    # Map residue number → secondary structure class
    ss_map = {}
    for (chain_id, resid) in dssp.keys():
        resnum = resid[1]
        ss_code = dssp[(chain_id, resid)][2]
        ss_map[resnum] = classify_ss(ss_code)

    prs_values = prs_df.loc[conformer]

    helix_vals, sheet_vals, other_vals = [], [], []
    residue_rows = []

    for col in prs_df.columns:
        try:
            resnum = int(col)
        except ValueError:
            raise ValueError(
                f"PRS TSV column '{col}' must be an integer residue number."
            )

        if col_resnums_set is not None and resnum not in col_resnums_set:
            continue

        prs_val = prs_values[col]
        ss = ss_map.get(resnum, "Other")

        if ss == "Helix":
            helix_vals.append(prs_val)
        elif ss == "Sheet":
            sheet_vals.append(prs_val)
        else:
            other_vals.append(prs_val)

        residue_rows.append({
            "Conformer": conformer,
            "Residue": resnum,
            "Secondary": ss,
            "PRS": prs_val,
        })

    # Mean PRS per residue
    helix_mean_per_conformer = sum(helix_vals) / len(helix_vals) if helix_vals else float("nan")
    sheet_mean_per_conformer = sum(sheet_vals) / len(sheet_vals) if sheet_vals else float("nan")
    other_mean_per_conformer = sum(other_vals) / len(other_vals) if other_vals else float("nan")

    # PRS contribution ratio
    helix_sum = sum(helix_vals)
    sheet_sum = sum(sheet_vals)
    other_sum = sum(other_vals)
    total_sum = helix_sum + sheet_sum + other_sum

    helix_ratio_prs = helix_sum / total_sum if total_sum > 0 else float("nan")
    sheet_ratio_prs = sheet_sum / total_sum if total_sum > 0 else float("nan")
    other_ratio_prs = other_sum / total_sum if total_sum > 0 else float("nan")

    row = {
        "Conformer": conformer,
        "Helix_mean_per_conformer": helix_mean_per_conformer,
        "Sheet_mean_per_conformer": sheet_mean_per_conformer,
        "Other_mean_per_conformer": other_mean_per_conformer,
        "Helix_ratio_prs": helix_ratio_prs,
        "Sheet_ratio_prs": sheet_ratio_prs,
        "Other_ratio_prs": other_ratio_prs,
    }

    return row, residue_rows


def save_outputs(df, residue_df, outdir, prefix):
    """
    Save TSVs and plots.
    prefix examples:
      - "all"
      - "r10_88"
    """
    ORDER = ["Other", "Sheet", "Helix"]

    # Conformer-level TSV
    tsv_path = os.path.join(outdir, f"secondary_prs_mean_and_ratio_{prefix}.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)

    """
    # Residue-level TSV (optional)
    residue_tsv = os.path.join(outdir, f"secondary_prs_residuelevel_{prefix}.tsv")
    residue_df.to_csv(residue_tsv, sep="\t", index=False)
    """

    # Mean PRS per residue (boxplot)
    mean_df = df.melt(
        id_vars="Conformer",
        value_vars=["Helix_mean_per_conformer", "Sheet_mean_per_conformer", "Other_mean_per_conformer"],
        var_name="Secondary",
        value_name="PRS"
    )
    mean_df["Secondary"] = mean_df["Secondary"].str.replace("_mean_per_conformer", "", regex=False)

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=mean_df, x="Secondary", y="PRS", order=ORDER)
    sns.stripplot(data=mean_df, x="Secondary", y="PRS", order=ORDER,
                  color="black", alpha=0.5)
    plt.ylim(bottom=0)
    plt.title(f"Mean PRS per Conformer [{prefix}]")
    plt.tight_layout()

    mean_png = os.path.join(outdir, f"secondary_prs_mean_per_conformer_boxplot_{prefix}.png")
    plt.savefig(mean_png, dpi=300)
    plt.close()

    """
    # PRS ratio plot (optional)
    ratio_df = df.melt(
        id_vars="Conformer",
        value_vars=["Helix_ratio_prs", "Sheet_ratio_prs", "Other_ratio_prs"],
        var_name="Secondary",
        value_name="PRS_ratio"
    )
    ratio_df["Secondary"] = ratio_df["Secondary"].str.replace("_ratio_prs", "", regex=False)

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=ratio_df, x="Secondary", y="PRS_ratio", order=ORDER)
    sns.stripplot(data=ratio_df, x="Secondary", y="PRS_ratio", order=ORDER,
                  color="black", alpha=0.5)
    plt.ylim(bottom=0)
    plt.title(f"PRS Contribution Ratio by Secondary Structure [{prefix}]")
    plt.tight_layout()

    ratio_png = os.path.join(outdir, f"secondary_prs_ratio_prs_boxplot_{prefix}.png")
    plt.savefig(ratio_png, dpi=300)
    plt.close()
    """

    # Residue-level pooled plot
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=residue_df, x="Secondary", y="PRS", order=ORDER)
    sns.stripplot(
        data=residue_df, x="Secondary", y="PRS", order=ORDER,
        color="black", alpha=0.15, jitter=0.25, size=2
    )
    plt.ylim(bottom=0)
    plt.title(f"PRS per Residue (All Conformers) [{prefix}]")
    plt.tight_layout()

    residue_png = os.path.join(outdir, f"secondary_prs_residuelevel_boxplot_{prefix}.png")
    plt.savefig(residue_png, dpi=300)
    plt.close()

    print(f"TSV saved : {tsv_path}")
    print(f"Figure saved: {mean_png}")


def main(pdb_dir, prs_tsv, output_root, range_ab=None):

    outdir = os.path.join(output_root, "result_7_secondaryPRS")
    os.makedirs(outdir, exist_ok=True)

    prs_df = pd.read_csv(prs_tsv, sep="\t", index_col=0)
    parser = PDBParser(QUIET=True)

    pdb_files = sorted(glob(os.path.join(pdb_dir, "*.pdb")))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in: {pdb_dir}")

    # ---------- all residues ----------
    rows_all = []
    residue_rows_all = []

    # ---------- optional range ----------
    rows_r = []
    residue_rows_r = []
    resnums_in_range_set = None
    prefix_r = None

    if range_ab is not None:
        a, b = range_ab

        if a <= 0 or b <= 0:
            raise ValueError("Range values must be positive integers.")
        if a > b:
            raise ValueError("Range must satisfy a <= b.")

        available = set(int(c) for c in prs_df.columns)
        resnums_in_range_set = set(r for r in range(a, b + 1) if r in available)
        prefix_r = f"r{a}_{b}"

    for pdb_file in pdb_files:
        conformer = os.path.basename(pdb_file).replace(".pdb", "")

        row_all, rr_all = aggregate_for_one_conformer(
            conformer, pdb_file, prs_df, parser, col_resnums_set=None
        )
        if row_all is not None:
            rows_all.append(row_all)
            residue_rows_all.extend(rr_all)

        if resnums_in_range_set is not None:
            row_r, rr_r = aggregate_for_one_conformer(
                conformer, pdb_file, prs_df, parser,
                col_resnums_set=resnums_in_range_set
            )
            if row_r is not None:
                rows_r.append(row_r)
                residue_rows_r.extend(rr_r)

    df_all = pd.DataFrame(rows_all)
    residue_df_all = pd.DataFrame(residue_rows_all)
    save_outputs(df_all, residue_df_all, outdir, prefix="all")

    if resnums_in_range_set is not None:
        df_r = pd.DataFrame(rows_r)
        residue_df_r = pd.DataFrame(residue_rows_r)
        save_outputs(df_r, residue_df_r, outdir, prefix=prefix_r)

    print("Done.")


def parse_range(s: str):
    """Parse '-r a,b' into (a, b)."""
    parts = s.split(",")

    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            "Range must be specified as 'a,b' (e.g., 10,88)."
        )

    try:
        a = int(parts[0])
        b = int(parts[1])
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Range values must be integers (e.g., 10,88)."
        )

    return a, b


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Secondary-structure stratified PRS analysis with optional residue range."
    )

    ap.add_argument("pdb_dir",
                    help="Directory containing conformer PDB files (*.pdb)")
    ap.add_argument("prs_tsv",
                    help="PRS TSV (rows=conformer, cols=residue numbers)")
    ap.add_argument("output_dir",
                    help="Output root directory")

    ap.add_argument(
        "-r", "--range",
        type=parse_range,
        default=None,
        help="Optional residue range 'a,b' (natural numbers)"
    )

    args = ap.parse_args()

    main(args.pdb_dir, args.prs_tsv, args.output_dir, range_ab=args.range)