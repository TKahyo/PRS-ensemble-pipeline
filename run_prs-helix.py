#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarize secondary-structure PRS results and generate:

- Summary TSV per protein
- Scatter plots (Helix vs Sheet)
- Scatter point tables
- Sorted bar chart (Helix/Sheet ratio)

Input list TSV (no header):
    <label> <Standalone|Domain|Tandem> <path_to_secondary_summary_tsv>

Each secondary summary TSV must contain:
  - Helix_mean_per_conformer, Sheet_mean_per_conformer
  - Helix_ratio_prs, Sheet_ratio_prs

Usage:
  python run_prs-helix.py protein_list.tsv out_dir
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REQUIRED_COLS = [
    "Helix_mean_per_conformer", "Sheet_mean_per_conformer",
    "Helix_ratio_prs", "Sheet_ratio_prs",
]


# ==========================================
# Input parsing
# ==========================================
def read_list_tsv(list_path: str):
    """Read list TSV lines: <label> <type> <tsv_path>."""
    items = []

    with open(list_path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, start=1):
            s = line.strip()

            if not s or s.startswith("#"):
                continue

            parts = s.split()

            if len(parts) < 3:
                raise ValueError(
                    f"Line {ln}: expected '<label> <type> <path>', got: {s}"
                )

            label = parts[0]
            raw_type = parts[1].strip()
            ptype = raw_type.lower()
            tsv_path = parts[2]

            if ptype not in ("standalone", "domain", "tandem"):
                raise ValueError(
                    f"Line {ln}: type must be Standalone / Domain / Tandem"
                )

            items.append((label, ptype.capitalize(), tsv_path))

    return items


def load_secondary_tsv(tsv_path: str) -> pd.DataFrame:
    """Load secondary PRS summary TSV."""
    df = pd.read_csv(tsv_path, sep="\t")

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {tsv_path}: {missing}")

    return df


# ==========================================
# Protein-level summarization
# ==========================================
def summarize_one(label: str, ptype: str, tsv_path: str) -> dict:
    """
    Aggregate conformer-level TSV into a single protein-level summary.
    """
    df = load_secondary_tsv(tsv_path)

    return {
        "label": label,
        "type": ptype,
        "tsv_path": tsv_path,
        "n_conformers": int(len(df)),
        "Helix_mean_per_conformer": float(df["Helix_mean_per_conformer"].mean()),
        "Sheet_mean_per_conformer": float(df["Sheet_mean_per_conformer"].mean()),
        "Helix_mean_fraction": float(df["Helix_ratio_prs"].mean()),
        "Sheet_mean_fraction": float(df["Sheet_ratio_prs"].mean()),
    }


# ==========================================
# Plot styles
# ==========================================
def _point_style(ptype: str):
    """
    Point style by protein type (grayscale):
      - Domain     : gray
      - Standalone : white
      - Tandem     : black
    """
    p = ptype.lower()

    if p == "domain":
        face = "0.7"
    elif p == "standalone":
        face = "1.0"
    elif p == "tandem":
        face = "0.0"
    else:
        face = "1.0"

    return dict(
        marker="o",
        s=60,
        facecolors=face,
        edgecolors="black",
        linewidths=1.0,
    )


# ==========================================
# Scatter plot
# ==========================================
def plot_scatter(df, xcol, ycol, out_png, xlabel, ylabel, title):

    fig = plt.figure(figsize=(6.0, 6.0))
    ax = fig.add_subplot(111)

    xs = pd.to_numeric(df[xcol], errors="coerce").to_numpy()
    ys = pd.to_numeric(df[ycol], errors="coerce").to_numpy()

    for _, r in df.iterrows():
        ax.scatter(
            float(r[xcol]),
            float(r[ycol]),
            **_point_style(str(r["type"]))
        )

    finite = np.isfinite(xs) & np.isfinite(ys)

    if finite.any():
        mn = float(np.min(np.concatenate([xs[finite], ys[finite]])))
        mx = float(np.max(np.concatenate([xs[finite], ys[finite]])))
        pad = 0.10 * (mx - mn if mx > mn else 1.0)

        ax.set_xlim(mn - pad, mx + pad)
        ax.set_ylim(mn - pad, mx + pad)

        lo = max(ax.get_xlim()[0], ax.get_ylim()[0])
        hi = min(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0, color="black")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.14, top=0.92)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def save_scatter_points_tsv(df, xcol, ycol, out_tsv):
    """Save scatter plot points as TSV."""
    cols = ["label", "type", "n_conformers", "tsv_path", xcol, ycol]

    out_df = df[cols].copy()
    out_df = out_df.rename(columns={
        xcol: f"x({xcol})",
        ycol: f"y({ycol})"
    })

    out_df.to_csv(out_tsv, sep="\t", index=False)


# ==========================================
# Sorted ratio bar chart
# ==========================================
def plot_sorted_ratio_barh(df, num_col, den_col, label_col, out_png, title, xlabel):

    d = df[[label_col, num_col, den_col]].copy()

    d[num_col] = pd.to_numeric(d[num_col], errors="coerce")
    d[den_col] = pd.to_numeric(d[den_col], errors="coerce")

    d["value"] = d[num_col] / d[den_col]
    d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=["value"])
    d = d.sort_values("value", ascending=False)

    if len(d) == 0:
        print(f"[Skip] No finite values for {num_col}/{den_col}.")
        return

    fig_h = max(4.0, 0.35 * len(d) + 1.0)

    fig = plt.figure(figsize=(9.0, fig_h))
    ax = fig.add_subplot(111)

    ax.barh(d[label_col].astype(str), d["value"].to_numpy())
    ax.invert_yaxis()

    ax.set_xlabel(xlabel)
    ax.set_title(title)

    # Optional annotation
    #for i, v in enumerate(d["value"].to_numpy()):
    #    ax.text(v, i, f" {v:.3f}", va="center", ha="left", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# ==========================================
# Main
# ==========================================
def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("list_tsv",
                    help="TSV: <label> <type> <secondary_summary_tsv>")
    ap.add_argument("out_dir",
                    help="Output directory")

    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    items = read_list_tsv(args.list_tsv)

    if len(items) == 0:
        print("No valid lines found in list file.")
        sys.exit(1)

    rows = []

    for label, ptype, tsv_path in items:
        tsv_path = os.path.expanduser(tsv_path)

        try:
            row = summarize_one(label, ptype, tsv_path)
            rows.append(row)
            print(f"[OK] {label} ({ptype}) n_conformers={row['n_conformers']}")
        except Exception as e:
            print(f"[Error] {label} ({ptype}) {tsv_path}: {e}")

    if len(rows) == 0:
        print("No datasets processed successfully.")
        sys.exit(1)

    df = pd.DataFrame(rows)

    # Summary table
    out_tsv = os.path.join(args.out_dir, "protein_secondary_summary.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[Saved] {out_tsv}")

    # Scatter point tables
    out_points_mean = os.path.join(args.out_dir, "helix_vs_sheet_mean_points.tsv")
    save_scatter_points_tsv(df, "Sheet_mean_per_conformer", "Helix_mean_per_conformer", out_points_mean)

    out_points_frac = os.path.join(args.out_dir, "helix_vs_sheet_mean_fraction_points.tsv")
    save_scatter_points_tsv(df, "Sheet_mean_fraction", "Helix_mean_fraction", out_points_frac)

    print(f"[Saved] {out_points_mean}")
    print(f"[Saved] {out_points_frac}")

    # Scatter plots
    plot_scatter(
        df,
        "Sheet_mean_per_conformer",
        "Helix_mean_per_conformer",
        os.path.join(args.out_dir, "helix_vs_sheet_mean.png"),
        "Sheet mean (per residue; across conformers)",
        "Helix mean (per residue; across conformers)",
        "Helix vs Sheet (mean per residue)"
    )

    plot_scatter(
        df,
        "Sheet_mean_fraction",
        "Helix_mean_fraction",
        os.path.join(args.out_dir, "helix_vs_sheet_mean_fraction.png"),
        "Sheet mean fraction (across conformers)",
        "Helix mean fraction (across conformers)",
        "Helix vs Sheet (mean fraction)"
    )

    print("[Saved] scatter plots")

    # Helix/Sheet ratio bar chart
    plot_sorted_ratio_barh(
        df,
        num_col="Helix_mean_fraction",
        den_col="Sheet_mean_fraction",
        label_col="label",
        out_png=os.path.join(args.out_dir, "helix_over_sheet_mean_fraction_sorted.png"),
        title="Helix / Sheet (mean fraction; sorted)",
        xlabel="Helix_mean_fraction / Sheet_mean_fraction"
    )

    print("[Saved] Helix/Sheet ratio plot")


if __name__ == "__main__":
    main()