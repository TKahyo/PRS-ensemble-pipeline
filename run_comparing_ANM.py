#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute RCk profile (C-terminal enrichment over last k residues)
based on ANM fluctuations.

For k = 1..K (default: 10), compute:
    RCk = mean(fluctuation of last k residues) / mean(fluctuation of all residues)

Output:
    - TSV table of RCk values per protein
    - Line plot of log2(RCk) vs k
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from prody import parsePDB, ANM, calcSqFlucts


# ==========================================
# Input parsing
# ==========================================
def read_list_file(tsv_path):
    items = []

    with open(tsv_path) as f:
        for line in f:
            s = line.strip()

            if not s or s.startswith("#"):
                continue

            parts = s.split()

            if len(parts) < 3:
                raise ValueError("Each line must be: <name> <type> <pdb_dir>")

            name, utype, pdb_dir = parts[0], parts[1].lower(), parts[2]

            if utype not in ("standalone", "domain", "tandem", "unknown"):
                utype = "unknown"

            items.append((name, utype, pdb_dir))

    return items


def list_pdbs(pdb_dir):
    exts = (".pdb", ".ent", ".pdb.gz", ".ent.gz")

    return sorted(
        os.path.join(pdb_dir, f)
        for f in os.listdir(pdb_dir)
        if f.lower().endswith(exts)
    )


# ==========================================
# CA extraction
# ==========================================
def select_ca(structure, pdb_file):
    ca = structure.select("name CA")

    if ca is None:
        raise ValueError(f"No CA atoms in {pdb_file}")

    return ca


def load_ca_coordsets(pdb_path):
    st = parsePDB(pdb_path)
    ca = select_ca(st, pdb_path)

    ag = ca.getAtomGroup()
    idx = ca.getIndices()

    coords = ag.getCoordsets()

    if coords is None:
        coords = ag.getCoords()[None, :, :]

    coords = np.asarray(coords)

    if coords.ndim == 2:
        coords = coords[None, :, :]

    return ca, coords[:, idx, :], np.asarray(ca.getResnums(), int)


# ==========================================
# ANM fluctuation calculation
# ==========================================
def compute_anm_mean(pdb_files, n_modes):

    flucs = []
    resnums = None

    # Single PDB with multiple models
    if len(pdb_files) == 1:
        ca_ref, coordsets, resnums = load_ca_coordsets(pdb_files[0])

        for cs in coordsets:
            ca = ca_ref.copy()
            ca.setCoords(cs)

            anm = ANM("anm")
            anm.buildHessian(ca)
            anm.calcModes(n_modes=n_modes)

            flucs.append(calcSqFlucts(anm[:n_modes]))

        return np.mean(flucs, axis=0), resnums

    # Multiple PDB files
    for pf in pdb_files:
        st = parsePDB(pf)
        ca = select_ca(st, pf)

        if resnums is None:
            resnums = ca.getResnums()

        anm = ANM("anm")
        anm.buildHessian(ca)
        anm.calcModes(n_modes=n_modes)

        flucs.append(calcSqFlucts(anm[:n_modes]))

    return np.mean(flucs, axis=0), np.asarray(resnums, int)


# ==========================================
# RCk calculation
# ==========================================
def rc_tail_ratio(values, k):
    k = max(1, min(k, len(values)))
    return np.mean(values[-k:]) / np.mean(values)


def _line_marker_face(ptype: str):
    p = ptype.lower()

    if p == "domain":
        return "0.7"
    elif p == "standalone":
        return "1.0"
    elif p == "tandem":
        return "0.0"
    else:
        return "1.0"


# ==========================================
# Main
# ==========================================
def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("list_tsv")
    ap.add_argument("out_dir")
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--modes", type=int, default=20)

    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    items = read_list_file(args.list_tsv)
    K = args.kmax

    rows = []

    for name, utype, pdb_dir in items:

        pdb_files = list_pdbs(pdb_dir)

        if len(pdb_files) == 0:
            continue

        anmfluc, _ = compute_anm_mean(pdb_files, args.modes)

        row = {
            "name": name,
            "type": utype,
            "pdb_dir": pdb_dir,
            "n_pdb_files": len(pdb_files),
            "n_residues": len(anmfluc),
        }

        for k in range(1, K + 1):
            row[f"RC{k}_ANM"] = rc_tail_ratio(anmfluc, k)

        rows.append(row)

    df = pd.DataFrame(rows)

    # Save TSV
    tsv_path = os.path.join(args.out_dir, "rck_anm_table.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ks = np.arange(1, K + 1)

    for _, r in df.iterrows():
        y = np.array([r[f"RC{k}_ANM"] for k in ks], float)
        y = np.log2(y)

        face = _line_marker_face(r["type"])

        ax.plot(
            ks, y,
            color="black",
            linewidth=1.2,
            marker="o",
            markersize=4.5,
            markerfacecolor=face,
            markeredgecolor="black",
        )

    ax.axhline(0.0, color="black", linewidth=1.0)

    ax.set_xlim(0, K + 0.5)
    ax.set_xticks(ks)
    ax.set_xlabel("Number of C-terminal residues (k)")
    ax.set_ylabel("log2(RCk) (ANM)")
    ax.set_title("C-terminal enrichment profile (log2 RCk, ANM)")

    fig.tight_layout()

    png_path = os.path.join(args.out_dir, "rck_anm_profile_log2.png")
    fig.savefig(png_path, dpi=300)
    plt.close(fig)

    print(f"[Saved] {tsv_path}")
    print(f"[Saved] {png_path}")


if __name__ == "__main__":
    main()