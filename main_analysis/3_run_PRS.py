#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from prody import *
from glob import glob
import pandas as pd
from multiprocessing import Pool, cpu_count

# ==========================================
# PRS computation
# ==========================================
def compute_one_prs(model):
    ca = model.select("calpha")

    anm = ANM("anm")
    anm.buildHessian(ca)
    anm.calcModes(n_modes=20)

    _, effectiveness, _ = calcPerturbResponse(anm)
    return effectiveness

def compute_prs(conformers, n_jobs=1):
    n_conf = len(conformers)

    ca0 = conformers[0].select("calpha")
    N = ca0.numAtoms()
    resnums = np.array(ca0.getResnums(), dtype=int)

    n_jobs = max(1, min(n_jobs, n_conf))

    if n_jobs == 1:
        eff_all = np.zeros((n_conf, N))
        for i, model in enumerate(conformers):
            eff_all[i, :] = compute_one_prs(model)
    else:
        with Pool(n_jobs) as p:
            results = p.map(compute_one_prs, conformers)
        eff_all = np.vstack(results)

    prs_mean = np.mean(eff_all, axis=0)
    prs_var  = np.var(eff_all, axis=0)

    return eff_all, prs_mean, prs_var, resnums

# ==========================================
# Save per-conformer PRS (TSV)
# ==========================================
def save_per_conformer_tsv(eff_all, resnums, pdb_files, outdir):
    conf_names = [os.path.splitext(os.path.basename(p))[0] for p in pdb_files]

    df = pd.DataFrame(eff_all, index=conf_names, columns=resnums)
    out_tsv = os.path.join(outdir, "PRS_per_conformer.tsv")

    df.to_csv(out_tsv, sep="\t")
    print(f"[per-conformer PRS TSV] → {out_tsv}")

# ==========================================
# Heatmap (with clustering)
# ==========================================
def plot_heatmap_clustermap(eff_all, resnums, conf_names, outdir):

    Nres = len(resnums)

    if Nres < 200:
        step, fontsize = 5, 7
    elif Nres < 400:
        step, fontsize = 10, 6
    else:
        step, fontsize = max(1, Nres // 50), 5

    tick_positions = np.arange(0, Nres, step)
    tick_labels = [str(resnums[i]) for i in tick_positions]

    cg = sns.clustermap(
        eff_all,
        cmap="viridis",
        metric="euclidean",
        method="ward",
        row_cluster=True,
        col_cluster=False,
        figsize=(24, 6),
        yticklabels=conf_names,
        xticklabels=False,
        cbar_kws={"shrink": 0.6},
    )

    ax = cg.ax_heatmap

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=90, fontsize=fontsize)
    ax.set_xlabel("Residue number (PDB numbering)")
    ax.set_ylabel("Conformer name")

    plt.savefig(
        os.path.join(outdir, "PRS_heatmap_clustermap.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

# ==========================================
# Heatmap (no clustering)
# ==========================================
def plot_heatmap_nocluster(eff_all, resnums, conf_names, outdir):

    fig = plt.figure(figsize=(24, 6))

    gs = fig.add_gridspec(1, 2, width_ratios=[0.03, 0.97], wspace=0.02)
    cax = fig.add_subplot(gs[0, 0])
    ax  = fig.add_subplot(gs[0, 1])

    sns.heatmap(
        eff_all,
        ax=ax,
        cmap="viridis",
        yticklabels=conf_names,
        xticklabels=resnums,
        cbar=True,
        cbar_ax=cax,
    )

    ax.set_xticklabels(resnums, rotation=90, fontsize=6)
    ax.set_xlabel("Residue number (PDB numbering)")

    ax.yaxis.tick_right()
    ax.set_ylabel("")

    ax.tick_params(axis="y", labelrotation=0)
    for t in ax.get_yticklabels():
        t.set_horizontalalignment("left")

    fig.subplots_adjust(left=0.02, right=0.95, bottom=0.18, top=0.98)

    fig.savefig(
        os.path.join(outdir, "PRS_heatmap_nocluster.png"),
        dpi=300
    )
    plt.close(fig)

# ==========================================
# PRS mean barplot
# ==========================================
def plot_prs(prs_mean, resnums, outdir):
    plt.figure(figsize=(16,4))
    plt.bar(resnums, prs_mean, color="steelblue")
    plt.xlabel("Residue number (PDB numbering)")
    plt.ylabel("PRS mean effectiveness")
    plt.title("PRS per-residue mean score")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "PRS_mean_barplot.png"), dpi=300)
    plt.close()

# ==========================================
# PyMOL coloring script
# ==========================================
def write_pymol_script(resnums, prs_mean, outdir, top_n=10):

    idx_sorted = np.argsort(prs_mean)[::-1]
    top_residues = resnums[idx_sorted[:top_n]]

    script_path = os.path.join(outdir, "PRS_color.pml")

    with open(script_path, "w") as f:
        f.write("reinitialize\n")
        f.write("load model.pdb\n")
        f.write("color gray80, all\n\n")

        for r in top_residues:
            f.write(f"select prs_{r}, resi {r}\n")
            f.write(f"color red, prs_{r}\n")

        f.write("show cartoon\n")
        f.write("set cartoon_transparency, 0.2\n")

    print(f"[PyMol script] → {script_path}")

# ==========================================
# Select base PDB for coloring
# ==========================================
def select_base_pdb(conformers, prs_mean):
    max_res = np.argmax(prs_mean)

    for model in conformers:
        ca = model.select("calpha")
        if ca is None:
            continue

        resnums = np.array(ca.getResnums(), dtype=int)

        if (max_res+1) in resnums:
            return model

    return conformers[0]

# ==========================================
# Write PRS-colored PDB (B-factor)
# ==========================================
def write_prs_pdb_allatoms(base_model, prs_mean, out_pdb):

    ca = base_model.select("calpha")
    if ca is None:
        raise RuntimeError("No C-alpha atoms in base_model")
    
    resnums = np.array(ca.getResnums(), dtype=int)
    resindices = ca.getResindices()

    if len(resnums) != len(prs_mean):
        raise ValueError("Residue count mismatch between PDB and PRS values")

    for i, residx in enumerate(resindices):
        sel_atoms = base_model.select(f"resindex {residx}")

        if sel_atoms is None or sel_atoms.numAtoms() == 0:
            continue

        sel_atoms.setBetas(float(prs_mean[i]))

    writePDB(out_pdb, base_model)
    print(f"[PRS-colored PDB] → {out_pdb}")

# ==========================================
# Main
# ==========================================
def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("conf_dir", help="directory of conformers (PDB)")
    parser.add_argument("out_dir", help="output directory")
    parser.add_argument("--top", type=int, default=10,
                        help="Top PRS residues for PyMol coloring")
    parser.add_argument("--n_jobs", type=int, default=1, help="number of parallel jobs")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    pdb_files = sorted(glob(os.path.join(args.conf_dir, "*.pdb")))
    if len(pdb_files) == 0:
        raise RuntimeError("No PDB files found.")

    print(f"Found {len(pdb_files)} conformers")

    conformers = [parsePDB(pdb) for pdb in pdb_files]
    conf_names = [os.path.splitext(os.path.basename(p))[0] for p in pdb_files]

    eff_all, prs_mean, prs_var, resnums = compute_prs(conformers, n_jobs=args.n_jobs)

    save_per_conformer_tsv(eff_all, resnums, pdb_files, args.out_dir)

    # Save mean and variance
    np.savetxt(
        os.path.join(args.out_dir, "PRS_mean.tsv"),
        np.vstack([resnums, prs_mean]).T,
        fmt="%d\t%.6f",
        header="resnum\tprs_mean"
    )

    np.savetxt(
        os.path.join(args.out_dir, "PRS_var.tsv"),
        np.vstack([resnums, prs_var]).T,
        fmt="%d\t%.6f",
        header="resnum\tprs_var"
    )

    # Visualization
    plot_prs(prs_mean, resnums, args.out_dir)

    if eff_all.shape[0] > 1:
        plot_heatmap_clustermap(eff_all, resnums, conf_names, args.out_dir)
    else:
        print("[Info] Only one conformer detected. Skipping clustermap.")

    plot_heatmap_nocluster(eff_all, resnums, conf_names, args.out_dir)

    # PyMOL script
    write_pymol_script(resnums, prs_mean, args.out_dir, top_n=args.top)

    # PRS-colored PDB
    base_model = select_base_pdb(conformers, prs_mean)
    out_pdb = os.path.join(args.out_dir, "PRS_avg_colored_allatoms.pdb")
    write_prs_pdb_allatoms(base_model, prs_mean, out_pdb)

    print("=== PRS analysis complete ===")

if __name__ == "__main__":
    main()