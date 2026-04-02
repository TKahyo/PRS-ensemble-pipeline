#!/usr/bin/env python3

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from prody import *
from multiprocessing import Pool

def process_one_pdb(pdb_file):
    structure = parsePDB(pdb_file)

    # ---------- CA atoms ----------
    ca = structure.select("name CA")
    if ca is None:
        return None

    resnums_ca = ca.getResnums()

    anm_ca = ANM("anm_ca")
    anm_ca.buildHessian(ca)
    anm_ca.calcModes(n_modes=20)
    fluc_ca = calcSqFlucts(anm_ca[:20])

    # ---------- All atoms ----------
    atoms_all = structure.select("all")
    if atoms_all is None:
        return None

    resnums_all = np.array([
        res.getResnum()
        for res in atoms_all.getHierView().iterResidues()
    ])

    anm_all = ANM("anm_all")
    anm_all.buildHessian(atoms_all)
    anm_all.calcModes(n_modes=20)

    fluc_all = calcSqFlucts(anm_all[:20])

    rmsf_per_res = []
    for res in atoms_all.getHierView().iterResidues():
        inds = res.getIndices()
        rmsf_per_res.append(np.mean(fluc_all[inds]))

    return fluc_ca, rmsf_per_res, resnums_ca, resnums_all

# =========================================================
# Arguments
# =========================================================
if len(sys.argv) < 3:
    print("Usage: python 2_run_ANM.py <pdb_dir> <out_dir>")
    sys.exit(1)

pdb_dir = sys.argv[1]
out_dir = sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

# =========================================================
# Collect PDB files
# =========================================================
pdb_files = [
    os.path.join(pdb_dir, f)
    for f in os.listdir(pdb_dir)
    if f.endswith(".pdb")
]

if not pdb_files:
    print("No PDB files found in", pdb_dir)
    sys.exit(1)

pdb_files.sort()
conf_names = [os.path.basename(f) for f in pdb_files]

print("Found PDB files:")
for f in pdb_files:
    print(f)


# ==========================================
# Parallel execution
# ==========================================
n_jobs = int(os.environ.get("N_JOBS", 1))
n_jobs = max(1, min(n_jobs, len(pdb_files)))

with Pool(n_jobs) as p:
    results = p.map(process_one_pdb, pdb_files)

valid_results = [r for r in results if r is not None]
valid_files = [f for f, r in zip(pdb_files, results) if r is not None]

if len(valid_results) == 0:
    raise RuntimeError("No valid ANM results computed.")

rmsf_matrix_ca = np.vstack([r[0] for r in valid_results])
rmsf_matrix_all = np.vstack([r[1] for r in valid_results])

resnums_ca = valid_results[0][2]
resnums_all = valid_results[0][3]

conf_names = [os.path.basename(f) for f in valid_files]


# =========================================================
# Save matrices (CSV / TSV)
# =========================================================
def save_matrix(matrix, resnums, names, out_csv, out_tsv):

    header = ["PDB"] + [f"Res{r}" for r in resnums]

    # CSV
    with open(out_csv, "w") as f:
        f.write(",".join(header) + "\n")
        for name, row in zip(names, matrix):
            f.write(name + "," + ",".join([f"{v:.6f}" for v in row]) + "\n")

    # TSV
    with open(out_tsv, "w") as f:
        f.write("\t".join(header) + "\n")
        for name, row in zip(names, matrix):
            f.write(name + "\t" + "\t".join([f"{v:.6f}" for v in row]) + "\n")


# CA RMSF
save_matrix(
    rmsf_matrix_ca,
    resnums_ca,
    conf_names,
    os.path.join(out_dir, "anm_rmsf_ca_matrix.csv"),
    os.path.join(out_dir, "anm_rmsf_ca_matrix.tsv")
)

# All-atom RMSF
save_matrix(
    rmsf_matrix_all,
    resnums_all,
    conf_names,
    os.path.join(out_dir, "anm_rmsf_allatom_matrix.csv"),
    os.path.join(out_dir, "anm_rmsf_allatom_matrix.tsv")
)

# =========================================================
# Heatmap (with clustering)
# =========================================================
def plot_heatmap_clustermap(matrix, resnums, conf_names, outpath):

    figsize = (
        max(10, matrix.shape[1] / 4),
        max(6,  matrix.shape[0] / 2.5)
    )

    cg = sns.clustermap(
        matrix,
        cmap="viridis",
        metric="euclidean",
        method="ward",
        figsize=figsize,
        row_cluster=True,
        col_cluster=False,
        yticklabels=conf_names,
        xticklabels=resnums,
        cbar_kws={"shrink": 0.6},
    )

    ax = cg.ax_heatmap

    ax.set_xticklabels(resnums, rotation=90, fontsize=6)
    ax.set_xlabel("Residue number")

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Conformer")

    reordered = cg.dendrogram_row.reordered_ind
    ax.set_yticklabels(
        [conf_names[i] for i in reordered],
        rotation=0,
        fontsize=8
    )

    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()

# =========================================================
# Heatmap (no clustering)
# =========================================================
def plot_heatmap_nocluster(matrix, resnums, conf_names, outpath):

    figsize = (
        max(10, matrix.shape[1] / 4),
        max(6,  matrix.shape[0] / 2.5)
    )

    plt.figure(figsize=figsize)

    ax = sns.heatmap(
        matrix,
        cmap="viridis",
        yticklabels=conf_names,
        xticklabels=resnums,
        cbar_kws={"shrink": 0.6},
    )

    ax.set_xticklabels(resnums, rotation=90, fontsize=6)
    ax.set_xlabel("Residue number")

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Conformer (original order)")

    ax.set_yticklabels(conf_names, rotation=0, fontsize=8)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()

# =========================================================
# Mean RMSF barplot
# =========================================================
def plot_mean_barplot(mean_values, resnums, ylabel, title, outpath):

    plt.figure(figsize=(16, 4))
    plt.bar(resnums, mean_values, color="steelblue")
    plt.xlabel("Residue number (PDB numbering)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

# =========================================================
# Generate outputs
# =========================================================
plot_heatmap_clustermap(
    rmsf_matrix_ca, resnums_ca, conf_names,
    os.path.join(out_dir, "anm_rmsf_ca_heatmap_clust.png")
)

plot_heatmap_clustermap(
    rmsf_matrix_all, resnums_all, conf_names,
    os.path.join(out_dir, "anm_rmsf_allatom_heatmap_clust.png")
)

plot_heatmap_nocluster(
    rmsf_matrix_ca, resnums_ca, conf_names,
    os.path.join(out_dir, "anm_rmsf_ca_heatmap_noclust.png")
)

plot_heatmap_nocluster(
    rmsf_matrix_all, resnums_all, conf_names,
    os.path.join(out_dir, "anm_rmsf_allatom_heatmap_noclust.png")
)

# =========================================================
# Mean RMSF plots
# =========================================================
if rmsf_matrix_ca.size > 0:
    mean_ca = np.mean(rmsf_matrix_ca, axis=0)
    plot_mean_barplot(
        mean_ca, resnums_ca,
        ylabel="ANM RMSF mean (CA)",
        title="ANM RMSF per-residue mean (CA)",
        outpath=os.path.join(out_dir, "anm_rmsf_ca_mean_barplot.png")
    )

if rmsf_matrix_all.size > 0:
    mean_all = np.mean(rmsf_matrix_all, axis=0)
    plot_mean_barplot(
        mean_all, resnums_all,
        ylabel="ANM RMSF mean (all atoms)",
        title="ANM RMSF per-residue mean (all atoms)",
        outpath=os.path.join(out_dir, "anm_rmsf_allatom_mean_barplot.png")
    )

print("=== ANM RMSF analysis complete ===")