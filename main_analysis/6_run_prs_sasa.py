#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

# ===============================
# Heatmap (PRS per conformer with optional clustering)
# ===============================
def draw_heatmap_with_highlight(prs_df, pocket_df, outdir):
    os.makedirs(outdir, exist_ok=True)

    # Map conformer → pocket residues
    pocket_dict = {}
    for conf_name in pocket_df.iloc[:, 0].unique():
        conf_residues = pocket_df[
            pocket_df.iloc[:, 0] == conf_name
        ]['resid'].astype(int).tolist()
        pocket_dict[conf_name] = set(conf_residues)

    # Generate both clustered and non-clustered heatmaps
    for row_cluster in [True, False]:

        figsize = (
            max(6, prs_df.shape[1] // 3),
            max(6, prs_df.shape[0] // 2)
        )

        cg = sns.clustermap(
            prs_df,
            cmap='viridis',
            row_cluster=row_cluster,
            col_cluster=False,
            figsize=figsize,
            xticklabels=True,
            yticklabels=True
        )

        ax = cg.ax_heatmap

        # Adjust row order after clustering
        if row_cluster:
            reordered_rows = cg.dendrogram_row.reordered_ind
            prs_index = prs_df.index[reordered_rows]
        else:
            prs_index = prs_df.index

        """
        # Highlight pocket residues (optional)
        for i, conf in enumerate(prs_index):
            residues = pocket_dict.get(conf, set())
            for j, resid in enumerate(prs_df.columns.astype(int)):
                if resid in residues:
                    rect = Rectangle(
                        (j, i), 1, 1,
                        fill=False,
                        edgecolor='white',
                        lw=1.5
                    )
                    ax.add_patch(rect)
        """

        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

        fname = f"heatmap_rowcluster_{row_cluster}.png"
        plt.savefig(
            os.path.join(outdir, fname),
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()

        print(f"[Heatmap saved] {fname}")

# ===============================
# PRS vs ASA scatter plot
# ===============================
def draw_prssasa_correlation(prs_df, pocket_df, outdir):

    rows = []

    # Merge PRS and ASA for pocket residues
    for conf_name in prs_df.index:
        conf_prs = prs_df.loc[conf_name]
        conf_pocket = pocket_df[pocket_df.iloc[:, 0] == conf_name]

        for _, row in conf_pocket.iterrows():
            resid = int(row['resid'])
            asa = float(row['ASA'])

            if resid in conf_prs.index.astype(int):
                prs_val = conf_prs[str(resid)]
                rows.append({
                    'PRS': prs_val,
                    'ASA': asa,
                    'conformer': conf_name,
                    'resid': resid
                })

    if len(rows) == 0:
        print("No matching pocket residues found for PRS vs ASA plot.")
        return

    corr_df = pd.DataFrame(rows)

    # Save TSV
    corr_df.to_csv(
        os.path.join(outdir, "PRS_ASA_values.tsv"),
        sep='\t',
        index=False
    )
    print("[TSV saved] PRS_ASA_values.tsv")

    # Scatter plot
    plt.figure(figsize=(6, 6))
    plt.scatter(
        corr_df['ASA'],
        corr_df['PRS'],
        s=20,
        color='lightgrey',
        edgecolor='black',
        alpha=0.5
    )

    plt.xlabel("ASA")
    plt.ylabel("PRS")
    plt.title("PRS vs ASA for pocket residues")
    plt.tight_layout()

    plt.savefig(
        os.path.join(outdir, "PRS_vs_ASA.png"),
        dpi=300
    )
    plt.close()

    print("[Correlation plot saved] PRS_vs_ASA.png")

# ===============================
# PRS comparison (pocket vs non-pocket)
# ===============================
def plot_prs_pocket_groups_boxplot(prs_df, pocket_df, outdir):

    pocket_dict = {}
    for conf_name in pocket_df.iloc[:, 0].unique():
        conf_residues = pocket_df[
            pocket_df.iloc[:, 0] == conf_name
        ]['resid'].astype(int).tolist()
        pocket_dict[conf_name] = set(conf_residues)

    rows = []

    for conf_name in prs_df.index:
        conf_prs = prs_df.loc[conf_name]
        pocket_res = pocket_dict.get(conf_name, set())

        for resid in prs_df.columns.astype(int):
            prs_val = conf_prs[str(resid)]
            group = 'pocket' if resid in pocket_res else 'non-pocket'

            rows.append({
                'PRS': prs_val,
                'group': group
            })

    plot_df = pd.DataFrame(rows)

    # Save TSV
    plot_df.to_csv(
        os.path.join(outdir, "PRS_pocket_groups.tsv"),
        sep='\t',
        index=False
    )
    print("[TSV saved] PRS_pocket_groups.tsv")

    # Boxplot + stripplot
    plt.figure(figsize=(6, 6))

    sns.boxplot(
        x='group',
        y='PRS',
        data=plot_df,
        showfliers=False,
        boxprops=dict(
            facecolor='lightgrey',
            edgecolor='black',
            linewidth=1
        )
    )

    sns.stripplot(
        x='group',
        y='PRS',
        data=plot_df,
        size=5,
        color='black',
        edgecolor='white',
        alpha=0.7
    )

    plt.ylabel("PRS")
    plt.xlabel("")
    plt.title("PRS values for pocket vs non-pocket residues")
    plt.tight_layout()

    plt.savefig(
        os.path.join(outdir, "PRS_pocket_groups_boxplot.png"),
        dpi=300
    )
    plt.close()

    print("[Plot saved] PRS_pocket_groups_boxplot.png")

# ===============================
# Main
# ===============================
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("prs_tsv",
                        help="PRS values per conformer TSV")
    parser.add_argument("pocket_tsv",
                        help="Pocket residues TSV with ASA")
    parser.add_argument("--outdir",
                        default="out",
                        help="Output directory")

    args = parser.parse_args()

    pocket_base = os.path.splitext(
        os.path.basename(args.pocket_tsv)
    )[0]
    pocket_base = pocket_base.replace("freesasa", "prssa")

    outdir = os.path.join(args.outdir, "result_6_prsasa")
    os.makedirs(outdir, exist_ok=True)

    prs_df = pd.read_csv(args.prs_tsv, sep='\t', index_col=0)
    pocket_df = pd.read_csv(args.pocket_tsv, sep='\t')

    # Skip if only one conformer
    if prs_df.shape[0] < 2:
        print("[Info] Only one conformer detected. Skipping PRSSA heatmap.")
        return

    draw_heatmap_with_highlight(prs_df, pocket_df, outdir)
    draw_prssasa_correlation(prs_df, pocket_df, outdir)
    plot_prs_pocket_groups_boxplot(prs_df, pocket_df, outdir)

if __name__ == "__main__":
    main()