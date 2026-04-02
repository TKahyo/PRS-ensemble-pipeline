#!/usr/bin/env python3
# Usage: python run_structure_png.py input.pdb output.png [chainID]

import sys
from Bio.PDB import PDBParser, DSSP
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow


# ==========================================
# DSSP parsing
# ==========================================
def run_dssp(pdb_file, chain_id=None):
    """
    Run DSSP and return a mapping:
        residue number → Helix / Sheet / Other
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]

    dssp = DSSP(model, pdb_file)

    ss_map = {}

    for (chain, res_id), dssp_data in dssp.items():

        if chain_id is not None and chain != chain_id:
            continue

        resnum = res_id[1]
        ss = dssp_data[2]

        if ss in ("H", "G", "I"):
            ss_simple = "Helix"
        elif ss in ("E", "B"):
            ss_simple = "Sheet"
        else:
            ss_simple = "Other"

        ss_map[resnum] = ss_simple

    return ss_map


# ==========================================
# Convert residue annotation to segments
# ==========================================
def make_segments(ss_map):
    """
    Convert residue-level annotation into contiguous segments.
    """
    segments = []

    residues = sorted(ss_map.keys())

    start = residues[0]
    prev = start
    prev_ss = ss_map[start]

    for r in residues[1:]:

        if ss_map[r] == prev_ss and r == prev + 1:
            prev = r
        else:
            segments.append((start, prev, prev_ss))
            start = r
            prev = r
            prev_ss = ss_map[r]

    segments.append((start, prev, prev_ss))

    return segments


# ==========================================
# Plot linear secondary structure diagram
# ==========================================
def plot_linear_secondary_structure(
    segments,
    out_png,
    figsize=(14, 2),
    helix_color="#d62728",
    sheet_color="#f1c40f",
    other_color="gray"
):
    fig, ax = plt.subplots(figsize=figsize)

    y = 0.5
    height = 0.35

    for start, end, ss in segments:
        width = end - start + 1

        if ss == "Helix":
            ax.add_patch(
                Rectangle(
                    (start, y),
                    width,
                    height,
                    facecolor=helix_color,
                    edgecolor="black",
                    linewidth=0.8
                )
            )

        elif ss == "Sheet":
            ax.add_patch(
                FancyArrow(
                    start,
                    y + height / 2,
                    width,
                    0,
                    width=height * 0.9,
                    length_includes_head=True,
                    facecolor=sheet_color,
                    edgecolor="black",
                    linewidth=0.8
                )
            )

        else:
            ax.plot(
                [start, end + 1],
                [y + height / 2, y + height / 2],
                color=other_color,
                linewidth=1.2
            )

    all_res = [r for seg in segments for r in (seg[0], seg[1])]

    ax.set_xlim(min(all_res) - 2, max(all_res) + 2)
    ax.set_ylim(0, 1.4)

    ax.set_yticks([])
    ax.set_xlabel("Residue number")

    for spine in ["left", "right", "top"]:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# ==========================================
# Main
# ==========================================
def main():

    if len(sys.argv) not in (3, 4):
        print("Usage:")
        print("  python run_structure_png.py input.pdb output.png [chainID]")
        sys.exit(1)

    pdb_file = sys.argv[1]
    out_png = sys.argv[2]
    chain_id = sys.argv[3] if len(sys.argv) == 4 else None

    ss_map = run_dssp(pdb_file, chain_id=chain_id)
    segments = make_segments(ss_map)

    plot_linear_secondary_structure(segments, out_png)

    print("Done.")
    print(f"Output figure: {out_png}")


if __name__ == "__main__":
    main()