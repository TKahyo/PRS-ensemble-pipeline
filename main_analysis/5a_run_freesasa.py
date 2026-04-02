#!/usr/bin/env python3
import freesasa
from Bio.PDB import PDBParser
import numpy as np
import argparse
import sys
import os
import glob

# ==========================================
# Argument parser
# ==========================================
parser = argparse.ArgumentParser(
    description="Calculate ASA of pocket-contact residues using FreeSASA"
)

parser.add_argument("pdb_file", type=str)
parser.add_argument("pocket_path", type=str, nargs="?",
                    help="LVpocket mol2 file (LVpocket mode)")
parser.add_argument("--cutoff", type=float, default=5.0)
parser.add_argument("--precutoff", type=float, default=5.0)
parser.add_argument("--fpocket", action="store_true")
parser.add_argument("--fpocket_out", type=str,
                    help="fpocket <basename>_out directory")
parser.add_argument("--outdir", type=str,
                    help="output directory (unused in stdout mode)")

args = parser.parse_args()

pdb_file = args.pdb_file
distance_cutoff = args.cutoff
pre_distance_cutoff = args.precutoff

# ==========================================
# Load pocket coordinates
# ==========================================
def load_mol2_coordinates(mol2_file):
    coords = []

    with open(mol2_file) as f:
        read_atoms = False

        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                read_atoms = True
                continue

            if line.startswith("@<TRIPOS>") and read_atoms:
                break

            if read_atoms:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        coords.append(list(map(float, parts[2:5])))
                    except:
                        pass

    return np.array(coords)


def load_fpocket_vertices(fpocket_outdir):
    coords = []

    vert_files = glob.glob(
        os.path.join(fpocket_outdir, "pockets", "*_vert.pqr")
    )

    for vf in vert_files:
        with open(vf) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                    except:
                        pass

    return np.array(coords)

# ==========================================
# Select pocket source (fpocket or LVpocket)
# ==========================================
if args.fpocket:
    if args.fpocket_out is None:
        print("Error: --fpocket_out must be specified in fpocket mode",
              file=sys.stderr)
        sys.exit(1)

    pocket_points = load_fpocket_vertices(args.fpocket_out)

else:
    if args.pocket_path is None:
        print("Error: LVpocket mol2 file must be specified",
              file=sys.stderr)
        sys.exit(1)

    pocket_points = load_mol2_coordinates(args.pocket_path)

if len(pocket_points) == 0:
    print("Error: No pocket points found", file=sys.stderr)
    sys.exit(1)

# ==========================================
# Load PDB structure
# ==========================================
parser_pdb = PDBParser(QUIET=True)
structure = parser_pdb.get_structure("protein", pdb_file)

bio_atoms = []
bio_atom_to_residue = []

for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] not in (" ", "H"):
                continue

            for atom in residue:
                if not atom.get_name().startswith("H"):
                    bio_atoms.append(atom)
                    bio_atom_to_residue.append(
                        (chain.id, residue.id[1], residue.resname)
                    )

# ==========================================
# FreeSASA calculation (per residue)
# ==========================================
fs_structure = freesasa.Structure(pdb_file)
result = freesasa.calc(fs_structure)

residue_asa = {}
limit = min(fs_structure.nAtoms(), len(bio_atoms))

for i in range(limit):
    key = bio_atom_to_residue[i]
    residue_asa.setdefault(key, 0.0)
    residue_asa[key] += result.atomArea(i)

# ==========================================
# Identify pocket-contact residues
# ==========================================
candidate_residues = set()

for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue

            coords = np.array([atom.coord for atom in residue])
            if coords.size == 0:
                continue

            dist = np.linalg.norm(
                coords[:, None, :] - pocket_points[None, :, :],
                axis=2
            )

            if np.min(dist) < pre_distance_cutoff:
                key = (chain.id, residue.id[1], residue.resname)

                if residue_asa.get(key, 0.0) > 0.0:
                    candidate_residues.add(key)

# ==========================================
# Output TSV to stdout
# ==========================================
basename = os.path.splitext(os.path.basename(pdb_file))[0]

if args.fpocket_out:
    basename = os.path.basename(os.path.normpath(args.fpocket_out))

basename = basename.replace("_out", "")

for (chain, resi, resn) in sorted(candidate_residues):
    asa = residue_asa.get((chain, resi, resn), 0.0)
    print(f"{basename}\t{chain}\t{resi}\t{resn}\t{asa:.2f}")