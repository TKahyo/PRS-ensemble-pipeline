# PRS-ensemble-pipeline
A pipeline for protein structural ensemble analysis integrating:

ANM-based dynamics
Perturbation Response Scanning (PRS) effectiveness
Pocket detection (LVpocket)
Solvent accessibility (FreeSASA)

This pipeline enables identification of functionally important and potentially druggable regions based on structural dynamics and residue-level perturbation response.

# Overview
[bash_analysis.sh]
Input PDB (multi-model or multiple structures)
        ↓
Split PDB models
        ↓
Structural superimposition
        ↓
Aligned structural ensemble
        ↓
ANM / PRS effectiveness analysis
        ↓
Pocket detection + SASA calculation
        ↓
Integrated analysis (PRS × surface accessibility)

[run_comparing_ANM.py]
Cross-protein ANM comparison (RCk analysis)

# Features
Structural ensemble-based analysis
ANM (Anisotropic Network Model)
PRS effectiveness (residue-level dynamic influence)
Automated pocket detection (LVpocket)
Solvent accessible surface area (FreeSASA)
Integrated analysis of dynamics and surface exposure
Batch processing pipeline
Cross-protein comparison of ANM-derived flexibility (RCk)

# Requirements
Languages
Python ≥ 3.8
R
Python packages
prody
numpy
pandas
matplotlib
seaborn
R packages
bio3d
ggplot2
pheatmap
umap

# External Tools
The following tools must be installed separately:
>LVpocket

>FreeSASA

These tools are not included in this repository and are distributed under their own licenses.

# Installation
Clone the repository:
```bash
git clone https://github.com/your-username/PRS-ensemble-pipeline.git
cd PRS-ensemble-pipeline
```
Install Python dependencies:
```bash
pip install prody numpy pandas matplotlib seaborn
```
Install R packages:
```bash
install.packages(c("bio3d","ggplot2","pheatmap","umap"))
```

# Usage
## 1. Prepare input files
Place your PDB files in the working directory.
- Supported formats:
  - Multi-model PDB (e.g., NMR ensembles)
  - Multiple single PDB files

## 2. Prepare configuration file
Create a configuration file (e.g., `config.sh`) specifying required paths:

```bash
# R scripts
SPLIT_PDBS_SCRIPT=/path/to/split_pdbs_bio3d.R
SUPERIMPOSE_SCRIPT=/path/to/superimpose_bio3d.R

# fpocket binary
FPOCKET_BIN=/path/to/fpocket
```
⚠️ This file is required. The pipeline will fail without it.

## 3. Run Example
```bash
# 1. Clone repository
git clone https://github.com/your-username/PRS-ensemble-pipeline.git
cd PRS-ensemble-pipeline

# 2. Prepare working directory
mkdir test_run
cd test_run
cp /path/to/your/*.pdb .

# 3. Create config file
nano config.sh

# 4. Run pipeline
bash ../bash_all_analysis.sh config.sh

# Optional: with residue range
bash ../bash_all_analysis.sh config.sh 10,88
```

---
## 4. Prepare ANM-derived flexibility profiles across multiple proteins
In addition to per-structure analysis, the pipeline provides a tool for
comparing ANM-derived flexibility profiles across multiple proteins.
This analysis is performed using:
>run_comparing_ANM.py

Prepare a .tsv file describing multiple protein datasets:
><name> <type> <pdb_dir e.g. out/superimposed/>
- Example
```
UBL3 standalone /UBL3/out/superimposed/
LC3B domain /LC3B/out/superimposed/
ATG8 tandem /ATG8/out/superimposed/
```

## 5. 



## Output Structure
```
- bash_analysis.sh
out/
├── split_pdbs/                # split structures (multi-model input)
├── superimposed/              # aligned structures
├── result_1_conformer/        # RMSD / PCA / UMAP / RMSF
├── result_2_ANM/              # ANM fluctuation analysis
├── result_3_PRS/              # PRS analysis
├── result_4_fpockets/         # pocket detection
├── result_5_freesasa.tsv      # solvent accessibility
├── result_6_prsasa/           # PRS × ASA integration
├── result_7_secondaryPRS/     # secondary structure PRS
└── auxiliary files (aln.fa, superimp_core.pdb, ...)
```

- run_comparing_ANM.py
./rck_anm_profile_log2.png



# Methodological Notes
All structures are aligned prior to analysis
ANM captures intrinsic collective motions
PRS effectiveness quantifies residue-level dynamic influence
Pocket and surface accessibility are integrated to identify functional regions

This implementation follows the formulation of PRS effectiveness as described in the associated manuscript.

# License
This repository is released under the MIT License.

🔗 Third-party Software
This pipeline relies on external tools:
LVpocket
FreeSASA
ProDy
These tools are distributed under their own respective licenses.
Users are responsible for complying with those licenses.

This repository does not redistribute third-party software.

# Citation
If you use this pipeline, please cite:
"UBL3 UBL domain exhibits distinct helix-centered dynamic control among ubiquitin-like proteins" (under handling)
