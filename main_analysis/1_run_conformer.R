#!/usr/bin/env Rscript
# Usage:
# Rscript 1_run_conformer.R <pdb_dir> <out_dir> [core_pdb_dir]

library(bio3d)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(umap)
library(ggdendro)
library(reshape2)
library(gridExtra)

# ------------------------------
# Plot theme (clean background, no grid, keep border)
# ------------------------------
theme_clean <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  )

# ------------------------------
# Arguments
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript 1_run_conformer.R <pdb_dir> <out_dir> [core_pdb_dir]")
}

pdb_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

core_pdb_dir <- ifelse(length(args) >= 3, args[3], pdb_dir)
core_file <- file.path(core_pdb_dir, "superimp_core_labelled_on_ref.pdb")

# ------------------------------
# Load PDB files
# ------------------------------
pdb_files <- list.files(pdb_dir, pattern="\\.pdb$", full.names = TRUE)
cat("Loading PDB files:\n")
print(pdb_files)

# Backbone alignment
pdbs_backbone <- pdbaln(pdb_files)

# All-atom loading
pdbs_allatom <- lapply(pdb_files, read.pdb)

# ------------------------------
# Backbone atom selection (core region if available)
# ------------------------------
if(file.exists(core_file)){
  core_pdb <- read.pdb(core_file)
  core_resno <- unique(core_pdb$atom$resno)
  cat("Using core PDB atoms from:", core_file, "\n")
} else {
  core_resno <- NULL
  cat("No core PDB found. Using all backbone atoms.\n")
}

backbone_inds <- atom.select(
  pdbs_backbone,
  resno = core_resno,
  elety = c("N","CA","C","O")
)$xyz

# ------------------------------
# RMSD calculation
# ------------------------------
cat("Calculating RMSD matrices...\n")

pdb_names <- basename(pdb_files)

# Backbone RMSD
rmsd_backbone <- rmsd(pdbs_backbone$xyz, a.inds=backbone_inds, fit=TRUE)
rownames(rmsd_backbone) <- colnames(rmsd_backbone) <- pdb_names

# All-atom RMSD
xyz_allatom <- do.call(rbind, lapply(pdbs_allatom, function(p) p$xyz))
allatom_inds <- 1:ncol(xyz_allatom)
rmsd_allatom <- rmsd(xyz_allatom, a.inds=allatom_inds, fit=TRUE)
rownames(rmsd_allatom) <- colnames(rmsd_allatom) <- pdb_names

# Save RMSD matrices
write.csv(rmsd_backbone, file=file.path(out_dir,"1_rmsd_backbone.csv"))
write.csv(rmsd_allatom, file=file.path(out_dir,"1_rmsd_allatom.csv"))

# ------------------------------
# RMSD triangular heatmap
# ------------------------------
rmsd_min <- min(rmsd_backbone, rmsd_allatom)
rmsd_max <- max(rmsd_backbone, rmsd_allatom)

plot_triangular_rmsd <- function(mat, hc, title, outfile){

  ord <- hc$order
  mat_ord <- mat[ord, ord]
  labels <- rownames(mat_ord)

  df <- melt(mat_ord, varnames = c("Row","Col"), value.name = "RMSD")
  df$Row <- factor(df$Row, levels = labels)
  df$Col <- factor(df$Col, levels = labels)

  # Lower triangle only
  df <- df[as.numeric(df$Row) > as.numeric(df$Col), ]

  p <- ggplot(df, aes(x = Col, y = Row, fill = RMSD)) +
    geom_tile() +
    coord_fixed() +
    scale_y_discrete(limits = rev(labels)) +
    scale_fill_gradient(
      low = "blue",
      high = "red",
      limits = c(rmsd_min, rmsd_max),
      na.value = "white"
    ) +
    labs(title = title, x = NULL, y = NULL, fill = "RMSD") +
    theme_clean +
    theme(
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  ggsave(outfile, p, width = 8, height = 8, dpi = 300, bg = "white")
}

hc_rmsd_backbone <- hclust(as.dist(rmsd_backbone))
hc_rmsd_allatom  <- hclust(as.dist(rmsd_allatom))

plot_triangular_rmsd(
  rmsd_backbone,
  hc_rmsd_backbone,
  "Backbone RMSD (Lower Triangle)",
  file.path(out_dir, "1_backbone_rmsd_triangle.png")
)

plot_triangular_rmsd(
  rmsd_allatom,
  hc_rmsd_allatom,
  "All-atom RMSD (Lower Triangle)",
  file.path(out_dir, "1_allatom_rmsd_triangle.png")
)

# ------------------------------
# PCA
# ------------------------------
cat("Running PCA...\n")

pc_backbone <- pca(pdbs_backbone$xyz, a.inds=backbone_inds, fit=TRUE)
scores_backbone <- as.data.frame(pc_backbone$z)
colnames(scores_backbone) <- paste0("PC",1:ncol(scores_backbone))
write.csv(scores_backbone, file=file.path(out_dir,"2_pca_scores_backbone.csv"), row.names=FALSE)

pc_allatom <- pca(xyz_allatom, a.inds=allatom_inds, fit=TRUE, use.svd=TRUE)
scores_allatom <- as.data.frame(pc_allatom$z)
colnames(scores_allatom) <- paste0("PC",1:ncol(scores_allatom))
write.csv(scores_allatom, file=file.path(out_dir,"2_pca_scores_allatom.csv"), row.names=FALSE)

# ------------------------------
# PCA plots (top 3 PCs)
# ------------------------------
pc_pairs <- list(c(1,2), c(1,3), c(2,3))

for(p in pc_pairs){

  ggplot(scores_backbone,
         aes(x=.data[[paste0("PC",p[1])]], y=.data[[paste0("PC",p[2])]])) +
    geom_point(color="red") +
    ggtitle(paste0("Backbone PCA: PC",p[1]," vs PC",p[2])) +
    theme_clean ->
    pplot

  ggsave(file.path(out_dir,
                   paste0("2_pca_backbone_PC",p[1],"_PC",p[2],".png")),
         pplot, width=8, height=6, bg="white")

  ggplot(scores_allatom,
         aes(x=.data[[paste0("PC",p[1])]], y=.data[[paste0("PC",p[2])]])) +
    geom_point(color="blue") +
    ggtitle(paste0("All-atom PCA: PC",p[1]," vs PC",p[2])) +
    theme_clean ->
    pplot

  ggsave(file.path(out_dir,
                   paste0("2_pca_allatom_PC",p[1],"_PC",p[2],".png")),
         pplot, width=8, height=6, bg="white")
}

# ------------------------------
# Residue-wise PC1 contribution
# ------------------------------
pc1_backbone <- apply(pc_backbone$au[,1, drop=FALSE], 1, function(x) x^2)

pc1_allatom <- if(!is.null(pc_allatom$au) && ncol(pc_allatom$au) >= 1){
  apply(pc_allatom$au[,1, drop=FALSE], 1, function(x) x^2)
} else {
  rep(0, nrow(xyz_allatom))
}

pc1_df_backbone <- data.frame(Residue=1:length(pc1_backbone),
                              PC1_Contribution=pc1_backbone)

ggplot(pc1_df_backbone, aes(x=Residue, y=PC1_Contribution)) +
  geom_bar(stat="identity", fill="steelblue") +
  ggtitle("Backbone PC1 Residue Contribution") +
  xlab("Residue") + ylab("Contribution") +
  theme_clean ->
  p_pc1b

ggsave(file.path(out_dir,"2_pc1_residue_contribution_backbone.png"),
       p_pc1b, width=10, height=4, bg="white")

# ------------------------------
# RMSF
# ------------------------------
ca_inds <- atom.select(pdbs_backbone, elety="CA")$xyz
rmsf_backbone <- rmsf(pdbs_backbone$xyz[, ca_inds])

pdb_first <- pdbs_allatom[[1]]
natoms <- length(pdb_first$atom$eleno)
rmsf_atom <- rmsf(xyz_allatom[, 1:(natoms*3)])
resnos <- pdb_first$atom$resno
rmsf_allatom <- tapply(rmsf_atom, resnos, mean)

rmsf_df_backbone <- data.frame(Residue=1:length(rmsf_backbone),
                               RMSF=rmsf_backbone)

rmsf_df_allatom <- data.frame(Residue=as.numeric(names(rmsf_allatom)),
                              RMSF=as.numeric(rmsf_allatom))

ggplot(rmsf_df_backbone, aes(x=Residue, y=RMSF)) +
  geom_bar(stat="identity", fill="coral") +
  ggtitle("Backbone RMSF") +
  xlab("Residue") + ylab("RMSF (Å)") +
  theme_clean ->
  p_rmsfb

ggsave(file.path(out_dir,"4_rmsf_backbone.png"),
       p_rmsfb, width=10, height=4, bg="white")

ggplot(rmsf_df_allatom, aes(x=Residue, y=RMSF)) +
  geom_bar(stat="identity", fill="orange") +
  ggtitle("All-atom RMSF (Residue-averaged)") +
  xlab("Residue") + ylab("RMSF (Å)") +
  theme_clean ->
  p_rmsfa

ggsave(file.path(out_dir,"4_rmsf_allatom.png"),
       p_rmsfa, width=10, height=4, bg="white")

# ------------------------------
# UMAP
# ------------------------------
umap_config <- umap::umap.defaults
umap_config$n_components <- 3
umap_config$random_state <- 42

n_samples <- nrow(scores_backbone)
umap_config$n_neighbors <- min(15, n_samples - 1)

cat("UMAP n_neighbors set to:", umap_config$n_neighbors, "\n")

umap_backbone <- umap(scores_backbone, config = umap_config)
umap_allatom  <- umap(scores_allatom,  config = umap_config)

df_umap_backbone <- as.data.frame(umap_backbone$layout)
colnames(df_umap_backbone) <- paste0("UMAP", seq_len(ncol(df_umap_backbone)))

df_umap_allatom <- as.data.frame(umap_allatom$layout)
colnames(df_umap_allatom) <- paste0("UMAP", seq_len(ncol(df_umap_allatom)))

write.csv(df_umap_backbone, file=file.path(out_dir,"3_umap_backbone.csv"), row.names=FALSE)
write.csv(df_umap_allatom, file=file.path(out_dir,"3_umap_allatom.csv"), row.names=FALSE)

# ------------------------------
# Clustering
# ------------------------------
cluster_list <- list(
  backbone_RMSD = cutree(hclust(as.dist(rmsd_backbone)), k=3),
  allatom_RMSD = cutree(hclust(as.dist(rmsd_allatom)), k=3),
  backbone_PCA = cutree(hclust(dist(scores_backbone)), k=3),
  allatom_PCA = cutree(hclust(dist(scores_allatom)), k=3),
  backbone_UMAP = cutree(hclust(dist(df_umap_backbone)), k=3),
  allatom_UMAP = cutree(hclust(dist(df_umap_allatom)), k=3)
)

cluster_df <- as.data.frame(cluster_list)
rownames(cluster_df) <- pdb_names

write.csv(cluster_df,
          file=file.path(out_dir,"5_cluster_attributions.csv"),
          row.names=TRUE)

pheatmap(
  as.matrix(cluster_df),
  color=colorRampPalette(c("blue","red"))(100),
  main="Cluster Attributions Heatmap",
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  filename=file.path(out_dir,"5_cluster_attributions_heatmap.png"),
  width=8,
  height=6,
  bg="white"
)

cat("Conformer analysis completed. Outputs saved to:", out_dir, "\n")