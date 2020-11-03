#### EVALUATION OF RESULTS ####

setwd("D:/Work (Yr 2 Sem 1)/Thesis/Scripts")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/"

library(Seurat)
library(synthvisium)
library(dplyr)
source("helperFunctions.R")
library(patchwork)

possible_dataset_types = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")

synthetic_visium_data <- readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds"))
seurat_obj_visium <- createAndPPSeuratFromVisium(synthetic_visium_data$counts)
spotlight_deconv = readRDS(paste0("rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_spotlight.rds"))

res = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(spotlight_deconv[[2]][,1:23]))
print(mean(diag(res), na.rm=TRUE))


###### PLOT PREDICTIONS ON UMAP #######

# Do some initialization of column names

for (dataset_type in possible_dataset_types){
  
  # Load reference data and deconvolution results
  synthetic_visium_data <- readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createAndPPSeuratFromVisium(synthetic_visium_data$counts)
  known_props <- synthetic_visium_data$relative_spot_composition[,1:23]
  spotlight_deconv <- readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_spotlight.rds"))
  music_deconv <- readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_music.rds"))
  
  celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:23]
  colnames(spotlight_deconv[[2]])[1:23] = celltypes
  colnames(music_deconv$Est.prop.weighted)[1:23] = celltypes
  
  # Correlation (by spots and by cell type)
  spotlight_corr_spots <- cor(t(known_props), t(spotlight_deconv[[2]][,1:23]))
  spotlight_corr_celltypes <- cor(known_props[,1:23], spotlight_deconv[[2]][,1:23], use="complete.obs")
  music_corr_spots <- cor(t(known_props), t(music_deconv$Est.prop.weighted[,1:23]))
  music_corr_celltypes <- cor(known_props[,1:23], music_deconv$Est.prop.weighted[,1:23], use="complete.obs")
  
  for (celltype in celltypes){
    # Add deconv result to visium metadata
    seurat_obj_visium@meta.data[celltype] = synthetic_visium_data$relative_spot_composition[,celltype]
    seurat_obj_visium@meta.data[paste0(celltype, "_spotlight")] = spotlight_deconv[[2]][,celltype]
    seurat_obj_visium@meta.data[paste0(celltype, "_music")] = music_deconv$Est.prop.weighted[,celltype]
    
    plot_dir <- paste0(path, "plots/allen_cortex_dwn/", dataset_type, "/")
    if (!dir.exists(plot_dir)){ dir.create(plot_dir) }
    
    plots <- FeaturePlot(seurat_obj_visium, c(celltype, paste0(celltype, "_spotlight"), paste0(celltype, "_music")), combine=FALSE)
    plots[[2]] <- plots[[2]] + ggtitle("SPOTlight", paste0("Corr=",
                                                           round(spotlight_corr_celltypes[celltype, celltype], 3)))
    plots[[3]] <- plots[[3]] + ggtitle("MuSiC", paste0("Corr=",
                                                       round(music_corr_celltypes[celltype, celltype], 3)))
    plots <- plots[[1]] + plot_spacer() + plots[[2]]+plots[[3]]
    
    png(paste0(plot_dir, str_replace(celltype, "/", "."), ".png"), width=1000, height=500)
    print(plots)
    dev.off()
  }
}