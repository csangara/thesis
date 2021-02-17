#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Please provide dataset to evaluate", call.=FALSE)
}

setwd("~/thesis/")
path <- "~/data/"
library(Seurat)
library(synthvisium)
library(dplyr)
library(stringr)
source("Scripts/helperFunctions.R")

possible_dataset_types = c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap",
                              "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                              "artificial_partially_dominant_celltype_diverse", "artificial_dominant_rare_celltype_diverse",
                              "artificial_regional_rare_celltype_diverse")

  datasets <- c('allen_cortex_dwn', 'brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
                'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation')

######## GENERATING NEEDED DATA ############
# Create synthetic data from scRNA data and save as RDS
# for (dataset_type in possible_dataset_types){
#   print(dataset_type)
#   createSynthvisiumRDS("rds/allen_cortex_dwn_original.rds", dataset_type,
#                        output_folder="Data/synthetic_datasets/allen_cortex_dwn/")
# }

dataset <- datasets[as.integer(inputs[1])]
scrna_dir <- "/group/irc/shared/synthetic_visium/test/"
scrna_path <- paste0(scrna_dir, str_remove(dataset, "_generation"), "_test.rds")
dir.create(paste0("results/", dataset))

######### SPOTLIGHT ######### 

library(SPOTlight)
dir.create(paste0("results/", dataset, "/spotlight/"))
# Extract the top marker genes from each cluster
seurat_obj_scRNA <- readRDS(scrna_path)
# seurat_obj_scRNA <- preprocessSeurat(seurat_obj_scRNA)
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA@meta.data$subclass
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA, assay = "SCT", slot = "data", verbose = TRUE,
                                      only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)

for (dataset_type in possible_dataset_types){
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts)
  spotlight_deconv <- spotlight_deconvolution(se_sc = seurat_obj_scRNA, counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                              clust_vr = "celltype", cluster_markers = cluster_markers_all, cl_n = 50,
                                              hvg = 3000, ntop = NULL, transf = "uv", method = "nsNMF", min_cont = 0.09)

  saveRDS(spotlight_deconv, paste0("results/", dataset, "/spotlight/", dataset, "_", dataset_type, "_spotlight.rds"))
}
