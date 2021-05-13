#!/usr/bin/env Rscript
inputs = commandArgs(trailingOnly=TRUE)
if (length(inputs) < 3) {
    stop("usage: < dataset index > < replication no. > < run no. >", call.=FALSE)
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

dataset <- datasets[as.integer(inputs[1])]
repl <- inputs[2]
run <- inputs[3]

# SCENARIO1
# results_path <- paste0("results/", dataset, "_s1/", repl, "_", run, "/")
# scrna_dir <- "/group/irc/shared/synthetic_visium/generation/"
# scrna_path <- paste0(scrna_dir, dataset, ".rds")

# SCENARIO2
results_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
scrna_dir <- "/group/irc/shared/synthetic_visium/test/"
scrna_path <- paste0(scrna_dir, str_remove(dataset, "_generation"), "_test.rds")

# SCENARIO3
# results_path <- paste0("results/", dataset, "_s4/", repl, "_", run, "/")
# scrna_dir <- "/group/irc/shared/synthetic_visium/test/"
# ref_dataset <- ifelse(dataset == "cerebellum_cell_generation", "cerebellum_nucleus_test.rds", "cerebellum_cell_test.rds")
# scrna_path <- paste0(scrna_dir, ref_dataset)
# print(paste("Dataset is: ", dataset))
# print(paste("Reference dataset is:", ref_dataset))

######### SPOTLIGHT ######### 

library(SPOTlight)
dir.create(paste0(results_path, "spotlight/"), recursive=TRUE)
# Extract the top marker genes from each cluster
seurat_obj_scRNA <- readRDS(scrna_path)
seurat_obj_scRNA <- preprocessSeurat(seurat_obj_scRNA)
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA@meta.data$celltype
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA, assay = "SCT", slot = "data", verbose = TRUE,
                                      only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)

for (dataset_type in possible_dataset_types){
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts)
  spotlight_deconv <- spotlight_deconvolution(se_sc = seurat_obj_scRNA, counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                              clust_vr = "celltype", cluster_markers = cluster_markers_all, cl_n = 50,
                                              hvg = 3000, ntop = NULL, transf = "uv", method = "nsNMF", min_cont = 0.09)

  saveRDS(spotlight_deconv, paste0(results_path, "spotlight/", dataset, "_", dataset_type, "_spotlight.rds"))
}
