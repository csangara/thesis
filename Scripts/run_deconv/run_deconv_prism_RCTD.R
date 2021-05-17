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
# results_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
# scrna_dir <- "/group/irc/shared/synthetic_visium/test/"
# scrna_path <- paste0(scrna_dir, str_remove(dataset, "_generation"), "_test.rds")

# SCENARIO3
results_path <- paste0("results/", dataset, "_s4/", repl, "_", run, "/")
scrna_dir <- "/group/irc/shared/synthetic_visium/test/"
ref_dataset <- ifelse(dataset == "cerebellum_cell_generation", "cerebellum_nucleus_test.rds", "cerebellum_cell_test.rds")
scrna_path <- paste0(scrna_dir, ref_dataset)
print(paste("Dataset is: ", dataset))
print(paste("Reference dataset is:", ref_dataset))

######### RCTD ##########

library(RCTD)
library(Matrix)
dir.create(paste0(results_path, "RCTD/"), recursive=TRUE)
seurat_obj_scRNA <- readRDS(scrna_path)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
seurat_obj_scRNA@meta.data$liger_ident_coarse <- factor(seurat_obj_scRNA@meta.data$celltype)
seurat_obj_scRNA@meta.data$nUMI <- colSums(seurat_obj_scRNA@assays$RNA@counts)

for (dataset_type in possible_dataset_types){

  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  spatialRNA_obj_visium <- RCTD:::SpatialRNA(counts=as(as(synthetic_visium_data$counts,"matrix"),"dgCMatrix"))

  RCTD_deconv <- create.RCTD(spatialRNA_obj_visium, seurat_obj_scRNA, max_cores = 4, CELL_MIN_INSTANCE=5)
  RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = FALSE)
  res = as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))

  saveRDS(res, paste0(results_path, "RCTD/", dataset, "_", dataset_type, "_RCTD.rds"))
}
