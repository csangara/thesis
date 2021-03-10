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
results_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
scrna_dir <- "/group/irc/shared/synthetic_visium/test/"
scrna_path <- paste0(scrna_dir, str_remove(dataset, "_generation"), "_test.rds")
dir.create(paste0("results/", dataset))

######### MuSiC ######### 

library(MuSiC)
library(xbioc)
dir.create(paste0(results_path, "music/"), recursive=TRUE)

# Load reference scRNA-seq data and convert to ExprSet
seurat_obj_scRNA = readRDS(scrna_path)
seurat_obj_scRNA@meta.data$celltype = seurat_obj_scRNA@meta.data$celltype
seurat_obj_scRNA = seurat_obj_scRNA %>% SetIdent(value = "celltype")
DefaultAssay(seurat_obj_scRNA) <- "RNA"
eset_obj_scRNA <- SeuratToExprSet(seurat_obj_scRNA)
rm(seurat_obj_scRNA)
res = list()

for (dataset_type in possible_dataset_types){
  # Load synthetic visium data and convert to Expreset
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts, PP=FALSE)
  eset_obj_visium <- SeuratToExprSet(seurat_obj_visium)
  rm(synthetic_visium_data, seurat_obj_visium)
  
  # Deconvolution
  music_deconv = music_prop(bulk.eset = eset_obj_visium, sc.eset = eset_obj_scRNA, clusters = 'celltype', samples='samples')
  # res[dataset_type] = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(music_deconv$Est.prop.weighted[,1:23]))
  saveRDS(music_deconv$Est.prop.weighted, paste0(results_path, "music/", dataset, "_", dataset_type, "_music.rds"))
  rm(music_deconv, eset_obj_visium)
  gc()
}