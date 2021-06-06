setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
source("Scripts/helperFunctions.R")
library(Seurat)
library(synthvisium)
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)
library(reshape2)
library(tidyr)

dataset_types_list <- list(
  original = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
               "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
               "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
               "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium"),
  recommended = c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
)
possible_dataset_types <- dataset_types_list[["recommended"]]

datasets <- c('allen_cortex_dwn', 'brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
              'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation')

methods <- c("spotlight", "music", "cell2location", "RCTD", "stereoscope")

dataset <- datasets[2]
repl <- "rep1"
run <- ""
result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
dataset_type = possible_dataset_types[1]