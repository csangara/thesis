setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"

library(Seurat)
library(synthvisium)
library(dplyr)
source("Scripts/helperFunctions.R")
library(patchwork)
library(ggplot2)
library(stringr)

artificial_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct")
datasets <- c('allen_cortex_dwn', 'brain_cortex_generation')

#### SAMPLE CODE FOR LOADING DATA ####
dataset_type <- possible_dataset_types[1]
synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
seurat_obj_visium <- createAndPPSeuratFromCounts(synthetic_visium_data$counts, PP=FALSE)

deconv = readRDS(paste0("results/", dataset, "/", "cell2location/", dataset, "_", 
                                  dataset_type, "_spotlight.rds"))

# Correlation
res = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(spotlight_deconv[[2]][,1:23]))
print(mean(diag(res), na.rm=TRUE))


#### COMPARE
methods <- c("cell2location")
dataset <- datasets[1]
all_results <- list()
n_celltypes <- 18

for (dataset_type in artificial_dataset_types){
  
  # Load reference data and deconvolution results
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_",
                                          dataset_type, "_synthvisium.rds"))
  
  # Initialization of column names
  celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:n_celltypes]
  celltypes <- str_replace(celltypes, "/", ".")
  colnames(synthetic_visium_data$relative_spot_composition)[1:n_celltypes] <- celltypes
  known_props <- synthetic_visium_data$relative_spot_composition[,1:n_celltypes]
  
  # Load deconvolution results
  deconv_list <- createDeconvResultList(methods, celltypes, paste0("results/", dataset, "/"), dataset)
  
  # Correlation and RMSE
  corr_spots <- sapply(deconv_list, function(k) mean(diag(cor(t(known_props), t(k[,1:n_celltypes]))), na.rm=TRUE))
  RMSE <- sapply(deconv_list, function(k) sqrt(sum((known_props-k[,1:n_celltypes])**2, na.rm=TRUE)/n_celltypes))
  
  # Classification metrics
  conf_matrices <- lapply(deconv_list, function(k) getConfusionMatrix(known_props, k))
  accuracy <- sapply(conf_matrices, function(k) round((k$tp + k$tn) / (k$tp + k$tn + k$fp + k$fn), 2))
  sensitivity <- sapply(conf_matrices, function(k) round(k$tp / (k$tp + k$fn), 2))
  specificity <- sapply(conf_matrices, function(k) round(k$tn / (k$tn + k$fp), 2))
  precision <- sapply(conf_matrices, function(k) round(k$tp / (k$tp + k$fp), 2))
  F1 <- sapply(methods, function(k) round(2 * ((precision[k] * sensitivity[k]) /
                                                 (precision[k] + sensitivity[k])), 2))
  # Get them into dataframe
  metrics <- data.frame(row.names=methods, "corr"=corr_spots, "RMSE"=RMSE,
                        "accuracy"=accuracy, "sensitivity"=sensitivity,
                        "specificity"=specificity, "precision"=precision, "F1"=F1)
  
  all_results[[dataset_type]] <- metrics
}

test <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered.rds")
