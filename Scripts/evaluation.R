#### EVALUATION OF RESULTS ####

setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"

library(Seurat)
library(synthvisium)
library(dplyr)
source("Scripts/helperFunctions.R")
library(patchwork)
library(ggplot2)
library(stringr)

possible_dataset_types <- c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")

artificial_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('allen_cortex_dwn', 'brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
              'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation')

#### SAMPLE CODE FOR LOADING DATA ####
dataset <- datasets[1]
dataset_type <- possible_dataset_types[1]
synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
seurat_obj_visium <- createAndPPSeuratFromCounts(synthetic_visium_data$counts, PP=FALSE)

spotlight_deconv = readRDS(paste0("results/", dataset, "/", "spotlight/", dataset, "_", 
                                  dataset_type, "_spotlight.rds"))

# Correlation
res = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(spotlight_deconv[[2]][,1:23]))
print(mean(diag(res), na.rm=TRUE))

###### PLOT PREDICTIONS ON UMAP #######

methods <- c("spotlight", "music", "cell2location", "RCTD", "stereoscope")
for (dataset_type in possible_dataset_types[1:2]){
  
  # Load reference data
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_",
                                          dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createAndPPSeuratFromCounts(synthetic_visium_data$counts)
  
  # Initialization of column names
  celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:23]
  celltypes <- str_replace(celltypes, "/", ".")
  colnames(synthetic_visium_data$relative_spot_composition)[1:23] <- celltypes
  known_props <- synthetic_visium_data$relative_spot_composition[,1:23]
  
  # Load deconvolution results
  deconv_list <- createDeconvResultList(methods, celltypes, paste0("results/", dataset, "/"), dataset)
  
  # Correlation (by spots and by cell type)
  corr_list <- lapply(deconv_list, function(k) cor(known_props[,1:23], k[,1:23], use="complete.obs"))
  
  for (celltype in celltypes){
    # Add deconv result to visium metadata
    seurat_obj_visium@meta.data[celltype] = synthetic_visium_data$relative_spot_composition[,celltype]
    
    for (method in methods){
      seurat_obj_visium@meta.data[paste0(celltype, "_", method)] = deconv_list[[method]][,celltype]
    }
    
    plot_dir <- paste0("plots/", dataset, "/", dataset_type, "/")
    if (!dir.exists(plot_dir)){ dir.create(plot_dir) }
    
    plots <- FeaturePlot(seurat_obj_visium, c(celltype, paste0(celltype, "_", methods)),
                                              combine=FALSE)
    
    for (i in 1:length(methods)){
      plots[[i+1]] <- plots[[i+1]] + ggtitle(methods[i], paste0("Corr=", round(corr_list[[methods[i]]][celltype, celltype], 3)))
      
    }

    plots <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]]
    #png(paste0(plot_dir, celltype, ".png"), width=1600, height=800)
    print(plots)
    #dev.off()
  }
}

#### CALCULATE PERFORMANCE METRICS ####

methods <- c("spotlight", "music", "cell2location", "RCTD", "stereoscope")
all_results <- list()
n_celltypes <- 23

for (dataset_type in possible_dataset_types){
  
  # Load reference data and deconvolution results
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_",
                                          dataset_type, "_synthvisium.rds"))
  
  # Initialization of column names
  celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:23]
  celltypes <- str_replace(celltypes, "/", ".")
  colnames(synthetic_visium_data$relative_spot_composition)[1:23] <- celltypes
  known_props <- synthetic_visium_data$relative_spot_composition[,1:23]
  
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

#### PLOT EACH METRIC ####
# corr, RMSE, accuracy, sensitivity, specificity, precision, F1
metric <- "F1" 
df <- data.frame(dataset_type = rep(names(all_results), length(methods)),
                 index = rep(1:13, length(methods)),
                 methods = rep(methods, each=13))
df$vals <- c(sapply(methods, function(u) sapply(possible_dataset_types, function(k) all_results[[k]][u,][metric])))

ggplot(data=df, aes(x=factor(index), y=as.numeric(vals), color=methods)) + geom_jitter(width=0.1) +
 labs(title=paste0(metric, " of different methods on all 13 dataset types")) + ylab(metric) + xlab("datasets")
