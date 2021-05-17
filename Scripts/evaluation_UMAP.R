inputs = commandArgs(trailingOnly=TRUE)
if (length(inputs) != 1) {
    stop("usage: < dataset index >", call.=FALSE)
}


setwd("~/thesis/")
path <- "~/data/"

library(Seurat)
library(synthvisium)
library(dplyr)
source("Scripts/helperFunctions.R")
library(patchwork)
library(ggplot2)
library(stringr)

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

dataset <- datasets[as.integer(inputs[1])]

###### PLOT PREDICTIONS ON UMAP #######
for (repl in c("rep1", "rep2", "rep3")[3]){
  for (run in c("run1", "run2", "run3")){
    result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
    if (repl == "rep1" && run == "run1") { next; }
    for (dataset_type in possible_dataset_types){
      
      # Load reference data
      synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                              dataset_type, "_synthvisium.rds"))
      seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts)
      
      # Initialization of column names
      ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
      celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:ncells]
      celltypes <- str_replace(celltypes, "/", ".")
      colnames(synthetic_visium_data$relative_spot_composition)[1:ncells] <- celltypes
      known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
      
      # Load deconvolution results
      deconv_list <- createDeconvResultList(methods, celltypes, result_path, dataset)
      
      # Correlation by cell type
      corr_list <- lapply(deconv_list, function(k) cor(known_props[,1:ncells],
                                                       k[,1:ncells], use="complete.obs"))
      
      for (celltype in celltypes){
        # Add deconv result to visium metadata
        seurat_obj_visium@meta.data[celltype] = synthetic_visium_data$relative_spot_composition[,celltype]
        
        for (method in methods){
          seurat_obj_visium@meta.data[paste0(celltype, "_", method)] = deconv_list[[method]][,celltype]
        }
        
        plot_dir <- paste0(result_path, "plots/", dataset_type, "/")
        if (!dir.exists(plot_dir)){ dir.create(plot_dir, recursive=TRUE) }
        
        # Plot proportion of each cell type for all methods
        plots <- FeaturePlot(seurat_obj_visium, c(celltype, paste0(celltype, "_", methods)),
                                                  combine=FALSE)
        # Add title for each method (first plot is the ground truth)
        for (i in 1:length(methods)){
          plots[[i+1]] <- plots[[i+1]] + ggtitle(methods[i], paste0("Corr=", round(corr_list[[methods[i]]][celltype, celltype], 3)))
          
        }
    
        plots <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]]
        png(paste0(plot_dir, celltype, ".png"), width=1600, height=800)
        print(plots)
        dev.off()
      }
    }
  }
}
