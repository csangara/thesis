#### DETAILED PLOTS OF EACH REPLICATE ####
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")

###### PLOT PREDICTIONS ON UMAP #######
for (dataset in datasets[2:length(datasets)]){
  for (repl in paste0('rep', 4:10)){
    result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
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

#### PLOT DISTRIBUTION OF CORRELATION ####
for (dataset in datasets[2:length(datasets)]){
  for (repl in paste0('rep', 4:10)){
    result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
    for (dataset_type in possible_dataset_types){
      # Load reference data and deconvolution results
      synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                              dataset_type, "_synthvisium.rds"))
      ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
      # Initialization of column names
      celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:ncells]
      celltypes <- str_replace(celltypes, "/", ".")
      colnames(synthetic_visium_data$relative_spot_composition)[1:ncells] <- celltypes
      known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
      
      # Load deconvolution results
      deconv_list <- createDeconvResultList(methods, celltypes, result_path, dataset)
      
      # Correlation and RMSE
      corr_spots <- lapply(deconv_list, function(k) diag(cor(t(known_props), t(k[,1:ncells]))))
      plot_dir <- paste0(result_path, "plots/corr_distribution/")
      dir.create(plot_dir, showWarnings = FALSE)
      df = data.frame(x=unlist(corr_spots),
                      method=rep(names(corr_spots), each=length(corr_spots[[1]])))
      png(paste0(plot_dir, dataset_type, ".png"), width=769, height=442)
      print(ggplot(df, aes(x=x, color=method)) + geom_density() +
              labs(title=paste0("Correlation distribution of spots; ", dataset_type, ", ", dataset)))
      dev.off()
      
    }
  }
}


#### PLOT EACH METRIC AS FACET ####
metrics <- c("RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")

for (dataset in datasets[2:length(datasets)]){
  for (repl in paste0("rep", 1:10)){
    result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
    
    if (!dir.exists(paste0(result_path, "plots/"))){ dir.create(paste0(result_path, "plots/"), recursive=TRUE) }
    
    all_results <- readRDS(paste0(result_path, "all_metrics_", dataset, ".rds"))
    all_results2 <- melt(all_results) %>% mutate("dataset_type"=str_replace(L1, "artificial_", "") %>%
                                                   factor(., levels=unique(.))) %>%
                                          filter(variable %in% metrics)
    all_results2$methods <- methods
    
    p <- ggplot(all_results2, aes(x=methods, y=value, shape=dataset_type, color=methods)) +
      scale_shape_manual(values=1:length(possible_dataset_types)) + geom_jitter(width=0.1) +
      facet_wrap(vars(variable), scales="free_y")
    
    png(paste0(result_path, "plots/all_metrics.png"), width=800, height=450)
    print(p)
    dev.off()

  }
}