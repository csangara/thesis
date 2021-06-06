
###### PLOT PREDICTIONS ON UMAP #######
for (dataset in datasets[2:length(datasets)]){
  for (repl in paste0('rep', 1:10)){
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
  for (repl in paste0('rep', 1:10)){
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

#### PRINT OUT RESULTS IN A TABLE ####
library(xlsx)
for (repl in paste0("rep", 1:10)){
  wb = createWorkbook()
  for (ds in datasets[-1]){
    res_path <- paste0("results/", ds, "/", repl, "_", run, "/")
    
    all_results <- readRDS(paste0(res_path, "all_metrics_", ds, ".rds"))
    sheet = createSheet(wb, ds)
    i = 1
    for (dataset_type in names(all_results)){
      addDataFrame(data.frame(x=dataset_type), sheet=sheet, startRow=i,
                   row.names=FALSE, col.names=FALSE)
      temp_data = t(all_results[[dataset_type]])
      addDataFrame(temp_data, sheet=sheet, startRow=i+1)
      i = i + nrow(temp_data) + 2
    }
  }
  dir.create(paste0("results/summary files/", repl, "_", run, "/"), recursive=TRUE, showWarnings = FALSE)
  saveWorkbook(wb, paste0("results/summary files/", repl, "_", run, "/all_results.xlsx"))
}

#### PLOT RESULTS BY DATASET ####
library(stringr)

for (repl in paste0("rep", 1:10)){
  df <- data.frame()
  for (ds in datasets[-1]){
    res_path <- paste0("results/", ds, "/", repl, "_", run, "/")
    all_results <- readRDS(paste0(res_path, "all_metrics_", ds, ".rds"))
    metric <- "RMSE"
    methods <- rownames(all_results[[1]])
    temp_df <- data.frame(sapply(all_results, function (k) k[metric]))
    temp_df <- data.frame(y=apply(temp_df, 1, mean),
                          ymin=apply(temp_df, 1, min),
                          ymax=apply(temp_df, 1, max),
                          method=methods,
                          dataset=rep(str_replace(ds, "_generation", ""), 
                                      length(methods)))
    df <- rbind(df, temp_df)
  }
  df$dataset <- str_replace(df$dataset, "cerebellum_", "cer_")
  df$dataset <- str_replace(df$dataset, "hippocampus", "hipp")
  dir.create(paste0("results/summary files/", repl, "_", run, "/"), recursive=TRUE, showWarnings = FALSE)
  png(paste0("results/summary files/", repl, "_", run, "/all_datasets_", metric, ".png"), width=800, height=500)
  print(ggplot(data=df, aes(x=dataset, y=y, color=method)) +
          geom_point(size=2, position=position_dodge(.5)) +
          geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2, position=position_dodge(.5)) +
          labs(title=paste0(metric, " of different methods on all datasets")) + ylab(metric) + xlab("datasets"))
  dev.off()
}

#### PLOT RESULTS BY DATASET TYPE ####
for (repl in paste0("rep", 1:10)){
  df <- data.frame()
  metric <- "RMSE"
  for (dataset_type in possible_dataset_types){
    all_results <- lapply(datasets[-1], function(k) {
      readRDS(paste0("results/", k, "/", repl, "_", run, "/all_metrics_", k, ".rds"))})
    temp_df <- data.frame(sapply(all_results, function(k) k[[dataset_type]][[metric]]))
    colnames(temp_df) <- datasets[-1]
    
    methods <- rownames(all_results[[1]][[1]])
    
    temp_df <- data.frame(y=apply(temp_df, 1, mean),
                          ymin=apply(temp_df, 1, min),
                          ymax=apply(temp_df, 1, max),
                          method=methods)
    df <- rbind(df, temp_df)
  }
  df$dataset_type = rep(1:length(possible_dataset_types), each=length(methods))
  
  png(paste0("results/summary files/", repl, "_", run, "/all_dataset_types_", metric, ".png"), width=800, height=500)
  print(ggplot(data=df, aes(x=factor(dataset_type), y=y, color=method)) +
          geom_point(size=2, position=position_dodge(.5)) +
          geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2, position=position_dodge(.5)) +
          labs(title=paste0(metric, " of different methods on all dataset types")) + ylab(metric) + xlab("dataset type"))
  dev.off()
}