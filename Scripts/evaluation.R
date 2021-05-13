#### EVALUATION OF RESULTS ####
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")

#### SAMPLE CODE FOR LOADING DATA ####
dataset_type <- possible_dataset_types[1]
synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                        dataset_type, "_synthvisium.rds"))
seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts, PP=FALSE)

method <- methods[1]
temp_deconv = readRDS(paste0(result_path, method, "/", dataset, "_", 
                                  dataset_type, "_", method, ".rds"))

# Correlation
ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
res = cor(t(synthetic_visium_data$relative_spot_composition[,1:ncells]), t(temp_deconv[[2]][,1:ncells]))
print(mean(diag(res), na.rm=TRUE))

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
#### CALCULATE PERFORMANCE METRICS ####
for (dataset in datasets[2:length(datasets)]){
  for (repl in paste0('rep', 1:10)){
    all_results <- list()
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
      corr_spots <- sapply(deconv_list, function(k) mean(diag(cor(t(known_props), t(k[,1:ncells]))), na.rm=TRUE))
      RMSE <- sapply(deconv_list, function(k) mean(sqrt(rowSums((known_props-k[,1:ncells])**2, na.rm=TRUE)/ncells)))
      reference_RMSE <- mean(sqrt(rowSums((known_props-(1/ncells))**2)/ncells))
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
    saveRDS(all_results, paste0(result_path, "all_metrics_", dataset, ".rds"))
  }
}


#### PLOT EACH METRIC ####
# corr, RMSE, accuracy, sensitivity, specificity, precision, F1
for (dataset in datasets[2:length(datasets)]){
  for (repl in paste0('rep', 1:10)){
    
    result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
    ntypes <- length(possible_dataset_types)
    all_results <- readRDS(paste0(result_path, "all_metrics_", dataset, ".rds"))
    metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1")
    if (!dir.exists(paste0(result_path, "plots/"))){ dir.create(paste0(result_path, "plots/"), recursive=TRUE) }
    for (metric in metrics){
      df <- data.frame(dataset_type = rep(names(all_results), length(methods)),
                       index = rep(1:ntypes, length(methods)),
                       methods = rep(methods, each=ntypes))
      df$vals <- c(sapply(methods, function(u) sapply(possible_dataset_types, function(k) all_results[[k]][u,][metric])))
      png(paste0(result_path, "plots/", metric, ".png"), width=769, height=442)
      print(ggplot(data=df, aes(x=factor(index), y=as.numeric(vals), color=methods)) + geom_jitter(width=0.1) +
        labs(title=paste0(metric, " of different methods on all ", ntypes, " dataset types; ", dataset)) + ylab(metric) + xlab("datasets"))
      dev.off()
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

#### DEVIATION BETWEEN REPLICATES ####
df <- data.frame()
for (dataset in datasets[-1]){
  all_reps <- lapply(paste0('rep', 1:10), function(repl){
      all_results <- readRDS(paste0(paste0("results/", dataset, "/", repl, "_/"),
                                    "all_metrics_", dataset, ".rds"))
      sapply(possible_dataset_types, function (u) all_results[[u]][["RMSE"]] )
  })
  
  mean_reps <- sapply(possible_dataset_types, function(u) {
    apply(sapply(1:10, function(k) all_reps[[k]][,u]), 1, mean)})
  sd_reps <- sapply(possible_dataset_types, function(u) {
    apply(sapply(1:10, function(k) all_reps[[k]][,u]), 1, sd)})
  
  #temp_df <- data.frame(mean=c(mean_reps),
  #                      sd=c(sd_reps),
  #                      method=rep(methods, length(possible_dataset_types)),
  #                      dataset_type=rep(possible_dataset_types, each=length(methods)),
  #                      dataset=rep(dataset,length(possible_dataset_types)*length(methods)))
  
  temp_df <- data.frame(all_RMSE=unlist(all_reps),
                        reps=rep(1:10, each=length(methods)*length(possible_dataset_types)),
                        method=c("SPOT", "MuSiC", "c2l", "RCTD", "stereo"),
                        dataset_type=rep(possible_dataset_types, each=length(methods)),
                        dataset=dataset)
  
  df <- rbind(df, temp_df)
}

df$dataset_type_index <-as.numeric(factor(df$dataset_type, levels=unique(possible_dataset_types)))
n_cells_list = c(18, 8, 8, 6, 16, 9, 14)
df$method <- sapply(df$method, str_replace, "music", "MuSiC")
df$method <- sapply(df$method, str_replace, "spotlight", "SPOTlight")
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)")
proper_dataset_names <- paste0(proper_dataset_names, " [", n_cells_list, "]")
df$proper_dataset_names <- rep(proper_dataset_names, each = length(possible_dataset_types)*length(methods))
df$dataset_type <- sapply(df$dataset_type, str_replace, "artificial_", "")

# SD ?
ggplot(df, aes(x=method, y=sd, shape=factor(str_replace(dataset_type, "artificial_", "")), color=str_replace(dataset, "_generation", ""))) +
  geom_jitter(width=0.2, alpha=0.7) +
  labs(title="SD of each method between reps", color="dataset", shape="dataset_type") + 
  ylab("Deviation of RMSE between 10 reps") + scale_shape_manual(values=1:length(possible_dataset_types)) +
  guides(color = guide_legend(order=1), shape=guide_legend(order=2)) + ylim(0, 0.1)

# MEAN
# by dataset type
ggplot(df, aes(shape= proper_dataset_names,
               y=mean, x=dataset_type_index, color=method, group=method)) +
  geom_point(size=2, position=position_dodge(0.6)) +
  labs(title="Mean of each method between reps", color="method", shape="dataset") + 
  ylab("RMSE") + scale_shape_manual(values=1:length(datasets)) +
  guides(color = guide_legend(order=1), shape=guide_legend(order=2)) + 
  scale_x_discrete("dataset type", limits=factor(1:8))

# by dataset
png("plots/mean_between_reps_by_dataset.png", width=1200, height=500)
ggplot(df, aes(shape=factor(dataset_type, levels=unique(dataset_type)),
               y=mean, x=as.numeric(as.factor(dataset)), color=method, group=method)) +
  geom_point(size=2, position=position_dodge(0.6)) +
  labs(title="Mean RMSE of each method by dataset", color="method", shape="dataset type") + 
  ylab("RMSE") + scale_shape_manual(values=1:length(possible_dataset_types)) +
  guides(color = guide_legend(order=1), shape=guide_legend(order=2)) + ylim(0, 0.275) +
  xlab("dataset") + scale_x_discrete("Dataset", limits=factor(1:7))
dev.off()

# by method
png("plots/mean_between_10reps_groupbydataset.png", width=800, height=400)
ggplot(df, aes(x=method, y=mean, color=proper_dataset_names,
               shape=factor(dataset_type, levels=unique(dataset_type)),
                            group=proper_dataset_names)) +
  geom_point(size=2, position=position_dodge(0.6)) +
  labs(title="Mean RMSE between reps", color="Dataset", shape="Dataset type") +
  ylab("RMSE") + xlab ("Method") + scale_shape_manual(values=1:length(possible_dataset_types)) +
  ylim(0, 0.25) +
  guides(color = guide_legend(order=1), shape=guide_legend(order=2))
dev.off()

# FACET GRID
df$dataset_type <- factor(df$dataset_type, levels=unique(df$dataset_type))
png("plots/all_RMSE_facet_annotate.png", width=1500, height=1000)
p <- ggplot(df, aes(x=method, y=all_RMSE, color=method)) + geom_boxplot() +
  ylab("Mean RMSE over 10 iterations") +stat_summary(geom="text", fun=median,
                                                      aes(label=sprintf("%1.3f", ..y..), color=method),
                                                      position=position_nudge(x=0.33), size=3.5)
print(p + facet_grid(rows=vars(proper_dataset_names), cols=vars(dataset_type)))
dev.off()

### WILCOXON TEST ###
pvalues <- list()
mean_RMSEs <- list()
median_RMSEs <- list()
# 10 = n*(n-1)/2, where n is the number of methods
n_comparisons <- 10*length(possible_dataset_types)*length(datasets)
for (dataset in datasets[2:length(datasets)]){
  for (dataset_type in possible_dataset_types){
    stacked_RMSE <- data.frame()
    for (repl in paste0('rep', 1:10)){
      result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
      
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
      RMSE <- sapply(deconv_list, function(k) sqrt(rowSums((known_props-k[,1:ncells])**2, na.rm=TRUE)/ncells))
      stacked_RMSE <- rbind(stacked_RMSE, stack(data.frame(RMSE)))
    }
    pval <- pairwise.wilcox.test(stacked_RMSE[,1], stacked_RMSE[,2],
                                 p.adjust.method = "none",
                                 paired=TRUE) 
    pval <- apply(pval$p.value, 2, function(k)
      p.adjust(k, method="bonferroni", n=n_comparisons))
    # Adjust pvalue for all comparisons
    pvalues[[dataset]][[dataset_type]] <- pval
    mean_RMSEs[[dataset]][[dataset_type]] <- apply(RMSE, 2, mean)
    median_RMSEs[[dataset]][[dataset_type]] <- apply(RMSE, 2, median)
  }
}

library(xlsx)
wb = createWorkbook()
sheet = createSheet(wb, "sheet1")
j = 1
for (dataset in datasets[2:length(datasets)]){
  addDataFrame(data.frame(x=dataset), sheet=sheet, startRow=1,
               startColumn=j, row.names=FALSE, col.names=FALSE)
  i = 2
  for (dataset_type in possible_dataset_types){
    addDataFrame(data.frame(x=dataset_type), sheet=sheet, startRow=i,
                 startColumn=j, row.names=FALSE, col.names=FALSE)
    temp_data = median_RMSEs[[dataset]][[dataset_type]]
    temp_data = as.matrix(temp_data) # have this for the means
    addDataFrame(temp_data, sheet=sheet, startRow=i+1, startColumn=j, col.names=FALSE)
    i = i + nrow(temp_data) + 2
  }
  j = j+4
}

saveWorkbook(wb, paste0("results/summary files/all_pvalues_wilcoxpaired_pairwise.xlsx"))
saveWorkbook(wb, paste0("results/summary files/all_meanRMSEs.xlsx"))
saveWorkbook(wb, paste0("results/summary files/all_medianRMSEs.xlsx"))

all_sums <- 0
for (dataset in datasets[2:length(datasets)]){
  print(paste0("Dataset: ", dataset))
  for (dataset_type in possible_dataset_types){
    pvalues[[dataset]][[dataset_type]]
    all_sums <- all_sums + sum(pvalues[[dataset]][[dataset_type]] > 0.05, na.rm=TRUE)
    nonsig <- which(pvalues[[dataset]][[dataset_type]] > 0.05, arr.ind=TRUE)
    if (length(nonsig) > 0){
      print(paste0("type: ", dataset_type))
      print(paste0(methods[1:4][nonsig[,2]], "/", methods[2:5][nonsig[,1]], collapse="; "))
      
    }
  }
  print("")
}
