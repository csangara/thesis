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

dataset <- datasets[6]
repl <- "rep1"
run <- "run1"
result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")

#### SAMPLE CODE FOR LOADING DATA ####
dataset_type <- possible_dataset_types[1]
synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                        dataset_type, "_synthvisium.rds"))
seurat_obj_visium <- createAndPPSeuratFromCounts(synthetic_visium_data$counts, PP=FALSE)

method <- methods[1]
temp_deconv = readRDS(paste0(result_path, method, "/", dataset, "_", 
                                  dataset_type, "_", method, ".rds"))

# Correlation
ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
res = cor(t(synthetic_visium_data$relative_spot_composition[,1:ncells]), t(temp_deconv[[2]][,1:ncells]))
print(mean(diag(res), na.rm=TRUE))

###### PLOT PREDICTIONS ON UMAP #######

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

#### CALCULATE PERFORMANCE METRICS ####

all_results <- list()

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
  RMSE <- sapply(deconv_list, function(k) sqrt(sum((known_props-k[,1:ncells])**2, na.rm=TRUE)/ncells))
  
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

#### PLOT EACH METRIC ####
# corr, RMSE, accuracy, sensitivity, specificity, precision, F1
ntypes <- length(possible_dataset_types)
all_results <- readRDS(paste0(result_path, "all_metrics_", dataset, ".rds"))
metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1")
metric <- metrics[1]
  
for (metric in colnames(metrics)){
  df <- data.frame(dataset_type = rep(names(all_results), length(methods)),
                   index = rep(1:ntypes, length(methods)),
                   methods = rep(methods, each=ntypes))
  df$vals <- c(sapply(methods, function(u) sapply(possible_dataset_types, function(k) all_results[[k]][u,][metric])))
  png(paste0(result_path, "plots/", metric, ".png"), width=769, height=442)
  print(ggplot(data=df, aes(x=factor(index), y=as.numeric(vals), color=methods)) + geom_jitter(width=0.1) +
    labs(title=paste0(metric, " of different methods on all ", ntypes, " dataset types; ", dataset)) + ylab(metric) + xlab("datasets"))
  dev.off()
}

#### PLOT DISTRIBUTION OF CORRELATION ####

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

#### PRINT OUT RESULTS IN A TABLE ####
library(xlsx)

wb = createWorkbook()
for (ds in datasets[-1]){
  res_path <- paste0("results/", ds, "/", repl, "_", run, "/")
  all_results <- readRDS(paste0(result_path, "all_metrics_", ds, ".rds"))
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
saveWorkbook(wb, paste0("results/summary_files/", repl, "_", run, "/all_results.xlsx"))

#### PLOT RESULTS BY DATASET ####
library(stringr)

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

png(paste0("results/summary_files/", repl, "_", run, "/all_datasets_", metric, ".png"), width=800, height=500)
ggplot(data=df, aes(x=dataset, y=y, color=method)) +
  geom_point(size=2, position=position_dodge(.5)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2, position=position_dodge(.5)) +
  labs(title=paste0(metric, " of different methods on all datasets")) + ylab(metric) + xlab("datasets")
dev.off()

#### PLOT RESULTS BY DATASET TYPE ####
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

png(paste0("results/summary_files/", repl, "_", run, "/all_dataset_types_", metric, ".png"), width=800, height=500)
ggplot(data=df, aes(x=factor(dataset_type), y=y, color=method)) +
  geom_point(size=2, position=position_dodge(.5)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2, position=position_dodge(.5)) +
  labs(title=paste0(metric, " of different methods on all dataset types")) + ylab(metric) + xlab("dataset type")
dev.off()
