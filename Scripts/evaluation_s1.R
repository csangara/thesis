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

dataset <- datasets[2]
repl <- "rep1"
run <- ""
result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
dataset_type = possible_dataset_types[1]

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
dataset <- datasets[2]

for (repl in paste0('rep', 1:10)){
  result_path <- paste0("results/", dataset, "_s1/", repl, "_", run, "/")
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

#### CALCULATE PERFORMANCE METRICS OF SCENARIO 1####
dataset <- datasets[2]
for (repl in paste0('rep', 1:10)){
    all_results <- list()
    result_path <- paste0("results/", dataset, "_s1/", repl, "_", run, "/")
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
    saveRDS(all_results, paste0(result_path, "all_metrics_", dataset, ".rds"))
  }
}


#### PLOT EACH METRIC ####
# corr, RMSE, accuracy, sensitivity, specificity, precision, F1
dataset <- datasets[2]
for (repl in paste0("rep", 1:10)){
  
  result_path <- paste0("results/", dataset, "_s1/", repl, "_", run, "/")
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

#### PLOT DISTRIBUTION OF CORRELATION ####
dataset <- datasets[2]
for (repl in paste0("rep", 1:10)){
  result_path <- paste0("results/", dataset, "_s1/", repl, "_", run, "/")
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

# Deviation between scenarios
df <- data.frame()
dataset <- datasets[2]
for (i in 1:2){
  ext <- c("_s1", "")[i]
  all_reps <- lapply(paste0('rep', 1:10), function(repl){
    
      all_results <- readRDS(paste0(paste0("results/", dataset, ext,"/", repl, "_/"),
                                    "all_metrics_", dataset, ".rds"))
      sapply(possible_dataset_types, function (u) all_results[[u]][["RMSE"]] )
    }
    
  )
  mean_reps <- sapply(possible_dataset_types, function(u) {
    apply(sapply(1:10, function(k) all_reps[[k]][,u]), 1, mean)})
  sd_reps <- sapply(possible_dataset_types, function(u) {
    apply(sapply(1:10, function(k) all_reps[[k]][,u]), 1, sd)})
  
  temp_df <- data.frame(mean=c(mean_reps),
                        sd=c(sd_reps),
                        method=rep(methods, length(possible_dataset_types)),
                        dataset_type=rep(possible_dataset_types, each=length(methods)),
                        dataset=rep(dataset,length(possible_dataset_types)*length(methods)))
  temp_df$scenario <- paste("Scenario", i)
  df <- rbind(df, temp_df)
}

df$dataset_type_index <-as.numeric(as.factor(df$dataset_type))
df$method <- sapply(df$method, str_replace, "music", "MuSiC")
df$method <- sapply(df$method, str_replace, "spotlight", "SPOTlight")
df$dataset_type <- sapply(df$dataset_type, str_replace, "artificial_", "")
df$dataset_type <- factor(df$dataset_type, levels=unique(df$dataset_type))
# mean
png("plots/mean_10reps_s1s2_braincortex.png", width=680, height=400)
ggplot(df, aes(x=method, y=mean, shape=dataset_type,
               color=scenario, group=dataset_type)) +
  geom_point(size=3, position=position_dodge(0.6)) +
  labs(title="", color="Dataset", shape="Dataset type") + xlab("Method") +
  ylab("RMSE") + scale_shape_manual(values=1:length(possible_dataset_types)) +
  ylim(0, 0.13) + guides(color = guide_legend(order=1), shape=guide_legend(order=2))
dev.off()

# Difference
dfs1 <- df[df$scenario=="Scenario 1",]
dfs2 <- df[df$scenario=="Scenario 2",]
for (method in methods){
  tempdfs1 <- dfs1[dfs1$method==method,]
  tempdfs2 <- dfs2[dfs2$method==method,]
  print(method)
  meandiff <- mean(tempdfs2$mean - tempdfs1$mean)
  scaled <- meandiff/mean(tempdfs1$mean)
  print(meandiff)
  print(scaled)
}

# SD
ggplot(df, aes(x=method, y=sd, shape=factor(str_replace(dataset_type, "artificial_", "")), color=str_replace(dataset, "_generation", ""))) +
  geom_jitter(width=0.2, alpha=0.7) +
  labs(title="SD of each method between reps", color="dataset", shape="dataset_type") + 
  ylab("Deviation of RMSE between 10 reps") + scale_shape_manual(values=1:length(possible_dataset_types)) +
  guides(color = guide_legend(order=1), shape=guide_legend(order=2)) + ylim(0, 0.1)