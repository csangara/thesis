#### EVALUATION OF RESULTS ####

setwd("D:/Work (Yr 2 Sem 1)/Thesis/Scripts")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/"

library(Seurat)
library(synthvisium)
library(dplyr)
source("helperFunctions.R")
library(patchwork)
library(ggplot2)
library(stringr)

possible_dataset_types = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")

#### SAMPLE CODE FOR LOADING DATA ####
synthetic_visium_data <- readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds"))
seurat_obj_visium <- createAndPPSeuratFromVisium(synthetic_visium_data$counts)
spotlight_deconv = readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_spotlight.rds"))
music_deconv = readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_music.rds"))
c2l_deconv = read.csv(paste0(path, "Misc/cell2location/W_cell_density_q05.csv"), row.names=1)
c2l_deconv = c2l_deconv/rowSums(c2l_deconv) # So everything sums to one

res = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(spotlight_deconv[[2]][,1:23]))
print(mean(diag(res), na.rm=TRUE))


###### PLOT PREDICTIONS ON UMAP #######

for (dataset_type in possible_dataset_types){
  
  # Load reference data and deconvolution results
  result_path <- paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_")
  synthetic_visium_data <- readRDS(paste0(result_path, "synthvisium.rds"))
  seurat_obj_visium <- createAndPPSeuratFromVisium(synthetic_visium_data$counts)
  
  spotlight_deconv <- readRDS(paste0(result_path, "spotlight.rds"))
  music_deconv <- readRDS(paste0(result_path, "music.rds"))
  # cibersort_deconv = read.table(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_cibersort.txt"),
  #                               sep="\t", row.names=1, header=TRUE)
  c2l_deconv <- read.csv(paste0(path, "Misc/cell2location/W_cell_density_q05.csv"), row.names=1)
  c2l_deconv <- c2l_deconv/rowSums(c2l_deconv) # So everything sums to one
  colnames(c2l_deconv) <- str_replace(colnames(c2l_deconv), "q05_spot_factors", "")
  
  # Initialization of column names
  celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:23]
  celltypes <- str_replace(celltypes, "/", ".")
  colnames(synthetic_visium_data$relative_spot_composition)[1:23] <- celltypes
  colnames(spotlight_deconv[[2]])[1:23] <- celltypes
  colnames(music_deconv$Est.prop.weighted)[1:23] <- celltypes
  # colnames(cibersort_deconv)[1:23] = celltypes
  c2l_deconv <- c2l_deconv[, match(celltypes, colnames(c2l_deconv))]
  
  # Correlation (by spots and by cell type)
  known_props <- synthetic_visium_data$relative_spot_composition[,1:23]
  spotlight_corr_celltypes <- cor(known_props[,1:23], spotlight_deconv[[2]][,1:23], use="complete.obs")
  music_corr_celltypes <- cor(known_props[,1:23], music_deconv$Est.prop.weighted[,1:23], use="complete.obs")
  # cibersort_corr_celltypes <- cor(known_props[,1:23], cibersort_deconv[,1:23], use="complete.obs")
  c2l_corr_celltypes <- cor(known_props[,1:23], c2l_deconv[,1:23], use="complete.obs")
  
  for (celltype in celltypes){
    # Add deconv result to visium metadata
    seurat_obj_visium@meta.data[celltype] = synthetic_visium_data$relative_spot_composition[,celltype]
    seurat_obj_visium@meta.data[paste0(celltype, "_spotlight")] = spotlight_deconv[[2]][,celltype]
    seurat_obj_visium@meta.data[paste0(celltype, "_music")] = music_deconv$Est.prop.weighted[,celltype]
    # seurat_obj_visium@meta.data[paste0(celltype, "_cibersort")] = cibersort_deconv[,celltype]
    seurat_obj_visium@meta.data[paste0(celltype, "_c2l")] = c2l_deconv[,celltype]
    
    plot_dir <- paste0(path, "plots/allen_cortex_dwn/", dataset_type, "/")
    if (!dir.exists(plot_dir)){ dir.create(plot_dir) }
    
    plots <- FeaturePlot(seurat_obj_visium, c(celltype,
                                              paste0(celltype, "_spotlight"),
                                              paste0(celltype, "_music"),
    #                                         paste0(celltype, "_cibersort")),
                                              paste0(celltype, "_c2l"),
                                              combine=FALSE))
    
    plots[[2]] <- plots[[2]] + ggtitle("SPOTlight", paste0("Corr=", round(spotlight_corr_celltypes[celltype, celltype], 3)))
    plots[[3]] <- plots[[3]] + ggtitle("MuSiC", paste0("Corr=", round(music_corr_celltypes[celltype, celltype], 3)))
    # plots[[4]] <- plots[[4]] + ggtitle("CIBERSORT", paste0("Corr=", round(cibersort_corr_celltypes[celltype, celltype], 3)))
    plots[[4]] <- plots[[4]] + ggtitle("Cell2Location", paste0("Corr=", round(c2l_corr_celltypes[celltype, celltype], 3)))
    
    #plots <- plots[[1]] + plot_spacer() + plots[[2]] + plots[[3]]
    plots <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]]
    png(paste0(plot_dir, celltype, "_c2l.png"), width=1000, height=500)
    print(plots)
    dev.off()
  }
}

#### CALCULATE PERFORMANCE METRICS ####

# Get the confusion matrix given ground truth and deconvoluted result
getConfusionMatrix <- function(known_props, test_props){
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  missing_rows <- which(rowSums(is.na(known_props)) > 0)
  
  if (length(missing_rows) > 0){
    test_props <- test_props[-missing_rows,]
    known_props <- known_props[-missing_rows,]
  }
  for (i in 1:nrow(known_props)){
    for (j in 1:ncol(known_props)){
      if (known_props[i, j] > 0 & test_props[i, j] > 0){
        tp <- tp + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] == 0){
        tn <- tn + 1
      } else if (known_props[i, j] > 0 & test_props[i, j] == 0){
        fn <- fn + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] > 0){
        fp <- fp + 1
      }
    }
  }
  return(list(tp=tp, tn=tn, fn=fn, fp=fp))
}

deconv_methods <- c("spotlight", "music")
all_results <- list()
for (dataset_type in possible_dataset_types){
  
  # Load reference data and deconvolution results
  result_path <- paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_")
  
  synthetic_visium_data <- readRDS(paste0(result_path, "synthvisium.rds"))
  known_props <- synthetic_visium_data$relative_spot_composition[,1:23]
  n_celltypes <- 23
  
  metrics <- data.frame(row.names=c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1"))
  
  for (deconv_method in deconv_methods){
    deconv_matrix <- readRDS(paste0(result_path, deconv_method, ".rds"))
    if (deconv_method == "spotlight"){
      deconv_matrix <- deconv_matrix[[2]]
    } else if (deconv_method == "music"){
      deconv_matrix <- deconv_matrix$Est.prop.weighted
    }
    
    # Correlation and RMSE
    corr_spots <- mean(diag(cor(t(known_props), t(deconv_matrix[,1:n_celltypes]))), na.rm=TRUE)
    RMSE <- sqrt(sum((known_props-deconv_matrix[,1:n_celltypes])**2, na.rm=TRUE)/n_celltypes)
    
    # Classification scores
    conf <- getConfusionMatrix(known_props, deconv_matrix)
    tp <- conf$tp; tn <- conf$tn; fp <- conf$fp; fn <- conf$fn
    accuracy <- round((tp + tn) / (tp + tn + fp + fn), 2)
    sensitivity <- round(tp / (tp + fn), 2)
    specificity <- round(tn / (tn + fp), 2)
    precision <- round(tp / (tp + fp), 2)
    F1 <- round(2 * ((precision * sensitivity) / (precision + sensitivity)), 2)
    
    # Get them into dataframe
    metrics[deconv_method] <- c(corr_spots, RMSE, accuracy, sensitivity, specificity, precision, F1)
  }
  all_results[[dataset_type]] <- metrics
}

#### PLOT EACH METRIC ####
i <- 7
metric <- "F1" 
df <- data.frame(x = rep(names(all_results), length(deconv_methods)),
                      index = rep(1:13, length(deconv_methods)),
                      methods = rep(deconv_methods, each=13))
df$y <- c(sapply(deconv_methods, function(u) sapply(possible_dataset_types, function(k) all_results[[k]][[u]][i])))
ggplot(df, aes(x=factor(index), y=y, color=methods)) + geom_point(size=2) +
 labs(title=paste0(metric, " of MuSiC and SPOTlight on all 13 dataset types")) + ylab(metric) + xlab("datasets")

#### PLOT TIME OF CIBERSORT ####

df <- data.frame(x = c(10, 25, 50, 75, 100),
                 y = c(9.15, 31.4, 49.9, 81.1, 108.3))
ggplot(df, aes(x=x, y=y)) + geom_line() + geom_point() + xlab("Number of spots") + ylab("Time (min)") +
  labs(title="CIBERSORT runtime in function of spots")
