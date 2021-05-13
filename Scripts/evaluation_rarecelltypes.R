### DIRICHLET DISTRIBUTION ###
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")

## CALCULATING METRICS ##
dataset_type <- possible_dataset_types[3]

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
      
      if (dataset_type == "artificial_dominant_celltype_diverse" | 
          dataset_type == "artificial_partially_dominant_celltype_diverse"){
        # Determine dominant cell type
        dominant_cell_type <- which.max(colSums(known_props))
        known_props <- known_props[,-dominant_cell_type]
        deconv_list <- lapply(deconv_list, function(k) k[,-dominant_cell_type])
        ncells <- ncells-1
      } else {
        priors <- synthetic_visium_data$gold_standard_priorregion
        priors <- priors[priors$freq > 0,]
        rare_celltype <- priors$celltype[which.min(priors$freq)]
        known_props <- as.matrix(known_props[,which(colnames(known_props) == rare_celltype)])
        deconv_list <- lapply(deconv_list, function(k) as.matrix(k[,which(colnames(k) == rare_celltype)]))
        ncells <- 1
      }

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
    saveRDS(all_results, paste0(result_path, "all_metrics_rare", dataset, ".rds"))
  }
}

df <- data.frame()
for (dataset in datasets[-1]){
  all_reps <- lapply(paste0('rep', 1:10), function(repl){
    all_results <- readRDS(paste0(paste0("results/", dataset, "/", repl, "_/"),
                                  "all_metrics_rare", dataset, ".rds"))
    sapply(possible_dataset_types, function (u) all_results[[u]][["F1"]] )
  })
  
  mean_reps <- sapply(possible_dataset_types, function(u) {
    apply(sapply(1:10, function(k) all_reps[[k]][,u]), 1, mean)})
  sd_reps <- sapply(possible_dataset_types, function(u) {
    apply(sapply(1:10, function(k) all_reps[[k]][,u]), 1, sd)})
  
  temp_df <- data.frame(mean=c(mean_reps),
                        sd=c(sd_reps),
                        method=rep(methods, length(possible_dataset_types)),
                        dataset_type=rep(possible_dataset_types, each=length(methods)),
                        dataset=rep(dataset,length(possible_dataset_types)*length(methods)))
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

# MEAN
# by dataset type
ggplot(df, aes(shape= proper_dataset_names,
               y=mean, x=dataset_type_index, color=method, group=method)) +
  geom_point(size=2, position=position_dodge(0.6)) +
  labs(title="Mean F1 of each method between reps", color="method", shape="dataset") + 
  ylab("F1") + scale_shape_manual(values=1:length(datasets)) +
  guides(color = guide_legend(order=1), shape=guide_legend(order=2)) + 
  scale_x_discrete("dataset type", limits=factor(1:length(possible_dataset_types)))

# by dataset
png("plots/mean_between_reps_by_dataset.png", width=1200, height=500)
ggplot(df, aes(shape=factor(dataset_type, levels=unique(dataset_type)),
               y=mean, x=as.numeric(as.factor(dataset)), color=method, group=method)) +
  geom_point(size=2, position=position_dodge(0.6)) +
  labs(title="Mean RMSE of each method by dataset", color="method", shape="dataset type") + 
  ylab("RMSE") + scale_shape_manual(values=1:length(possible_dataset_types)) +
  guides(color = guide_legend(order=1), shape=guide_legend(order=2)) +
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
    temp <- df[df$dataset==dataset& df$dataset_type==dataset_type,]
    rownames(temp) <- temp$method
    temp <- temp[,"mean", drop=FALSE]
    addDataFrame(temp, sheet=sheet, startRow=i+1, startColumn=j, col.names=FALSE)
    
    i = i + nrow(temp) + 2
  }
  j = j + 4
}
saveWorkbook(wb, paste0("results/summary files/all_meanF1s.xlsx"))
