#### EVALUATION OF RESULTS ####
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")
dataset <- datasets[3]

#### CALCULATE PERFORMANCE METRICS OF SCENARIO 4####
for (dataset in datasets[3:4]){
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
      saveRDS(all_results, paste0(result_path, "all_metrics_", dataset, ".rds"))
    }
  }
}

# Deviation between scenarios
library(patchwork)
for (dataset in datasets[3:4]){
  df <- data.frame()
  for (i in 1:2){
    ext <- c("", "_s4")[i]
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
    temp_df$scenario <- c("Scenario 2", "Scenario 3")[i]
    df <- rbind(df, temp_df)
  
  }
  
  df$dataset_type_index <-as.numeric(as.factor(df$dataset_type))
  df$method <- factor(df$method, levels=sort(methods))
  proper_method_names <- c("c2l", "MuSiC", "RCTD", "SPOT", "stereo")
  names(proper_method_names) <- sort(methods)
  df$dataset_type <- sapply(df$dataset_type, str_replace, "artificial_", "")
  df$dataset_type <- factor(df$dataset_type, levels=unique(df$dataset_type))
  title <- ifelse(dataset == "cerebellum_cell_generation", "Cerebellum (sc)", "Cerebellum (sn)")
  # mean
  # png("plots/scenario3_cell.png", width=680, height=400) # or scenario3_nucleus.png
  p_temp <- ggplot(df, aes(x=method, y=mean, shape=dataset_type,
                 color=scenario, group=dataset_type)) +
    geom_point(size=2, stroke=1, position=position_dodge(0.6)) +
    labs(title=title, color="Scenario", shape="Dataset type") + xlab("Method") +
    ylab("RMSE") + scale_shape_manual(values=1:length(possible_dataset_types)) +
    ylim(0, 0.3) + guides(color = guide_legend(order=1), shape=guide_legend(order=2)) +
    scale_x_discrete(labels=proper_method_names)
  # dev.off()
  if (dataset == "cerebellum_cell_generation"){
    p <- p_temp
  } else {
    p <- p + p_temp
  }
}
p <- p + plot_layout(guides = "collect")
png("plots/scenario3_both.png", width=210, height=100, units="mm", res=200)
print(p)
dev.off()