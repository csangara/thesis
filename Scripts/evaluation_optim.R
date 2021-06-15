#### EVALUATION OF (OPTIMIZED) RESULTS ####
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")

#### CALCULATE PERFORMANCE METRICS OF SCENARIO 1####
dataset <- datasets[2]
methods_optim <- c(methods, c("spotlight_optim", "stereoscope_optim"))

for (repl in paste0('rep', 9:10)){
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
    deconv_list <- createDeconvResultList(methods_optim, celltypes, result_path, dataset)
    
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
    F1 <- sapply(methods_optim, function(k) round(2 * ((precision[k] * sensitivity[k]) /
                                                   (precision[k] + sensitivity[k])), 2))
    # Get them into dataframe
    metrics <- data.frame(row.names=methods_optim, "corr"=corr_spots, "RMSE"=RMSE,
                          "accuracy"=accuracy, "sensitivity"=sensitivity,
                          "specificity"=specificity, "precision"=precision, "F1"=F1)
    
    all_results[[dataset_type]] <- metrics
    saveRDS(all_results, paste0(result_path, "all_metrics_", dataset, "_withoptim.rds"))
  }
}

### STABILITY PLOTS ###
for (dataset in datasets[3:4]){
  result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
  file_name <- paste0(result_path, "spotlight/", dataset, "_", dataset_type, "_spotlight")
  #file_name_opt <- paste0(result_path, "spotlight_optim/", dataset, "_", dataset_type, "_spotlight")
  deconv <- readRDS(paste0(file_name, ".rds"))
  #deconv_opt <- readRDS(paste0(file_name_opt, ".rds"))
  nmf_mod <- deconv[[1]]
  h <- NMF::coef(nmf_mod[[1]])
  rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
  topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
    h = h,
    train_cell_clust = nmf_mod[[2]])
  
  print(topic_profile_plts[[2]] + ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90), 
    axis.text = ggplot2::element_text(size = 12)) +
    ggtitle(dataset))
}

possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)
all_colors <- c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#00b0f6", "#e76bf3", "#e76bf3") %>%
  setNames(c("c2l", "MuSiC", "RCTD", "SPOT", "SPOT_opt", "stereo", "stereo_opt"))
metric <- "RMSE"
proper_dataset_names <- c("brain_cortex_generation"= "Brain cortex")

df <- data.frame()
for (metric in possible_metrics[c(2, 4 ,7)]){
  all_reps <- lapply(paste0('rep', 1:10), function(repl){
    all_results <- readRDS(paste0(paste0("results/", dataset, "/", repl, "_/"),
                                  "all_metrics_", dataset, "_withoptim.rds"))
    sapply(possible_dataset_types, function (u) all_results[[u]][[metric]] )
  })
  temp_df <- data.frame(all_values=unlist(all_reps),
                        reps=rep(1:10, each=length(methods_optim)*length(possible_dataset_types)),
                        method=c("SPOT", "MuSiC", "c2l", "RCTD", "stereo", "SPOT_opt", "stereo_opt"),
                        dataset_type=rep(possible_dataset_types, each=length(methods_optim)),
                        dataset=dataset,
                        metric=metric)
  
  df <- rbind(df, temp_df)
}

df$dataset_type <- sapply(df$dataset_type, str_replace, "artificial_", "")
df$method <- factor(df$method, levels=sort(unique(df$method)))
df <- df %>% mutate(optim=grepl("opt", df$method),
                    method_group=sapply(df$method, function(k) str_split(k, "_")[[1]][1]))
df$metric <- factor(df$metric, levels=c("RMSE", "F1", "sensitivity"))
df$width <- ifelse(grepl("c2l|MuSiC|RCTD", df$method), 10, 90)
# FACET GRID
df <- mutate(df, dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20))
df$dt_linebreak <- factor(df$dt_linebreak, levels=unique(df$dt_linebreak))

p <- ggplot(df, aes(x=method_group, y=all_values, color=method, fill=optim)) +
  geom_boxplot(position=position_dodge2()) +
  ylab(paste0("Average ", proper_metric_names[metric])) + labs(color="Method", fill="Optimized") +
  scale_color_manual(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight",
                                "SPOTlight_opt",  "stereoscope", "stereoscope_opt"),
                       values=all_colors) +
  scale_fill_manual(values=c("white", "yellow"), labels=c("No", "Yes")) +
  theme(legend.position="bottom", legend.direction = "horizontal", legend.box="vertical",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(cols=vars(metric), scales="free_y", #cols=vars(dt_linebreak), 
             labeller=labeller(metric=proper_metric_names[unique(as.character(df$metric))]))

png(paste0("plots/metrics_facet/optim_facet.png"), width=297, height=190, units="mm", res=200)
print(p)
dev.off()

