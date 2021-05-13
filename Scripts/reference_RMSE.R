### CALCULATE REFERENCE RMSE USING DIRICHLET DISTRIBUTION ###
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")

dataset <- datasets[8]
RMSE_avg_list = list()
RMSE_dir_list = list()
for (dataset_type in possible_dataset_types){
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                          dataset_type, "_synthvisium.rds"))
  ncells <- ncol(synthetic_visium_data$spot_composition)-2
  nspots <- nrow(synthetic_visium_data$spot_composition)
  known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
  
  all_dir_dist = c()
  for (i in 1:1000){
    dir_dist <- rdirichlet(nspots, rep(1.0, ncells))
    RMSE_dir <- mean(sqrt(rowSums((known_props-dir_dist)**2)/ncells))
    all_dir_dist[i] = RMSE_dir
  }    
  RMSE_dir_list[[dataset_type]] <- mean(all_dir_dist)
  RMSE_avg_list[[dataset_type]] <- mean(sqrt(rowSums((known_props-(1/ncells))**2)/ncells))
}
sapply(names(RMSE_dir_list), function(k) print(c(RMSE_dir_list[[k]], RMSE_avg_list[[k]])))


### PLOTTING ###

# REFERENCE RMSE
ref_RMSE_df <- data.frame()
for (dataset in datasets[2:length(datasets)]){
  all_reps <- lapply(c('rep1', 'rep2', 'rep3'), function(repl){
    all_results <- lapply(possible_dataset_types, function(dataset_type) {
      # Load reference data and deconvolution results
      synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                              dataset_type, "_synthvisium.rds"))
      ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
      known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
      
      mean(sqrt(rowSums((known_props-(1/ncells))**2)/ncells))
    })
  })
  mean_reps <- sapply(1:length(possible_dataset_types), function(u) {
    mean(sapply(c(1,2,3), function(k) all_reps[[k]][[u]]))})
  sd_reps <- sapply(1:length(possible_dataset_types), function(u) {
    sd(sapply(c(1,2,3), function(k) all_reps[[k]][[u]]))})
  
  
  temp_df <- data.frame(mean=mean_reps,
                        sd=sd_reps,
                        dataset_type=possible_dataset_types,
                        dataset=rep(dataset,length(possible_dataset_types)))
  ref_RMSE_df <- rbind(ref_RMSE_df, temp_df)
}


# Plot just the lines
colors <- c("#F8766D", "#C49A00", "#53B400", "#00C094",
            "#00B6EB", "#A58AFF", "#FB61D7")

for (i in 1:7){
  png(paste0("plots/ref_RMSE_", datasets[i+1], ".png"), width=1200, height=500)
  print(ggplot(ref_RMSE_df[ref_RMSE_df$dataset==datasets[i+1],]) +
          geom_hline(aes(yintercept=mean), color=colors[i]) +
          geom_point(aes(x=1, y=mean, shape=factor(dataset_type)), color=colors[i], size=3) +
          ylim(0, 0.275) + scale_shape_manual(values=1:length(possible_dataset_types)) +
          theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_blank()))
  dev.off()
}
