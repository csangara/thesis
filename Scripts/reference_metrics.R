### CALCULATE REFERENCE RMSE USING DIRICHLET DISTRIBUTION ###
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")
library(DirichletReg)

metric_dir_list = list()
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
for (dataset in datasets[2:length(datasets)]){
  for (dataset_type in possible_dataset_types){
    synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                            dataset_type, "_synthvisium.rds"))
    ncells <- ncol(synthetic_visium_data$spot_composition)-2
    nspots <- nrow(synthetic_visium_data$spot_composition)
    known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
    known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% melt() %>% select(value)
    
    iters = 100
    all_metrics <- matrix(, nrow=iters, ncol=length(possible_metrics))
    for (i in 1:iters){
      dir_dist <- rdirichlet(nspots, rep(1.0, ncells))
      
      RMSE <- mean(sqrt(rowSums((known_props-dir_dist)**2)/ncells))
      corr <- mean(diag(cor(t(known_props), t(dir_dist[,1:ncells]))), na.rm=TRUE)

      conf <- getConfusionMatrix(known_props, dir_dist)
      accuracy <- round((conf$tp + conf$tn) / (conf$tp + conf$tn + conf$fp + conf$fn), 2)
      sensitivity <- round(conf$tp / (conf$tp + conf$fn), 2)
      specificity <- round(conf$tn / (conf$tn + conf$fp), 2)
      precision <- round(conf$tp / (conf$tp + conf$fp), 2)
      F1 <- round(2 * ((precision * sensitivity) / (precision + sensitivity)), 2)

      model <- mmdata(c(dir_dist), known_binary_all)
      curve <- evalmod(model)
      prcs <- subset(auc(curve), curvetypes == "PRC")
      
      all_metrics[i,] <- sapply(possible_metrics, get) %>% setNames(possible_metrics)
    }
    
    colnames(all_metrics) <- possible_metrics
    metric_dir_list[[dataset]][[dataset_type]] <- colMeans(all_metrics)
    
  }
}

saveRDS(metric_dir_list, "rds/ref_all_metrics.rds")
#metric_dir_list <- readRDS("rds/ref_all_metrics.rds")

## EXPORTING TABLE
ref_metric_list <- readRDS("rds/ref_all_metrics.rds")
ref_metric_df <- reshape2::melt(ref_metric_list) %>% mutate("metric"=rep(possible_metrics, 56)) %>%
  `colnames<-`(c("value", "dataset_type", "dataset", "metric"))
ref_metric_df$dataset_type <- factor(ref_metric_df$dataset_type, levels=possible_dataset_types)

for (metric in possible_metrics[8]){
  temp_metric_df <- ref_metric_df[ref_metric_df$metric==metric,]
  
  temp_metric_print <- t(dcast(temp_metric_df, dataset_type~dataset, value.var="value"))
  colnames(temp_metric_print) <- temp_metric_print[1,]
  temp_metric_print <- temp_metric_print[-1,]
  write.table(temp_metric_print, paste0("Misc/reference_metrics/ref_", metric, ".tsv"), quote=FALSE, sep="\t")
}
