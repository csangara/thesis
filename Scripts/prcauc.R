#### EVALUATION OF PRECISION-RECALL ####
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")
library(precrec)
method <- methods[1]

#### CALCULATE PERFORMANCE METRICS ####
for (dataset in datasets[2:length(datasets)]){
  for (dataset_type in possible_dataset_types){
    all_matrices <- list()
    all_known_matrices <- list()
    nspots_vec <- c()
    for (repl in paste0('rep', 1:10)){
      result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
      
      # Load reference data and deconvolution results
      synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                              dataset_type, "_synthvisium.rds"))
      ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
      nspots_vec <- c(nspots_vec, nrow(synthetic_visium_data$spot_composition))
      
      # Initialization of column names
      celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:ncells]
      celltypes <- str_replace(celltypes, "/", ".")
      colnames(synthetic_visium_data$relative_spot_composition)[1:ncells] <- celltypes
      known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
      known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% as.matrix()
      all_known_matrices[[repl]] <- known_binary_all
      
      # Load deconvolution results
      deconv_list <- createDeconvResultList(methods, celltypes, result_path, dataset)
      all_matrices[[repl]] <- lapply(deconv_list, as.matrix)
    }
    
    # set.seed(10)
    # sampled_rows <- sample.int(min(nspots_vec))
    # all_matrices_ds <- lapply(1:10, function (k) lapply(all_matrices[[k]],
    #                                   function (u) c(u[sampled_rows,])))
    # all_known_matrices_ds <- lapply(1:10, function(k) c(all_known_matrices[[k]][sampled_rows,]))
    # scores <- join_scores(all_matrices_ds)
    # labels <- join_labels(rep(all_known_matrices_ds, each=5))
    
    scores <- join_scores(lapply(1:10, function (k) lapply(all_matrices[[k]], c)), chklen=FALSE)
    labels <- join_labels(rep(lapply(all_known_matrices, c), each=5), chklen=FALSE)
    
    # Make mode
    model <- mmdata(scores, labels, modnames=rep(methods, 10), dsids = rep(1:10, each=5))
    curve <- evalmod(model)
    autoplot(curve, "PRC", show_cb = TRUE)
    prcs <- subset(auc(curve), curvetypes == "PRC")
    
    # Get them into dataframe
    metrics <- data.frame(row.names=methods, "corr"=corr_spots, "RMSE"=RMSE,
                          "accuracy"=accuracy, "sensitivity"=sensitivity,
                          "specificity"=specificity, "precision"=precision, "F1"=F1)
      
      
    
    
  }
}
