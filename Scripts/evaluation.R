#### EVALUATION OF RESULTS ####
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")
method <- methods[1]

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
      
      # Area under precision-recall curve
      known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% melt() %>% select(value)
      deconv_unlist <- lapply(deconv_list, function (k) c(as.matrix(k)))
      scores <- join_scores(deconv_unlist)
      model <- mmdata(scores, known_binary_all, modnames=methods) # make model
      curve <- evalmod(model)
      prcs <- subset(auc(curve), curvetypes == "PRC")
      
      # Get them into dataframe
      metrics <- data.frame(row.names=methods, "corr"=corr_spots, "RMSE"=RMSE,
                            "accuracy"=accuracy, "sensitivity"=sensitivity,
                            "specificity"=specificity, "precision"=precision, "F1"=F1,
                            "prc"=prcs$aucs)
      
      all_results[[dataset_type]] <- metrics
    }
    saveRDS(all_results, paste0(result_path, "all_metrics_", dataset, ".rds"))
  }
}

#### PLOTTING ####
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)
metric <- "RMSE"

for (metric in possible_metrics[8]){
  df <- data.frame()
  for (dataset in datasets[-1]){
    all_reps <- lapply(paste0('rep', 1:10), function(repl){
        all_results <- readRDS(paste0(paste0("results/", dataset, "/", repl, "_/"),
                                      "all_metrics_", dataset, ".rds"))
        sapply(possible_dataset_types, function (u) all_results[[u]][[metric]] )
    })
    temp_df <- data.frame(all_values=unlist(all_reps),
                          reps=rep(1:10, each=length(methods)*length(possible_dataset_types)),
                          method=c("SPOT", "MuSiC", "c2l", "RCTD", "stereo"),
                          dataset_type=rep(possible_dataset_types, each=length(methods)),
                          dataset=dataset)
    df <- rbind(df, temp_df)
  }
  
  proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                            "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)")
  names(proper_dataset_names) <- datasets[2:length(datasets)]
  df$dataset_type <- sapply(df$dataset_type, str_replace, "artificial_", "")
  
  # FACET GRID
  df <- mutate(df, dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20))
  df$dt_linebreak <- factor(df$dt_linebreak, levels=unique(df$dt_linebreak))
  
  p <- ggplot(df, aes(x=method, y=all_values, color=method)) + geom_boxplot(width=0.75) +
    ylab(paste0("Average ", proper_metric_names[metric])) + labs(color="Method") +
    scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
    theme(legend.position="bottom", legend.direction = "horizontal",
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(rows=vars(dataset), cols=vars(dt_linebreak), scales="free_y",
               labeller=labeller(dataset=proper_dataset_names))
  
  png(paste0("plots/metrics_facet/all_", metric, "_facet.png"), width=297, height=190, units="mm", res=200)
  print(p)
  dev.off()
  
  # ANNOTATING MEDIAN VALUE
  # png("plots/all_RMSE_facet_annote.png", width=297, height=190, units="mm", res=200)
  # p + stat_summary(geom="text", fun=median,
  #              aes(label=sprintf("%1.3f", ..y..), color=method),
  #              position=position_nudge(x=0.33), size=2)
  # dev.off()
  
  # PLOTTING REFERENCE LINE
  if (metric == "RMSE"){
    RMSE_ref_list <- melt(readRDS("rds/ref_RMSE_all.rds")) %>% group_by(L2, L1) %>% summarise(ref_RMSE=mean(value), .groups="drop")
    RMSE_dir_df2 <- RMSE_ref_list %>% `colnames<-`(c("dataset_type", "dataset", "value"))
    RMSE_dir_df2$dataset_type <- sapply(RMSE_dir_df2$dataset_type, str_replace, "artificial_", "")
    RMSE_dir_df2 <- mutate(RMSE_dir_df2, dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20))
    RMSE_dir_df2$dt_linebreak <- factor(RMSE_dir_df2$dt_linebreak, levels=unique(RMSE_dir_df2$dt_linebreak))
    p <- p + geom_hline(data=RMSE_dir_df2, aes(yintercept=value), linetype="dashed", color="gray")
  } else {
    ref_metric_list <- readRDS("rds/ref_all_metrics.rds")
    ref_metric_df <- reshape2::melt(ref_metric_list) %>% mutate("metric"=rep(possible_metrics, 56)) %>%
      `colnames<-`(c("value", "dataset_type", "dataset", "metric"))
    ref_metric_df <- ref_metric_df[ref_metric_df$metric==metric,]
    ref_metric_df$dataset_type <- sapply(ref_metric_df$dataset_type, str_replace, "artificial_", "")
    ref_metric_df <- mutate(ref_metric_df, dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20))
    ref_metric_df$dt_linebreak <- factor(ref_metric_df$dt_linebreak, levels=unique(ref_metric_df$dt_linebreak))
    p <- p + geom_hline(data=ref_metric_df, aes(yintercept=value), linetype="dashed", color="gray")
  }
  
  png(paste0("plots/metrics_facet/all_", metric, "_facet_ref.png"), width=297, height=190, units="mm", res=200)
  print(p)
  dev.off()
}

## BEST VALUES ##
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")

for (metric in possible_metrics[8]){
  df <- data.frame()
  for (dataset in datasets[-1]){
    all_reps <- lapply(paste0('rep', 1:10), function(repl){
      all_results <- readRDS(paste0(paste0("results/", dataset, "/", repl, "_/"),
                                    "all_metrics_", dataset, ".rds"))
      sapply(possible_dataset_types, function (u) all_results[[u]][[metric]] )
    })
    temp_df <- data.frame(all_RMSE=unlist(all_reps),
                          reps=rep(1:10, each=length(methods)*length(possible_dataset_types)),
                          method=c("SPOT", "MuSiC", "c2l", "RCTD", "stereo"),
                          dataset_type=rep(possible_dataset_types, each=length(methods)),
                          dataset=dataset)
    
    df <- rbind(df, temp_df)
  }
  df$dataset_type <- factor(df$dataset_type, levels=possible_dataset_types)
  best <- df %>% group_by(dataset, dataset_type, method) %>% summarise(medi=median(all_RMSE)) %>%
    slice_max(order_by=medi, n=1)
  write.table(best, paste0("Misc/best_values/best_values_", metric, ".tsv"), sep="\t", row.names=FALSE, quote=FALSE)
}

## Read in file from python
best_df <- read.csv("Misc/best_values/best_values_count.csv") %>% melt(id.var="X") %>% setNames(c("method", "metric", "count"))
best_df$metric <- best_df$metric %>% str_replace("prc", "PRC AUC") %>% R.utils::capitalize() %>%
  factor(., levels=c("Accuracy", "Specificity", "Sensitivity", "Precision", "F1", "PRC AUC"))
p <- ggplot(best_df, aes(x=metric, y=count, fill=method)) + geom_bar(width=0.75, position="fill", stat="identity") +
  ylab("% Best performing") + xlab("Metric") + labs(fill="Method") +
  scale_fill_manual(values = c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#e76bf3", "#a1a1a1"), 
                        labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope", "Tie"))
best_df$pct <- best_df %>% group_by(metric) %>% mutate(pct = count/sum(count)) %>% ungroup %>% select(pct)
png("plots/class_metrics_barplot.png", width=210, height=100, units="mm", res=200)
print(p)
dev.off()

#### WILCOXON TEST ####
### OLD - concatenate all values ###
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
      
      # RMSE_all spots
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

### NEW - mean of 10 iterations ###

all_reps <- lapply(paste0('rep', 1:10), function(repl){
  all_results <- readRDS(paste0(paste0("results/", dataset, "/", repl, "_/"),
                                "all_metrics_", dataset, ".rds"))
  sapply(possible_dataset_types, function (u) all_results[[u]][["RMSE"]] )
})
melt_all_reps <- melt(all_reps)
melt_all_reps$methods <- methods

pvalues <- list()
median_RMSEs <- list()
ref_RMSE_list <- readRDS("rds/ref_RMSE_all.rds")
# 15 = n*(n-1)/2, where n is the number of methods+reference
n_comparisons <- 15*length(possible_dataset_types)*(length(datasets)-1)
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
      
      # Mean RMSE
      RMSE <- sapply(deconv_list, function(k) mean(sqrt(rowSums((known_props-k[,1:ncells])**2, na.rm=TRUE)/ncells)))
      ref_RMSE <- mean(unlist(ref_RMSE_list[[dataset]][[dataset_type]])) %>% setNames("ref")
      stacked_RMSE <- rbind(stacked_RMSE, stack(c(ref_RMSE, RMSE)))
    }
    pval <- pairwise.wilcox.test(stacked_RMSE[,1], stacked_RMSE[,2],
                                 p.adjust.method = "none",
                                 paired=TRUE) 

    pvalues[[dataset]][[dataset_type]] <- melt(pval$p.value, na.rm=TRUE)
    median_RMSEs[[dataset]][[dataset_type]] <- stacked_RMSE %>% group_by(ind) %>%
      summarise(median_RMSE=median(values), .groups="drop")
  }
}
pvalues_melt <- melt(pvalues, value.name=c("pvals"), id.vars=c("Var1", "Var2"))
pvalues_melt$padj <- p.adjust(pvalues_melt$pvals, method="BY")

temp %>% dcast(Var1~Var2, value.var="padj")

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
    #temp_data <- pvalues_melt[pvalues_melt$L2==dataset_type & pvalues_melt$L1==dataset, c("Var1", "Var2", "padj")]
    #temp_data <- temp_data %>% dcast(Var1~Var2, value.var="padj")
    temp_data <- as.matrix(median_RMSEs[[dataset]][[dataset_type]])
    addDataFrame(temp_data, sheet=sheet, startRow=i+1, startColumn=j, row.names=FALSE, col.names=FALSE)
    i = i + nrow(temp_data) + 2
  }
  j = j+4
}

#saveWorkbook(wb, paste0("results/summary files/all_pvalues_wilcoxpaired_pairwise.xlsx"))
#saveWorkbook(wb, paste0("results/summary files/all_meanRMSEs.xlsx"))
saveWorkbook(wb, paste0("results/summary files/all_medianRMSEs.xlsx"))
