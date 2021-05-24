source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")
library(caret)
library(reshape2)
library(patchwork)

## EVALUATING "DOMINANT CELL TYPE (5)" AND "DOMINANT RARE CELL TYPE (7)" ##
dataset_type <- possible_dataset_types[7] # 5 or 7

all_results <- list()
for (dataset in datasets[2:length(datasets)]){
  for (repl in paste0('rep', 1:10)){
    print(paste(dataset, repl))
    result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
    metrics <- list()
    
    # Load reference data
    synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                            dataset_type, "_synthvisium.rds"))
    ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
    n_spots <- ncol(synthetic_visium_data$counts)
    
    # Initialization of column names
    celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:ncells]
    celltypes <- str_replace(celltypes, "/", ".")
    colnames(synthetic_visium_data$relative_spot_composition)[1:ncells] <- celltypes
    known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]

    # Load deconvolution results
    deconv_list <- createDeconvResultList(methods, celltypes, result_path, dataset)
    
    # Convert to binary vectors (one vs all) and put in a list of length ncells
    known_binary <- lapply(celltypes, function (celltype) {
        factor(ifelse(known_props[[celltype]] > 0, celltype, "Others"),
               levels=c(celltype, "Others"))}) %>% setNames(., celltypes)
    
    # Do this for each method, but round the results (helps c2l, RCTD and stereoscope do better)
    pred_binary <- lapply(methods, function (method) {
      lapply(celltypes, function (celltype) {
        factor(ifelse(round(deconv_list[[method]][,celltype],2) > 0, celltype, "Others"),
        levels=c(celltype, "Others"))  }) %>% setNames(., celltypes)
      }) %>% setNames(., methods)
    
    # Get confusion matrix for each cell type
    all_confusion_mats <- lapply(methods, function(method) {
      lapply(celltypes, function(celltype) {
        confusionMatrix(pred_binary[[method]][[celltype]], known_binary[[celltype]])
      }) %>% setNames (., celltypes)
      }) %>% setNames(., methods)
      
    metric_names <- c("Sensitivity", "Specificity", "Precision", "F1", "Balanced Accuracy")
    
    # Create a list of metrics, divided by each cell type
    metrics_by_group <- lapply(metric_names, function(metric_name) {
      sapply(methods, function(method) {
          sapply(celltypes, function(celltype) {
            all_confusion_mats[[method]][[celltype]]$byClass[metric_name]
          })}) %>% `rownames<-`(celltypes) }) %>% setNames(., metric_names)
    
    # Add accuracy and RMSE
    metrics_by_group[["Accuracy"]] <- sapply(methods, function(method) {
      sapply(celltypes, function(celltype) {
        all_confusion_mats[[method]][[celltype]]$overall["Accuracy"]})
      }) %>% `rownames<-`(celltypes)
    
    metrics_by_group[["RMSE"]] <-  sapply(deconv_list, function(k) 
      sqrt(colSums((known_props-k[,1:ncells])**2, na.rm=TRUE)/n_spots))
    
    metrics_avg <- sapply(metrics_by_group, function(metrics) apply(metrics, 2, mean))
    
    if (dataset_type == "artificial_dominant_celltype_diverse" | 
        dataset_type == "artificial_partially_dominant_celltype_diverse"){
      # Determine dominant cell type
      dominant_cell_type <- which.max(colSums(known_props))
      print(dominant_cell_type)
      metrics[["dominant"]] <- sapply(metrics_by_group, function(k) k[dominant_cell_type,])
      metrics[["others"]] <- sapply(metrics_by_group, function(k) apply(k[-dominant_cell_type,], 2, mean))
      
      # Do this if you want to plot values of all groups instead of averaging them
      #metrics[["others"]] <- melt(lapply(metrics_by_group, function(k) k[-dominant_cell_type,])) %>%
      #  `colnames<-`(c("celltype", "method", "value","metrics"))
      
    } else {
      # Rare cell type
      priors <- synthetic_visium_data$gold_standard_priorregion
      priors <- priors[priors$freq > 0,]
      rare_celltype <- priors[which.min(priors$freq),]$celltype
      print(rare_celltype)
      metrics[["rare"]] <- sapply(metrics_by_group, function(k) k[rare_celltype,])
      metrics[["others"]] <- sapply(metrics_by_group, function(k) apply(k[rownames(k) != rare_celltype,], 2, mean))
      
      # Also added code for computing metrics per region, but this isn't used in the paper
      if (dataset_type == "artificial_dominant_rare_celltype_diverse"){
        region_label=synthetic_visium_data$relative_spot_composition$region
        confusion_mats_by_region <- lapply(methods, function(method) {
          lapply(paste0("priorregion", 1:5), function(region) {
            confusionMatrix(pred_binary[[method]][[rare_celltype]][region_label==region],
                            known_binary[[rare_celltype]][region_label==region])
          }) %>% setNames (., paste0("priorregion", 1:5))
        }) %>% setNames(., methods)
        
        metrics_by_region <- lapply(metric_names, function(metric_name) {
          sapply(methods, function(method) {
            sapply(paste0("priorregion", 1:5), function(region) {
              confusion_mats_by_region[[method]][[region]]$byClass[metric_name]
            }) })  }) %>% setNames(., metric_names)
        #print(metrics_by_region[["F1"]])
        #region_comp <- getregionComp(synthetic_visium_data, verbose=F)
        #print(region_comp[grepl(rare_celltype, region_comp)])
      }
    }
    all_results[[dataset]][[repl]] <- metrics
  }
}

## PLOTTING TIME ##
df <- melt(all_results)
colnames(df) <- c("method", "metrics", "value", "type", "rep", "dataset")
df <- df[!grepl("Accuracy", df$metrics),]
chosen <- c("dominant", "others","rare", "both")[4] # for dominant, choose 1,2,4, for rare, choose 1,3,4 
title <- c("dominant", "rare")[sum(grepl("dominant", dataset_type), grepl("rare", dataset_type))]

# Note: NA values occur because there are no true negatives in some cases
if (chosen != "both"){ 
  temp_df <- df[df$type==chosen,]
  temp_df$method <- factor(temp_df$method, levels=sort(methods))
  p <- ggplot(temp_df[temp_df$metrics!="RMSE",], aes(x=method, color=method, y=value)) +
    geom_boxplot(width=0.5) +   ylab("Score") + labs(color="Method") +
    scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
    theme(legend.position="bottom", legend.direction = "horizontal",
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(cols=vars(metrics))
  png(paste0("plots/rare_celltypes_evaluation/", title, "_cell_type_facet_", chosen, "_together.png"), width=210, height=80, units="mm", res=200)
  print(p)
  dev.off()
} else {
  p <- ggplot(df[df$metrics!="RMSE",], aes(x=factor(method, levels=sort(methods)), color=type, y=value)) +
    geom_boxplot(width=0.5) + ylab("Score") + labs(color="Type") + xlab("Methods") +
    facet_grid(cols=vars(metrics))
  png(paste0("plots/rare_celltypes_evaluation/", title, "_cell_type_facet_", chosen, "_together.png"), width=210, height=80, units="mm", res=200)
  print(p)
  dev.off()
}

## Combine results for all datasets together (Appendix)
chosen <- c("dominant", "others","rare", "both")[4] # for dominant, choose 1,2,4, for rare, choose 1,3,4 
title <- c("dominant", "rare")[sum(grepl("dominant", dataset_type), grepl("rare", dataset_type))]
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)")
proper_method_names <- c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")
names(proper_dataset_names) <- datasets[2:length(datasets)]

if (chosen != "both"){
  temp_df <- df[df$type==chosen,]
  temp_df$method <- factor(temp_df$method, levels=sort(methods))
  p1 <- ggplot(temp_df[temp_df$metrics!="RMSE",], aes(x=method, color=method, y=value)) +
    geom_boxplot(width=0.5) + ylab("Score") + labs(color="Method") +
    scale_color_discrete(labels=proper_method_names) +
    facet_grid(cols=vars(metrics), rows=vars(dataset)) +
    theme(strip.background.y = element_blank(), strip.text.y = element_blank())
  
  p2 <- ggplot(temp_df[temp_df$metrics=="RMSE",], aes(x=method, color=method, y=value)) +
    geom_boxplot(width=0.5) + ylab("Score") + labs(color="Method") +
    scale_color_discrete(labels=proper_method_names) +
    facet_grid(cols=vars(metrics), rows=vars(dataset), labeller =labeller(dataset=proper_dataset_names)) +
    theme(axis.title.y=element_blank())
  
  combined <- p1 + p2 & theme(legend.position = "bottom", legend.direction = "horizontal",
                              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  
  png(paste0("plots/rare_celltypes_evaluation/", title, "_cell_type_facet_", chosen, "_all.png"), width=297, height=210, units="mm", res=200)
  print(combined + plot_layout(widths = c(4, 1), guides="collect"))
  dev.off()
} else {
  
  p1 <- ggplot(df[df$metrics!="RMSE",], aes(x=factor(method, levels=sort(methods)), color=type, y=value)) +
    geom_boxplot(width=0.5) + ylab("Score") + labs(color="Type") + xlab("Method") +
    facet_grid(cols=vars(metrics), rows=vars(dataset)) +
    theme(strip.background.y = element_blank(), strip.text.y = element_blank())
  
  p2 <- ggplot(df[df$metrics=="RMSE",], aes(x=factor(method, levels=sort(methods)), color=type, y=value)) +
    geom_boxplot(width=0.5) + ylab("Score") + labs(color="Type") + xlab("Method") +
    facet_grid(cols=vars(metrics), rows=vars(dataset), labeller=labeller(dataset=proper_dataset_names)) +
    theme(axis.title.y=element_blank())
  
  combined <- p1 + p2 & theme(legend.position = "bottom", legend.direction = "horizontal",
                              axis.title.x=element_blank())
  
  png(paste0("plots/rare_celltypes_evaluation/", title, "_cell_type_facet_", chosen, "_all.png"), width=297, height=210, units="mm", res=200)
  print(combined + plot_layout(widths = c(4, 1), guides="collect"))
  dev.off()
}

### PRINTING THE BEST PERFORMERS ###
summaries <- df %>% group_by(metrics, dataset, method) %>% summarise(med = median(value))
for (metric in unique(df$metrics)){
  for (dataset in datasets[2:length(datasets)]){
    print(paste(dataset, metric))
    if (metric == "RMSE"){
      # print(summaries[summaries$dataset==dataset & summaries$metrics == metric, ])
      print(methods[which.min(summaries[summaries$dataset==dataset & summaries$metrics == metric, ]$med)])
    } else {
      print(methods[which.max(summaries[summaries$dataset==dataset & summaries$metrics == metric, ]$med)])
    }
  }
}


### PLOTTING PRECISION-RECALL CURVE ("micro") ###
library(MLeval)
deconv_melt <- melt(lapply(deconv_list, data.frame))
deconv_melt$absent <- 1-deconv_melt$value
known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% melt() %>% select(value)
deconv_melt <- cbind(deconv_melt, known_binary_all)
colnames(deconv_melt) <- c("celltype", "present", "Group", "absent", "obs")
deconv_melt <- deconv_melt[, c("absent", "present", "obs", "Group")]
test <- evalm(deconv_melt)

### BARPLOTS ####
library(RColorBrewer)
dataset <- datasets[2]
repl <- 'rep4'
result_path <- paste0("results/", dataset, "/", repl, "_", run, "/")
dataset_type <- possible_dataset_types[7]
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

dominant_cell_type <- which.max(colSums(known_props))
deconv_temp <- deconv_list %>% purrr::list_modify("known"=known_props)
mean_props <- melt(sapply(deconv_temp, colMeans))
colnames(mean_props) <- c("celltype", "method", "mean_prop")
mean_props_without_dom <- mean_props[mean_props$celltype != names(dominant_cell_type),]
my_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(ncells)
ggplot(mean_props, aes(fill=factor(celltype, levels=celltypes), y=mean_prop, x=method)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  #scale_fill_brewer(palette="Dark2")
  scale_fill_manual(values=my_colors) +
  ylab("Proportions") + xlab("Method") +labs(fill="Cell types")
  
### CODE DUMP ###
#geom_point(position=position_jitterdodge(), size=0.01, alpha=0.5) # adding dots to the graph
#df <- melt(all_results, id.vars=c("celltype", "method", "metrics"))
#colnames(df) <- c("celltype", "method", "metrics", "var", "value", "type", "rep", "dataset")
#colnames(df)[9:10] <- c("repl", "dataset")
# gridExtra::grid.arrange(p1, p2, ncol=2, widths=c(3.5, 1.5)) + theme(legend.position="horizontal")

  