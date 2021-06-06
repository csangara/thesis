#### PLOTTING PRECISION-RECALL CURVEL ####
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")
library(precrec)
method <- methods[1]

i = 1
all_plots = list()
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
      known_binary_all <- ifelse(known_props > 0, "present", "absent") 
      all_known_matrices[[repl]] <- c(as.matrix(known_binary_all))
      
      # Load deconvolution results
      deconv_list <- createDeconvResultList(methods, celltypes, result_path, dataset)
      all_matrices[[repl]] <- lapply(sort(methods), function(k) c(as.matrix(deconv_list[[k]])))
    }
    
    # Create score and labels
    scores <- join_scores(all_matrices, chklen=FALSE)
    labels <- join_labels(rep(all_known_matrices, each=5), chklen=FALSE)
    
    # Make model
    proper_method_names <- c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")
    model <- mmdata(scores, labels, modnames=rep(proper_method_names, 10), dsids = rep(1:10, each=5))
    curve <- evalmod(model)
    p <- autoplot(curve, "PRC", show_cb = TRUE)
    
    all_plots[[i]] <- p
    i = i + 1
    
  }
}
# Best performing (PRCAUC) for the border
best_prc <- read.table("Misc/best_values/best_values_prc.tsv", sep="\t", header=TRUE)
colors <- c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#e76bf3") %>%
  setNames(c("c2l", "MuSiC", "RCTD", "SPOT", "stereo"))

# For each plot, turn off legend, axis labels, and color the border
all_plots_mod <- all_plots
for (i in 1:56){
  temp_p <- all_plots_mod[[i]]
  
  temp_p <- temp_p + theme(plot.title = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      legend.position="none",
                      axis.text.y = element_blank(),
                      axis.text.x = element_blank(),
                      plot.margin=unit(c(1,0,1,0), "mm"),
                      panel.border = element_rect(colour = colors[best_prc$method[i]],
                                                  size = 1, linetype = "solid"))
  
  all_plots_mod[[i]] <- temp_p
}

# Arrnage in a grid
png(paste0("plots/metrics_facet/all_prcurves_colored.png"), width=267, height=170, units="mm", res=200)
gridExtra::grid.arrange(grobs=all_plots_mod, nrow=7, ncol=8)
dev.off()

# This is to plot the legend by itself
legend <- cowplot::get_legend(p + theme(legend.direction="horizontal"))
grid::grid.newpage()
png(paste0("plots/metrics_facet/all_prcurves_onlylegend.png"), width=150, height=20, units="mm", res=200)
grid::grid.draw(legend)
dev.off()
