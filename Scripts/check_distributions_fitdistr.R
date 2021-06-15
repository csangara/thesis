library(Seurat)
library(fitdistrplus)
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")

##### CHECK DISTRIBUTIONS OF SYNTHETIC DATA WITH FITDISTR ######
#### HELPER FUNCTIONS ####
getDistributionInfo <- function(data, path, dataset_name, dist_type="discrete"){
  png(paste0(path, dataset_name, "_hist_cumdist.png"), width=800, height=600)
  plotdist(data, histo = TRUE, demp = TRUE)
  mtext(dataset_name, line=3)
  dev.off()
  # descdist(data, boot=1000) # another plot
  
  distributions <- c("norm", "pois", "nbinom")
  plot_legend <- c("normal", "poisson", "negative binomial")
  if (dist_type == "continuous") { 
    distributions <- c("norm", "beta", "gamma")
    data <- normalize(data)
    plot_legend <- c("normal", "beta", "gamma")
  }
  
  all_models <- list()
  for (distribution in distributions) {
    all_models[[distribution]] = fitdist(data, distribution)
  }
  
  png(paste0(path, dataset_name, "_QQplot_", dist_type, ".png"))
  # also denscomp(), cdfcomp(), ppcomp() but qq seems the most useful
  qqcomp(all_models, legendtext = plot_legend,
         main=paste("QQ-plot for", dataset_name))
  dev.off()
  
  gof <- gofstat(all_models, fitnames=distributions)
  gof_df <- data.frame(gof[c("chisq", "chisqpvalue", "ks", "cvm", "ad", "bic", "aic")])
  write.table(matrix(paste("GOF statistics for", dataset_name)),
              file=paste0(path, dataset_name, "_", dist_type, "_GOF.txt"),
              col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  write.table(round(gof_df, 4), file=paste0(path, dataset_name,  "_", dist_type, "_GOF.txt"),
              quote=FALSE, sep="\t", append=TRUE)
  write.table(matrix(c("ks = Kolomogrov-Smirnov statistic, cvm=Cramer-von Mises statistic, ad=Anderson-Darling statistic,
                       bic=Bayesian Information Criterion, aic=Akaike's Information criterion",
                       "ks, cvm, and ad are recommended for continuous data, while chisq is for discrete ones")),
              file=paste0(path, dataset_name,  "_", dist_type, "_GOF.txt"),
              col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)  
  
}

# Convert data from 0 to 1
normalize <- function(x) {
  x_scaled <- (x - min(x) + 0.001) / (max(x) - min(x) + 0.002)
  return(x_scaled)
}


###### CHECK DISTRIBUTION OF VISIUM DATA ###### 
# Histograms for brain sections "stxBrain"
library(SeuratData)
sections = c("posterior1", "posterior2", "anterior1", "anterior2")
stats = data.frame("mean", "sd")

for (i in 1:4){
  section = sections[i]
  brain <- LoadData("stxBrain", type = section)
  brain <- colSums(as.matrix(GetAssayData(brain)))
  png(file=paste0("plots/data_distribution_fitdistr/stxBrain_", section, "_hist.png"))
  print(hist(brain, main=paste0("UMI distribution per spot of stxBrain, ", section)))
  mtext(paste0("Mean: ", round(mean(brain),2),
               "; SD: ", round(sd(brain), 2)), side=3)
  dev.off()
  getDistributionInfo(brain, "plots/data_distribution_fitdistr/", paste0("stxBrain_", section), "continuous")
  getDistributionInfo(brain, "plots/data_distribution_fitdistr/", paste0("stxBrain_", section), "discrete")
  stats[i,] = c(mean(brain), sd(brain))
}

# Histogram for "stxKidney"
kidney <- LoadData("stxKidney")
kidney <- colSums(as.matrix(GetAssayData(kidney)))
png(file="plots/data_distribution_fitdistr/stxKidney_hist.png")
print(hist(kidney, main="UMI distribution per spot of stxKidney"))
mtext(paste0("Mean: ", round(mean(kidney),2),
             "; SD: ", round(sd(kidney), 2)), side=3)
dev.off()
print(c(mean(kidney), sd(kidney)))
getDistributionInfo(kidney, "plots/data_distribution_fitdistr/", "stxKidney", "continuous")
getDistributionInfo(kidney, "plots/data_distribution_fitdistr/", "stxKidney", "discrete")

###### CHECK DISTRIBUTION OF SCRNA-SEQ DATA ###### 
# SMART-Seq mouse brain data
allen_cortex <- readRDS("rds/allen_cortex.rds")
allen_cortex <- colSums(as.matrix(GetAssayData(allen_cortex)))
print(hist(allen_cortex, main="UMI distribution per cell of SMART-Seq scRNA brain"))
print(c(mean(allen_cortex), sd(allen_cortex)))

# 10x scRNA-seq mouse brain data (Robin's file)
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered.rds")
DefaultAssay(brain_sc_10x) <- "RNA"
brain_sc_10x <- colSums(as.matrix(GetAssayData(brain_sc_10x)))
print(hist(brain_sc_10x, main="UMI distribution per cell of 10x scRNA brain"))
print(c(mean(brain_sc_10x), sd(brain_sc_10x)))

# Create new object with rawcounts (to reduce file size)
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered.rds")
new_brain <- CreateSeuratObject(counts=brain_sc_10x@assays$RNA@counts, meta.data = brain_sc_10x@meta.data)
saveRDS(new_brain, "rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")

###### SYNTHVISIUM ######
synthetic_visium_data <- readRDS("Data/synthetic_datasets/for_comparing_dists/synthvisium_sc_10x_braincortex_robin_rawcounts_3224spots.rds")
png(file="plots/data_distribution_fitdistr/synthvisium_mean30000,sd8000_hist.png")
print(hist(colSums(synthetic_visium_data$counts), main="UMI distribution per spot of synthvisium synthetic data"))
mtext(paste0("Mean: ", round(mean(colSums(synthetic_visium_data$counts)),2),
             "; SD: ", round(sd(colSums(synthetic_visium_data$counts)), 2)), side=3)
dev.off()
getDistributionInfo(colSums(synthetic_visium_data$counts), "plots/data_distribution_fitdistr/",
                    "synthvisium_mean30000,sd8000", "continuous")
getDistributionInfo(colSums(synthetic_visium_data$counts), "plots/data_distribution_fitdistr/",
                    "synthvisium_mean30000,sd8000", "discrete")

#### SPOTLIGHT ####
spotlight_synth <- readRDS("Data/synthetic_datasets/for_comparing_dists/spotlight_sc_10x_braincortex_robin_rawcounts_2500spots.rds")
hist(colSums(spotlight_synth$topic_profiles), breaks = c(16000, 18000, 20000, 22000),
     main="UMI distribution per spot of SPOTlight synthetic data")

###### STEREOSCOPE ######
# There are three result files
expression_data <- read.table("Data/synthetic_datasets/for_comparing_dists/stereoscope_sc_10x_braincortex_robin_rawcounts_default/counts.st_synth.tsv",
                              header=TRUE, sep="\t", row.names=1)
png(file="plots/data_distribution_fitdistr/stereoscope_default_hist.png")
print(hist(rowSums(expression_data), main="UMI distribution per spot of stereoscope synthetic data"))
mtext(paste0("Mean: ", round(mean(rowSums(expression_data)),2),
             "; SD: ", round(sd(rowSums(expression_data)), 2)), side=3)
dev.off()
print(paste0("Mean: ", round(mean(rowSums(expression_data)),2),
             "; SD: ", round(sd(rowSums(expression_data)), 2)))
getDistributionInfo(rowSums(expression_data),
                    "plots/data_distribution_fitdistr/", "stereoscope_synthetic_data", "continuous")
getDistributionInfo(rowSums(expression_data),
                    "plots/data_distribution_fitdistr/", "stereoscope_synthetic_data", "discrete")

####### CELL2LOCATION ########
# Read in results
library(stringr)
# Get seeds
file_names <- list.files("Data/synthetic_datasets/for_comparing_dists/cell2location_sc_10x_braincortex_robin_rawcounts/")
seeds <- str_extract(gsub(".p", "", gsub(".*_", "", file_names)), regex("[0-9]+"))
seeds <- unique(seeds[!is.na(seeds)])

for (seed in seeds){
  expression_data <- read.table(paste0("Data/synthetic_datasets/for_comparing_dists/cell2location_sc_10x_braincortex_robin_rawcounts/2000spots/synthetic_ST_seed",
                                       seed, "_1_counts.csv"), header=TRUE, sep=",", row.names=1)
  png(file=paste0("plots/data_distribution_fitdistr/cell2location_seed", seed, "_2000spots_hist.png"))
  print(hist(rowSums(expression_data), main="UMI distribution per spot of cell2location synthetic data"))
  mtext(paste0("Mean: ", round(mean(rowSums(expression_data)),2),
               "; SD: ", round(sd(rowSums(expression_data)), 2)), side=3)
  dev.off()
  print(paste0("Mean: ", round(mean(rowSums(expression_data)),2),
               "; SD: ", round(sd(rowSums(expression_data)), 2)))
  
  # had to remove lognormal distribution from the function
  getDistributionInfo(rowSums(expression_data), "plots/data_distribution_fitdistr/",
                      paste0("cell2location_synthetic_data_seed", seed), "continuous")
  getDistributionInfo(rowSums(expression_data), "plots/data_distribution_fitdistr/",
                      paste0("cell2location_synthetic_data_seed", seed), "discrete")
}