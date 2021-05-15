#!/usr/bin/env Rscript
inputs = commandArgs(trailingOnly=TRUE)
if (length(inputs) < 2) {
  stop("usage: < dataset index > < replicate >", call.=FALSE)
}
# setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
setwd("~/thesis/")

library(Seurat)
library(synthvisium)
library(dplyr)
library(stringr)
source("Scripts/helperFunctions.R")

datasets <- c('allen_cortex_dwn', 'brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'pbmc', 'scc_p5')

dataset_type <- "artificial_uniform_distinct"
dataset <- datasets[as.integer(inputs[1])]
repl <- inputs[2]
k <- Sys.getenv("SGE_TASK_ID")

# SERVER
path <- paste0("~/thesis/downsampling/", dataset, "/", repl, "/")
results_path <- paste0(path, "results/")
scrna_path <- paste0("/group/irc/shared/synthetic_visium/test/", dataset, "_test.rds")

# LOCAL
# path <- paste0("D:/Work (Yr 2 Sem 1)/Thesis/downsampling/", dataset, "/", repl, "/")
# results_path <- paste0(path, "results/")
# scrna_path <- paste0("Data/synthetic_datasets/test_set/", dataset, "_test.rds")

######### SPOTLIGHT ######### 

library(SPOTlight)
dir.create(paste0(results_path, "spotlight/"), recursive=TRUE)

# Extract the top marker genes from each cluster
seurat_obj_scRNA <- readRDS(scrna_path)
seurat_obj_scRNA <- preprocessSeurat(seurat_obj_scRNA)
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA@meta.data$celltype
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA, assay = "SCT", slot = "data", verbose = TRUE,
                                      only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)

# Load synthetic visium data and downsample
synthetic_visium_data <- readRDS(paste0(path, list.files(path, pattern="rds")[1]))
spots <- unlist(read.table(paste0(path, "synthvisium_downsample_info/cells", k, ".txt")))
genes <- unlist(read.table(paste0(path, "synthvisium_downsample_info/genes", k, ".txt")))
synthetic_visium_data <- synthetic_visium_data$counts[genes, spots]
seurat_obj_visium <- CreateSeuratObject(counts = synthetic_visium_data, assay = "Spatial")

start_time <- Sys.time()
spotlight_deconv <- spotlight_deconvolution(se_sc = seurat_obj_scRNA, counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                            clust_vr = "celltype", cluster_markers = cluster_markers_all, cl_n = 50,
                                            hvg = 3000, ntop = NULL, transf = "uv", method = "nsNMF", min_cont = 0.09)
end_time <- Sys.time()
saveRDS(spotlight_deconv, paste0(results_path, "spotlight/", dataset, "_", k, "_spotlight.rds"))

job_info <- data.frame(JOB_ID=Sys.getenv("JOB_ID"),
                       SGE_TASK_ID=k,
                       START_TIME=start_time,
                       END_TIME=end_time,
                       TIME=end_time-start_time)
write.table(t(job_info), paste0(results_path, "spotlight/", dataset, "_", k, "_spotlight_info.txt"),
            quote = FALSE, sep="\t", row.names = TRUE, col.names = FALSE)

