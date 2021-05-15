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

######### RCTD ##########

library(RCTD)
library(Matrix)
dir.create(paste0(results_path, "RCTD/"), recursive=TRUE)

# Load reference scRNA-seq
seurat_obj_scRNA <- readRDS(scrna_path)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
seurat_obj_scRNA@meta.data$liger_ident_coarse <- factor(seurat_obj_scRNA@meta.data$celltype)
seurat_obj_scRNA@meta.data$nUMI <- colSums(seurat_obj_scRNA@assays$RNA@counts)

# Load synthetic visium data and downsample
synthetic_visium_data <- readRDS(paste0(path, list.files(path, pattern="rds")[1]))
spots <- unlist(read.table(paste0(path, "synthvisium_downsample_info/cells", k, ".txt")))
genes <- unlist(read.table(paste0(path, "synthvisium_downsample_info/genes", k, ".txt")))

synthetic_visium_data <- synthetic_visium_data$counts[genes, spots]
spatialRNA_obj_visium <- RCTD:::SpatialRNA(counts=as(as(synthetic_visium_data,"matrix"),"dgCMatrix"))

start_time <- Sys.time()
RCTD_deconv <- create.RCTD(spatialRNA_obj_visium, seurat_obj_scRNA, max_cores = 8, CELL_MIN_INSTANCE=5)
RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = FALSE)
end_time <- Sys.time()
res = as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))

saveRDS(res, paste0(results_path, "RCTD/", dataset, "_", k, "_RCTD.rds"))

job_info <- data.frame(JOB_ID=Sys.getenv("JOB_ID"),
                       SGE_TASK_ID=k,
                       START_TIME=start_time,
                       END_TIME=end_time,
                       TIME=end_time-start_time)
write.table(t(job_info), paste0(results_path, "RCTD/", dataset, "_", k, "_RCTD_info.txt"),
            quote = FALSE, sep="\t", row.names = TRUE, col.names = FALSE)
