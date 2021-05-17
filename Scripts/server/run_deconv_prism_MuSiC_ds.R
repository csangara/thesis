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
results_path <- paste0(path, "results_sc/")
# scrna_path <- paste0("/group/irc/shared/synthetic_visium/test/", dataset, "_test.rds")
scrna_path <- "/group/irc/shared/synthetic_visium/raw_data/brain_cortex/scRNAseq/seurat_obj_scrnaseq_cortex_filtered.rds"

# LOCAL
# path <- paste0("D:/Work (Yr 2 Sem 1)/Thesis/downsampling/", dataset, "/", repl, "/")
# results_path <- paste0(path, "results_sc/")
# scrna_path <- paste0("Data/synthetic_datasets/test_set/", dataset, "_test.rds")

######### MuSiC ######### 

library(MuSiC)
library(xbioc)
dir.create(paste0(results_path, "music/"), recursive=TRUE)

# Load reference scRNA-seq data and convert to ExprSet
seurat_obj_scRNA = readRDS(scrna_path)

### DOWNSAMPLING SCRNA-SEQ
cells <- unlist(read.table(paste0(path, "scref_downsample_info/ref_cells", k, ".txt")))
genes <- unlist(read.table(paste0(path, "scref_downsample_info/ref_genes", k, ".txt")))
seurat_obj_scRNA <- subset(seurat_obj_scRNA, cells=cells, features=genes)

seurat_obj_scRNA = seurat_obj_scRNA %>% SetIdent(value = "subclass")
DefaultAssay(seurat_obj_scRNA) <- "RNA"
eset_obj_scRNA <- SeuratToExprSet(seurat_obj_scRNA)
rm(seurat_obj_scRNA)

# Read in synthvisium
synthetic_visium_data <- readRDS(paste0(path, list.files(path, pattern="rds")[1]))
seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts, PP=FALSE)

### DOWNSAMPLING SYNTHVISIUM
# synthetic_visium_data <- readRDS(paste0(path, list.files(path, pattern="rds")[1]))
# spots <- unlist(read.table(paste0(path, "synthvisium_downsample_info/cells", k, ".txt")))
# genes <- unlist(read.table(paste0(path, "synthvisium_downsample_info/genes", k, ".txt")))
# synthetic_visium_data <- synthetic_visium_data$counts[genes, spots]
# seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data, PP=FALSE)

# Convert to ExprSet
eset_obj_visium <- SeuratToExprSet(seurat_obj_visium)
rm(synthetic_visium_data, seurat_obj_visium)

start_time <- Sys.time()
# Deconvolution
music_deconv = music_prop(bulk.eset = eset_obj_visium, sc.eset = eset_obj_scRNA, clusters = 'subclass', samples='samples')
end_time <- Sys.time()
# res[dataset_type] = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(music_deconv$Est.prop.weighted[,1:23]))
saveRDS(music_deconv$Est.prop.weighted, paste0(results_path, "music/", dataset, "_", k, "_music.rds"))
rm(music_deconv, eset_obj_visium)
gc()

job_info <- data.frame(JOB_ID=Sys.getenv("JOB_ID"),
           SGE_TASK_ID=k,
           START_TIME=start_time,
           END_TIME=end_time,
           TIME=end_time-start_time)
write.table(t(job_info), paste0(results_path, "music/", dataset, "_", k, "_music_info.txt"),
                                quote = FALSE, sep="\t", row.names = TRUE, col.names = FALSE)
