inputs = commandArgs(trailingOnly=TRUE)
if (length(inputs) < 1) {
    stop("Please enter whether you want sc, sn, or both references", call.=FALSE)
}

#!/usr/bin/env Rscript

setwd("~/thesis/")
path <- "~/data/"
library(Seurat)
library(dplyr)
library(stringr)
source("Scripts/helperFunctions.R")

# Load in file using Seurat, and convert to SpatialRNA (do this locally)
# See: https://rdrr.io/github/dmcable/RCTD/src/R/SpatialRNA.R (read.VisiumSpatialRNA)
# seurat_obj_visium <- Load10X_Spatial("Data/Liver/JBO01")
# coords <- seurat_obj_visium@images$slice1@coordinates
# coords <- coords[,c("row", "col")]
# colnames(coords) <- c("x", "y")
# 
# GetAssayData(seurat_obj_visium)
# spatialRNA_obj_visium <- RCTD:::SpatialRNA(coords = coords,
#                                        counts = GetAssayData(seurat_obj_visium))
# saveRDS(spatialRNA_obj_visium, "Data/Liver/JBO01_spatialRNA.rds")

results_path <- paste0("~/data/Liver/results/")
scrna_path <- "~/data/Liver/anndatafineAnnot_new.rds"

######### RCTD ##########

library(RCTD)
library(Matrix)
dir.create(paste0(results_path, "RCTD/"), recursive=TRUE)

seurat_obj_scRNA <- readRDS(scrna_path)
techtype <- inputs[1]
if (techtype == 'sc'){
    seurat_obj_scRNA <- seurat_obj_scRNA[,seurat_obj_scRNA$type=="RnaSeq"]
} else if (techtype == 'sn'){
    seurat_obj_scRNA <- seurat_obj_scRNA[,seurat_obj_scRNA$type=="nucSeq"]
}
seurat_obj_scRNA <- seurat_obj_scRNA[,seurat_obj_scRNA$fine_annot != "Capsular Fibroblasts"]

DefaultAssay(seurat_obj_scRNA) <- "RNA"
seurat_obj_scRNA@meta.data$liger_ident_coarse <- factor(seurat_obj_scRNA@meta.data$fine_annot)
seurat_obj_scRNA@meta.data$nUMI <- colSums(seurat_obj_scRNA@assays$RNA@counts)

spatialRNA_obj_visium <- readRDS("~/data/Liver/JBO01_spatialRNA.rds")

RCTD_deconv <- create.RCTD(spatialRNA_obj_visium, seurat_obj_scRNA, max_cores = 8)
RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = FALSE)
res = as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))

saveRDS(res, paste0(results_path, "RCTD/liver_deconv_RCTD_", techtype, ".rds"))
