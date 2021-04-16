library(Seurat)
library(SeuratData)
path <- "D:/Work (Yr 2 Sem 1)/Thesis/"
cortex <- readRDS(paste0(path, "rds/cortex.rds"))
allen_reference <- readRDS(paste0(path, "rds/allen_cortex.rds"))
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE)
allen_reference <-  RunPCA(allen_reference, assay = "SCT", verbose = FALSE)
allen_reference <-  RunUMAP(allen_reference, reduction = "pca", dims = 1:30)