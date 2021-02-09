library(Seurat)
library(SeuratData)
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")


# Histograms for brain sections
sections = c("posterior1", "posterior2", "anterior1", "anterior2")
stats = data.frame("mean", "sd")

for (i in 1:4){
  section = sections[i]
  brain <- LoadData("stxBrain", type = section)
  brain <- colSums(as.matrix(GetAssayData(brain)))
  png(file=paste0("plots/data_distribution/stxBrain_", section, "_hist.png"))
  print(hist(brain, main=paste0("UMI distribution per spot of stxBrain, ", section)))
  mtext(paste0("Mean: ", round(mean(brain),2),
               "; SD: ", round(sd(brain), 2)), side=3)
  dev.off()
  stats[i,] = c(mean(brain), sd(brain))
}

# Histogram for kidney section
kidney <- LoadData("stxKidney")
kidney <- colSums(as.matrix(GetAssayData(kidney)))
png(file="plots/data_distribution/stxKidney_hist.png")
print(hist(kidney, main="UMI distribution per spot of stxKidney"))
mtext(paste0("Mean: ", round(mean(kidney),2),
             "; SD: ", round(sd(kidney), 2)), side=3)
dev.off()
print(c(mean(kidney), sd(kidney)))

# SMART-Seq mouse brain data
allen_cortex <- readRDS("rds/allen_cortex.rds")
allen_cortex <- colSums(as.matrix(GetAssayData(allen_cortex)))
print(hist(allen_cortex, main="UMI distribution per cell of SMART-Seq scRNA brain"))
print(c(mean(allen_cortex), sd(allen_cortex)))

# Load 10x mouse brain data
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered.rds")
DefaultAssay(brain_sc_10x) <- "RNA"
brain_sc_10x@assays$SCT = NULL
brain_sc_10x <- colSums(as.matrix(GetAssayData(brain_sc_10x)))
print(hist(brain_sc_10x, main="UMI distribution per cell of 10x brain"))
print(c(mean(brain_sc_10x), sd(brain_sc_10x)))

# Generating SPOTlight synthetic data from this
source("Scripts/synthetic_data_generation/spotlight_test_spot_fun.R")
spotlight_synth <- test_spot_fun(brain_sc_10x, clust_vr="subclass", n=2500)
# saveRDS(spotlight_synth, "Data/synthetic_datasets/spotlight_2500spots_brain10x_robin.rds")
hist(colSums(spotlight_synth$topic_profiles), breaks = c(16000, 18000, 20000, 22000))

# Creating h5ad file
# First, update object
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered.rds")
new_brain <- CreateSeuratObject(counts=brain_sc_10x@assays$RNA@counts, meta.data = brain_sc_10x@meta.data)
saveRDS(new_brain, "rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")

source("Scripts/helperFunctions.R")
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")
convertSeuratRDSToh5ad("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds") # Doesn't work, fails at Convert() 
convertSeuratRDSToLoom("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")
