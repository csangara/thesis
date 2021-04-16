### MISCELLANEOUS (PRE)PROCESSING SCRIPTS ###

setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/"
source("Scripts/helperFunctions.R")

possible_dataset_types = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")

### PREPROCESSING DATA ###
# Preprocess scRNA reference data
seurat_obj_scRNA =  readRDS("rds/allen_cortex_dwn_original.rds")
seurat_obj_scRNA <- SCTransform(seurat_obj_scRNA, assay="RNA", verbose = FALSE)
seurat_obj_scRNA <- RunPCA(seurat_obj_scRNA, verbose = FALSE)
seurat_obj_scRNA <- RunUMAP(seurat_obj_scRNA, dims = 1:30, verbose = FALSE)
seurat_obj_scRNA <- FindNeighbors(seurat_obj_scRNA, dims = 1:30, verbose = FALSE)
seurat_obj_scRNA <- FindClusters(seurat_obj_scRNA, verbose = FALSE)
# Set cell type and object identity
seurat_obj_scRNA@meta.data$celltype = seurat_obj_scRNA@meta.data$subclass
seurat_obj_scRNA = seurat_obj_scRNA %>% SetIdent(value = "celltype")
DimPlot(seurat_obj_scRNA, reduction = "umap",pt.size = 0.5, label = T)

saveRDS(seurat_obj_scRNA, paste0(path,"rds/allen_cortex_dwn.rds"))
seurat_obj_scRNA <- readRDS(paste0(path, "rds/allen_cortex_dwn.rds"))
DimPlot(seurat_obj_scRNA, reduction = "umap", label = TRUE, group.by = "celltype")

### EXPLORING SYNTHETIC VISIUM DATA ###
dataset_type <- "artificial_diverse_distinct"
synthetic_visium_data <- readRDS("rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds")

# Explore synthetic visium data
synthetic_visium_data$counts %>% as.matrix() %>% .[1:5,1:5] #  Gene counts for each spot
# Cell type composition of each spot
synthetic_visium_data$spot_composition %>% .[1:10,] # absolute
synthetic_visium_data$relative_spot_composition %>% .[1:10,] # relative
# Which cell type is present in which region, and with which prior frequency
synthetic_visium_data$gold_standard_priorregion %>% head()

# To create seurat object from synthetic visium data, use
seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts)
# Visualize a priori defined regions vs clusters from gene expression
p_priorregion = DimPlot(seurat_obj_visium, reduction = "umap", label = TRUE, group.by = "orig.ident") # a priori defined regions
p_exprs_clusters = DimPlot(seurat_obj_visium, reduction = "umap", label = TRUE) #
patchwork::wrap_plots(list(p_priorregion, p_exprs_clusters), nrow = 1)

### DOWNSAMPLING SPOTS FOR CIBERSORT ###
input_path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/CIBERSORTx/cibersort_mixture_ArtUniDistinct.txt"
reduceSpotsCS(input_path, 75)

### CONVERTING SEURAT OBJECT TO LOOM ###
for (dataset_type in possible_dataset_types){
  print(dataset_type)
  convertSeuratRDSToLoom(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type,
                                "_synthvisium.rds"), isSeurat = FALSE)
}

## SAVING NON-NORMALIZED DATA AS H5AD ##
# Reference file
convertSeuratRDSToh5ad("D:/Work (Yr 2 Sem 1)/Thesis/allen_cortex_dwn_original.rds")

# Synthvisium files
for (dataset_type in possible_dataset_types[10:13]){
  print(dataset_type)
  input_path <- paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds")
  convertSeuratRDSToh5ad(input_path, isSeurat = FALSE, PP= FALSE, raw=TRUE)
}

#### PLOT TIME OF CIBERSORT ####

df <- data.frame(x = c(10, 25, 50, 75, 100),
                 y = c(9.15, 31.4, 49.9, 81.1, 108.3))
ggplot(df, aes(x=x, y=y)) + geom_line() + geom_point() + xlab("Number of spots") + ylab("Time (min)") +
  labs(title="CIBERSORT runtime in function of spots")

######## PLOT PROPS #########
# library(ggplot2)
# 
# n = nrow(decon_mtrx)
# knownP = synthetic_visium_data$relative_spot_composition[1:n,1:23]
# predP = decon_mtrx[1:n,1:23]
# new_pred = data.frame("pred" = c(t(predP)))
# new_pred$known = as.numeric(c(t(knownP)))
# new_pred$celltype = rep(colnames(knownP), n)
# 
# ggplot(new_pred, aes(x=known, y=pred, shape=celltype)) +
#   geom_abline(slope=1, intercept=0, linetype=2, colour="gray20") +
#   scale_shape_manual(values = 1:23) + geom_point(size=1) +
#   labs(x="Known Proportions", y="Predicted Proportions")

#### UMAP of generation data ####
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
datasets <- c('brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
              'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation')
for (dataset in datasets){
  seurat_obj <- readRDS(paste0(path, "generation_set/", dataset, ".rds"))
  png(paste0("D:/Work (Yr 2 Sem 1)/Thesis/plots/UMAP_generation/", dataset, ".png"))
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "celltype"))
  dev.off()
}

# UMAP of test
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
test_set <- c('brain_cortex_test.rds', 'cerebellum_cell_test.rds', 'cerebellum_nucleus_test.rds',
              'hippocampus_test.rds', 'kidney_test.rds',
              'pbmc_test.rds', 'scc_p5_test.rds')
for (dataset in test_set){
  seurat_obj <- readRDS(paste0(path, "test_set/", dataset))
  png(paste0("D:/Work (Yr 2 Sem 1)/Thesis/plots/UMAP_test/", dataset, ".png"))
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "celltype"))
  dev.off()
}

#### VIOLIN PLOT nCount_RNA ####
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/"
library(patchwork)
library(ggplot2)
ncount_list <- readRDS(paste0(path, "rds/nCounts_list.rds"))
datasets <- c("brain_cortex", "cerebellum_cell", "cerebellum_nucleus",
              "hippocampus", "kidney", "pbmc", "scc_p5")
df <- data.frame()
foi <- "nFeature_RNA"
for (i in 1:length(datasets)){
  temp_df <- data.frame(counts = ncount_list[[i]][,foi],
                        logcounts = log10(ncount_list[[i]][,foi]),
                        dataset = rep(datasets[i], length(ncount_list[[i]][,foi])))
  
  df <- rbind(df, temp_df)
}
# sp <- split(df$counts, df$dataset)
# a <- lapply(seq_along(sp), function(i){
#   d <- density(sp[[i]])
#   k <- which.max(d$y)
#   data.frame(dataset = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
# })
# a <- do.call(rbind, a)
ggplot(df, aes(x=counts, color=dataset, fill=dataset)) + geom_density(alpha=0.2) +
  #geom_text(data = a, aes(x = xmax, y = ymax, label = dataset)) + xlim(0, 25000)
  scale_x_continuous(trans='log10') + labs(title=foi)
