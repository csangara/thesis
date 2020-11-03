setwd("D:/Work (Yr 2 Sem 1)/Thesis/Scripts")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/"

library(Seurat)
library(synthvisium)
library(dplyr)
library(SPOTlight)
library(MuSiC)
library(xbioc)
source("helperFunctions.R")

######## GENERATING NEEDED DATA ############

# Preprocess scRNA reference data
seurat_obj_scRNA =  readRDS(paste0(path, "allen_cortex_dwn_original.rds"))
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

# Create synthetic data from scRNA data and save as RDS
possible_dataset_types = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")

for (dataset_type in possible_dataset_types){
  print(dataset_type)
  createSynthvisiumRDS(paste0(path,"rds/allen_cortex_dwn.rds"), dataset_type)
}

# Explore synthetic visium data
synthetic_visium_data$counts %>% as.matrix() %>% .[1:5,1:5] #  Gene counts for each spot
# Cell type composition of each spot
synthetic_visium_data$spot_composition %>% .[1:10,] # absolute
synthetic_visium_data$relative_spot_composition %>% .[1:10,] # relative
# Which cell type is present in which region, and with which prior frequency
synthetic_visium_data$gold_standard_priorregion %>% head()

# To create seurat object from synthetic visium data, use
seurat_obj_visium <- createAndPPSeuratFromVisium(synthetic_visium_data$counts)
# Visualize a priori defined regions vs clusters from gene expression
p_priorregion = DimPlot(seurat_obj_visium, reduction = "umap", label = TRUE, group.by = "orig.ident") # a priori defined regions
p_exprs_clusters = DimPlot(seurat_obj_visium, reduction = "umap", label = TRUE) #
patchwork::wrap_plots(list(p_priorregion, p_exprs_clusters), nrow = 1)

######### RUNNING SPOTLIGHT ######### 

# Extract the top marker genes from each cluster
seurat_obj_scRNA <- readRDS(paste0(path,"rds/allen_cortex_dwn.rds"))
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA@meta.data$subclass
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA, assay = "SCT", slot = "data", verbose = TRUE, 
                                              only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)

for (dataset_type in possible_dataset_types){
  set.seed(123)
  synthetic_visium_data <- readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createAndPPSeuratFromVisium(synthetic_visium_data$counts)
  spotlight_deconv <- spotlight_deconvolution(se_sc = seurat_obj_scRNA, counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                          clust_vr = "subclass", cluster_markers = cluster_markers_all, cl_n = 50,
                                          hvg = 3000, ntop = NULL, transf = "uv", method = "nsNMF", min_cont = 0.09)
  
  decon_mtrx <- spotlight_deconv[[2]]
  
  res = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(decon_mtrx[,1:23]))
  print(mean(diag(res), na.rm=TRUE))
  
  #res2 = cor(synthetic_visium_data$relative_spot_composition[,1:23], decon_mtrx[,1:23], use="complete.obs")
  #mean(diag(res2), na.rm=TRUE)
  
  saveRDS(spotlight_deconv, paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_spotlight.rds"))
}


######### MuSiC ######### 

# Load reference scRNA-seq data and convert to ExprSet
seurat_obj_scRNA = readRDS(paste0(path, "rds/allen_cortex_dwn.rds"))
eset_obj_scRNA <- SeuratToExprSet(seurat_obj_scRNA)
res = list()

for (dataset_type in possible_dataset_types){
  # Load synthetic visium data and convert to Expreset
  synthetic_visium_data <- readRDS(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createAndPPSeuratFromVisium(synthetic_visium_data$counts)
  eset_obj_visium <- SeuratToExprSet(seurat_obj_visium)
  
  # Deconvolution
  music_deconv = music_prop(bulk.eset = eset_obj_visium, sc.eset = eset_obj_scRNA, clusters = 'celltype', samples='samples')
  res[dataset_type] = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(music_deconv$Est.prop.weighted[,1:23]))
  saveRDS(music_deconv, paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_music.rds"))
}


