setwd("D:/Work (Yr 2 Sem 1)/Thesis")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
library(Seurat)
library(synthvisium)
library(dplyr)
source("Scripts/helperFunctions.R")

possible_dataset_types = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")

recommended_dataset_types = c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap",
                              "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                              "artificial_partially_dominant_celltype_diverse", "artificial_dominant_rare_celltype_diverse",
                              "artificial_regional_rare_celltype_diverse")

datasets <- c('allen_cortex_dwn', 'brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
              'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation')

######## GENERATING NEEDED DATA ############
# Create synthetic data from scRNA data and save as RDS
for (dataset_type in possible_dataset_types){
  print(dataset_type)
  createSynthvisiumRDS("rds/allen_cortex_dwn_original.rds", dataset_type,
                       output_folder="Data/synthetic_datasets/allen_cortex_dwn/")
}

dataset <- datasets[1]
scrna_path <- "rds/allen_cortex_dwn_original.rds"
  
######### SPOTLIGHT ######### 

library(SPOTlight)

# Extract the top marker genes from each cluster
seurat_obj_scRNA <- readRDS(scrna_path)
seurat_obj_scRNA <- preprocessSeurat(seurat_obj_scRNA)
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA@meta.data$subclass
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA, assay = "SCT", slot = "data", verbose = TRUE, 
                                              only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)

for (dataset_type in possible_dataset_types){
  set.seed(123)
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts)
  start_time <- Sys.time()
  spotlight_deconv <- spotlight_deconvolution(se_sc = seurat_obj_scRNA, counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                          clust_vr = "subclass", cluster_markers = cluster_markers_all, cl_n = 50,
                                          hvg = 3000, ntop = NULL, transf = "uv", method = "nsNMF", min_cont = 0.09)
  end_time <- Sys.time()
  decon_mtrx <- spotlight_deconv[[2]]
  
  res = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(decon_mtrx[,1:23]))
  print(mean(diag(res), na.rm=TRUE))
  
  #res2 = cor(synthetic_visium_data$relative_spot_composition[,1:23], decon_mtrx[,1:23], use="complete.obs")
  #mean(diag(res2), na.rm=TRUE)
  
  saveRDS(spotlight_deconv, paste0("results/", dataset, "/spotlight/", dataset, "_", dataset_type, "_spotlight.rds"))
}

######### MuSiC ######### 

library(MuSiC)
library(xbioc)

# Load reference scRNA-seq data and convert to ExprSet
seurat_obj_scRNA = readRDS(scrna_path)
seurat_obj_scRNA@meta.data$celltype = seurat_obj_scRNA@meta.data$subclass
seurat_obj_scRNA = seurat_obj_scRNA %>% SetIdent(value = "celltype")
eset_obj_scRNA <- SeuratToExprSet(seurat_obj_scRNA)
res = list()

for (dataset_type in possible_dataset_types){
  # Load synthetic visium data and convert to Expreset
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  seurat_obj_visium <- createSeuratFromCounts(synthetic_visium_data$counts, PP=FALSE)
  eset_obj_visium <- SeuratToExprSet(seurat_obj_visium)
  
  # Deconvolution
  music_deconv = music_prop(bulk.eset = eset_obj_visium, sc.eset = eset_obj_scRNA, clusters = 'celltype', samples='samples')
  res[dataset_type] = cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(music_deconv$Est.prop.weighted[,1:23]))
  saveRDS(music_deconv$Est.prop.weighted, paste0("results/", dataset, "/music/", dataset, "_", dataset_type, "_music.rds"))
}

######### RCTD ##########

library(RCTD)
library(Matrix)

seurat_obj_scRNA <- readRDS(scrna_path)
seurat_obj_scRNA@meta.data$liger_ident_coarse <- factor(seurat_obj_scRNA@meta.data$subclass)
seurat_obj_scRNA@meta.data$nUMI <- colSums(seurat_obj_scRNA@assays$RNA@counts)

for (dataset_type in possible_dataset_types[2:length(possible_dataset_types)]){
  
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", dataset, "_", dataset_type, "_synthvisium.rds"))
  spatialRNA_obj_visium <- RCTD:::SpatialRNA(counts=as(as(synthetic_visium_data$counts,"matrix"),"dgCMatrix"))
  
  start_time <- Sys.time()
  RCTD_deconv <- create.RCTD(spatialRNA_obj_visium, seurat_obj_scRNA, max_cores = 4, CELL_MIN_INSTANCE = 5)
  RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = FALSE)
  end_time <- Sys.time()
  res = as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))
  
  saveRDS(res, paste0("results/", dataset, "/RCTD/", dataset, "_", dataset_type, "_RCTD.rds"))
}