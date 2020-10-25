library(Seurat)
library(synthvisium)
library(dplyr)
library(SPOTlight)

seurat_obj_scRNA =  readRDS("../rds/allen_cortex_dwn.rds")

DimPlot(seurat_obj_scRNA, group.by = "subclass", label = T)
DimPlot(seurat_obj, group.by = "brain_subregion")

# Preprocess scRNA data
seurat_obj_scRNA <- SCTransform(seurat_obj_scRNA, assay="RNA", verbose = FALSE)
seurat_obj_scRNA <- RunPCA(seurat_obj_scRNA, verbose = FALSE)
seurat_obj_scRNA <- RunUMAP(seurat_obj_scRNA, dims = 1:30, verbose = FALSE)
seurat_obj_scRNA <- FindNeighbors(seurat_obj_scRNA, dims = 1:30, verbose = FALSE)
seurat_obj_scRNA <- FindClusters(seurat_obj_scRNA, verbose = FALSE)

# Set celltype
seurat_obj_scRNA@meta.data$celltype = seurat_obj_scRNA@meta.data$subclass
seurat_obj_scRNA = seurat_obj_scRNA %>% SetIdent(value = "celltype")

DimPlot(seurat_obj_scRNA, reduction = "umap",pt.size = 0.5, label = T)

# Create synthetic visium data from scRNA data, then preprocess
synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj_scRNA, dataset_type = "real", 
                                            clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                            n_spots_min = 50, n_spots_max = 200, visium_mean = 20000, visium_sd = 5000)

# Explore synthetic visium data
synthetic_visium_data$counts %>% as.matrix() %>% .[1:5,1:5] #  Gene counts for each spot
# Cell type composition of each spot
synthetic_visium_data$spot_composition %>% .[1:10,] # absolute
synthetic_visium_data$relative_spot_composition %>% .[1:10,] # relative
# Which cell type is present in which region, and with which prior frequency
synthetic_visium_data$gold_standard_priorregion %>% head() 

seurat_obj_visium = CreateSeuratObject(counts = synthetic_visium_data$counts, min.cells = 2, min.features = 200, assay = "Spatial")
seurat_obj_visium = SCTransform(seurat_obj_visium, assay = "Spatial", verbose = FALSE)
seurat_obj_visium = RunPCA(seurat_obj_visium, assay = "SCT", verbose = FALSE)
seurat_obj_visium = RunTSNE(seurat_obj_visium, reduction = "pca", dims = 1:30)
seurat_obj_visium = RunUMAP(seurat_obj_visium, reduction = "pca", dims = 1:30)
seurat_obj_visium = FindNeighbors(seurat_obj_visium, reduction = "pca", dims = 1:30)
seurat_obj_visium = FindClusters(seurat_obj_visium, verbose = FALSE, resolution = 0.5)

# Visualize a priori defined regions vs clusters from gene expression
p_priorregion = DimPlot(seurat_obj_visium, reduction = "umap", label = TRUE, group.by = "orig.ident") # a priori defined regions
p_exprs_clusters = DimPlot(seurat_obj_visium, reduction = "umap", label = TRUE) # 
patchwork::wrap_plots(list(p_priorregion, p_exprs_clusters), nrow = 1)



######### SPOTLIGHT ######### 

# Extract the top marker genes from each cluster
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA@meta.data$subclass
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA, assay = "SCT", slot = "data", verbose = TRUE, 
                                              only.pos = TRUE, logfc.threshold = 1, min.pct = 0.9)

set.seed(123)
spotlight_deconv <- spotlight_deconvolution(se_sc = seurat_obj_scRNA, counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                        clust_vr = "subclass", cluster_markers = cluster_markers_all, cl_n = 50,
                                        hvg = 3000, ntop = NULL, transf = "uv", method = "nsNMF", min_cont = 0.09)

decon_mtrx <- spotlight_deconv[[2]]

res = diag(cor(t(synthetic_visium_data$relative_spot_composition[,1:23]), t(decon_mtrx[,1:23])))
mean(diag(res), na.rm=TRUE)

res2 = cor(synthetic_visium_data$relative_spot_composition[,1:23], decon_mtrx[,1:23], use="complete.obs")
mean(diag(res2), na.rm=TRUE)


######## PLOT PROPS #########
library(ggplot2)

n = nrow(decon_mtrx)
knownP = synthetic_visium_data$relative_spot_composition[1:n,1:23]
predP = decon_mtrx[1:n,1:23]
new_pred = data.frame("pred" = c(t(predP)))
new_pred$known = as.numeric(c(t(knownP)))
new_pred$celltype = rep(colnames(knownP), n)

ggplot(new_pred, aes(x=known, y=pred, shape=celltype)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="gray20") +
  scale_shape_manual(values = 1:23) + geom_point(size=1) +
  labs(x="Known Proportions", y="Predicted Proportions")
  
