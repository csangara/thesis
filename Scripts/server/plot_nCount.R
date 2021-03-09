data_path <- "/group/irc/shared/synthetic_visium/raw_data/"

paths <- c("brain_cortex/scRNAseq/seurat_obj_scrnaseq_cortex_filtered.rds",
           "cerebellum/seurat_obj_sn_sc_cerebellum.rds", 
           "real/seurat_obj_hippocampus_filtered.rds",
           "scRNAseq/seurat_obj_scrnaseq_kidney_filtered.rds",
           "pbmc_data/seurat_obj_pbmc_filtered.rds",
           "scRNAseq/seurat_obj_scrnaseq_scc_p5_filtered.rds")

datasets <- c("brain_cortex", "cerebellum_cell", "cerebellum_nucleus",
              "hippocampus", "kidney", "pbmc", "scc_p5")
p <- list()
for (i in 1:length(datasets)){
  seuratObj <- readRDS(paste0(data_path, paths[i]))
  Idents(seuratObj) <- datasets[i]
  p[[i]] <- VlnPlot(seuratObj, features="nCount_RNA")
  png(paste0("~/thesis/plots/", datasets[i], "_nCount_RNA.png"))
  print(VlnPlot(seuratObj, features="nCount_RNA"))
  dev.off()
}