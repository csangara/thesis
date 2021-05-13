source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/init.R")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"

###### SCRIPT FOR DOWNSAMPLING SYNTHVISIUM AND SC-REF DATA ######
#### DOWNSAMPLING SYNTHVISIUM SPOTS ####
downsampleSynthvisium <- function(counts, n_cells, n_features, seed=10){
  set.seed(seed)
  cell_index <- sample.int(ncol(counts), n_cells)
  cell_names <- colnames(counts)[cell_index]
  feature_index <- sample.int(nrow(counts), n_features)
  feature_names <- rownames(counts)[feature_index]
  
  return (counts[feature_names, cell_names])
  
}

# Generate synthetic data with at least 10000 spots
#synthetic_visium_data <- generate_synthetic_visium(seurat_obj_scRNA, dataset_type,
#                                                   "celltype", n_regions=8, n_spots_min=2000, n_spots_max=2500)
synthvisium_filepath <- paste0("downsampling/", str_replace(dataset, "_generation", ""),
                               "/brain_cortex_art_uni_distinct_18594spots_17538genes.rds")
#saveRDS(synthetic_visium_data, synthvisium_filepath)
synthetic_visium_data <- readRDS(synthvisium_filepath)

# Convert to h5ad
convertSeuratRDSToh5ad(synthvisium_filepath, output_path=NULL, isSeurat=FALSE, raw=FALSE, update=FALSE)

# Make dataframe of combinations (10 cells & 10 genes, 10 cells & 100 genes, etc...)
all_scales <- gtools::permutations(4,2,10^(1:4), repeats.allowed=TRUE)
colnames(all_scales) <- c("n_cells", "n_genes")

# Write conversion file (only run this once)
# write.table(all_scales, file = paste0("downsampling/", str_replace(dataset, "_generation", ""),
#                                       "/synthvisium_downsample_info/conversion.txt"),
#                                       sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Generate list of cells and genes for each case
for (k in 1:nrow(all_scales)){
  n_cells <- all_scales[k,1]
  n_features <- all_scales[k,2]
  print(c(n_cells, n_features))
  ds_count <- downsampleSynthvisium(synthetic_visium_data$counts, n_cells, n_features)

  write.table(colnames(ds_count), file = paste0("downsampling/", str_replace(dataset, "_generation", ""),
              "/synthvisium_downsample_info/cells", k, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(rownames(ds_count), file = paste0("downsampling/", str_replace(dataset, "_generation", ""),
              "/synthvisium_downsample_info/genes", k, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#### DOWNSAMPLING SINGLE-CELL REFERENCE DATA ####
downsampleSCref <- function(cells_df, genes, n_cells, n_features, seed=10){
  set.seed(seed)
  cell_pct <- n_cells/nrow(cells_df)
  # This is so the sampling is stratified
  chosen_cells <- cells_df %>% group_by(celltype) %>% sample_frac(cell_pct)
  
  feature_index <- sample.int(length(genes), n_features)
  feature_names <- genes[feature_index]
  
  return (list(cells=chosen_cells$cells, genes=feature_names))
  
}

all_scales <- gtools::permutations(4,2,sort(c(10^(2:4), 5000)), repeats.allowed=TRUE)
colnames(all_scales) <- c("n_cells", "n_genes")
# Write conversion file (only run this once)
# write.table(all_scales, file = paste0("downsampling/", str_replace(dataset, "_generation", ""),
#                                       "/scref_downsample_info/conversion.txt"),
#                                        sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Only get the names of cells and genes (and metadata) to save space
seurat_obj_scRNA <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered.rds")
cells_df <- data.frame(cells=colnames(seurat_obj_scRNA), celltype = seurat_obj_scRNA$subclass)
genes <- rownames(seurat_obj_scRNA)

# Generate list of cells and genes for each case
for (k in 1:nrow(all_scales)){
  n_cells <- all_scales[k,1]
  n_features <- all_scales[k,2]
  print(c(n_cells, n_features))
  ds_object <- downsampleSCref(cells_df, genes, n_cells, n_features)
  
  write.table(ds_object$cells, file = paste0("downsampling/", str_replace(dataset, "_generation", ""),
                                                "/scref_downsample_info/ref_cells", k, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(ds_object$genes, file = paste0("downsampling/", str_replace(dataset, "_generation", ""),
                                                "/scref_downsample_info/ref_genes", k, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Check equal proportions
library(ggplot2)
df <- data.frame(vals = c(table(ds_object$cells), table(seurat_obj_scRNA$celltype)),
                 celltypes = names(table(ds_object$cells)),
                 source=rep(c("downsampled", "full"), each=length(unique(seurat_obj_scRNA$celltype))))
ggplot(df, aes(fill=celltypes, y=vals, x=source)) + 
  geom_bar(position="fill", stat="identity")
