library(Seurat)
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")

###### GENERATING SYNTHETIC SPOTS OF DIFFERENT METHODS ######
## SYNTHVISIUM ##
library(synthvisium)
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")

set.seed(10)
synthetic_visium_data = generate_synthetic_visium(seurat_obj = brain_sc_10x, dataset_type = "artificial_diverse_distinct",
                                                  clust_var = "subclass", n_regions = 5, n_spots_min = 600, n_spots_max = 700,
                                                  visium_mean = 20000, visium_sd = 6000)
saveRDS(synthetic_visium_data, "Data/synthetic_datasets/for_comparing_dists/synthvisium_sc_10x_braincortex_robin_rawcounts_3224spots.rds")

#kidney
kidney <- readRDS("Scripts/robin/raw_data/scRNAseq/seurat_obj_scrnaseq_kidney_filtered.rds")
kidney <- CreateSeuratObject(GetAssayData(kidney, slot="counts"), meta.data = data.frame(kidney$celltype))

synthetic_visium_data = generate_synthetic_visium(seurat_obj = kidney, dataset_type = "artificial_diverse_distinct",
                                                  clust_var = "kidney.celltype", n_regions = 3, n_spots_min = 400, n_spots_max = 500,
                                                  visium_mean = 35000, visium_sd = 8000)
saveRDS(synthetic_visium_data, "Data/synthetic_datasets/for_comparing_dists/synthvisium_sc_10x_kidney_robin_rawcounts_1314spots.rds")

## SPOTLIGHT ##
source("Scripts/synthetic_data_generation/spotlight_test_spot_fun.R")
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")
spotlight_synth <- test_spot_fun(brain_sc_10x, clust_vr="subclass", n=3500)
saveRDS(spotlight_synth, "Data/synthetic_datasets/for_comparing_dists/spotlight_sc_10x_braincortex_robin_rawcounts_3500spots.rds")

kidney <- readRDS("Scripts/robin/raw_data/scRNAseq/seurat_obj_scrnaseq_kidney_filtered.rds")
spotlight_synth <- test_spot_fun(kidney, clust_vr="celltype", n=1500)
saveRDS(spotlight_synth, "Data/synthetic_datasets/for_comparing_dists/spotlight_sc_10x_kidney_robin_rawcounts_1500spots.rds")

## STEREOSCOPE ##
# need two .tsv files, also at least 30 cells required per cell type
selected_celltypes <- names(table(brain_sc_10x$subclass)[table(brain_sc_10x$subclass) >= 30])
Idents(brain_sc_10x) <- brain_sc_10x$subclass
brain_sc_10x <- subset(brain_sc_10x, idents = selected_celltypes) # Select cell types >= 30 cells

# Counts file, need cell x gene matrix
write.table(t(as.matrix(brain_sc_10x[["RNA"]]@counts)),
            file="rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts_counts.tsv",
            quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

# Metadata file
metadata <- data.frame("cell" = colnames(brain_sc_10x),
                       "bio_celltype" = brain_sc_10x$subclass)
write.table(metadata, file="rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts_metadata.tsv",
            quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

# kidney
kidney <- readRDS("Scripts/robin/raw_data/scRNAseq/seurat_obj_scrnaseq_kidney_filtered.rds")
selected_celltypes <- names(table(kidney$celltype)[table(kidney$celltype) >= 30])
Idents(kidney) <- kidney$celltype
kidney <- subset(kidney, idents = selected_celltypes) # Select cell types >= 30 cells
write.table(t(as.matrix(kidney[["RNA"]]@counts)),
            file="rds/seurat_obj_scrnaseq_kidney_filtered_rawcounts_counts.tsv",
            quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
metadata <- data.frame("cell" = colnames(kidney),
                       "bio_celltype" = kidney$celltype)
write.table(metadata, file="rds/seurat_obj_scrnaseq_kidney_filtered_rawcounts_metadata.tsv",
            quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

# Run the make_st_set.py file using these two files!
# Read expression data
expression_data <- read.table("Data/synthetic_datasets/for_comparing_dists/stereoscope_sc_10x_braincortex_robin_rawcounts_default/counts.st_synth.tsv",
                              header=TRUE, sep="\t", row.names=1)

## CELL2LOCATION ##
# need a .h5ad file and an annotation (.csv) file
# Creating h5ad file
# Convert to loom/h5ad (h5ad doesn't work as of 09/02/2021)
source("Scripts/helperFunctions.R")
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")
convertSeuratRDSToh5ad("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds") # Doesn't work, fails at Convert() 
convertSeuratRDSToLoom("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")

# Metadata file
metadata <- data.frame("cell" = colnames(brain_sc_10x),
                       "bio_celltype" = brain_sc_10x$subclass)
write.table(metadata, file="rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts_metadata.csv",
            quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)

# kidney
kidney <- readRDS("Scripts/robin/raw_data/scRNAseq/seurat_obj_scrnaseq_kidney_filtered.rds")
metadata <- data.frame("cell" = colnames(kidney),
                       "bio_celltype" = kidney$celltype)
write.table(metadata, file="rds/seurat_obj_scrnaseq_kidney_filtered_rawcounts_metadata.csv",
            quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)


# Run the pipeline (create_synthetic_data_cell2location.sh) using those two files
# Read in results as follows
library(stringr)
# Get seeds
file_names <- list.files("Data/synthetic_datasets/for_comparing_dists/cell2location_sc_10x_braincortex_robin_rawcounts/")
seeds <- str_extract(gsub(".p", "", gsub(".*_", "", file_names)), regex("[0-9]+"))
seeds <- unique(seeds[!is.na(seeds)])
expression_data <- read.table(paste0("Data/synthetic_datasets/for_comparing_dists/cell2location_sc_10x_braincortex_robin_rawcounts/2000spots/synthetic_ST_seed",
                                       seeds[1], "_1_counts.csv"), header=TRUE, sep=",", row.names=1)
  