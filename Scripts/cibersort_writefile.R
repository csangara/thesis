## REFERENCE ##
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
seurat_obj_scRNA =  readRDS("rds/allen_cortex_dwn.rds")

# Extract transformed counts
SCT_matrix <- as.matrix(GetAssayData(seurat_obj_scRNA)) # SCTransformed
colnames(SCT_matrix) <- seurat_obj_scRNA@meta.data$subclass

# Can also do raw counts
rawcounts <- as.matrix(GetAssayData(seurat_obj_scRNA[["RNA"]])) # Raw counts
rawcounts_CPM <- rawcounts/sum(rawcounts)*10^6
colnames(rawcounts_CPM) <- seurat_obj_scRNA@meta.data$subclass

# Write to CIBERSORT-compatible format
write.table("Gene", file = "Data/CIBERSORTx/cibersort_scRNAreference_SCT.txt", sep = "\t",eol="\t", append = FALSE, quote = FALSE, row.names=FALSE, col.names = FALSE)
write.table(SCT_matrix, file = "Data/CIBERSORTx/cibersort_scRNAreference_SCT.txt", sep = "\t", append = TRUE, quote = FALSE, row.names = TRUE, col.names = TRUE)

## MIXTURE (SPATIAL) ##
mixture = readRDS("rds/synthvisium_spatial/allencortexdwn_artificial_uniform_distinct_synthvisium.rds")
seurat_obj_mixture = CreateSeuratObject(counts = mixture$counts, min.cells = 2, min.features = 200, assay = "Spatial")
seurat_obj_mixture = SCTransform(seurat_obj_mixture, assay = "Spatial", verbose = FALSE)

SCT_matrix_mixture <- as.matrix(GetAssayData(seurat_obj_mixture))
write.table("GeneSymbol", file = "Data/CIBERSORTx/cibersort_mixture_ArtUniDistinct.txt", sep = "\t",eol="\t", append = FALSE, quote = FALSE, row.names=FALSE, col.names = FALSE)
write.table(SCT_matrix_mixture, file = "Data/CIBERSORTx/cibersort_mixture_ArtUniDistinct.txt", sep = "\t", append = TRUE, quote = FALSE, row.names = TRUE, col.names = TRUE)

## COMPARING RESULTS
results <- read.table("Data/CIBERSORTx/output/", sep="\t")
res = cor(t(mixture$relative_spot_composition[,1:23]), t(decon_mtrx[,1:23]))
mean(diag(res), na.rm=TRUE)