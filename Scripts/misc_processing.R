### MISCELLANEOUS (PRE)PROCESSING SCRIPTS ###

setwd("D:/Work (Yr 2 Sem 1)/Thesis/Scripts")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/"
source("helperFunctions.R")

possible_dataset_types = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium")

### DOWNSAMPLING SPOTS FOR CIBERSORT ###

input_path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/CIBERSORTx/cibersort_mixture_ArtUniDistinct.txt"
reduceSpotsCS(input_path, 75)

### CONVERTING SEURAT OBJECT TO LOOM ###

for (dataset_type in possible_dataset_types){
  print(dataset_type)
  convertSeuratRDSToLoom(paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type,
                                "_synthvisium.rds"), TRUE)
}

## SAVING NON-NORMALIZED DATA AS H5AD ##
# Reference file
convertSeuratRDSToh5ad("D:/Work (Yr 2 Sem 1)/Thesis/allen_cortex_dwn_original.rds")

# Synthvisium files
for (dataset_type in possible_dataset_types[10:13]){
  print(dataset_type)
  input_path <- paste0(path, "rds/synthvisium_spatial/allen_cortex_dwn_", dataset_type, "_synthvisium.rds")
  convertSeuratRDSToh5ad(input_path, createSeuratFromRDS = TRUE, PP=FALSE)
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

