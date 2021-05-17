library(Seurat)
library(SeuratData)
# library(synthvisium)

DownsampleMatrix <- function(mat){
  while (nrow(mat) > 1500) {
    rowsamp <- sample(1:nrow(mat), nrow(mat) %/% 2)
    mat <- mat[rowsamp, ]
  }
  
  while (ncol(mat) > 500) {
    colsamp <- sample(1:ncol(mat), ncol(mat) %/% 2)
    mat <- mat[,colsamp]
  }
  return(mat)
}

#### PART 1: COMPARING DIFFERENT METHODS ####
## BRAIN ##
# brain <- LoadData("stxBrain.SeuratData", type = "posterior1")
# 
# # Synthvisium - brain
# synthvisium_synth <- readRDS("data/synthetic_datasets/synthvisium/brain_3224spots.rds")
# 
# # SPOTLIGHT
# spotlight_synth <- readRDS("data/synthetic_datasets/spotlight/braincortex/3500spots.rds")
# 
# # STEREOSCOPE
# stereoscope_synth <- read.table("data/synthetic_datasets/stereoscope/braincortex/3500spots/counts.st_synth.tsv",
#                               header=TRUE, sep="\t", row.names=1)
# 
# # CELL2LOCATION
# cell2location_synth <- read.table(paste0("data/synthetic_datasets/cell2location/braincortex/3500spots/synthetic_ST_seed307_3500_counts.csv"),
#                                   header=TRUE, sep=",", row.names=1)
# cell2location_synth <- cell2location_synth[rowSums(cell2location_synth)>0,] # Remove spots with zero counts
# 
# ddsList <- list(
#   realBrain = as.matrix(GetAssayData(brain)),
#   spotlight_synth = as.matrix(spotlight_synth$topic_profiles),
#   synthvisium = as.matrix(synthvisium_synth$counts),
#   stereoscope = as.matrix(t(stereoscope_synth)),
#   cell2location = as.matrix(t(cell2location_synth))
# )
# 
# outputFileName <- "four_methods_brain_withstats.html"
# description <- "comparing all four methods of generating synth data"

## KIDNEY ##
kidney <- LoadData("stxKidney.SeuratData")

# Synthvisium - brain
synthvisium_synth <- readRDS("data/synthetic_datasets/synthvisium/kidney_1314spots.rds")
kidney <- subset(kidney, features = rownames(synthvisium_synth$counts))

# SPOTLIGHT
spotlight_synth <- readRDS("data/synthetic_datasets/spotlight/kidney/1500spots.rds")

# STEREOSCOPE
stereoscope_synth <- read.table("data/synthetic_datasets/stereoscope/kidney/1500spots/counts.st_synth.tsv",
                                header=TRUE, sep="\t", row.names=1)

# CELL2LOCATION
cell2location_synth <- read.table(paste0("data/synthetic_datasets/cell2location/kidney/1500spots/synthetic_ST_seed371_1500_counts.csv"),
                                  header=TRUE, sep=",", row.names=1)
cell2location_synth <- cell2location_synth[rowSums(cell2location_synth)>0,] # Remove spots with zero counts

ddsList <- list(
  realKidney = as.matrix(GetAssayData(kidney)),
  spotlight_synth = as.matrix(spotlight_synth$topic_profiles),
  synthvisium = as.matrix(synthvisium_synth$counts),
  stereoscope = as.matrix(t(stereoscope_synth)),
  cell2location = as.matrix(t(cell2location_synth))
)

outputFileName <- "four_methods_kidney_oridwn_withstats.html"
description <- "comparing all four methods of generating synth data, original data is downsampled"


#### PART 2: COMPARING DIFFERENT GENERATION SCHEMES ####
# 
# dataset_types_list <- list(
#   original = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
#                "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
#                "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
#                "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium"),
#   recommended = c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
#                   "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
#                   "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
# )
# possible_dataset_types <- dataset_types_list[["recommended"]]
# 
# # brain <- LoadData("stxBrain.SeuratData", type = "posterior1")
# # ddsList <- list(realBrain = as.matrix(GetAssayData(brain)))
# 
# kidney <- LoadData("stxKidney.SeuratData")
# ddsList <- list(realKidney = as.matrix(GetAssayData(kidney)))
# 
# path <- "~/data/"
# repl <- "rep1"
# dataset <- "kidney_generation"
# for (dataset_type in possible_dataset_types){
#   temp_synth_data <- readRDS(paste0(path, dataset, "/", repl, "/",
#                                     dataset, "_", dataset_type, "_synthvisium.rds"))
#   ddsList[[stringr::str_remove(dataset_type, "artificial_")]] <- as.matrix(temp_synth_data$counts)
# }
# 
# outputFileName <- "kidney_allArtificial_withstats.html"
# description <- "comparing real kidney with all artificial types"

#### PART 3: COMPARING DIFFERENT SPATIAL DATASETS ####

# ddsList <- list()
# 
# for (section in c("posterior1", "posterior2", "anterior1", "anterior2")){
#   brain <- LoadData("stxBrain.SeuratData", type = section)
#   ddsList[[paste0("brain_", section)]] <- as.matrix(GetAssayData(brain))
# }
# 
# # 10x Kidney
# kidney <- LoadData("stxKidney.SeuratData")
# ddsList[["kidney"]] <- as.matrix(GetAssayData(kidney))
# 
# outputFileName <- "stxDatasets_withstats.html"
# description <- "comparing brain and kidney datasets"

#### RUNNING COUNTSIMQC ####
print(Sys.getenv("RSTUDIO_PANDOC"))
print(outputFileName)
countsimQC::countsimQCReport(ddsList = ddsList, 
                             outputFile = outputFileName,
                             outputDir = "countsimQC_results/",
                             description = description,
                             forceOverwrite = TRUE,
                             calculateStatistics = TRUE)