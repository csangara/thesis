setwd("D:/Work (Yr 2 Sem 1)/Thesis")

library(Seurat)
library(SeuratData)
library(synthvisium)

dataset_types_list <- list(
  original = c("real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
               "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
               "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
               "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium"),
  recommended = c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
)

possible_dataset_types <- dataset_types_list[["recommended"]]

#### COMPARING ARTIFICIAL SYNTHVISIUM TYPES TO REAL ####

brain <- LoadData("stxBrain", type = "posterior1")
ddsList <- list(realBrain = as.matrix(GetAssayData(brain)))

path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
repl <- "rep1"
dataset <- "brain_cortex_generation"
for (dataset_type in possible_dataset_types){
  temp_synth_data <- readRDS(paste0(path, dataset, "/", repl, "/", 
                                    dataset, "_", dataset_type, "_synthvisium.rds"))
  ddsList[[dataset_type]] <- as.matrix(temp_synth_data$counts)
}


countsimQC::countsimQCReport(ddsList = ddsList, 
                             outputFile = "brain_allArtificial.html",
                             outputDir = "thesis/rds/countsimQC",
                             outputFormat = "html_document", 
                             showCode = FALSE, forceOverwrite = TRUE,
                             savePlots = TRUE, description = "This is my test report.", 
                             maxNForCorr = 25, maxNForDisp = Inf, 
                             calculateStatistics = TRUE, subsampleSize = 25,
                             kfrac = 0.01, kmin = 5, ignorePandoc = TRUE,
                             permutationPvalues = FALSE, nPermutations = NULL)


#### COMPARING DIFFERENT METHODS OF DATA GENERATION ####

# Synthvisium - brain
path <- "/group/irc/shared/synthetic_visium/generation/"
brain_sc_10x <- readRDS(paste0(path, "brain_cortex_generation.rds"))
DefaultAssay(brain_sc_10x) <- "RNA"

path <- "D:/Work (Yr 2 Sem 1)/Thesis/rds/"
brain_sc_10x <- readRDS(paste0(path, "seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds"))


set.seed(10)
synthetic_visium_data = generate_synthetic_visium(seurat_obj = brain_sc_10x, dataset_type = "artificial_diverse_distinct",
                                                  clust_var = "subclass", n_regions = 5, n_spots_min = 50, n_spots_max = 500,
                                                  visium_mean = 30000, visium_sd = 8000)


# STEREOSCOPE
stereoscope_synth <- read.table("Data/synthetic_datasets/stereoscope_sc_10x_braincortex_robin_rawcounts_default/counts.st_synth.tsv",
                                header=TRUE, sep="\t", row.names=1)
# CELL2LOCATION
library(stringr)
# Get seeds
file_names <- list.files("Data/synthetic_datasets/cell2location_sc_10x_braincortex_robin_rawcounts_2000spots/")
seeds <- str_extract(gsub(".p", "", gsub(".*_", "", file_names)), regex("[0-9]+"))
seeds <- unique(seeds[!is.na(seeds)])

for (seed in seeds[1]){
  cell2location_synth <- read.table(paste0("data/synthetic_datasets/cell2location_sc_10x_braincortex_robin_rawcounts_2000spots/synthetic_ST_seed",
                                           seed, "_1_counts.csv"), header=TRUE, sep=",", row.names=1)
}

cell2location_synth <- cell2location_synth[rowSums(cell2location_synth)>0,]

## COUNTSIMQC ##
ddsList <- list(
  realBrain = as.matrix(DownsampleMatrix(GetAssayData(brain))),
  synthvisium = as.matrix(DownsampleMatrix(synthetic_visium_data$counts))
)

countsimQC::countsimQCReport(ddsList = ddsList, 
                             outputFile = "countsim_report_brain&synth.html",
                             outputDir = "thesis/", outputFormat = "html_document", 
                             showCode = FALSE, forceOverwrite = TRUE,
                             savePlots = TRUE, description = "This is my test report.", 
                             maxNForCorr = 25, maxNForDisp = Inf, 
                             calculateStatistics = TRUE, subsampleSize = 25,
                             kfrac = 0.01, kmin = 5, ignorePandoc = TRUE,
                             permutationPvalues = FALSE, nPermutations = NULL)


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
