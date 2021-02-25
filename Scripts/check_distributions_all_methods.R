library(Seurat)
library(fitdistrplus)
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")

##### LOAD DIFFERENT DATASETS #####
# 10x Brain
library(SeuratData)
sections = c("posterior1", "posterior2", "anterior1", "anterior2")

for (i in 1:4){
  section = sections[i]
  brain <- LoadData("stxBrain", type = section)
}

# 10x Kidney
kidney <- LoadData("stxKidney")

# Synthvisium - brain
library(synthvisium)
brain_sc_10x <- readRDS("rds/seurat_obj_scrnaseq_cortex_filtered_rawcounts.rds")

set.seed(10)
synthetic_visium_data = generate_synthetic_visium(seurat_obj = brain_sc_10x, dataset_type = "artificial_diverse_distinct",
                                                  clust_var = "subclass", n_regions = 5, n_spots_min = 50, n_spots_max = 500,
                                                  visium_mean = 30000, visium_sd = 8000)

# SPOTLIGHT
spotlight_synth <- readRDS("Data/synthetic_datasets/spotlight_sc_10x_braincortex_robin_rawcounts_2500spots.rds")

# STEREOSCOPE
stereoscope_synth <- read.table("Data/synthetic_datasets/stereoscope_sc_10x_braincortex_robin_rawcounts_2000spots_10000genes/counts.st_synth.tsv",
                              header=TRUE, sep="\t", row.names=1)
# CELL2LOCATION 
library(stringr)
# Get seeds
file_names <- list.files("Data/synthetic_datasets/cell2location_sc_10x_braincortex_robin_rawcounts_2000spots/")
seeds <- str_extract(gsub(".p", "", gsub(".*_", "", file_names)), regex("[0-9]+"))
seeds <- unique(seeds[!is.na(seeds)])

for (seed in seeds[1]){
  cell2location_synth <- read.table(paste0("Data/synthetic_datasets/cell2location_sc_10x_braincortex_robin_rawcounts_2000spots/synthetic_ST_seed",
                                       seed, "_1_counts.csv"), header=TRUE, sep=",", row.names=1)
}

## COUNTSIMQC ##
ddsList <- list(
  realBrain = as.matrix(DownsampleMatrix(GetAssayData(brain))),
  synthvisium = as.matrix(DownsampleMatrix(synthetic_visium_data$counts)),
  stereoscope = as.matrix(DownsampleMatrix(t(stereoscope_synth))),
  cell2location = as.matrix(DownsampleMatrix(t(cell2location_synth)))
)


countsimQC::countsimQCReport(ddsList = ddsList, 
                 outputFile = "countsim_report_dwn.html",
                 outputDir = "D:/Work (Yr 2 Sem 1)/Thesis/plots/", outputFormat = "html_document", 
                 showCode = FALSE, forceOverwrite = TRUE,
                 savePlots = TRUE, description = "This is my test report.", 
                 maxNForCorr = 25, maxNForDisp = Inf, 
                 calculateStatistics = TRUE, subsampleSize = 25,
                 kfrac = 0.01, kmin = 5, permutationPvalues = FALSE, nPermutations = NULL)


DownsampleMatrix <- function(mat){
  while (nrow(mat) > 2000) {
    rowsamp <- sample(1:nrow(mat), nrow(mat) %/% 2)
    mat <- mat[rowsamp, ]
  }
  return(mat)
}
