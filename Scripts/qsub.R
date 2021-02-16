library(qsub)
config <- create_qsub_config(
  remote = "chananchidas@prism.psb.ugent.be:7777", 
  local_tmp_path = "D:/Work (Yr 2 Sem 1)/Thesis/.r2gridengine",
  remote_tmp_path = "/scratch/irc/personal/chananchidas/.r2gridengine")

set_default_qsub_config(config)

test <- qsub_lapply(seq_len(10), function(x) x + 1)

# Testing the generation of synthetic data
createSynthvisiumRDS <- function(inputscRNA_rds, dataset_type, output_folder="synthvisium_spatial/"){
  seurat_obj_scRNA =  readRDS(inputscRNA_rds)
  
  # Create synthetic visium data from scRNA data
  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj_scRNA, dataset_type = dataset_type, 
                                                    clust_var = "subclass", region_var = "brain_subregion" , n_regions = NULL,
                                                    n_spots_min = 50, n_spots_max = 200, visium_mean = 20000, visium_sd = 5000)
  
  directory = dirname(inputscRNA_rds)
  inputscRNA_name = str_split(basename(inputscRNA_rds), "\\.")[[1]][1]
  
  saveRDS(synthetic_visium_data, paste0(directory, "/", output_folder, inputscRNA_name, "_", dataset_type, "_synthvisium.rds"))
  print(paste0("Dataset saved at ", directory, "/", output_folder, inputscRNA_name, "_", dataset_type, "_synthvisium.rds"))
}

library(synthvisium)
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
allen_cortex_dwn = readRDS("rds/allen_cortex_dwn_original.rds")
test <- qsub_lapply(list(allen_cortex_dwn), generate_synthetic_visium, dataset_type = "real", clust_var="subclass", n_regions=NULL,
               region_var="brain_subregion")

# Getting metadata of datasets on PRISM
path <- "/group/irc/shared/synthetic_visium/generation/"
datasets <- c("brain_cortex_generation.rds", "cerebellum_cell_generation.rds",
              "cerebellum_nucleus_generation.rds", "hippocampus_generation.rds",
              "kidney_generation.rds", "pbmc_generation.rds", "scc_p5_generation.rds")

# Works locally
metadata <- qsub_lapply(paste0(path, datasets), function (k) colnames(readRDS(k)@meta.data))
metadata <- readRDS("Scripts/server/dataset_metadata.rds")
