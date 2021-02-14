path <- "/group/irc/shared/synthetic_visium/generation/"
datasets <- c("brain_cortex_generation.rds", "cerebellum_cell_generation.rds",
              "cerebellum_nucleus_generation.rds", "hippocampus_generation.rds",
              "kidney_generation.rds", "pbmc_generation.rds", "scc_p5_generation.rds")
metadata <- lapply(paste0(path, datasets), function (k) colnames(readRDS(k)@meta.data))

dataset_types = c("artificial_uniform_distinct", "artificial_diverse_distinct", 
                  "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse",
                  "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse",
                  "artificial_regional_rare_celltype_diverse")

### Actual code ###
createSynthvisiumRDS <- function(inputscRNA_rds, dataset_type, output_path="", celltype_var = "celltype"){
  seurat_obj_scRNA =readRDS(inputscRNA_rds)
  
  # Create synthetic visium data from scRNA data
  synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj_scRNA, dataset_type = dataset_type, 
                                                    clust_var = "celltype", n_regions = 5,
                                                    n_spots_min = 100, n_spots_max = 200,
                                                    visium_mean = 20000, visium_sd = 5000)
  
  directory = dirname(inputscRNA_rds)
  inputscRNA_name = str_split(basename(inputscRNA_rds), "\\.")[[1]][1]
  dir.create((paste0(output_folder, inputscRNA_name)))
  output_folder <- paste0(output_folder, inputscRNA_name, "/")
  saveRDS(synthetic_visium_data, paste0(output_folder, inputscRNA_name, "_", dataset_type, "_synthvisium.rds"))
  print(paste0("Dataset saved at ", output_folder, inputscRNA_name, "_", dataset_type, "_synthvisium.rds"))
}

dataset_type = dataset_types[1]
dataset <- datasets[5]
output_path = "/home/chananchidas/data/"
createSynthvisiumRDS(paste0(path, dataset), dataset_type, output_path=output_path)