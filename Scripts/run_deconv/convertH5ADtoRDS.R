#!/usr/bin/env Rscript
inputs = commandArgs(trailingOnly=TRUE)
if (length(inputs) < 1) {
  stop("usage: < input file > < output path >", call.=FALSE)
}

library(SeuratDisk)
library(Seurat)

input_path <- inputs[1]
Convert(input_path, dest = "h5seurat", overwrite = TRUE)
file_name_with_dir <- tools::file_path_sans_ext(input_path)
seurat_obj <- LoadH5Seurat(paste0(file_name_with_dir, ".h5seurat"))

output_path <- ifelse(length(inputs)==2, inputs[2], dirname(input_path))
file_name <- stringr::str_split(basename(input_path), "\\.")[[1]][1]
saveRDS(seurat_obj, paste0(output_path, "/", file_name, ".rds"))
