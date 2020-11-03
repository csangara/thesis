setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
mixture = readRDS("rds/synthvisium_spatial/allencortexdwn_artificial_uniform_distinct_synthvisium.rds")
seurat_obj_mixture = CreateSeuratObject(counts = mixture$counts, min.cells = 2, min.features = 200, assay = "Spatial")
seurat_obj_mixture = SCTransform(seurat_obj_mixture, assay = "Spatial", verbose = FALSE)


sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=names(Idents(seurat_obj_mixture)),
                       cellType=Idents(seurat_obj_mixture))
sc.pdata <- new("AnnotatedDataFrame",
                         data=sc.pheno)
sc.data <- as.matrix(GetAssayData(seurat_obj_mixture)[,names(Idents(seurat_obj_mixture)),drop=F])
sc.eset <- ExpressionSet(assayData=sc.data, phenoData=sc.pdata)
