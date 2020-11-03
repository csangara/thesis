library(MuSiC)
library(xbioc)
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path = "Data/MuSiC/"

#############################
# Bulk Tissue Cell Type Estimation
# Bulk data human pancreas
GSE50244.bulk.eset = readRDS(paste0(path, 'GSE50244bulkeset.rds'))
GSE50244.bulk.eset

# scRNA-seq data
EMTAB.eset = readRDS(paste0(path, 'EMTABesethealthy.rds'))
EMTAB.eset

# DECONVOLUTION (using GSE50244 bulk and EMTAB reference)
# Only 6 cell types
Est.prop.GSE50244 = music_prop(bulk.eset = GSE50244.bulk.eset, sc.eset = EMTAB.eset, clusters = 'cellType',
                               samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma',
                                                                   'acinar', 'ductal'), verbose = F)
names(Est.prop.GSE50244)

##############################
# Estimation of cell type proportions with pre-grouping of cell types
# Bulk data
Mouse.bulk.eset = readRDS(paste0(path, 'Mousebulkeset.rds'))
Mouse.bulk.eset

# ScRNA-seq reference
Mousesub.eset = readRDS(paste0(path, 'Mousesubeset.rds'))
Mousesub.eset

# DECONVOLUTION
# Clustering single cell data
Mousesub.basis = music_basis(Mousesub.eset, clusters = 'cellType', samples = 'sampleID', 
                             select.ct = c('Endo', 'Podo', 'PT', 'LOH', 'DCT', 'CD-PC', 'CD-IC', 'Fib',
                                           'Macro', 'Neutro','B lymph', 'T lymph', 'NK'))

# Manually specify the cluster and annotated single cell data with cluster information.
clusters.type = list(C1 = 'Neutro', C2 = 'Podo', C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'), C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))

cl.type = as.character(Mousesub.eset$cellType)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
pData(Mousesub.eset)$clusterType = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))

# 13 selected cell types
s.mouse = unlist(clusters.type)
s.mouse

load(paste0(path, 'IEmarkers.RData'))
Est.mouse.bulk = music_prop.cluster(bulk.eset = Mouse.bulk.eset, sc.eset = Mousesub.eset,
                                    group.markers = list(C1=NULL, C2=NULL, C3=Epith.marker, C4=Immune.marker),
                                    clusters = 'cellType', group = 'clusterType',
                                    samples = 'sampleID', clusters.type = clusters.type)

##############################
# Benchmark evaluation from pseudobulk
XinT2D.eset = readRDS(paste0(path, 'XinT2Deset.rds'))
XinT2D.eset

XinT2D.construct.full = bulk_construct(XinT2D.eset, clusters = 'cellType', samples = 'SubjectName')
names(XinT2D.construct.full)
XinT2D.construct.full$Bulk.counts
head(XinT2D.construct.full$num.real)

# calculate cell type proportions
XinT2D.construct.full$prop.real = relative.ab(XinT2D.construct.full$num.real, by.col = FALSE)
head(XinT2D.construct.full$prop.real)

# Estimate cell type proportions of artificial bulk data
Est.prop.Xin = music_prop(bulk.eset = XinT2D.construct.full$Bulk.counts, sc.eset = EMTAB.eset,
                          clusters = 'cellType', samples = 'sampleID', 
                          select.ct = c('alpha', 'beta', 'delta', 'gamma'))

Eval_multi(prop.real = data.matrix(XinT2D.construct.full$prop.real), 
           prop.est = list(data.matrix(Est.prop.Xin$Est.prop.weighted), 
                           data.matrix(Est.prop.Xin$Est.prop.allgene)), 
           method.name = c('MuSiC', 'NNLS'))
