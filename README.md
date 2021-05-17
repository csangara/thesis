# Cell Type Deconvolution in Spatial Transcriptomics

In this repository, you can find scripts to evaluate five deconvolution methods: cell2location, MuSiC, stereoscope, RCTD, and SPOTlight. I used seven scRNA-seq datasets to generate synthetic spatial data using the package **synthvisium** (not yet publicly available). The raw datasets along with the download links are listed below.

|      Dataset     |                                                                     Direct download link                                                                    |
|:----------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------:|
|   Brain cortex   |                                           [Link](https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1)                                           |
| Cerebellum (sc)  | [Link](https://singlecell.broadinstitute.org/single_cell/study/SCP948/robust-decomposition-of-cell-type-mixtures-in-spatial-transcriptomics#study-download) |
| Cerebellum (sn)  |                                                                                                                                                             |
|    Hippocampus   |                                        [Link](https://storage.googleapis.com/linnarsson-lab-loom/l1_hippocampus.loom)                                       |
|      Kidney      |                [Link](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107585/suppl/GSE107585_Mouse_kidney_single_cell_datamatrix.txt.gz)               |
|       PBMC       |                                [Link](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)                               |
| SCC (patient 5)  |                            [Link](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144236/suppl/GSE144236_cSCC_counts.txt.gz)                           |

(Both cerebellum datasets can be downloaded from the link.)

You can also follow this [Seurat vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) for a quick preprocessing pipeline of the PBMC data. Instead of synthvisium you can generate synthetic data using scripts from SPOTlight, stereoscope, or cell2location (recommended) as well. You can find the scripts for the first two in `Scripts/synthetic_data_generation` and for cell2location at [this repository](https://github.com/emdann/ST_simulation/). There is a tutorial in that folder on how to use those functions.
