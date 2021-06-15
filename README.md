# Cell Type Deconvolution in Spatial Transcriptomics

In this repository, you can find scripts to run and evaluate five deconvolution methods: cell2location, MuSiC, stereoscope, RCTD, and SPOTlight. I used seven scRNA-seq datasets to generate synthetic spatial data using the package **synthvisium** (not yet publicly available). The raw datasets along with the download links are listed below.

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

You can also follow this [Seurat vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) for a quick preprocessing pipeline of the PBMC data.

## Synthetic data generation
Instead of synthvisium you can generate synthetic data using scripts from SPOTlight, stereoscope, or cell2location as well. Some sample code for running these functions can be found at `Scripts/synthetic_data_generation`, although the cell2location functions have to be cloned from [here](https://github.com/emdann/ST_simulation).

## Running deconvolution methods
Scripts for running the deconvolution methods can be found at `Scripts/run_deconv` along with a description for using those files.

For scripts to generate downsampled data and get the runtime of each method, check out `Scripts/run_deconv_downsample`.

## Evaluation
Evaluation scripts are found at `Scripts/` with the prefix `evaluation_`.

## Plots
Along with a more high-resolution of the plots in the thesis, you can also find scripts that are used to generate the plots. The plots are in the directory `plots/` and there I try to make a link with the corresponding scripts.
