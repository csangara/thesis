# Cell Type Deconvolution in Spatial Transcriptomics

In this repository, you can find the analysis scripts and plots pertaining to the dissertation. There are scripts to run and evaluate the five deconvolution methods (cell2location, MuSiC, stereoscope, RCTD, and SPOTlight). Later, I also apply cell2location and RCTD on real data.

## Benchmarking

To perform benchmarking, synthetic data has to be created from a reference scRNA-seq dataset. Then, we run different deconvolution methods on the datasets and evaluate them. I made use of seven scRNA-seq datasets to generate synthetic spatial data using the package **synthvisium** (not yet publicly available). The raw datasets along with the download links are listed below.

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

I did not preprocess the scRNA-seq data myself so I cannot share the scripts here, but the procedure is described in section 5.1 of the text. You can also follow this [Seurat vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) for a quick preprocessing of the PBMC data.

### Synthetic data generation
As an alternative to synthvisium, you can generate synthetic data using scripts from SPOTlight, stereoscope, or cell2location. Some sample code for running these functions can be found at `Scripts/synthetic_data_generation`, although the cell2location functions have to be cloned from [here](https://github.com/emdann/ST_simulation).

The countsimQC reports between different synthetic data generation algorithms can be found in the folder `countsimQC/`.

### Running deconvolution methods
Scripts for running the deconvolution methods can be found at `Scripts/run_deconv` along with a description for using those files. The deconvolution results are compiled in the folder `results/`.

For scripts to generate downsampled data and get the runtime of each method, check out `Scripts/run_deconv_downsample`.

### Evaluation
Evaluation scripts are found at `Scripts/` with the prefix `evaluation_`. These make use of the deconvolution results saved in `results/`. The folder structure of `results/` is dataset → replicate → method output. Within each replicate folder (`rep` prefix), you will find the `all_metrics*.rds` file which has the computed metrics (RMSE and six classification metrics) for all 8 dataset types. In addition the method outputs, there is also the `plots/` subfolder that contains plots of the predictions on a UMAP for each cell type. The `corr_distribution/` subfolder contains density plots of the correlation across all spots. These plots were not shown in the dissertation.

### Plots
Along with high-resolution of the plots found in the thesis, you can also find scripts that are used to generate the plots. The plots are in the directory `plots/` and there I try to make a link with the corresponding scripts.

## Application on real data
You can find the scripts for  preprocessing, running deconvolution tools, and evaluating the liver dataset in the folder `Scripts/liver/`.
