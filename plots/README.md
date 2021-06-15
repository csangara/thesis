# Plots description

| Folder/File                     | Description                                                                              | Figure in text         | Script to create plot (`/Scripts/`)                |
|---------------------------------|------------------------------------------------------------------------------------------|------------------------|----------------------------------------------------|
| UMAP_generation/ and UMAP_test/ | UMAP plots of the seven reference datasets                                               | A2                     | misc_processing.R                                                |
| UMAP_synthvisium/               | UMAP plots of the synthetic dataset from the single-cell cerebellum dataset.             | 3-3                    | misc_processing.R                                  |
| data_distribution_fitdistr/     | Distributions of synthetic datasets by `fitdistr` instead of `countsimQC`                | N/A                    | check_distributions_fitdistr.R                     |
| downsampling/                   | Runtimes of scalability tests                                                            | 3-9                    | run_deconv_downsample/ downsample_plot_heatmaps.py |
| liver/                          | Plots pertaining to the liver dataset                                                    | 3-12 to 15,  A13-A17 | liver/*.pynb, liver/evaluate_liver_data.R          |
| metrics_facets/                 | Faceted plots of all metrics                                                             | 3-5, A6-A12            | evaluation.R                                       |
| nCount_RNA/                     | Density and violin plots of count and gene distributions of the seven reference datasets | 3-1, A1                | plot_nCount.R                                      |
| SD_between_reps_runs.png        | Stochasticity between different runs and replicates                                      | 3-8                    | evaluation_old_3reps3runs.R                        |
| classification_barplot.png      | Relative frequency bar chart of best-performing methods for classification metrics       | 3-6                    | evaluation.R                                       |
| scenario1.png                   | RMSE between using the same vs different sets of cells                                   | 3-4                    | evaluation_s1.R                                    |
| scenario3_*.png                 | RMSE between using the same vs different reference datasets                              | 3-7                    | evaluation_s3.R                                    |
