import argparse as arp
import sys
import os
import matplotlib as mpl
mpl.use('Agg')

def main():
    ##### PARSING COMMAND LINE ARGUMENTS #####
    prs = arp.ArgumentParser()
    
    prs.add_argument('sc_data_path',
                     type = str, help = 'path to single cell h5ad count data')
    
    prs.add_argument('cuda_device', type = str, help = "index of cuda device ID (from 0-7) or cpu")

    prs.add_argument('-o','--out_dir', default = None,
                     type = str, help = 'directory for regression model')

    prs.add_argument('-a','--annotation_column', default = 'celltype',
                 type = str, help = 'column name for covariate')
    
    args = prs.parse_args()
    cuda_device = args.cuda_device

    assert (cuda_device in ["0", "1", "2", "3", "4", "5", "6", "7"] or cuda_device == "cpu"), "invalid device id"
    
    ## scRNA reference (raw counts)
    sc_data_path = args.sc_data_path
    import anndata
    adata_scrna_raw = anndata.read_h5ad(sc_data_path)
    
    if args.out_dir is None:
        results_folder = os.path.dirname(sc_data_path) + '/cell2location_results/'
    else:
        results_folder = args.out_dir

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
        
    celltype = args.annotation_column
    
    ##### MAIN CODE #####
    if cuda_device.isdigit():
        os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

    import sys
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import cell2location
    import matplotlib as mpl
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    import seaborn as sns

    # silence scanpy that prints a lot of warnings
    import warnings
    warnings.filterwarnings('ignore')

    # Reduce the number of genes by discarding lowly expressed genes
    # remove cells and genes with 0 counts everywhere
    sc.pp.filter_cells(adata_scrna_raw, min_genes=1)
    sc.pp.filter_genes(adata_scrna_raw, min_cells=1)

    # calculate the mean of each gene across non-zero cells
    adata_scrna_raw.var['n_cells'] = (adata_scrna_raw.X.toarray() > 0).sum(0)
    adata_scrna_raw.var['nonz_mean'] = adata_scrna_raw.X.toarray().sum(0) / adata_scrna_raw.var['n_cells']
    nonz_mean_cutoff = 0.05
    cell_count_cutoff = np.log10(adata_scrna_raw.shape[0] * 0.0005)
    cell_count_cutoff2 = np.log10(adata_scrna_raw.shape[0] * 0.03)

    # select genes based on mean expression in non-zero cells
    adata_scrna_raw = adata_scrna_raw[:,(np.array(np.log10(adata_scrna_raw.var['nonz_mean']) > nonz_mean_cutoff)
             | np.array(np.log10(adata_scrna_raw.var['n_cells']) > cell_count_cutoff2))
          & np.array(np.log10(adata_scrna_raw.var['n_cells']) > cell_count_cutoff)]
    # add count matrix as raw
    adata_scrna_raw.raw = adata_scrna_raw

    # Run the pipeline
    from cell2location import run_regression
    r, adata_scrna_raw = run_regression(adata_scrna_raw, # input data object]

                       verbose=True, return_all=True,

                       train_args={
                        'covariate_col_names': [celltype], # column listing cell type annotation
                        #'sample_name_col': 'donor', # column listing sample ID for each cell

                        # column listing technology, e.g. 3' vs 5', 
                        # when integrating multiple single cell technologies corresponding 
                        # model is automatically selected
                        'tech_name_col': None, 

                        'stratify_cv': celltype, # stratify cross-validation by cell type annotation

                        'n_epochs': 100, 'minibatch_size': 1024, 'learning_rate': 0.01,

                        'use_cuda': cuda_device.isdigit(), # use GPU?

                        'train_proportion': 0.9, # proportion of cells in the training set (for cross-validation)
                        'l2_weight': True,  # uses defaults for the model

                        'use_raw': True},

                       model_kwargs={}, # keep defaults
                       posterior_args={}, # keep defaults

                       export_args={'path': results_folder + 'regression_model/', # where to save results
                                    'save_model': True, # save pytorch model?
                                    'run_name_suffix': ''})


if __name__ == '__main__':
    main()