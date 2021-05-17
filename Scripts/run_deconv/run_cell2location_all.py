import argparse as arp
import sys
import os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
mpl.use('Agg')

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

def main():
    ##### PARSING COMMAND LINE ARGUMENTS #####
    prs = arp.ArgumentParser()
    
    prs.add_argument('sc_data_path',
                     type = str, help = 'path to single cell h5ad count data')
    
    prs.add_argument('sp_data_prefix', type = str, help = 'path and file prefix of spatial data')
    # e.g., /srv/scratch/chananchidas/data/brain_cortex_generation/brain_cortex_generation_
    
    prs.add_argument('cuda_device', type = str, help = "index of cuda device ID, from 0-7")
    
    prs.add_argument('dataset_type', nargs='+',
                        help = 'multiple arguments of dataset types to be deconvolved, or "all" or "recommended"')

    prs.add_argument('-o','--out_dir', default = None,
                     type = str, help = 'directory for regression model')

    prs.add_argument('-a','--annotation_column', default = 'celltype',
                 type = str, help = 'column name for covariate')
    
    possible_dataset_types = ["real", "real_top1","real_top1_uniform","real_top2_overlap",
    "real_top2_overlap_uniform", "real_missing_celltypes_visium", "artificial_uniform_distinct", 
    "artificial_diverse_distinct",  "artificial_uniform_overlap", "artificial_diverse_overlap",
    "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
    "artificial_missing_celltypes_visium", "artificial_dominant_rare_celltype_diverse",
    "artificial_regional_rare_celltype_diverse", "artificial_diverse_distinct_missing_celltype_sc",
    "artificial_diverse_overlap_missing_celltype_sc"]

    recommended_dataset_types = possible_dataset_types[6:12] + possible_dataset_types[13:15]
    
    # Get command line arguments
    args = prs.parse_args()
    sc_data_path = args.sc_data_path
    cuda_device = args.cuda_device
    sp_data_prefix = args.sp_data_prefix
    dataset_list = args.dataset_type
    celltype = args.annotation_column
    
    # Create result directory
    if args.out_dir is None:
        results_folder = os.path.dirname(sc_data_path) + '/cell2location_results/'
    else:
        results_folder = args.out_dir

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)    
    
    assert os.path.isfile(sc_data_path), "invalid input sc file path"
    assert cuda_device in ["0", "1", "2", "3", "4", "5", "6", "7"], "invalid device id"
    
    if dataset_list[0] == "all":
        dataset_list = possible_dataset_types
    elif dataset_list[0] == "recommended":
        dataset_list = recommended_dataset_types    
    assert all(dataset_type in possible_dataset_types for dataset_type in dataset_list), "invalid dataset type"
    assert all(os.path.isfile(sp_data_prefix + dataset_type + "_synthvisium.h5ad") \
                              for dataset_type in dataset_list), "invalid file path"
    
    os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device
    os.environ["CPATH"]="/usr/local/cuda/include:$CPATH" #To use cuDNN
    data_type = 'float32'
    os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=' + data_type + ',force_device=True'  
    
    run_cell2location_regression(sc_data_path, results_folder, celltype)
    regression_model_output = os.listdir(results_folder + "/regression_model")[0]
    run_cell2location(sp_data_prefix, dataset_list, results_folder, regression_model_output,
                     celltype)
    
def run_cell2location_regression(sc_data_path, results_folder, celltype):

    import cell2location
    
    ## scRNA reference (raw counts)
    adata_scrna_raw = anndata.read_h5ad(sc_data_path)
    
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

                        'use_cuda': True, # use GPU?

                        'train_proportion': 0.9, # proportion of cells in the training set (for cross-validation)
                        'l2_weight': True,  # uses defaults for the model

                        'use_raw': True},

                       model_kwargs={}, # keep defaults
                       posterior_args={}, # keep defaults

                       export_args={'path': results_folder + 'regression_model/', # where to save results
                                    'save_model': True, # save pytorch model?
                                    'run_name_suffix': ''})


def run_cell2location(sp_data_prefix, dataset_list, results_folder, regression_model_output,
                     covariate_col_names):
    

    import cell2location
    
    reg_path = f'{results_folder}regression_model/{regression_model_output}/'

    if not os.path.exists(results_folder + "std_model/"):
        os.makedirs(results_folder + "std_model/")

    ## READ IN SPATIAL DATA ##
    for dataset_type in dataset_list:
        print("\nDataset type: " + dataset_type)
        adata = sc.read_h5ad(sp_data_prefix + dataset_type + "_synthvisium.h5ad")
        print("Read in file from " + sp_data_prefix + dataset_type + "_synthvisium.h5ad")
        adata.obs['sample'] = dataset_type
        adata.var['SYMBOL'] = adata.var_names

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        adata.var['mt'] = [gene.startswith('mt') for gene in adata.var['SYMBOL']]
        adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

        # mitochondria-encoded (MT) genes should be removed for spatial mapping
        adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
        adata = adata[:, ~adata.var['mt'].values]

        adata_vis = adata.copy()
        adata_vis.raw = adata_vis

        ## READ IN REFERENCE DATA
        adata_raw = sc.read(f'{reg_path}sc.h5ad')

        # Export cell type expression signatures:
        inf_aver = adata_raw.raw.var.copy()
        inf_aver = inf_aver.loc[:, [f'means_cov_effect_{covariate_col_names}_{i}' for i in adata_raw.obs[covariate_col_names].unique()]]
        from re import sub
        inf_aver.columns = [sub(f'means_cov_effect_{covariate_col_names}_{i}', '', i) for i in adata_raw.obs[covariate_col_names].unique()]
        inf_aver = inf_aver.iloc[:, inf_aver.columns.argsort()]

        # scale up by average sample scaling factor
        inf_aver = inf_aver * adata_raw.uns['regression_mod']['post_sample_means']['sample_scaling'].mean()

        ## RUN CELL2LOCATION ##
        r = cell2location.run_cell2location(

              # Single cell reference signatures as pd.DataFrame 
              # (could also be data as anndata object for estimating signatures analytically - `sc_data=adata_snrna_raw`)
              sc_data=inf_aver, 
              # Spatial data as anndata object
              sp_data=adata_vis,

              # the column in sc_data.obs that gives cluster idenitity of each cell
              summ_sc_data_args={'cluster_col': covariate_col_names},

              train_args={'use_raw': True, # By default uses raw slots in both of the input datasets.
                          'n_iter': 30000, # Increase the number of iterations if needed (see below)

                          # Whe analysing the data that contains multiple samples, 
                          # cell2location will select a model version which pools information across samples
                          # For details see https://cell2location.readthedocs.io/en/latest/cell2location.models.html#module-cell2location.models.CoLocationModelNB4E6V2
                          'sample_name_col': 'sample'}, # Column in sp_data.obs with Sample ID

              # Number of posterios samples to use for estimating parameters,
              # reduce if not enough GPU memory
              posterior_args={'n_samples': 1000}, 


              export_args={'path': results_folder + 'std_model/', # path where to save results
                           'run_name_suffix': '_' + dataset_type # optinal suffix to modify the name the run
                          },

              model_kwargs={ # Prior on the number of cells, cell types and co-located combinations

                            'cell_number_prior': {
                                # Use visual inspection of the tissue image to determine 
                                # the average number of cells per spot,
                                # an approximate count is good enough:
                                'cells_per_spot': 8, 
                                # Prior on the number of cell types (or factors) in each spot
                                'factors_per_spot': 7, 
                                # Prior on the number of correlated cell type combinations in each spot
                                'combs_per_spot': 2.5
                            },

                             # Prior on change in sensitivity between technologies
                            'gene_level_prior':{
                                # Prior on average change in expression level from scRNA-seq to spatial technology,
                                # this reflects your belief about the sensitivity of the technology in you experiment
                                'mean': 1/2, 
                                # Prior on how much individual genes differ from that average,
                                # a good choice of this value should be lower that the mean
                                'sd': 1/4
                            }
              }
        )

            
if __name__ == '__main__':
    main()