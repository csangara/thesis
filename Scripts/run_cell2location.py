##### COMMAND LINE PART #####
import sys

if len(sys.argv) < 3:
    sys.stderr.write("USAGE: python %s <CUDA_DEVICE_ID> <DATASET_TYPE_1> ... <DATASET_TYPE_N>\n" % sys.argv[0])
    sys.exit(1)

cuda_device = sys.argv[1]
dataset_list = sys.argv[2:]

possible_dataset_types = ["real", "real_top1","real_top1_uniform","real_top2_overlap","real_top2_overlap_uniform",
                           "real_missing_celltypes_visium", "artificial_uniform_distinct", "artificial_diverse_distinct",
                           "artificial_uniform_overlap", "artificial_diverse_overlap", "artificial_dominant_celltype_diverse",
                           "artificial_partially_dominant_celltype_diverse", "artificial_missing_celltypes_visium"]

assert cuda_device in ["0", "1", "2"], "invalid device id"
assert all(dataset_type in possible_dataset_types for dataset_type in dataset_list), "invalid dataset type"
##### MAIN PART #####

import os
os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np

data_type = 'float32'

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=' + data_type + ',force_device=True'

sys.path.insert(1, '/srv/scratch/chananchidas/data/')

import cell2location
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

sp_data_folder = '/srv/scratch/chananchidas/data/'
sc_data_folder = '/srv/scratch/chananchidas/data/'
results_folder = '/srv/scratch/chananchidas/results/'

## READ IN SPATIAL DATA ##
for dataset_type in dataset_list:
    print("\nDataset type: " + dataset_type)
    adata = sc.read_h5ad(sp_data_folder + "allen_cortex_dwn_" + dataset_type + "_synthvisium.h5ad")
    print("Read in file from " + sp_data_folder + "allen_cortex_dwn_" + dataset_type + "_synthvisium.h5ad")
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
    reg_mod_name = 'RegressionGeneBackgroundCoverageTorch_24covariates_1404cells_29093genes'
    reg_path = f'{results_folder}regression_model/{reg_mod_name}/'

    adata_raw = sc.read(f'{reg_path}sc.h5ad')

    # Export cell type expression signatures:
    covariate_col_names = 'celltype'

    inf_aver = adata_raw.raw.var.copy()
    inf_aver = inf_aver.loc[:, [f'means_cov_effect_{covariate_col_names}_{i}' for i in adata_raw.obs[covariate_col_names].unique()]]
    from re import sub
    inf_aver.columns = [sub(f'means_cov_effect_{covariate_col_names}_{i}', '', i) for i in adata_raw.obs[covariate_col_names].unique()]
    inf_aver = inf_aver.iloc[:, inf_aver.columns.argsort()]

    # scale up by average sample scaling factor
    inf_aver = inf_aver * adata_raw.uns['regression_mod']['post_sample_means']['sample_scaling'].mean()

    sc.settings.set_figure_params(dpi = 100, color_map = 'RdPu', dpi_save = 100,
                                  vector_friendly = True, format = 'pdf',
                                  facecolor='white')

    ## RUN CELL2LOCATION ##
    r = cell2location.run_cell2location(

          # Single cell reference signatures as pd.DataFrame 
          # (could also be data as anndata object for estimating signatures analytically - `sc_data=adata_snrna_raw`)
          sc_data=inf_aver, 
          # Spatial data as anndata object
          sp_data=adata_vis,

          # the column in sc_data.obs that gives cluster idenitity of each cell
          summ_sc_data_args={'cluster_col': "celltype"},

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
