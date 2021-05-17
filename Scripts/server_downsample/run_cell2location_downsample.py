##### PARSING COMMAND LINE ARGUMENTS #####
import argparse as arp
import os

def main():
    prs = arp.ArgumentParser()
    
    prs.add_argument('sp_data', type = str, help = 'path to full spatial dataset file')

    prs.add_argument('result_dir', type = str, help = 'directory to regression model and results')
    
    prs.add_argument('cuda_device', type = str, help = "index of cuda device ID (from 0-7) or cpu")
    
    prs.add_argument('ds_file_path', type = str, help = "path to files containing cell and genes names for downsampling")
                     
    prs.add_argument('task_id', type = str, help = "task id of the downsampling, from 1-16")
    
    prs.add_argument('-a','--annotation_column', default = 'celltype',
             type = str, help = 'column name for covariate')
    
    prs.add_argument('-r', '--regression_model_path', default = None,
                     type = str, help = 'path to regression model')
    
    args = prs.parse_args()
    
    cuda_device = args.cuda_device
    ds_file_path = args.ds_file_path
    task_id = args.task_id
    sp_data = args.sp_data
    results_folder = args.result_dir
    covariate_col_names = args.annotation_column
      
    if args.regression_model_path is None:
        regression_model_output = os.listdir(results_folder + "/regression_model")[0]
        reg_path = f'{results_folder}regression_model/{regression_model_output}/'
    else:
        reg_path = args.regression_model_path

    assert (cuda_device in ["0", "1", "2", "3", "4", "5", "6", "7"] or cuda_device == "cpu"), "invalid device id"
    assert os.path.isfile(sp_data), "invalid file path"
    
    ##### MAIN PART #####
    if cuda_device.isdigit():
        os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device
        os.environ["CPATH"]="/usr/local/cuda/include:$CPATH" #To use cuDNN
        device='cuda'
    else:
        device='cpu'
        
    import sys
    import scanpy as sc
    import anndata
    import pandas as pd
    import numpy as np

    # this line forces theano to use the GPU and should go before importing cell2location
    os.environ["THEANO_FLAGS"] = 'device=' + device + ',floatX=float32' + ',force_device=True'
    
    import cell2location
    import matplotlib as mpl
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    import seaborn as sns
    mpl.use('Agg')

    # silence scanpy that prints a lot of warnings
    import warnings
    warnings.filterwarnings('ignore')

    if not os.path.exists(results_folder + "std_model/"):
        os.makedirs(results_folder + "std_model/")
        
    ## READ IN SPATIAL DATA ##
    adata = sc.read_h5ad(sp_data)
    print("Read in file from " + sp_data)
    # Downsample
    with open(ds_file_path+"cells"+task_id+".txt", "r") as f:
        spots = [line.strip() for line in f]
    with open(ds_file_path+"genes"+task_id+".txt", "r") as f:
        genes = [line.strip() for line in f]
    adata = adata[spots,genes]
    print("Data now has {} spots and {} genes".format(*adata.shape))
    
    adata.obs['sample'] = "artificial_uniform_distinct"
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
          verbose=False,
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
                       'run_name_suffix': '_' + task_id # optinal suffix to modify the name the run
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