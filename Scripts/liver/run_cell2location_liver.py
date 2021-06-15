##### PARSING COMMAND LINE ARGUMENTS #####
import argparse as arp
import os

def main():
    prs = arp.ArgumentParser()
    
    prs.add_argument('sp_data_path', type = str, help = 'path to spatial data')
    
    prs.add_argument('result_dir', type = str, help = 'directory to regression model and results')
    
    prs.add_argument('cuda_device', type = str, help = "index of cuda device ID, from 0-7")
    
    prs.add_argument('-a','--annotation_column', default = 'celltype',
             type = str, help = 'column name for covariate')
    
    prs.add_argument('-r', '--regression_model_path', default = None,
                     type = str, help = 'path to regression model')
    
    prs.add_argument('-s', '--slide', default="1",
                     type = str, help = 'select slide 1-4, or all')

    args = prs.parse_args()
    
    cuda_device = args.cuda_device
    sp_data_path = args.sp_data_path
    results_folder = args.result_dir
    covariate_col_names = args.annotation_column
    slide = args.slide
      
    if args.regression_model_path is None:
        regression_model_output = os.listdir(results_folder + "/regression_model")[0]
        reg_path = f'{results_folder}regression_model/{regression_model_output}/'
    else:
        reg_path = args.regression_model_path
        
    assert cuda_device in ["0", "1", "2", "3", "4", "5", "6", "7"], "invalid device id"
    assert slide in ["1", "2", "3", "4"] or slide == "all", "slide does not exist"
    if slide.isdigit():
        assert 'filtered_feature_bc_matrix.h5' in os.listdir(sp_data_path + "/JBO0" + slide), "file path does not contain h5 feature matrix"
    else:
        assert all('filtered_feature_bc_matrix.h5' in os.listdir(sp_data_path + "/JBO0" + str(i)) for i in range(1,5)),\
        "one or more file path does not contain h5 feature matrix"
    
    ##### MAIN PART #####
    os.environ["CUDA_VISIBLE_DEVICES"]=cuda_device
    os.environ["CPATH"]="/usr/local/cuda/include:$CPATH" #To use cuDNN

    import sys
    import scanpy as sc
    import anndata
    import pandas as pd
    import numpy as np

    data_type = 'float32'

    # this line forces theano to use the GPU and should go before importing cell2location
    os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=' + data_type + ',force_device=True'

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
    if slide == "all":
        # We will merge all slides together in one adata object
        adata_list, sample_name = [], []
        for i in range(1,5):
            name = 'JBO0'+ str(i)
            temp_adata = sc.read_visium(sp_data_path + "/" + name)
            print("Read in file from " + sp_data_path + "/" + name)
            temp_adata.var_names_make_unique()
            temp_adata.var["mt"] = temp_adata.var_names.str.startswith("mt-")
            sc.pp.calculate_qc_metrics(temp_adata, qc_vars=["mt"], inplace=True)
            temp_adata.obs['sample'] = name
            sample_name.append(name)
            adata_list.append(temp_adata)

        adata = adata_list[0].concatenate(adata_list[1:], batch_key="sample", uns_merge="unique", \
                                      batch_categories = sample_name, index_unique=None)
    else:
        adata = sc.read_visium(sp_data_path + "/JBO0" + slide)     
        print("Read in file from " + sp_data_path)
        adata.var_names_make_unique()
        adata.obs['sample'] = "JBO0" + slide
        adata.var['mt'] = adata.var_names.str.startswith("mt-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    
    adata.obs_names_make_unique()    
    # Calculate QC metrics and filter    
    print("Before filtering: {} spots and {} genes".format(*adata.shape))
    adata.var['SYMBOL'] = adata.var_names
    sc.pp.filter_cells(adata, min_counts=11000)
    sc.pp.filter_cells(adata, max_counts=50000)
    adata = adata[adata.obs["pct_counts_mt"] < 20]
    sc.pp.filter_genes(adata, min_cells=10)

    # mitochondria-encoded (MT) genes should be removed for spatial mapping
    adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
    adata = adata[:, ~adata.var['mt'].values]
    print("After filtering: {} spots and {} genes".format(*adata.shape))
    
    adata_vis = adata.copy()
    adata_vis.raw = adata_vis

    ## READ IN REFERENCE DATA
    adata_raw = sc.read(f'{reg_path}sc.h5ad')

    # Export cell type expression signatures:
    inf_aver = adata_raw.raw.var.copy()
    inf_aver.index = adata_raw.raw.var['SYMBOL']
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
          verbose=True,
          # the column in sc_data.obs that gives cluster idenitity of each cell
          summ_sc_data_args={'cluster_col': covariate_col_names},

          train_args={'use_raw': True, # By default uses raw slots in both of the input datasets.
                      'n_iter': 15000, # Increase the number of iterations if needed (see below)

                      # Whe analysing the data that contains multiple samples, 
                      # cell2location will select a model version which pools information across samples
                      # For details see https://cell2location.readthedocs.io/en/latest/cell2location.models.html#module-cell2location.models.CoLocationModelNB4E6V2
                      'sample_name_col': 'sample'}, # Column in sp_data.obs with Sample ID

          # Number of posterios samples to use for estimating parameters,
          # reduce if not enough GPU memory
          posterior_args={'n_samples': 1000}, 


          export_args={'path': results_folder + 'std_model/', # path where to save results
                       'run_name_suffix': '' # optinal suffix to modify the name the run
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