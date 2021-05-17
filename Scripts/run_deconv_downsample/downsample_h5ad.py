import scanpy as sc
import sys
import os

if len(sys.argv) < 4:
    print("usage: downsample_h5ad.py < data_path > < downsample_file_path > < task_id > < opt: cells, genes or both > < opt: output_name >")
    exit()
    
data_path = sys.argv[1]
ds_file_path = sys.argv[2]
task_id = sys.argv[3]
if len(sys.argv) >= 5:
    ds_type = sys.argv[4]
else:
    ds_type = "both"
output_name = "temp_ds.h5ad"

if len(sys.argv) == 6:
    output_name = sys.argv[5]
    
assert ds_type in ['cells', 'genes', 'both'], "invalid downsampling option"

#sp_data = "/srv/scratch/chananchidas/thesis/downsampling/brain_cortex/rep1/brain_cortex_art_uni_distinct_18594spots_17538genes.h5ad"
#ds_file_path = "/srv/scratch/chananchidas/thesis/downsampling/brain_cortex/rep1/synthvisium_downsample_info/"
adata = sc.read_h5ad(data_path)

spots, genes = adata.obs_names, adata.var_names

if ds_type in ['cells', 'both']:
    with open(ds_file_path + "cells" + task_id + ".txt", "r") as f:
        spots = [line.strip() for line in f]
if ds_type in ['genes', 'both']:
    with open(ds_file_path + "genes" + task_id + ".txt", "r") as f:
        genes = [line.strip() for line in f]
    
adata = adata[spots,genes]
adata.raw = adata
print("Data now has {} spots and {} genes".format(*adata.shape))

# Save new file as temp_ds.h5ad
new_file_name = os.path.join(os.path.dirname(data_path), output_name)
adata.write_h5ad(new_file_name)