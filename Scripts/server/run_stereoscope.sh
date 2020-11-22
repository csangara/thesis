echo "Using GPU$1"
export CUDA_VISIBLE_DEVICES=$1
export LD_LIBRARY_PATH=/srv/scratch/chananchidas/anaconda3/envs/stereoscope/lib/
stereoscope test

possible_dataset_types="real real_top1 real_top1_uniform real_top2_overlap real_top2_overlap_uniform real_missing_celltypes_visium artificial_uniform_distinct artificial_diverse_distinct artificial_uniform_overlap artificial_diverse_overlap artificial_dominant_celltype_diverse artificial_partially_dominant_celltype_diverse artificial_missing_celltypes_visium"

for dataset_type in $possible_dataset_types
do
    echo "Reading data/raw/allen_cortex_dwn_${dataset_type}_synthvisium.h5ad"
    mkdir -p stereoscope_results/${dataset_type}
    
    # Running only the fitted model to the data
    stereoscope run --sc_fit stereoscope_test/R*.tsv stereoscope_test/logits*.tsv \
    --st_cnt data/raw/allen_cortex_dwn_${dataset_type}_synthvisium.h5ad \
    -o stereoscope_results/${dataset_type}/ -n 5000 --gpu
    
    # Running both the model and the data
    # stereoscope run --sc_cnt data/raw/allen_cortex_dwn_original.h5ad --label_colname subclass \
    # --st_cnt data/raw/allen_cortex_dwn_${dataset_type}_synthvisium.h5ad \
    # -o stereoscope_results/${dataset_type}/ -n 5000 --gpu
done