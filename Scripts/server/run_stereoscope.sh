#!/bin/bash
if [ "$#" -lt 4 ]; then
  echo "This script fits a model to spatial data"
  echo "Use 'all' or 'recommended' or type dataset types with a space in between."
  echo "USAGE: run_stereoscope.sh < path to sc fitted model > < prefix to spatial data > < output directory > < cuda device > < dataset 1 > ... < dataset n >"
  exit
fi

echo "Using GPU$4"
export CUDA_VISIBLE_DEVICES=$4
export LD_LIBRARY_PATH=/srv/scratch/chananchidas/anaconda3/envs/stereoscope/lib/

SC_MODEL_PATH=$1
SP_PREFIX=$2
OUT_DIR=$3

if [ "$5" == "all" ]
then
  possible_dataset_types="real real_top1 real_top1_uniform real_top2_overlap real_top2_overlap_uniform real_missing_celltypes_visium artificial_uniform_distinct artificial_diverse_distinct artificial_uniform_overlap artificial_diverse_overlap artificial_dominant_celltype_diverse artificial_partially_dominant_celltype_diverse artificial_missing_celltypes_visium"
elif [ "$5" == "recommended" ]
then
  possible_dataset_types="artificial_uniform_distinct artificial_diverse_distinct artificial_uniform_overlap artificial_diverse_overlap artificial_dominant_celltype_diverse artificial_partially_dominant_celltype_diverse artificial_dominant_rare_celltype_diverse artificial_regional_rare_celltype_diverse"
else
  possible_dataset_types=${@:5}
fi

for dataset_type in $possible_dataset_types
do
    mkdir -p $OUT_DIR/$dataset_type

    stereoscope run --sc_fit $SC_MODEL_PATH/R*.tsv $SC_MODEL_PATH/logits*.tsv \
    --st_cnt ${SP_PREFIX}${dataset_type}_synthvisium.h5ad \
    -o $OUT_DIR/$dataset_type/ -n 5000 --gpu
done