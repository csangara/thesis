#!/bin/bash
if [ "$#" -lt 7 ]; then
  echo "$# arguments"
  echo "This script creates a model using the single-cell data, then fits that model to spatial data"
  echo "Use 'all' or 'recommended' or type dataset types with a space in between."
  echo "USAGE: run_stereoscope_regression_all.sh < path to sc file > < model output directory > < prefix to spatial data > < results output directory > < cuda device > < covariate column name > < dataset 1 > ... < dataset n >"
  exit
fi

echo "Using GPU$5"
export CUDA_VISIBLE_DEVICES=$5
export LD_LIBRARY_PATH=/srv/scratch/chananchidas/anaconda3/envs/stereoscope/lib/

SC_PATH=$1
MODEL_OUT_DIR=$2
COV_NAME=$6

echo "Reading $SC_PATH"
mkdir -p $MODEL_OUT_DIR

stereoscope run --sc_cnt $SC_PATH --label_colname $COV_NAME -o $MODEL_OUT_DIR -n 5000 --gpu

if [ "$7" == "all" ]
then
  possible_dataset_types="real real_top1 real_top1_uniform real_top2_overlap real_top2_overlap_uniform real_missing_celltypes_visium artificial_uniform_distinct artificial_diverse_distinct artificial_uniform_overlap artificial_diverse_overlap artificial_dominant_celltype_diverse artificial_partially_dominant_celltype_diverse artificial_missing_celltypes_visium"
elif [ "$7" == "recommended" ]
then
  possible_dataset_types="artificial_uniform_distinct artificial_diverse_distinct artificial_uniform_overlap artificial_diverse_overlap artificial_dominant_celltype_diverse artificial_partially_dominant_celltype_diverse artificial_dominant_rare_celltype_diverse artificial_regional_rare_celltype_diverse"
else
  possible_dataset_types=${@:7}
fi

SP_PREFIX=$3
RESULT_OUT_DIR=$4

for dataset_type in $possible_dataset_types
do
    mkdir -p $RESULT_OUT_DIR/$dataset_type

    stereoscope run --sc_fit $MODEL_OUT_DIR/R*.tsv $MODEL_OUT_DIR/logits*.tsv \
    --st_cnt ${SP_PREFIX}${dataset_type}_synthvisium.h5ad \
    -o $RESULT_OUT_DIR/$dataset_type/ -n 5000 --gpu
done
