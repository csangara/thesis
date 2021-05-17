#!/bin/bash
if [ "$#" -lt 6 ]; then
  echo "This script fits a model to downsampled spatial data"
  echo "USAGE: run_stereoscope.sh < path to sc fitted model > < path to spatial data > < output directory > < path to downsample file > < cuda device > < task_id > "
  exit
fi

echo "Using GPU$5"
export CUDA_VISIBLE_DEVICES=$5
export LD_LIBRARY_PATH=/srv/scratch/chananchidas/anaconda3/envs/stereoscope/lib/

SC_MODEL_PATH=$1
SP_DATA_PATH=$2
SP_DATA_DIR=$(dirname $SP_DATA_PATH)
OUT_DIR=$3
DS_FILE_PATH=$4
i=$6

mkdir -p $OUT_DIR/$i
# Downsample data
python /srv/scratch/chananchidas/scripts/downsample_h5ad.py $SP_DATA_PATH $DS_FILE_PATH $i

SECONDS=0

stereoscope run --sc_fit $SC_MODEL_PATH/R*.tsv $SC_MODEL_PATH/logits*.tsv \
--st_cnt $SP_DATA_DIR/temp_ds.h5ad -o $OUT_DIR/$i/ --gpu

echo $SECONDS > $OUT_DIR/brain_cortex_${i}_stereoscope_info.txt

rm $SP_DATA_DIR/temp_ds.h5ad

  
