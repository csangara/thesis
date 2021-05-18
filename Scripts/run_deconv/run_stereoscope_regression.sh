#!/bin/bash
if [ "$#" -lt 4 ]; then
  echo "This script creates a model using the single-cell data"
  echo "USAGE: run_stereoscope_regression.sh < path to sc file > < covariate column name > < output directory > < cuda device > < save_time >"
  exit
fi

echo "Using GPU$4"
export CUDA_VISIBLE_DEVICES=$4
export LD_LIBRARY_PATH=/srv/scratch/chananchidas/anaconda3/envs/stereoscope/lib/

SC_PATH=$1
COV_NAME=$2
OUT_DIR=$3

echo "Reading $SC_PATH"
mkdir -p $OUT_DIR

SECONDS=0
stereoscope run --sc_cnt $SC_PATH --label_colname $COV_NAME -o $OUT_DIR -n 5000 --gpu

if [ "${@: -1}" = "true" ]; then
  echo $SECONDS > $OUT_DIR/../model_stereoscope_info.txt
fi