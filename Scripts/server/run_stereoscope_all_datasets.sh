#!/bin/bash
REP=$1
cuda_device=$2
RUN=""
# datasets: brain_cortex cerebellum_cell cerebellum_nucleus hippocampus kidney pbmc scc_p5
#for RUN in run3
#do
for DATASET in brain_cortex cerebellum_cell cerebellum_nucleus hippocampus kidney pbmc scc_p5
do
  ./scripts/run_stereoscope_all.sh  \
  data/test_set/${DATASET}_test.h5ad  data/${DATASET}_generation/$REP/$RUN/stereoscope_results/model/  \
  data/${DATASET}_generation/$REP/${DATASET}_generation_  data/${DATASET}_generation/$REP/$RUN/stereoscope_results/ \
   $cuda_device  celltype  recommended
done
#done