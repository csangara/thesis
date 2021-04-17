REP=$1
cuda_device=$2

# datasets: brain_cortex cerebellum_cell cerebellum_nucleus hippocampus kidney pbmc scc_p5

# for RUN in run1 run2 run3
#do
for DATASET in brain_cortex cerebellum_cell cerebellum_nucleus hippocampus kidney pbmc scc_p5
do
  python scripts/run_cell2location_all.py data/test_set/${DATASET}_test.h5ad data/${DATASET}_generation/$REP/${DATASET}_generation_ $cuda_device recommended -o data/${DATASET}_generation/$REP/$RUN/cell2location_results/
done
#done