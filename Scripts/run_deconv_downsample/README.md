# Running deconvolution tools on downsampled inputs

This is the script for evaluating the runtime and scalability of each method. First, run `downsample_dataset.R` to generate a list of txt files containing the names of cells/spots and genes for the downsampled matrices. By default the combinations tested are as follows, for synthetic data

|         |   1  |   2  |   3   |   4   |   5  |   6  |   7   |   8   |   9  |  10  |   11  |   12  |   13  |   14  |   15  |   16  |
|:-------:|:----:|:----:|:-----:|:-----:|:----:|:----:|:-----:|:-----:|:----:|:----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
| n_spots |  100 |  100 |  100  |  100  | 1000 | 1000 |  1000 |  1000 | 5000 | 5000 |  5000 |  5000 | 10000 | 10000 | 10000 | 10000 |
| n_genes | 1000 | 5000 | 10000 | 15000 | 1000 | 5000 | 10000 | 15000 | 1000 | 5000 | 10000 | 15000 |  1000 |  5000 | 10000 | 15000 |

And for the reference data, [1000, 5000, 10000] cells are tested and the genes are not downsampled. So you should have 16 x 2 files for synthetic data and 3 files for single-cell data.

## Methods in R

The scripts for MuSiC, RCTD and SPOTlight have been modified to work with array jobs in a grid engine. The scripts retrieve the `SGE_TASK_ID` variable from the environment to choose the correct downsampling file. So if you can run the file on the cluster, simply add the parameter `-t 1-16` (for synthetic data) and `-t 1-3` (for single-cell data) to `qsub`. You will have to comment and uncomment the respective code chunks for downsampling synthetic and/or single-cell data. The output is a text file (and the proportions matrix) containing the runtime of the process.

## Methods in Python

Things are a little more complicated with cell2location and stereoscope. We will divide it into two cases, depending on whether you want to downsample the synthetic data or the single-cell data.

### Downsampling synthetic data

In this case, you can just build one reference model for all 16 cases. For **cell2location**, the script `run_cell2location_regression_time.sh` keeps track of the time to build the model. This is actually just the wrapper for the function `run_cell2location_regression.py`.

`./run_cell2location_regression_time.sh data/brain_cortex_test.h5ad $cuda_device -o $OUT_DIR`

To fit the model, the script `run_cell2location_downsample.py` is a modified version of `run_cell2location.py` which takes the downsample files path and task_id as extra arguments.

```
run_cell2location_downsample.py [-h] [-a ANNOTATION_COLUMN]
                                       [-r REGRESSION_MODEL_PATH]
                                       sp_data result_dir cuda_device
                                       ds_file_path task_id

positional arguments:
  sp_data               path to full spatial dataset file
  result_dir            directory to regression model and results
  cuda_device           index of cuda device ID, from 0-7
  ds_file_path          path to files containing cell and genes names for
                        downsampling
  task_id               task id of the downsampling, from 1-16
```

A wrapper script for this has been provided called the `run_cell2location_downsample_time.sh`, i.e.,

```
DIR=thesis/downsampling/brain_cortex/rep1
./scripts/run_cell2location_downsample_time.sh $DIR/brain_cortex_synthvisium_full.h5ad \
$DIR/results/cell2location/ $cuda_device $DIR/synthvisium_downsample_info/ \
-r $DIR/results/cell2location/regression_model/$REGRESSION_MODEL_NAME/
```

For **stereoscope**, you can use the original script `../run_deconv/run_stereoscope_regression.sh` to build the model. If you add the `true` argument this will also keep the time for model building:

`./run_stereoscope_regression.sh data/brain_cortex_test.h5ad celltype $DIR/results/stereoscope/model/ $cuda_device true`

Then, run `run_stereoscope_downsample.sh`. This has the arguments `< path to sc fitted model > < path to spatial data > < output directory > < path to downsample file > < cuda device > < task_id >`. You can put it in a for loop as follows

```
for i in {1..16}; do
./run_stereoscope_downsample.sh $DIR/results/stereoscope/model/ $DIR/ brain_cortex_synthvisium_full.h5ad \
$DIR/results/stereoscope/ $DIR/synthvisium_downsample_info/ $cuda_device $i
done
```

### Downsampling single-cell reference

In this case you will have to build multiple models to be used with the same synthetic data. First, we have to downsample the reference data with `downsample_h5ad.py`.

```
for i in {1..3}; do python downsample_h5ad.py scrnaseq_reference_full.h5ad \
$DIR/scref_downsample_info/ref_ $i cells brain_temp_${i}.h5ad; done
```

To build the model in **cell2location**, I wrote a for loop wrapping the `run_cell2location_regression.py` function:

```
for i in {1..3}; do
SECONDS=0
python run_cell2location_regression.py brain_temp_${i}.h5ad $cuda_device -o $OUT_DIR
echo $SECONDS > $OUT_DIR/model_${i}_cell2location_info.txt
done
```

To fit the model, you can use `run_cell2location_downsample.py` and enter the path to each regression model separately. I also used a downsampled version of the synthetic data with 1000 spots and 15000 genes (case 8):

```
python run_cell2location_downsample.py $DIR/brain_cortex_synthvisium_full.h5ad \
$OUT_DIR $cuda_device $DIR/synthvisium_downsample_info/ 8 \
-r $DIR/results_sc/cell2location/regression_model/$REGRESSION_MODEL_NAME/
```

A similar process can be followed using **stereoscope**, starting from model building:

```
for i in {1..3}; do
SECONDS=0
stereoscope run --sc_cnt brain_temp_${i}.h5ad --label_colname celltype -o $OUT_DIR -n 5000 --gpu
echo $SECONDS > $OUT_DIR/../model_${i}_stereoscope_info.txt
done
```

To fit the model, first downsample the synthetic data

`python downsample_h5ad.py  $DIR/brain_cortex_synthvisium_full.h5ad  $DIR/synthvisium_downsample_info/ 8`

This temporary downsampled file is saved as `temp_ds.h5ad`. We can then fit the model in a for loop `i in {1..3}`

```
SC_MODEL_PATH=$OUT_DIR/model_${i}/
stereoscope run --sc_fit $SC_MODEL_PATH/R*.tsv $SC_MODEL_PATH/logits*.tsv \
--st_cnt temp_ds.h5ad -o $OUT_DIR/$i/ --gpu
```
