# Running deconvolution methods

## Methods in R

MuSiC, RCTD, and SPOTlight can be run from the command line as follows:

``Rscript --vanilla run_deconv_prism_\*.R <dataset_index> <replication no.> <run no.>``

(This was just tweaked from the file `../run_deconv.R` to accept command line arguments so it can be run on a gridengine more easily.)

## Methods in python

The file `convertRDStoLoom_synthvisium.R` provides functions and sample code for converting a rds file to loom or h5ad. Both cell2location and stereoscope contain two steps, model building and model fitting.
I have provided three scripts for each method, two for each step and the third script is for running both together.

### cell2location

This can be run by either using the script `run_cell2location_regression.py` then `run_cell2location.py`, or by only running `run_cell2location_all.py`.

**Running them separately:**
```
python scripts/run_cell2location_regression.py data/brain_cortex_test.h5ad $cuda_device \
-o data/brain_cortex_generation/rep1/run1/cell2location_results/

python scripts/run_cell2location.py data/brain_cortex_generation/rep1/brain_cortex_generation_ \
data/brain_cortex_generation/rep1/run1/cell2location_results/ $cuda_device recommended \
-r RegressionGeneBackgroundCoverageTorch_9covariates_1693cells_9458genes
```

**Running both together:**
```
python scripts/run_cell2location_all.py data/brain_cortex_test.h5ad data/brain_cortex_generation/rep1/brain_cortex_generation_ \
$cuda_device recommended -o data/brain_cortex_generation/rep1/run1/cell2location_results/
```

**Docs**
```
run_cell2location_regression.py [-h] [-o OUT_DIR]
                                       [-a ANNOTATION_COLUMN]
                                       sc_data_path cuda_device

positional arguments:
  sc_data_path          path to single cell h5ad count data
  cuda_device           index of cuda device ID, from 0-7

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_DIR, --out_dir OUT_DIR
                        directory for regression model
  -a ANNOTATION_COLUMN, --annotation_column ANNOTATION_COLUMN
                        column name for covariate
```

```
run_cell2location.py [-h] [-a ANNOTATION_COLUMN] [-r REGRESSION_MODEL_PATH]
                            sp_data_prefix result_dir cuda_device dataset_type
                            [dataset_type ...]

positional arguments:
  sp_data_prefix        path and file prefix of spatial data
  result_dir            directory to results (and regression model)
  cuda_device           index of cuda device ID, from 0-7
  dataset_type          multiple arguments of dataset types to be deconvolved,
                        or "all" or "recommended"

optional arguments:
  -h, --help            show this help message and exit
  -a ANNOTATION_COLUMN, --annotation_column ANNOTATION_COLUMN
                        column name for covariate
  -r REGRESSION_MODEL_PATH, --regression_model_path REGRESSION_MODEL_PATH
                        path to regression model
```

```
run_cell2location_all.py [-h] [-o OUT_DIR] [-a ANNOTATION_COLUMN]
                                sc_data_path sp_data_prefix cuda_device
                                dataset_type [dataset_type ...]
```

### stereoscope

This can be run by either using the script `run_stereoscope_regression.sh` then `run_stereoscope.sh`, or by only running `run_stereoscope_all.sh`.

**Running them separately:**
```
./scripts/run_stereoscope_regression.sh data/brain_cortex_test.h5ad celltype \
data/brain_cortex_generation/rep1/run1/stereoscope_results/model/ $cuda_device

./scripts/run_stereoscope.sh data/brain_cortex_generation/rep1/run1/stereoscope_results/model/ \
data/brain_cortex_generation/rep1/brain_cortex_generation_ \
data/brain_cortex_generation/rep1/run1/stereoscope_results/ $cuda_device recommended
```

**Running both together:**
```
./scripts/run_stereoscope_all.sh data/brain_cortex_test.h5ad
data/brain_cortex_generation/rep1/run1/stereoscope_results/model/ \
data/brain_cortex_generation/rep1/brain_cortex_generation_ \
data/brain_cortex_generation/rep1/run1/stereoscope_results/ \
$cuda_device celltype recommended
```

**Docs**
```
run_stereoscope_regression.sh < path to sc file > < covariate column name > < output directory > < cuda device >

run_stereoscope.sh < path to sc fitted model > < prefix to spatial data > < output directory > < cuda device >
< dataset 1 > ... < dataset n >

run_stereoscope_all.sh < path to sc file > < model output directory > < prefix to spatial data > < results output directory >
< cuda device > < covariate column name > < dataset 1 > ... < dataset n >
```

### Compiling result files
Since the output of these two methods are stored in different folders (and also with ambiguous names), the `copy_files_` scripts can be used to compile results. Both have three arguments,
the input directory, output directory, and dataset name.

```
./scripts/copy_files_cell2location.sh data/brain_cortex_generation/rep1/run1/cell2location_results/
results/brain_cortex_generation/rep1_run1/cell2location brain_cortex_generation

./scripts/copy_files_stereoscope.sh data/brain_cortex_generation/rep1/run1/stereoscope_results/
/results/brain_cortex_generation/rep1_run1/stereoscope/ brain_cortex_generation
```
