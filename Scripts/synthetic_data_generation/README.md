# Generating synthetic spots

## SPOTlight

https://github.com/MarcElosua/SPOTlight/blob/master/R/test_spot_fun.R

**Input:** Seurat object

**Output:** A list of two elements, (1) dataframe where each column is the proportion of each spot (cell_composition), and (2) a sparse matrix with the expression of synthetic spots, where the rows are genes and columns are spots (topic_profiles)

**Parameters:**
* *clust_vr*: name of the variable containing the cell clustering
* *n*: number of spots to generate

`test_spot_fun(seurat_obj, clust_vr="celltype", n=1000)`

## Stereoscope

https://github.com/almaan/stereoscope/tree/master/comparison/synthetic_data_generation

**Input:**
* Count data: a .tsv file with cells as observations (rows) and genes as variables (columns), each cell (row) should have a unique label
* Annotation data: a .tsv file with the same rownames as the count data file, either with one single column listing the annotations, or multiple columns where the column containing the labels should be named *bio_celltype*, i.e.:

|           cell           | bio_celltype |
|:------------------------:|:------------:|
| 10X04_1_AAACATACCTTATC-1 |   Oligos_5   |

Requires at least 30 cells per cell type!
 
**Output:**
* counts.st_synth.tsv : expression data (spot x gene) [by default, 500 genes and 1000 spots]
* members.st_synth.tsv : number of cells (spot x number of cells in each spot)
* proportions.st_synth.tsv : proportion values (spot x proportions of each cell type)

```
make_st_set.py [-h] -c SC_COUNTS -l SC_LABELS [-ns N_SPOTS]
                      [-ng N_GENES] [-o OUT_DIR]
                      [-ncr N_CELL_RANGE N_CELL_RANGE] [-t TAG]
                      
optional arguments:
  -h, --help            show this help message and exit
  -c SC_COUNTS, --sc_counts SC_COUNTS
                        path to single cell count data
  -l SC_LABELS, --sc_labels SC_LABELS
                        path to single cell labels/annotation data
  -ns N_SPOTS, --n_spots N_SPOTS
                        number of spots
  -ng N_GENES, --n_genes N_GENES
                        number of genes
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory
  -ncr N_CELL_RANGE N_CELL_RANGE, --n_cell_range N_CELL_RANGE N_CELL_RANGE
                        lower bound (first argument) and upper bound (second
                        argument) for the number of cells at each spot
  -t TAG, --tag TAG     tag to mark data se with
```

**Sample code**

`
python make_st_set.py --sc_counts data/data_rawcounts.tsv
--sc_labels data/data_metadata.tsv  --n_spots 2000 --n_genes 10000
--n_cell_range 2 10 --out_dir data/synthetic_data_stereoscope/
`

## cell2location

https://github.com/emdann/ST_simulation

**Input:**
* Count data: h5ad file
* Annotation data: .csv file (make sure ANNO_COL argument below matches annotation column name in csv file)

**Output:**
* csv file containing the number of cells per cell type in each spot
* csv file containing the count matrix for the simulated ST spots
* csv file containing the number of UMIs per cell type in each spot, for benchmarking deconvolution methods that model number of UMIs

Sample code (step-by-step explanation below)
```
SCRIPT_PATH=thesis/Scripts/synthetic_data_generation/cell2location
OUT_DIR=data/synthetic_data_cell2location/
mkdir -p $OUT_DIR

N_SPOTS=1000

# Step 1 (returns three seeds)
python $SCRIPT_PATH/split_sc.py data/data_rawcounts.h5ad \
data/data_metadata.csv --out_dir $OUT_DIR --annotation_col bio_celltype
 
# Step 2
# This gives all three seeds that were returned
seed=$(ls $OUT_DIR/labels_generation* | sed 's/.*_//' | sed 's/.p//')
# Loop through each seed
for SEED in $seed
do
    python $SCRIPT_PATH/assemble_design.py --tot_spots $N_SPOTS --mean_high 3 --mean_low 1 \
    --out_dir $OUT_DIR --annotation_col bio_celltype $SEED

    # Step 3
    python $SCRIPT_PATH/assemble_composition.py --tot_spots $N_SPOTS \
    --annotation_col bio_celltype --out_dir $OUT_DIR $SEED
 
    # Step 4
    python $SCRIPT_PATH/assemble_st.py --out_dir $OUT_DIR --annotation_col bio_celltype $SEED
done
```

**Step 1:** Split single-cell dataset. we split the cells in the single-cell dataset in a 'generation set', that will be used to simulate the ST spots, and a 'validation' set, that will be used to train the deconvolution models that we want to benchmark.


```
split_sc.py [-h] [--annotation_col ANNO_COL] [--out_dir OUT_DIR]
                   h5ad_file annotation_file

positional arguments:
  h5ad_file             path to h5ad file with raw single-cell data
  annotation_file       path to csv file with cell annotations


optional arguments:
  -h, --help            show this help message and exit
  --annotation_col ANNO_COL
                        Name of column to use in annotation file (default:
                        annotation_1)
  --out_dir OUT_DIR     Output directory
```
  
**Output:** generation and validation count matrices and cell type annotations are saved as pickle files (.p), with a random seed identifying the split (by default, three seeds are generated). You will see 12 files (3 each of counts_generation_<seed>.p, labels_generation_<seed>.p, counts_validation_<seed>.p, and labels_validation_<seed>.p).

**Step 2:** Build design matrix: in this step we define which cell types are (A) low/high density and (B) Uniformly present in all the spots or localized in few spots (regional). To generate synthetic spots with ~10 cells per spot (as seen with nuclear segmentation on Visium spots) they recommend setting the mean number of cells per spot per cell type < 5.

```
assemble_design.py [-h] [--tot_spots TOT_SPOTS] [--mean_high MEAN_HIGH]
               [--mean_low MEAN_LOW] [--percent_uniform PERCENT_UNIFORM]
               [--percent_sparse PERCENT_SPARSE] [--annotation_col ANNO_COL]
               [--out_dir OUT_DIR] [--assemble_id ASSEMBLE_ID]
               seed

positional arguments:
  seed                  random seed of split

optional arguments:
  -h, --help            show this help message and exit
  --tot_spots TOT_SPOTS
                        Total number of spots to simulate
  --mean_high MEAN_HIGH
                        Mean cell density for high-density cell types
  --mean_low MEAN_LOW   Mean cell density for low-density cell types
  --percent_uniform PERCENT_UNIFORM
                        Sparsity of uniform cell types (% non-zero spots of
                        total spots)
  --percent_sparse PERCENT_SPARSE
                        Sparsity of sparse cell types (% non-zero spots of
                        total spots)
  --annotation_col ANNO_COL
                        Name of column to use in annotation file (default:
                        annotation_1)
  --out_dir OUT_DIR     Output directory
  --assemble_id ASSEMBLE_ID
                        ID of ST assembly
 ```
                        
**Output:** synthetic_ST_seed\${seed}_design.csv contains the design used for the simulation, e.g.,

|       | uniform | density | nspots | mean_ncells |
|-------|:-------:|:-------:|:------:|:-----------:|
|   Vip |    0    |    1    |    2   |   0.721145  |
| Lamp5 |    1    |    1    |   100  |    0.3146   |

**Step 3:** Assemble cell type composition per spot: based on the design matrix, we define the cell type composition of each spot i.e. how many cells per cell type are in each spot. An assemble ID is used to identify the assembly (we assemble many composition matrices with the same design).

```
assemble_composition.py [-h] [--tot_spots TOT_SPOTS]
                               [--out_dir OUT_DIR] [--assemble_id ASSEMBLE_ID]
                               [--annotation_col ANNO_COL]
                               seed

positional arguments:
  seed                  random seed of split

optional arguments:
  -h, --help            show this help message and exit
  --tot_spots TOT_SPOTS
                        Total number of spots to simulate
  --out_dir OUT_DIR     Output directory
  --assemble_id ASSEMBLE_ID
                        ID of ST assembly
  --annotation_col ANNO_COL
                        Name of column to use in annotation file (default:
                        annotation_1)
 ```
 
 **Output:** synthetic_ST_seed\${seed}_\${assemble_id}_composition.csv contains the number of cells per cell type in each spot, for benchmarking deconvolution models.
 
 **Step 4:** Assemble simulated ST spots
 ```
 assemble_st.py [-h] [--out_dir OUT_DIR] [--assemble_id ASSEMBLE_ID]
                      [--annotation_col ANNO_COL]
                      seed

positional arguments:
  seed                  random seed of split


optional arguments:
  -h, --help            show this help message and exit
  --out_dir OUT_DIR     Output directory
  --assemble_id ASSEMBLE_ID
                        ID of ST assembly
  --annotation_col ANNO_COL
                        Name of column to use in annotation file (default:
                        annotation_1)
```                        
**Output:**
* synthetic_ST_seed\${seed}_\${assemble_id}_counts.csv contains the count matrix for the simulated ST spots
* synthetic_ST_seed\${seed}_\${assemble_id}_umis.csv contains the number of UMIs per cell type in each spot, for benchmarking deconvolution methods that model number of UMIs
