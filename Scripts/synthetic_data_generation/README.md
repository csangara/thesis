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

