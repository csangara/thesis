cd /srv/scratch/chananchidas/
SCRIPT_PATH=thesis/Scripts/synthetic_data_generation/cell2location
OUT_DIR=data/cell2location_sc_10x_braincortex_robin_rawcounts_default/
mkdir -p $OUT_DIR

N_SPOTS=100
if [ "$1" != "" ]; then
    N_SPOTS=$1
fi

# Step 1 (returns three seeds)
python $SCRIPT_PATH/split_sc.py data/seurat_obj_scrnaseq_cortex_filtered_rawcounts.h5ad \
data/seurat_obj_scrnaseq_cortex_filtered_rawcounts_metadata.csv \
--out_dir $OUT_DIR --annotation_col bio_celltype
 
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