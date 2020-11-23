possible_dataset_types="real real_top1 real_top1_uniform real_top2_overlap real_top2_overlap_uniform real_missing_celltypes_visium artificial_uniform_distinct artificial_diverse_distinct artificial_uniform_overlap artificial_diverse_overlap artificial_dominant_celltype_diverse artificial_partially_dominant_celltype_diverse artificial_missing_celltypes_visium"

for dataset_type in $possible_dataset_types
do
    newname="allen_cortex_dwn_${dataset_type}_stereoscope.tsv"
    #echo $newname
    cp ${dataset_type}/allen_cortex_dwn_${dataset_type}_synthvisium/W*.tsv all_results/$newname
done
