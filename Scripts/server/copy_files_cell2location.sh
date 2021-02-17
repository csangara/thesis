DATASET=$1
for dir in data/$DATASET/cell2location_results/std_model/LocationModelLinearDependent*/
do
    term='genes_(\w+)/'
    [[ $dir =~ $term ]]
    newname="${DATASET}_${BASH_REMATCH[1]}_cell2location.csv"
    echo $dir
    # echo "${BASH_REMATCH[1]}"
    echo $newname
    cp $dir/W_cell_density_q05.csv data/$DATASET/cell2location_results_compiled/$newname
done
