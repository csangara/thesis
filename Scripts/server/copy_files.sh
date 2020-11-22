for dir in LocationModelLinearDependent*/
do
    term='genes_(\w+)/'
    [[ $dir =~ $term ]]
    newname="allen_cortex_dwn_${BASH_REMATCH[1]}_cell2location.csv"
    # echo $dir
    # echo "${BASH_REMATCH[1]}"
    echo $newname
    cp $dir/W_cell_density_q05.csv allen_cortex_dwn_results/$newname
done
