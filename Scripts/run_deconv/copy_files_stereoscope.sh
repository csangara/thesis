DATASET=$1
OUT_DIR=$2

mkdir -p $OUT_DIR

for dir in data/$DATASET/stereoscope_results/*
do
  if [ $(basename $dir) != "model" ]
  then
    DATASET_TYPE=$(basename $dir)
    newname="${DATASET}_${DATASET_TYPE}_stereoscope.tsv"
    echo $DATASET_TYPE
    echo $newname
    
    cp data/$DATASET/stereoscope_results/$DATASET_TYPE/${DATASET}_${DATASET_TYPE}_synthvisium/W*.tsv \
    $OUT_DIR/$newname
  fi
done
exit
