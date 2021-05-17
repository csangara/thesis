#!/usr/bin/env bash


for i in {1..16}
do
  echo $i
  SECONDS=0
  python /srv/scratch/chananchidas/scripts/run_cell2location_downsample.py "$@" $i
  OUT_DIR=$2
  echo $SECONDS > $OUT_DIR/brain_cortex_${i}_cell2location_info.txt  
done
