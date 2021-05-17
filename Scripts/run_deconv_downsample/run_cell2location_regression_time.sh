#!/bin/bash

SECONDS=0
python /srv/scratch/chananchidas/scripts/run_cell2location_regression.py "$@"
OUT_DIR=${@: -1}
echo $SECONDS > $OUT_DIR/model_cell2location_info.txt


