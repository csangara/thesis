# EXECUTE THIS AT THE CORRECT DIRECTORY!
: <<'ADDING_REP'
DATASET=$1
for DIR in $DATASET/*

do
	REP=$(basename $DIR)
	# echo $DIR
	# echo $REP
	for FILE in $DIR/*
	do
		FILENAME=$(basename $FILE .rds)
		# echo $FILENAME
		NEWNAME="${FILENAME}_${REP}.rds"
		mv $FILE $DATASET/$REP/$NEWNAME
	done
'
ADDING_REP

: <<'RENAMING_FILES'
DATASET=$1
for DIR in $DATASET/*
do
	find $DIR -exec rename 's/_generation//' {} + # Step 3 (comment out others)
	for FILE in $DIR/*
	do
		find $FILE -exec rename -v 's/_generation_/_/' *.rds {} + # Step 1
		find $FILE -exec rename -n 's/_synthvisium_/_/' *.rds {} + # Step 2
	done
	
done
'
RENAMING_FILES

# find . -maxdepth 1 -exec rename 's/_generation//' {} + # For renaming base directories

: <<'MOVING FILES'
DATASET_TYPES="artificial_uniform_distinct artificial_diverse_distinct artificial_uniform_overlap artificial_diverse_overlap artificial_dominant_celltype_diverse artificial_partially_dominant_celltype_diverse artificial_dominant_rare_celltype_diverse artificial_regional_rare_celltype_diverse"
DATASET=$1
for DATASET_TYPE in $DATASET_TYPES
do	
	mkdir -p $DATASET/$DATASET_TYPE
	echo $DATASET_TYPE
	#echo $DATASET/rep*/*$DATASET_TYPE*.rds
	mv $DATASET/rep*/*$DATASET_TYPE*.rds $DATASET/$DATASET_TYPE
done
rmdir $DATASET/rep*
'
MOVING_FILES
