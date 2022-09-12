#!/bin/bash
# Extract information about the lanes in the 02-raw folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

for i in $(ls -1 02-raw/*.fastq.gz)
do
    # TEST:
    # echo ${$i%.fastq}.gz
    basename $i | grep -v "_R2\.fastq\.gz" | perl -pe 's/\.fastq\.gz//'
done > 01-info_files/lane_info.txt

