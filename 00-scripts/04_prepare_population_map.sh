#!/bin/bash
# Create population map file
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
INFO_FILES="01-info_files"

grep -v "#" $INFO_FILES/sample_information.csv | \
    awk '{print $3 "_" $4 "\t" $5}' | \
    sort -u > $INFO_FILES/population_map.txt
