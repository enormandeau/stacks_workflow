#!/bin/bash
# Launch ustacks to treat all the samples individually
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Global variables
NUM_CPUS=16

# Copy script as it was run
SCRIPT=$0
NAME="ustacks_iteration.sh"
LOG_FOLDER="10-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

maxid=$(ls -1 04-all_samples/*.fq.gz | wc -l)
seq $maxid | parallel -j $NUM_CPUS 00-scripts/ustacks_iteration.sh
