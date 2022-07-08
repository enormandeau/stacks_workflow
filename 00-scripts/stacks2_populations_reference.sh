#!/bin/bash
# Launch populations
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
SAMPLE_FOLDER="04-all_samples"
STACKS_FOLDER="05-stacks"
LOG_FOLDER="10-log_files"
INFO_FILES_FOLDER="01-info_files"
POP_MAP="population_map.txt"

# Populations does not seem to use more than 1 CPU most of the time
NUM_CPU="10"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"
cp $INFO_FILES_FOLDER/$POP_MAP $LOG_FOLDER/"$TIMESTAMP"_"$POP_MAP"

# module load gcc/6.2.0
# module load stacks/2.3e

populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 2 -r 0.6 \
    --ordered-export --fasta-loci --vcf
