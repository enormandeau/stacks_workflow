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

# Number of CPUs
NUM_CPU="20"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"
cp $INFO_FILES_FOLDER/$POP_MAP $LOG_FOLDER/"$TIMESTAMP"_"$POP_MAP"

# module load gcc/6.2.0
# module load stacks/2.3e

# Running gstacks
# NOTE need --unpaired for single-end?
gstacks -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" -t "$NUM_CPU"
