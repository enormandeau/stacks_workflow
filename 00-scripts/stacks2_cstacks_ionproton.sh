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

cstacks -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/population_map_for_cstacks.txt -n 1 -p "$NUM_CPU" --disable-gapped
