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
NUM_CPU="60"

# Gnu Parallel version
ls -1 -S "$SAMPLE_FOLDER"/*.gz |
    parallel -j "$NUM_CPU" ustacks -f {} -o "$STACKS_FOLDER" -i {#} -p 1 -m 3 -M 3 -N 5 --disable-gapped # -H --deleverage
