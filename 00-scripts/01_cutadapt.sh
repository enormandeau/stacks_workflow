#!/bin/bash
# Removing adapters using cutadapt
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Create directory for untrimmed files
mkdir 02-raw/trimmed 2>/dev/null

for i in $(ls -1 02-raw/*.fastq.gz)
do
    name=$(basename $i)
    echo "--- Cutadapt ---"
    echo "Name: $name"
    echo
    cutadapt -a file:01-info_files/adapters.fasta \
        -o 02-raw/trimmed/"$name" \
        -e 0.2 \
        -m 50 \
        "$i" 2>&1 | tee 98-log_files/"$TIMESTAMP"_01_cutadapt"${i%.fastq.gz}".log
done

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

