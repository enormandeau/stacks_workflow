#!/bin/bash
# Removing adapters using cutadapt

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
        "$i" 2>&1 | tee 98-log_files/01_cutadapt"${i%.fastq.gz}".log
done

