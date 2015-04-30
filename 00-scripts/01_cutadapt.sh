#!/bin/bash
# Removing adapters using cutadapt

for i in $(ls -1 02-raw/*.fastq.gz)
do
    echo "--- Cutadapt ---"
    echo
    cutadapt -a file:01-info_files/adapters.fasta \
        -o "$i".trimmed.fastq.gz \
        -e 0.2 \
        -m 50 \
        "$i" 2>&1 | tee 98-log_files/01_cutadapt"${i%.fastq.gz}".log
done

