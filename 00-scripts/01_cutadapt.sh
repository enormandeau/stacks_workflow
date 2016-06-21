#!/bin/bash
# Removing adapters using cutadapt
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
NCPU=$1

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi

# Create directory for untrimmed files
mkdir 02-raw/trimmed 2>/dev/null

rm 10-log_files/"$TIMESTAMP"_01_cutadapt"${i%.fastq.gz}".log 2> /dev/null
ls -1 02-raw/*.fastq.gz |
parallel -j $NCPU cutadapt -a file:01-info_files/adapters.fasta \
        -o 02-raw/trimmed/{/} \
        -e 0.2 \
        -m 50 \
        {} '2>&1' '>>' 10-log_files/"$TIMESTAMP"_01_cutadapt"${i%.fastq.gz}".log

