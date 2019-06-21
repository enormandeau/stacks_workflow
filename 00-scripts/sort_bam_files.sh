#!/bin/bash
# Sort yet unsorted bam files
#
# srun -c 8 --mem 20G -p large --time 21-00:00

# Modules
module load samtools

# Samtools sort
ls -1 *.bam | grep -v sorted | parallel -j 8 samtools sort -o {.}.sorted.bam {}
