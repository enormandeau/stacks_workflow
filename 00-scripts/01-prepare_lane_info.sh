#!/bin/bash
# Automatically extract information about sequencing lanes

for i in $(ls -1 02-raw/*.fastq.gz)
do
    basename $i | perl -pe 's/\.fastq\.gz//'
done > 01-info_files/lane_info.txt

