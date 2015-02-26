#!/bin/bash
# Extract information about the lanes in the 02-raw folder

for i in $(ls -1 02-raw/*.fastq.gz)
do
    # TEST:
    # echo ${$i%.fastq}.gz
    basename $i | perl -pe 's/\.fastq\.gz//'
done > 01-info_files/lane_info.txt

