#!/bin/bash

for i in $(ls -1 04-all_samples/*.fq)
do
    name=$(basename $i)
    bwa aln -n 5 -k 3 -t 10 -f 04-all_samples/$name.aligned ./01-info_files/genome 04-all_samples/$name
    bwa samse -r "@RG\tID:"$name"\tSM:"$name"\tPL:Illumina" ./01-info_files/genome 04-all_samples/$name.aligned ./04-all_samples/$name > 04-all_samples/$name.sam

done

