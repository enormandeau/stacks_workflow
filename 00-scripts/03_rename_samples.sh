#!/bin/bash
# Renaming the extracted sample files

# Renaming files
cat 01-info_files/lane_info.txt |
    while read file
    do
        for s in $(ls -1 03-samples/$file/sample*)
        do
            sample=$(basename $s)
            mv 03-samples/$file/$sample 03-samples/$file/$file"_"$sample
        done
    done

# Linking files
cat 01-info_files/lane_info.txt |
    while read file
    do
        for s in $(ls -1 03-samples/$file/*.fq)
        do
            sample=$(basename $s)
            cp -l 03-samples/$file/$sample 04-all_samples
        done
    done

