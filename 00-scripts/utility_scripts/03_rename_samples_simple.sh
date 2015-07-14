#!/bin/bash
# Renaming the extracted sample files
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Global variables
INFO_FILES="01-info_files"
SAMPLES_FOLDER="03-samples"
ALL_SAMPLES_FOLDER="04-all_samples"

# Renaming files
cat $INFO_FILES/lane_info.txt |
    while read file
    do
        for s in $(ls -1 $SAMPLES_FOLDER/$file/sample* 2> /dev/null)
        do
            sample=$(basename $s)
            mv $SAMPLES_FOLDER/$file/$sample $SAMPLES_FOLDER/$file/$file"_"$sample
        done
    done

# Linking files from 03-samples in 04-all_samples
cat $INFO_FILES/lane_info.txt |
    while read lane
    do
        for sample_file in $(ls -1 $SAMPLES_FOLDER/$lane/*.fq.gz)
        do
            barcode=$(echo $sample_file | perl -pe 's/^.*_sample_//; s/\.fq\.gz//')
            sample_info=$(grep $lane $INFO_FILES/sample_information.csv | grep -E "[[:space:]]$barcode[[:space:]]")
            population=$(echo $sample_info | cut -d " " -f 3)
            sample_name=$(echo $sample_info | cut -d " " -f 4)
            new_name=$(echo "$population"_"$sample_name".fq.gz)

            #cp -l $sample_file $ALL_SAMPLES_FOLDER/$new_name
            ln $sample_file $ALL_SAMPLES_FOLDER/$new_name
        done
    done

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

