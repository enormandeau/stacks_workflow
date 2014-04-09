#!/bin/bash
# Renaming the extracted sample files

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

# Linking files
cat $INFO_FILES/lane_info.txt |
    while read lane
    do
        for sample_file in $(ls -1 $SAMPLES_FOLDER/$lane/*.fq)
        do
            barcode=$(echo $sample_file | perl -pe 's/^.*_sample_//; s/\.fq//')
            sample_info=$(grep $lane $INFO_FILES/sample_information.csv | grep -E "[[:space:]]$barcode[[:space:]]")
            population=$(echo $sample_info | cut -d " " -f 3)
            sample_name=$(echo $sample_info | cut -d " " -f 4)
            new_name=$(echo "$population"_"$sample_name".fq)

            #cp -l $sample_file $ALL_SAMPLES_FOLDER/$new_name
            ln $sample_file $ALL_SAMPLES_FOLDER/$new_name
        done
    done

