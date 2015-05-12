#!/bin/bash
# Renaming the extracted sample files

# Global variables
INFO_FILES="01-info_files"
SAMPLES_FOLDER="03-samples"
ALL_SAMPLES_FOLDER="04-all_samples"

# Renaming files
cat 01-info_files/sample_information.csv |
    grep -vE '^#' |
    cut -f 1-4 |
    perl -pe 's/\.fastq\.gz//' |
    perl -pe 's/\t/\/sample_/' |
    perl -pe 's/([ACTG]+)\t/\1.fq.gz\t/' |
    perl -pe 's/(\w+)\t(\S+)$/04-all_samples\/\1_\2.fq.gz/' > renaming_01.txt

cut -f 2 renaming_01.txt | sort -u > renaming_02.txt

cat renaming_02.txt |
    while read i
    do
        echo Treating: "$i"
        rm $ALL_SAMPLES_FOLDER/"$i" 2> /dev/null
        grep "$i" renaming_01.txt |
            cut -f 1 |
            while read j
            do
                echo "  Copying: 03-samples/""$j"
                cat 03-samples/"$j" >> "$i"
            done
    done | tee 98-log_files/03_rename_samples_complex.log

rm renaming_01.txt renaming_02.txt 2> /dev/null

