#!/bin/bash
# Launch 'process_radtags' on all the lanes

# Global variable
TRIM_LENGTH=$1 # Length to which the reads are cut after process_radtags
ENZYME=$2      # Name of the enzyme (run 'process_radtags', without options, for list)

# Prepare bacode_lengths.txt from barcodes.txt, 
perl -ne 'chomp; print length."\n"' 01-info_files/barcodes.txt |
    sort -un > 01-info_files/barcode_lengths.txt

# Create barcode files, eg: barcodes_4.txt, barcodes_5.txt ...
cat 01-info_files/barcode_lengths.txt |
    while read b
    do
        perl -sane 'chomp; if (length eq $b) {print $_."\n"}'\
            -- -b=$b 01-info_files/barcodes.txt |
        sort -u > 01-info_files/barcodes_$b".txt"
    done

# Comment 
cat 01-info_files/lane_info.txt |
    while read f
    do
        cat 01-info_files/barcode_lengths.txt |
        while read b
        do
            mkdir 03-samples/$f
            process_radtags -i gzfastq -f 02-raw/$f".fastq.gz" -o \
                03-samples/$f -b 01-info_files/barcodes_"$b".txt \
                -c -q -r -t $TRIM_LENGTH --barcode_dist\
                --filter_illumina 3 -E phred33 -e $ENZYME
        done
    done

