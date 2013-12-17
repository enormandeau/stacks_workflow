#!/bin/bash
# Launch 'process_radtags' on all the lanes

### Global variable
INFO_FILES="01-info_files"
TRIM_LENGTH=$1 # Length to cut reads after process_radtags
ENZYME=$2 # Name of the enzyme (run 'process_radtags' without options for list)

cat $INFO_FILES/lane_info.txt |
while read f
do
    grep -vE "^$" $INFO_FILES/sample_information.csv | \
        grep -v "Barcode" | \
        cut -f 2 > $INFO_FILES/barcodes.txt

    # Prepare bacode_lengths.txt from barcodes.txt 
    perl -ne 'chomp; print length."\n"' $INFO_FILES/barcodes.txt | \
        sort -un > $INFO_FILES/barcode_lengths.txt

    # Create barcode files, eg: barcodes_4.txt, barcodes_5.txt ...
    cat $INFO_FILES/barcode_lengths.txt | \
    while read b
    do
        perl -sane 'chomp; if (length eq $b) {print $_."\n"}'\
            -- -b=$b $INFO_FILES/barcodes.txt | \
        sort -u > $INFO_FILES/barcodes_$b".txt"
    done

    # Extract reads for each individuals using barcodes
    cat $INFO_FILES/barcode_lengths.txt | \
    while read b
    do
        mkdir 03-samples/$f
        process_radtags -i gzfastq -f 02-raw/$f".fastq.gz" -o \
            03-samples/$f -b $INFO_FILES/barcodes_"$b".txt \
            -c -q -r -t $TRIM_LENGTH --barcode_dist 2 \
            --filter_illumina 3 -E phred33 -e $ENZYME
    done

    #rm $INFO_FILES/barcode_lengths.txt 2> /dev/null
    #rm $INFO_FILES/barcodes*.txt 2> /dev/null
done

