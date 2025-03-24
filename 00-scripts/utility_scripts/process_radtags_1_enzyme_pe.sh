#!/bin/bash

### Global variable
INFO_FILES="01-info_files"
TRIM_LENGTH=$1 # Length to cut reads after process_radtags
ENZYME1=$2 # Name of the enzyme (run 'process_radtags' without options for list)
LANE=$3
TEMPBARCODES=.temp.${LANE}.barcodes

grep -vE "^$" "$INFO_FILES"/sample_information.csv  | \
    grep "$LANE" | \
    grep -v "Barcode" | \
    cut -f 2 > "$INFO_FILES"/"$TEMPBARCODES"

# Extract the barcodes of all different lengths at once, STACKS 1.23+
mkdir 03-samples/"$LANE" 2> /dev/null

process_radtags \
    -i gzfastq \
    -1 02-raw/trimmed/"$LANE"".fastq.gz" \
    -2 02-raw/trimmed/$(echo "$LANE" | perl -pe 's/_R1/_R2/')".fastq.gz" \
    -o 03-samples/"$LANE" \
    -b "$INFO_FILES"/"$TEMPBARCODES" \
    -c -r -t "$TRIM_LENGTH" \
    -q -s 0 \
    --barcode_dist_1 2 \
    -E phred33 \
    --renz_1 "$ENZYME1" 

# Copy log files to ./10-log_files/
cp 03-samples/"$LANE"/process_radtags.log 10-log_files/"$TIMESTAMP"_02_process_radtags_2_enzymes_"$LANE".log

# Cleanup temp files
#rm "$INFO_FILES"/"$TEMPBARCODES"
