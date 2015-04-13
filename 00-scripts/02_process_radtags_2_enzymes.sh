#!/bin/bash
# Launch process_radtags on all the lanes

### Global variable
INFO_FILES="01-info_files"
TRIM_LENGTH=$1 # Length to cut reads after process_radtags
ENZYME1=$2 # Name of the enzyme (run 'process_radtags' without options for list)
ENZYME2=$3 # Name of the enzyme (run 'process_radtags' without options for list)

# Write command to file
echo -e "process_radtags commande used:\n\n\
$(echo process_radtags_2_enzymes.sh $TRIM_LENGTH $ENZYME1 $ENZYME2)" \
> 98-log_files/process_radtags_2_enzymes_command.log

# Extract reads
cat $INFO_FILES/lane_info.txt |
while read f
do
    grep -vE "^$" $INFO_FILES/sample_information.csv  | \
        grep "$f" | \
        grep -v "Barcode" | \
        cut -f 2 > $INFO_FILES/barcodes.txt

    # Extract the barcodes of all different lengths at once, STACKS 1.23+
    mkdir 03-samples/$f 2> /dev/null
    process_radtags \
        -i gzfastq \
        -f 02-raw/$f".fastq.gz" \
        -o 03-samples/$f \
        -b $INFO_FILES/barcodes.txt \
        -c -q -r -t $TRIM_LENGTH \
        --barcode_dist_1 2 \
        -E phred33 \
        --renz_1 $ENZYME1 \
        --renz_2 $ENZYME2

    # Copy log files to ./98_log_files/
    cp 03-samples/$f/process_radtags.log 98-log_files/process_radtags_2_enzymes_"$f".log
done

