#!/bin/bash
# Launch process_radtags on all the lanes
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"
INFO_FILES_FOLDER="01-info_files"
SAMPLE_INFO="sample_information.csv"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"
cp $INFO_FILES_FOLDER/$SAMPLE_INFO $LOG_FOLDER/"$TIMESTAMP"_"$SAMPLE_INFO"

### Global variable
INFO_FILES="01-info_files"
TRIM_LENGTH=$1 # Length to cut reads after process_radtags
ENZYME1=$2 # Name of the enzyme (run 'process_radtags' without options for list)
ENZYME2=$3 # Name of the enzyme (run 'process_radtags' without options for list)

# Write command to file
echo -e "process_radtags commande used:\n\n\
    $(echo process_radtags_2_enzymes.sh $TRIM_LENGTH $ENZYME1 $ENZYME2)" \
> 10-log_files/process_radtags_2_enzymes_command.log

# Extract reads
cat $INFO_FILES/lane_info.txt |
while read f1
do
    f2=$(echo "$f1" | perl -pe 's/_R1$/_R2/')
    echo "$f1" "$f2"
    grep -vE "^$" $INFO_FILES/sample_information.csv  | \
        grep "$f1" | \
        grep -v "Barcode" | \
        cut -f 2 > $INFO_FILES/barcodes.txt

    # Extract the barcodes of all different lengths at once, STACKS 1.23+
    mkdir 03-samples/$f1 2> /dev/null
    process_radtags \
        -i gzfastq \
        -1 02-raw/trimmed/$f1".fastq.gz" \
        -2 02-raw/trimmed/$f2".fastq.gz" \
        -o 03-samples/$f1 \
        -b $INFO_FILES/barcodes.txt \
        -c -r -t $TRIM_LENGTH \
        -q -s 0 \
        -P \
        --barcode_dist_1 2 \
        -E phred33 \
        --renz_1 $ENZYME1 \
        --renz_2 $ENZYME2

    # Copy log files to ./10-log_files/
    cp 03-samples/$f1/process_radtags.log 10-log_files/"$TIMESTAMP"_02_process_radtags_2_enzymes_"$f1".log
done

