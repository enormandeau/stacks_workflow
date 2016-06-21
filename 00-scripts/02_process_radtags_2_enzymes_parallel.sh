#!/bin/bash
# Launch process_radtags on all the lanes
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
NCPU=$4

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"
INFO_FILES_FOLDER="01-info_files"
SAMPLE_INFO="sample_information.csv"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"
cp $INFO_FILES_FOLDER/$SAMPLE_INFO $LOG_FOLDER/"$TIMESTAMP"_"$SAMPLE_INFO"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi


### Global variable
INFO_FILES="01-info_files"
TRIM_LENGTH=$1 # Length to cut reads after process_radtags
ENZYME1=$2 # Name of the enzyme (run 'process_radtags' without options for list)
ENZYME2=$3 # Name of the enzyme (run 'process_radtags' without options for list)

# Write command to file
echo -e "process_radtags command used:\n\n\
    $(echo process_radtags_2_enzymes.sh $TRIM_LENGTH $ENZYME1 $ENZYME2)" \
> 10-log_files/process_radtags_2_enzymes_command.log

# Extract reads
cat $INFO_FILES/lane_info.txt |
    parallel -j $NCPU 00-scripts/utility_scripts/process_radtags_2_enzymes.sh $TRIM_LENGTH $ENZYME1 $ENZYME2 $LANE

