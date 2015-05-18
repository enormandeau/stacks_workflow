#!/bin/bash
# Create population map file
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Global variables
INFO_FILES="01-info_files"
ALL_SAMPLES_FOLDER="04-all_samples"

grep -v "#" $INFO_FILES/sample_information.csv | \
    awk '{print $3 "_" $4 "\t" $5}' | \
    sort -u > $INFO_FILES/population_map.txt

#ls -1 $ALL_SAMPLES_FOLDER/ |
#while read i
#do
#    echo -e "$i" |
#    perl -pe 's/\.fq.gz//'
#done > $INFO_FILES/population_map_temp.txt

#cat $INFO_FILES/population_map_temp.txt |
#while read i
#do
#    echo -ne $i"\t"
#    pop=$(echo $i | cut -d "_" -f 1)
#    ind=$(echo $i | cut -d "_" -f 2)
#    echo $(grep -E "[[:space:]]$pop[[:space:]]$ind[[:space:]]" \
#        $INFO_FILES/sample_information.csv | \
#        cut -f 5)
#done > $INFO_FILES/population_map.txt

rm $INFO_FILES/population_map_temp.txt 2> /dev/null

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

