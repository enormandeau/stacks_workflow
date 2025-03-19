#!/bin/bash
# Reduces population map to only include the samples that were actually present after filters
#  for example samples with too few reads or other quality issues.
#  This enables removal of samples that are in sample information file
#  without disrupting downstream functions
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
INFO_FILES="01-info_files"
SAMPLE_FILES="04-all_samples"

# Create file listing aligned samples that are present in sample folder
basename -a $SAMPLE_FILES/*.sorted.bam | \
    sed 's/\.sorted\.bam//g' - > $INFO_FILES/samplelist_present.txt

# Reporting
echo "Number of samples in the current pop map (full):"  
wc -l $INFO_FILES/population_map.txt

echo "Reducing pop map down to the number of samples present in $SAMPLE_FILES folder: "
wc -l $INFO_FILES/samplelist_present.txt

# Select only the present samples from the population map
grep -f $INFO_FILES/samplelist_present.txt $INFO_FILES/population_map.txt > $INFO_FILES/population_map_present_samples.txt 

# Write over full population map with reduced population map
mv $INFO_FILES/population_map_present_samples.txt $INFO_FILES/population_map.txt

# Clean workspace
rm $INFO_FILES/samplelist_present.txt

# Reporting
echo "Number of samples in the current pop map (reduced):"
wc -l $INFO_FILES/population_map.txt

