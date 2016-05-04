#!/bin/bash
# Use number of reads per sample to calculate the needed volume per sample to
# prepare second library in order to normalize read depth per sample

# Iterate over extracted sequencing chips
cat 01-info_files/lane_info.txt |
while read i
do
    # Variables
    numseq=03-samples/"$i".numseq
    infos=03-samples/"$i".infos

    # Remove count files
    rm "$numseq" 2>/dev/null
    rm "$infos" 2>/dev/null

    # Get informations from sample_information.csv for each chip
    grep -E "^$i(\.f.*q\.gz)?\s" 01-info_files/sample_information.csv |
        perl -pe 's/\.f(ast)?q\.gz//' > "$infos"

    # Iterate over samples in a chip
    for j in 03-samples/"$i"/*.f*q.gz
    do
        # Get number of reads per sample per chip
        echo "$j" $(echo $(gunzip -c "$j" | wc -l) / 4 | bc) |
            perl -pe 's/.*sample_//; s/\.f(ast)?q\.gz//' >> "$numseq"

    done
done
