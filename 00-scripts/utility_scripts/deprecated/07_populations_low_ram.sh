#!/bin/bash

# Run populations for each single pop

## Create list of populations to use based on present *tags* files
ls -1 05-stacks/*tags* |
    grep -v catalog |
    cut -d "/" -f 2 |
    cut -d "_" -f 1 |
    sort |
    uniq -c |
    sort -n |
    awk '{print $2,$1}' > 01-info_files/populations_for_low_ram.infos

## Backup population_map file
cp 01-info_files/population_map.txt 01-info_files/population_map.txt.bak.$(date +%Y%m%d_%H%M%S)

## Create output folder for each population's VCF
mkdir 09-vcfs_by_pop 2>/dev/null

## Create population maps for each pop
awk '{print $1}' 01-info_files/populations_for_low_ram.infos |
    while read pop
    do
        grep $pop"_" 01-info_files/population_map.txt > 09-vcfs_by_pop/population_map_pop_"$pop".txt
    done

## Run populations on each population

for i in $(ls -1 09-vcfs_by_pop/population_map_pop_*.txt)
    do
        # Get population name
        pop=$(echo $i | cut -d "/" -f 2 | cut -d "_" -f 4 | cut -d "." -f 1)
        echo "---"
        echo $i
        echo $pop
        echo "---"

        # Copy pop map
        cp $i 09-vcfs_by_pop/population_map.txt

        # Run populations
        ./00-scripts/utility_scripts/stacks_4_populations_by_pop.sh

        # Rename result vcf file to include pop name
        mv 05-stacks/batch_1.vcf 09-vcfs_by_pop/"$i".vcf
    done
