#!/bin/bash
# Testing the full stacks_workflow pipeline
# WARNING! Do not use this script to run your analyses!

echo "Mini test should take approximately 1m 15sec to run"

# Clean directories
rm 01-info_files/barcodes.txt 2> /dev/null
rm 01-info_files/lane_info.txt 2> /dev/null
rm 01-info_files/population_map.txt 2> /dev/null
rm 02-raw/*.fastq.gz 2> /dev/null
rm -r 02-raw/trimmed 2> /dev/null
rm -r 03-samples/* 2> /dev/null
rm -r 04-all_samples/* 2> /dev/null
rm -r 05-stacks/* 2> /dev/null
rm -r 06-stacks_rx/* 2> /dev/null
rm -r 10-log_files/* 2> /dev/null

# Get raw data, and sample information
cp -l ~/temp.backup/stacks_workflow_test_data/*_mini_*.fastq.gz 02-raw/data.fastq.gz
cp ~/temp.backup/stacks_workflow_test_data/sample_information.csv 01-info_files

# Preparatory scripts
./00-scripts/00_prepare_lane_info.sh
./00-scripts/01_cutadapt.sh

# Splitting the samples
./00-scripts/02_process_radtags_2_enzymes.sh 70 pstI mspI

# Renaming the samples
./00-scripts/03_rename_samples.sh

# Creating the population map
./00-scripts/04_prepare_population_map.sh

# Running STACKS
./00-scripts/stacks_1a_ustacks.sh
#./00-scripts/stacks_1b_pstacks.sh
./00-scripts/stacks_2_cstacks.sh
./00-scripts/stacks_3_sstacks.sh
./00-scripts/stacks_4_populations.sh
#../../00-scripts/stacks_5a_rxstacks_likelihoods.sh # Not needed for tests
./00-scripts/stacks_5b_rxstacks.sh
./00-scripts/stacks_6_cstacks_rx.sh
./00-scripts/stacks_7_sstacks_rx.sh
./00-scripts/stacks_8_populations_rx.sh

# Filtering
./00-scripts/05_filter_vcf.py \
    -i 06-stacks_rx/batch_1.vcf \
    -o filtered.vcf \
    -c 2 -m 7 -I 4 -l 10 -C 30 \
    -p 70 --use_percent \
    -a 0.05 -A 0.1 \
    -H 0.5 -y 1 \
    -f -0.3 -F 0.1 \
    -s 6 -q


