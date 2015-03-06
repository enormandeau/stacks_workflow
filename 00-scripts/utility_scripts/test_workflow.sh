#!/bin/bash
# Testing the full stacks_workflow pipeline
# WARNING! Do not use this script to run your analyses!

# TODO - Add timestamp to logs
# TODO - Add reporting

# Clean directories
rm 01-info_files/barcodes.txt 2> /dev/null
rm 01-info_files/lane_info.txt 2> /dev/null
rm 01-info_files/population_map.txt 2> /dev/null
rm 02-raw/*.fastq.gz 2> /dev/null
rm -r 03-samples/* 2> /dev/null
rm -r 04-all_samples/* 2> /dev/null
rm -r 05-stacks/* 2> /dev/null
rm -r 05-stacks_rx/* 2> /dev/null

# TODO remove when timestamps will be implemented
rm -r 98-log_files/* 2> /dev/null

# Get raw data and sample information
cp -l ~/temp.backup/stacks_workflow_test_data/*.fastq.gz 02-raw
cp ~/temp.backup/stacks_workflow_test_data/sample_information.csv 01-info_files

# Preparatory scripts
./00-scripts/01_prepare_lane_info.sh
./00-scripts/02_process_radtags_2_enzymes.sh 70 pstI mspI
./00-scripts/03_rename_samples.sh
./00-scripts/04_prepare_population_map.sh

# STACKS executables
./00-scripts/stacks_1a_ustacks.sh
./00-scripts/stacks_1b_pstacks.sh
./00-scripts/stacks_2_cstacks.sh
./00-scripts/stacks_3_sstacks.sh
./00-scripts/stacks_4_populations.sh
#../../00-scripts/stacks_5a_rxstacks_likelihoods.sh
./00-scripts/stacks_5b_rxstacks.sh
./00-scripts/stacks_6_cstacks_rx.sh
./00-scripts/stacks_7_sstacks_rx.sh
./00-scripts/stacks_8_populations_rx.sh

