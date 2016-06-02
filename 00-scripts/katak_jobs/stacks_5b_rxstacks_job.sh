#!/bin/bash

#SBATCH -D ./ 
#SBATCH --job-name="rxstacks"
#SBATCH -o log-rxstacks.out
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH -A ibismini
#SBATCH --mail-type=ALL
#SBATCH --mail-user=type_your_mail@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=50000

cd $SLURM_SUBMIT_DIR

#./00-scripts/00_prepare_lane_info.sh
#./00-scripts/01_cutadapt.sh
#./00-scripts/02_process_radtags_2_enzymes.sh 80 pstI mspI
#./00-scripts/02_process_radtags.sh
#./00-scripts/03_rename_samples_complex.sh
#./00-scripts/03_rename_samples_simple.sh
#./00-scripts/04_prepare_population_map.sh
#./00-scripts/05_filter_vcf.py
#./00-scripts/stacks_1a_ustacks.sh
#./00-scripts/stacks_1b_pstacks.sh
#./00-scripts/stacks_2_cstacks.sh
#./00-scripts/stacks_3_sstacks.sh
#./00-scripts/stacks_4_populations.sh
#./00-scripts/stacks_5a_rxstacks_likelihoods.sh
./00-scripts/stacks_5b_rxstacks.sh
#./00-scripts/stacks_6_cstacks_rx.sh
#./00-scripts/stacks_7_sstacks_rx.sh
#./00-scripts/stacks_8_populations_rx.sh

