#!/bin/bash

#SBATCH -D ./ 
#SBATCH --job-name="extract"
#SBATCH -o log-extract.out
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH -A ibismini
#SBATCH --mail-type=ALL
#SBATCH --mail-user=type_your_mail@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=50000

cd $SLURM_SUBMIT_DIR

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

#global variables
OUTPUT="./05-stacks/filtered"

FILE_REP="05-stacks"
INPUT="batch_1_filtered.vcf"
OUTPUT="05-stacks/batch_1_extracted.vcf"

#run script
./00-scripts/utility_scripts/extract_first_snp.py $FILE_REP/"$INPUT" $OUTPUT 2>&1 | tee 98-log_files/"$TIMESTAMP"_extract.log
