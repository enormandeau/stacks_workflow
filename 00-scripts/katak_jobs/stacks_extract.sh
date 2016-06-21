#!/bin/bash
#SBATCH -J "extract"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH -A ibismini
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=50000

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

#global variables
OUTPUT="./05-stacks/filtered"

FILE_REP="05-stacks"
INPUT="batch_1_filtered.vcf"
OUTPUT="05-stacks/batch_1_extracted.vcf"

# Extract first SNP
./00-scripts/utility_scripts/extract_first_snp.py \
    $FILE_REP/"$INPUT" $OUTPUT 2>&1 | tee 10-log_files/"$TIMESTAMP"_extract.log
