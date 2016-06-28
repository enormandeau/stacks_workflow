#!/bin/bash
#SBATCH -J "paralogs"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH -A ibismini
#SBATCH --mail-type=END,FAIL
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
input_vcf="extracted_temp.vcf"
paralog_info="blacklist.loci.paralogs.ind.filtered.txt"
threshold=6
output_vcf="filtered_no_paralog.vcf"

./00-scripts/utility_scripts/vcf_remove_paralogs.py $input_vcf $paralog_info \
    $threshold $output_vcf 2>&1 | tee 10-log_files/"$TIMESTAMP"_paralogs.log
