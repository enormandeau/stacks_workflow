#!/bin/bash
#$ -N filter_paralogs
#$ -M your.addresse@service.com
#$ -m beas
#$ -pe smp 8
#$ -l h_vmem=40G
#$ -l h_rt=30:00:00
#$ -cwd
#$ -S /bin/bash


TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"


#move work dir
cd $SGE_O_WORKDIR

#global variables
input_vcf="extracted_temp.vcf"
paralog_info="blacklist.loci.paralogs.ind.filtered.txt"
threshold=6
output_vcf="filtered_no_paralog.vcf"

./00-scripts/utility_scripts/vcf_remove_paralogs.py $input_vcf $paralog_info $threshold $output_vcf 2>&1 | tee 98-log_files/"$TIMESTAMP"_paralogs.log
