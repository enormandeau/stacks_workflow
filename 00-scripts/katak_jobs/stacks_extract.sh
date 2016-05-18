#!/bin/bash
#$ -N extract
#$ -M your.addresse@service.co
#$ -m beas
#$ -pe smp 1
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
OUTPUT="./05-stacks/filtered"

FILE_REP="05-stacks"
INPUT="batch_1_filtered.vcf"
OUTPUT="05-stacks/batch_1_extracted.vcf"

#run script
./00-scripts/utility_scripts/extract_first_snp.py $FILE_REP/"$INPUT" $OUTPUT 2>&1 | tee 98-log_files/"$TIMESTAMP"_extract.log
