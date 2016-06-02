#!/bin/bash

#SBATCH -D ./ 
#SBATCH --job-name="filter"
#SBATCH -o log-filter.out
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
MAF_GLOBAL="-a 0.01"	# global maf
MAF_POP="-A 0.05"		# population maf
MIN_DEPTH="-m 7"		# min depth
MIN_ALLELE_COVERAGE="-c 1"	# min allele
LOG_LIKELIHOOD="-l 6"		# log_likelihood
FIS="-f 0.4"			# fis
MIN_PRES="-p 70"			# min presence
JOKER_POP="-x 1"			# joker by pop
JOKER_HET="-y 1"			# joker for het
HET="-H 0.6"			# het
All_IMB="-I 4"			#allelic imbalance
MAX_NB_SNPS="-s 10"		# max number of snps by loci
MAX_FIS="-F 0.4"		#max fis
#JOKER_FIS="-z 1"		#joker fis

INPUT="-i filtered_no_paralog.vcf"
OUTPUT="-o batch_1_filtered.vcf"

./00-scripts/05_filter_vcf.py -q $INPUT $OUTPUT \
 		$MIN_PRES $MIN_ALLELE_COVERAGE $ALL_IMB $MIN_DEPTH $ALL_IMB $LOG_LIKELIHOOD \
		--use_percent $MAF_GLOBAL $MAF_POP $HET $JOKER_HET \
		$FIS $MAX_FIS $JOKER_POP $MAX_NB_SNPS $JOKER_FIS 2>&1 | tee 98-log_files/"$TIMESTAMP"_filter.log 
