#!/bin/bash
#SBATCH -J "filter"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismax
#SBATCH -A ibismax
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=10G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

#global variables
INPUT="-i filtered_no_paralog.vcf"
OUTPUT="-o batch_1_filtered.vcf"

MIN_ALLELE_COVERAGE="-c 1"	# min allele
MIN_DEPTH="-m 7"		    # min depth

LOG_LIKELIHOOD="-l 6"		# log_likelihood
All_IMB="-I 4"			    # allelic imbalance

MIN_PRES="-p 70"			# min presence
USE_PERCENT="--use_percent" # min presence in percents
JOKER_POP="-x 0"			# min presence joker num pops

MAF_GLOBAL="-a 0.01"	    # global maf
MAF_POP="-A 0.05"		    # population maf

HET="-H 0.6"			    # heterozygozity
JOKER_HET="-y 0"			# joker for heterozygozity num pops

FIS="-f 0.4"			    # fis
MAX_FIS="-F 0.4"		    # max fis
JOKER_FIS="-z 0"		    # joker for fis num pops

MAX_NB_SNPS="-s 10"		    # max number of snps by loci

# Filter
./00-scripts/05_filter_vcf.py -q $INPUT $OUTPUT \
    $MIN_ALLELE_COVERAGE $MIN_DEPTH $LOG_LIKELIHOOD $ALL_IMB \
    $MIN_PRES $USE_PERCENT $JOKER_POP $MAF_GLOBAL $MAF_POP \
    $HET $JOKER_HET $FIS $MAX_FIS $JOKER_FIS $MAX_NB_SNPS 2>&1 |
    tee 10-log_files/"$TIMESTAMP"_filter.log
