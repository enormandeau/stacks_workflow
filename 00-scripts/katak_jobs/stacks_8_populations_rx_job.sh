#!/bin/bash
#SBATCH -J "rxpopulations"
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

# Create population ouputs (vcf...)
./00-scripts/stacks_8_populations_rx.sh
