#!/bin/bash
#SBATCH -J "rxsstacks"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismax
#SBATCH -A ibismax
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=6-00:00
#SBATCH --mem=30G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Assign genotypes to individuals for loci in catalog
./00-scripts/stacks_7_sstacks_rx.sh
