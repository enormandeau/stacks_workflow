#!/bin/bash
#SBATCH -J "ustacks"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH -A ibismini
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=2-00:00
#SBATCH --mem=40000

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Find stacks for each sample
./00-scripts/stacks_1a_ustacks.sh
