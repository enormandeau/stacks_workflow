#!/bin/bash
#SBATCH -J "cstacks_rx"
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

# Create stacks catalog
./00-scripts/stacks_6_cstacks_rx.sh
