#!/bin/bash
#SBATCH -J "norm2"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismax
#SBATCH -A ibismax
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=0-01:00
#SBATCH --mem=4000

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Create stacks catalog
./00-scripts/sequencing_normalization_02.py 2000000 20000 480
