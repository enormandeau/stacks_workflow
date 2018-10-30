#!/bin/bash
#SBATCH -J "norm2"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=0-01:00
#SBATCH --mem=1G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Create stacks catalog
./00-scripts/sequencing_normalization_02.py 5000000 40000 480
