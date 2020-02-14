#!/bin/bash
#SBATCH -J "cstacks"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=8G

module load stacks/1.48

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Create stacks catalog
./00-scripts/stacks1_2_cstacks.sh
