#!/bin/bash
#SBATCH -J "pstacks"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=7-00:00
#SBATCH --mem=10G

module load stacks/1.48

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Find stacks for each sample
./00-scripts/stacks1_1b_pstacks.sh
