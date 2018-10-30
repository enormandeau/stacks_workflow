#!/bin/bash
#SBATCH -J "cutadapt"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=2G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Load cutadapt module
module load cutadapt

# Reduce number of CPUs here and above if you have less than 4 chips/lanes
./00-scripts/01_cutadapt.sh 1
