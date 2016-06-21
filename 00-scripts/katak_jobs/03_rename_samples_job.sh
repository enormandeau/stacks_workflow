#!/bin/bash
#SBATCH -J "rename"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH -A ibismini
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=2000

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Process radtags with one CPU or more (also change CPU number above)
./00-scripts/03_rename_samples.sh
