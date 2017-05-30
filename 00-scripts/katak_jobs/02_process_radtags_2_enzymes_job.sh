#!/bin/bash
#SBATCH -J "proc_rad"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismax
#SBATCH -A ibismax
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=2G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Process radtags with one CPU or more (also change CPU number above)
./00-scripts/02_process_radtags_2_enzymes_parallel.sh 80 pstI mspI 1
