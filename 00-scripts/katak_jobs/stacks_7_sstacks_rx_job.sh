#!/bin/bash
#SBATCH -J "sstacks_rx"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH -A ibismini
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@service.com
#SBATCH --time=1-00:00
#SBATCH --mem=40000

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Assign genotypes to individuals for loci in catalog
./00-scripts/stacks_7_sstacks_rx.sh
