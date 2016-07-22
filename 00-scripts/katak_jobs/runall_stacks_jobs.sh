#!/bin/bash
# Launch all katak jobs with dependencies
# WARNING! All scripts must be edited prior to submitting and not changed until run!
# HINT: Give meaningful job names in submission scripts to distinguish them

# Global variables
DEPENDS="--dependency=afterok:"
SCRIPTPATH="./00-scripts/katak_jobs"

# Processing the reads and normalization
p1=$(sbatch "$SCRIPTPATH"/01_cutadapt_job.sh                                | awk '{print $4}')
p2=$(sbatch "$DEPENDS"$p1 "$SCRIPTPATH"/02_process_radtags_2_enzymes_job.sh | awk '{print $4}')
p3=$(sbatch "$DEPENDS"$p2 "$SCRIPTPATH"/03_rename_samples_job.sh            | awk '{print $4}')
p4=$(sbatch "$DEPENDS"$p2 ./00-scripts/sequencing_normalization_01.sh       | awk '{print $4}')
p5=$(sbatch "$DEPENDS"$p4 ./00-scripts/sequencing_normalization_02.py       | awk '{print $4}')

# STACKS first half (before rxstacks)
s1=$(sbatch "$DEPENDS"$p3 "$SCRIPTPATH"/stacks_1a_ustacks_job.sh            | awk '{print $4}')
s2=$(sbatch "$DEPENDS"$s1 "$SCRIPTPATH"/stacks_2_cstacks_job.sh             | awk '{print $4}')
s3=$(sbatch "$DEPENDS"$s2 "$SCRIPTPATH"/stacks_3_sstacks_job.sh             | awk '{print $4}')
s4=$(sbatch "$DEPENDS"$s3 "$SCRIPTPATH"/stacks_4_populations_job.sh         | awk '{print $4}')

# STACKS second half (rxstacks)
s5=$(sbatch "$DEPENDS"$s3 "$SCRIPTPATH"/stacks_5b_rxstacks_job.sh           | awk '{print $4}')
s6=$(sbatch "$DEPENDS"$s5 "$SCRIPTPATH"/stacks_6_cstacks_rx_job.sh          | awk '{print $4}')
s7=$(sbatch "$DEPENDS"$s6 "$SCRIPTPATH"/stacks_7_sstacks_rx_job.sh          | awk '{print $4}')
s8=$(sbatch "$DEPENDS"$s7 "$SCRIPTPATH"/stacks_8_populations_rx_job.sh      | awk '{print $4}')
