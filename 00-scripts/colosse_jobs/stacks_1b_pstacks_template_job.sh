#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N pstacks__LIST__
#PBS -o log.pstacks__LIST__.out
#PBS -e log.pstacks__LIST__.err
#PBS -l walltime=20:00:00
#PBS -M your.address@here
###PBS -m ea 
#PBS -l nodes=1:ppn=8
#PBS -r n

# Launch pstacks to treat all the samples individually
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT="$0"
NAME="$(basename $0)"
LOG_FOLDER="10-log_files"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Move to job submission directory
cd "$PBS_O_WORKDIR"

# Export STACKS path
export PATH="/rap/ihv-653-ab/software/stacks/current/bin":$PATH

# Options
# Comment out options that you do not wish to use
t="-t bam"            # t: input file Type. Supported types: bowtie, sam, or bam
o="-o 05-stacks"      # o: output path to write results.
#i="-i 1"             # i: SQL ID to insert into the output to identify this
                      #   sample
m="-m 1"              # m: minimum depth of coverage to report a stack
                      #   (default 1)
p="-p 8"              # p: enable parallel execution with num_threads threads
model_type="--model_type snp"
alpha="--alpha 0.05"
#bound_low="--bound_low 0"
#bound_high="--bound_high 1"
#bc_err_freq="--bc_err_freq 1"
 
file_list=__LIST__
# Launch pstacks for all the individuals
id=__ID__
for file in $(cat "$file_list")
do
    echo -e "\n\n##### Treating individual $id: $file\n\n"
    pstacks $t $o $i $m $p $max_locus_stacks $model_type $alpha $bound_low \
        $bound_high $bc_err_freq -f $file -i $id
    id=$(echo $id + 1 | bc)
done 2>&1 #| tee 10-log_files/"$TIMESTAMP"_stacks_1b_pstacks.log
