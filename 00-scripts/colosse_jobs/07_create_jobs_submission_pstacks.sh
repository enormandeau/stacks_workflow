#!/bin/bash

# Global variables
NUMSAMPLES=50
TEMPLATE=./00-scripts/colosse_jobs/stacks_1b_pstacks_template_job.sh

# Cleanup first
rm temp.pstacks.list.* 2>/dev/null
rm 00-scripts/colosse_jobs/PSTACKS_temp.pstacks.list.* 2>/dev/null

# Create list of samples
ls 04-all_samples/*bam | split -a 4 -d -l "$NUMSAMPLES" - temp.pstacks.list.

#create list of jobs
for file_list in $(ls temp.pstacks.list.*)
do
    JOB=00-scripts/colosse_jobs/PSTACKS_"$file_list".sh
    toEval="cat "$TEMPLATE" | sed 's/__LIST__/"$file_list"/g'"
    eval "$toEval" > "$JOB"
    id=$(echo "$file_list" | cut -d "." -f 4)
    echo "$file_list"
    echo $id
    starting=$(echo "1 + $id * $NUMSAMPLES" | bc)
    echo $starting
    sed -i "s/id=__ID__/id=$starting/" "$JOB"
done

exit

#submit jobs
for i in $(ls 00-scripts/colosse_jobs/PSTACKS*sh)
do
    msub "$i"
done
