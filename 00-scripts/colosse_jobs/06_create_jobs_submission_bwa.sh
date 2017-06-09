#!/bin/bash

# Create list of samples
ls 04-all_samples/*f*q.gz | split -l 10 - temp.list.

#change directory in gsnap template job
i=$(pwd)
sed "s#__PWD__#$i#g" 00-scripts/colosse_jobs/bwa_stacks_template.sh > 00-scripts/colosse_jobs/bwa_stacks.sh

# Create list of jobs
for base in $(ls temp.list.*)
do
    toEval="cat 00-scripts/colosse_jobs/bwa_stacks.sh | sed 's/__LIST__/$base/g'"
    eval $toEval > 00-scripts/colosse_jobs/BWA_$base.sh
done

# Submit jobs
for i in $(ls 00-scripts/colosse_jobs/BWA_*.sh)
do
    msub $i
done
