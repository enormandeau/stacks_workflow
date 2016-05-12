#!/bin/bash

# Create list of samples
ls 04-all_samples/*f*q.gz | split -l 30 - list

#change drectory in gsnap template job
i=$(pwd)
sed "s#__PWD__#$i#g" 00-scripts/colosse_jobs/gsnap_stacks_template.sh > 00-scripts/colosse_jobs/gsnap_stacks.sh

#create list of jobs
for base in $(ls list*)
do
    toEval="cat 00-scripts/colosse_jobs/gsnap_stacks.sh | sed 's/__LIST__/$base/g'"
    eval $toEval > 00-scripts/colosse_jobs/GSNAP_$base.sh
done

exit

#submit jobs
for i in $(ls 00-scripts/colosse_jobs/GSNAP*sh)
do
    msub $i
done
