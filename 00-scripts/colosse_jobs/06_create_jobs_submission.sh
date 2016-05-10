#!/bin/bash



# Create list of samples
ls 04-all_samples/*f*q.gz|split -l 30 - list


#change drectory

for i in $(pwd); do sed -i "s#__PWD__#$i#g" 00-scripts/colosse_jobs/gsnap_stacks_template.sh; done

#create list of jobs
for base in $(ls list*); do toEval="cat 00-scripts/colosse_jobs/gsnap_stacks_template.sh| sed 's/__LIST__/$base/g'"; eval $toEval > 00-scripts/colosse_jobs/GSNAP_$base.sh;done

#submit jobs
for i in $(ls 00-scripts/colosse_jobs/GSNAP*sh); do msub $i; done
