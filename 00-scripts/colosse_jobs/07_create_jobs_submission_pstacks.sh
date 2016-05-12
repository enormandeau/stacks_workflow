#!/bin/bash

# Create list of samples
ls 04-all_samples/*bam | split -l 30 - list

#give proper direction
for i in $(pwd); do sed -i "s#__PWD__#$i#g" 00-scripts/colosse_jobs/stacks_1b_template_colosse.sh; done

#create list of jobs
for base in $(ls list*); do toEval="cat 00-scripts/colosse_jobs/stacks_1b_template_colosse.sh | sed 's/__LIST__/$base/g'"; eval $toEval > 00-scripts/colosse_jobs/PSTACKS_$base.sh;done

#give proper direction
sed -i 's#__PWD__#/rap/ihv-653-ab/jeremy_leluyer/jeremy_leluyer/04-stacks_workflow_2016-04-14#g' 00-scripts/colosse_jobs/stacks_1b_template_colosse.sh

#change SQL id
sed -i 's/id=1/id=31/g' 00-scripts/colosse_jobs/PSTACKS_listab.sh
sed -i 's/id=1/id=61/g' 00-scripts/colosse_jobs/PSTACKS_listac.sh
sed -i 's/id=1/id=91/g' 00-scripts/colosse_jobs/PSTACKS_listad.sh
sed -i 's/id=1/id=121/g' 00-scripts/colosse_jobs/PSTACKS_listae.sh
sed -i 's/id=1/id=151/g' 00-scripts/colosse_jobs/PSTACKS_listaf.sh
sed -i 's/id=1/id=181/g' 00-scripts/colosse_jobs/PSTACKS_listag.sh

exit
#submit jobs

for i in $(ls 00-scripts/colosse_jobs/PSTACKS*sh); do msub $i; done
