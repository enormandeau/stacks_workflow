#!/bin/bash
# Lauch cstacks on all the .tsv files

# Comment out options that you do not wish to use

b="-b 1"            # b: MySQL ID of this batch
#g="-g"             # g: base catalog matching on genomic location, not sequence identity
#m="-m"             # m: include tags in catalog that match more than one entry
n="-n 0"            # n: number of mismatches allowed between sample tags when generating the catalog (default 0)            
o="-o 05-stacks"    # o: output path to write results

p="-p 16"           # p: enable parallel execution with num_threads threads

#s                  # TSV file from which to load radtags into the catalog
#S                  # MySQL ID of the sample (for web interface)

#Advanced options:
#rm="--report_mmatches="                #report query loci that match more than one catalog locus


# DO NOT MODIFY THE FOLLWING OPTION UNLESS CSTACKS IS NOT RUN CONTINUOUSLY
# Automatically create the list of filenames for cstacks 's' and append a MySQL ID for that sample 'S'

id=1
s="$(for file in $(ls -1 /media/laura/stacks_workflow-map/05-stacks/*.tags.tsv | perl -pe 's/\.tags\.tsv//'); do echo -s $file - S $id; id=$(echo $id + 1 | bc); done)"

# Run cstacks
cstacks $b $s $o $g $m $n $p $rm

