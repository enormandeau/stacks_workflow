#!/bin/bash

# Works with version 0.99995

# Options
# Comment out options that you do not wish to use

b="-b 1"            #b: MySQL ID of this batch.
o="-o 05-stacks"    #o: output path to write results.
#g="-g"              #g: base catalog matching on genomic location, not sequence identity.
#m="-m"
n="-n 0"            #n: number of mismatches allowed between sample tags when
                    #   generating the catalog (default 0).
p="-p 16"            #p: enable parallel execution with num_threads threads.
#report_mmatches="--report_mmatches"


# s: filename prefix from which to load loci into the catalog.
# -=( DO NOT MODIFY THE FOLLWING OPTION! )=-
s="$(for file in $(ls -1 05-stacks/*.tags.tsv | perl -pe 's/\.tags\.tsv//'); do echo -s $file; done)"

cstacks $b $s $o $g $m $n $p

