#!/bin/bash
# Lauch cstacks on all the .tags.tsv files

# Options
# Comment out options that you do not wish to use

b="-b 1"            # b: MySQL ID of this batch
o="-o 05-stacks_rx" # o: output path to write results
#g="-g"             # g: base catalog matching on genomic location, not sequence
                    #   identity
#m="-m"             # m: include tags in catalog that match more than one entry
n="-n 1"            # n: number of mismatches allowed between sample tags when
                    #   generating the catalog (default 0)
p="-p 8"           # p: enable parallel execution with num_threads threads
#catalog="--catalog PATH"
#report_mmatches="--report_mmatches"    # --report_mmatches: report query loci
                    #   that match more than one catalog locus


# -=( DO NOT MODIFY THE FOLLWING OPTION! )=-
# This will automatically create the list of filenames for cstacks
# s: filename prefix from which to load loci into the catalog
s="$(for file in $(ls -1 05-stacks_rx/*.tags.tsv.gz | perl -pe 's/\.tags\.tsv\.gz//'); do echo -s $file; done)"

# Run cstacks
cstacks $b $s $o $g $m $n $p $catalog $report_mmatches 2>&1 | tee stacks_6_cstacks_rx.log

