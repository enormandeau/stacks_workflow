#!/bin/bash
# Launch sstacks to treat all the samples indivually
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# OPTIONS: Comment out options that you do not wish to use
p="-p 16"                  # p: enable parallel execution with num_threads
                           #   threads
b="-b 1"                   # b: MySQL ID of this batch
c="-c 05-stacks/batch_1"   # c: TSV file from which to load the catalog loci
o="-o 05-stacks"           # o: output path to write results
#g="-g"                    # base matching on genomic location, not sequence identity.
#x="-x"                    # don't verify haplotype of matching locus.
#v="-v"                    # print program version.
#h="-h"                    # display this help messsage.


# -=( DO NOT MODIFY THE FOLLWING OPTION! )=-
# This will automatically create the list of filenames for sstacks
# s: filename prefix from which to load sample loci.
s="$(for file in $(ls -1 05-stacks/*.tags.tsv.gz | grep -v catalog | perl -pe 's/\.tags\.tsv\.gz//'); do echo -s $file; done)"
sstacks $p $b $c $s $o $g $x $v $h 2>&1 | tee 98-log_files/"$TIMESTAMP"_stacks_3_sstacks.log

