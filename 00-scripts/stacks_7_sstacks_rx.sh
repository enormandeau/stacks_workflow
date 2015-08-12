#!/bin/bash
# Launch sstacks to treat all the samples indivually
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# OPTIONS: Comment out options that you do not wish to use
p="-p 16"                  # p: enable parallel execution with num_threads
                           #   threads
b="-b 1"                   # b: MySQL ID of this batch
c="-c 05-stacks_rx/batch_1"   # c: TSV file from which to load the catalog loci
o="-o 05-stacks_rx"           # o: output path to write results
#g="-g"                    # g: base matching on genomic location, not sequence
                           #   identity
#x="-x"                    # x: don't verify haplotype of matching locus

# Launch sstacks on all samples
for file in $(ls -1 05-stacks_rx/*.tags.tsv.gz | grep -v catalog | perl -pe 's/\.tags\.tsv\.gz//')
do
    sstacks $p $b $c $r $o $g $x $h -s $file
done 2>&1 | tee 98-log_files/"$TIMESTAMP"_stacks_7_sstacks_rx.log

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

