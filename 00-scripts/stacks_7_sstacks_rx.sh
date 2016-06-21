#!/bin/bash
# Launch sstacks to treat all the samples indivually
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# OPTIONS: Comment out options that you do not wish to use
p="-p 1"                  # p: enable parallel execution with num_threads threads
b="-b 1"                   # b: MySQL ID of this batch
c="-c 06-stacks_rx/batch_1"   # c: TSV file from which to load the catalog loci
o="-o 06-stacks_rx"           # o: output path to write results
#g="-g"                    # g: base matching on genomic location, not sequence
                           #   identity
#x="-x"                    # x: don't verify haplotype of matching locus

#Gapped assembly options:
#gap=" --gapped" 	#preform gapped alignments between stacks.

# Prepare list of samples to treat
s="$(for file in $(ls -1 06-stacks_rx/*.tags.tsv.gz | grep -v catalog | \
    perl -pe 's/\.tags\.tsv\.gz//'); do echo -s $file; done)"

# Launch sstacks on all samples
sstacks $p $b $c $gap $s $o $g $x $v $h 2>&1 | tee 10-log_files/"$TIMESTAMP"_stacks_7_sstacks_rx.log
