#!/bin/bash
# Launch sstacks to treat all the samples indivually

# OPTIONS: Comment out options that you do not wish to use
p="-p 8"                  # p: enable parallel execution with num_threads
                           #   threads
b="-b 1"                   # b: MySQL ID of this batch
c="-c 05-stacks/batch_1"   # c: TSV file from which to load the catalog loci
#r="-r"                    # r: Load the TSV file of a single sample instead of a catalog
o="-o 05-stacks"           # o: output path to write results
#g="-g"                    # g: base matching on genomic location, not sequence
                           #   identity
#x="-x"                    # x: don't verify haplotype of matching locus

# Launch sstacks on all samples
for file in $(ls -1 05-stacks/*.tags.tsv.gz | grep -v catalog | perl -pe 's/\.tags\.tsv\.gz//')
do
    sstacks $p $b $c $r $o $g $x $h -s $file
done 2>&1 | tee stacks_3_sstacks.log

