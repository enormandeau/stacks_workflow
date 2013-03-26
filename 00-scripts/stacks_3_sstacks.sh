#!/bin/bash

# Works with version 0.99995

# Options
# Comment out options that you do not wish to use

b="-b 1"                    #b: MySQL ID of this batch.
c="-c 05-stacks/batch_1"    #c: TSV file from which to load the catalog loci.
o="-o 05-stacks"            #o: output path to write results.
p="-p 16"                   #p: enable parallel execution with num_threads threads.
#g="-g"                      #g: base matching on genomic location, not sequence identity.
#x="-x"                      #x: don't verify haplotype of matching locus.

for file in $(ls -1 05-stacks/*sample*.tags.tsv | perl -pe 's/\.tags\.tsv//'); do sstacks $b $c $o $p $g $x $h -s $file; done  

