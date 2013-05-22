#!/bin/bash
# Launch sstacks to treat all the samples indivually
# Sets of stacks constructed by the ustacks or pstacks programs can be searched against a catalog produced by cstacks
# In the case of a genetic map, stacks from the progeny would be matched against the catalog to determine which progeny contain which parental alleles
# For large number of samples you can split the running of sstacks on different computers or cluster nodes

# Comment out options that you do not wish to use

p="-p 16"                  # p: enable parallel execution with num_threads threads
b="-b 1"                   # b: MySQL ID of this batch
c="-c 05-stacks/batch_1"   # c: TSV file from which to load the catalog loci
#r="-r"                    # r: Load the TSV file of a single sample instead of a catalog
#s                         # s: all type of TSV file from which to load sample RAD-Tags (included below)
o="-o 05-stacks"           # o: output path to write results
#g="-g"                    # g: base matching on genomic location, not sequence identity
#x="-x"                    # x: don't verify haplotype of matching locus
v="-v"                     # v: printe program version to keep track

#h                         # display help message

# Launch sstacks on all samples

for file in $(ls -1 05-stacks/*sample*.tags.tsv | perl -pe 's/\.tags\.tsv//'); do sstacks $p $b $c $r $o $g $x $v -s $file; done  

