#!/bin/bash
# Launch sstacks to treat all the samples indivually
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# OPTIONS: Comment out options that you do not wish to use
b="-b 1"            # b: Batch ID to examine when exporting from the catalog
#c="-c"             # c: make automated corrections to the data.
P="-P 05-stacks"    # P: output path to write results
t="-t F2"           # t: map type to write. 'CP', 'DH', 'F2', 'BC1' and 'GEN' are the currently supported map types
o="-o joinmap"      # o: output file type to write, 'joinmap', 'onemap', 'rqtl', and 'genomic' are currently supported.
m="-m 5"            # m: specify a minimum stack depth required before exporting a locus in a particular individual.
R="-r 2"            # r: minimum number of progeny required to print a marker.
#s="file SQL"       # s: output a file to import results into an SQL database.
#B="Blacklist"      # B: specify a file containing Blacklisted markers to be excluded from the export.
#W="Whitelist"      # W: specify a file containign Whitelisted markers to include in the export.
#e="renz"           # e: restriction enzyme, required if generating 'genomic' output

#filtering option
log_likelihood="--lnl_lim -10"    #: filter loci with log likelihood values below this threshold.

#Automated corrections options:
#min_homo="--min_hom_seqs 5"      #: minimum number of reads required at a stack to call a homozygous genotype (default 5).
#min_hetero="--min_het_seqs 0.05" #: below this minor allele frequency a stack is called a homozygote, above it (but below --max_het_seqs) it is called unknown (default 0.05).

#max_hetero="--max_het_seqs 0.1"  #: minimum frequency of minor allele to call a heterozygote (default 0.1).

#Manual corrections options:
#man_genotype="--cor_path <path>" #: path to file containing manual genotype corrections from a Stacks SQL database to incorporate into output.

# Launch genotypes on all samples
genotypes $b $P $t $m $R $o \
    $B $W \
    $c $s $e \
    $log_likelihood $min_hom $min_hetero $max_hetero $man_genotype 2>&1 | \
    tee 98-log_files/"$TIMESTAMP"_stacks_9_genotypes.log
