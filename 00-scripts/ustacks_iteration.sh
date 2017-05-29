#!/bin/bash

# Global variables
id=$[ $1 - 1 ]
files=(04-all_samples/*.fq.gz)
file=${files[$id]}

#file=$(ls -1 04-all_samples/*.fq.gz | head -n $id | tail -n 1)

# OPTIONS: Comment out options that you do not wish to use
t="-t gzfastq"    # t: input file Type. Supported types: fasta, fastq, gzfasta,
                  #   or gzfastq
o="-o 05-stacks"  # o: output path to write results.
#i="-i 1"         # i: SQL ID to insert into the output to identify this sample
m="-m 4"          # m: Minimum depth of coverage required to create a stack
                  #   (default 3).
M="-M 5"          # M: Maximum distance (in nucleotides) allowed between stacks
                  #   (default 2).
N="-N 7"          # N: Maximum distance allowed to align secondary reads to
                  #   primary stacks (default: M + 2).
#R="-R"           # R: retain unused reads.
H="-H"            # H: disable calling haplotypes from secondary reads.
p="-p 1"         # p: enable parallel execution with num_threads threads.
r="-r"            # r: enable the Removal algorithm, to drop highly-repetitive
                  #   stacks (and nearby errors) from the algorithm.
d="-d"            # d: enable the Deleveraging algorithm, used for resolving
                  #   over merged tags.
max_locus_stacks="--max_locus_stacks 2"
#k_len="--k_len 31" # --k_len: specify k-mer size for matching between
                    # alleles and loci (automatically calculated by default).
model_type="--model_type bounded" #--model_type: either 'snp' (default), 'bounded', or 'fixed'

#Gapped assembly options:
#gap="--gapped"  #preform gapped alignments between stacks.
#maxgap="--max_gaps 2"	# number of gaps allowed between stacks before merging (default: 2). 
#minallen="--min_aln_len 0.80"	#minimum length of aligned sequence in a gapped alignment (default: 0.80).

#alpha="--alpha 0.05"
bound_low="--bound_low 0"
bound_high="--bound_high 1"
#bc_err_freq="--bc_err_freq 1"

# Run ustacks
echo -e "\n\n##### Treating individual $id: $file\n\n"
sample_id=$[ $id + 1 ]
ustacks $t $o $i $m $M $N $R $H $p $r $d $max_locus_stacks $k_len \
    $gap $maxgap $minallen \
    $model_type $alpha $bound_low $bound_high $bc_err_freq -f $file -i $sample_id
