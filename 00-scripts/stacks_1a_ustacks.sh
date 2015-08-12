#!/bin/bash
# Launch ustacks to treat all the samples individually
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# OPTIONS: Comment out options that you do not wish to use
t="-t gzfastq"    # t: input file Type. Supported types: fasta, fastq, gzfasta,
                  #   or gzfastq
o="-o 05-stacks"  # o: output path to write results.
#i="-i 1"         # i: SQL ID to insert into the output to identify this sample
m="-m 4"          # m: Minimum depth of coverage required to create a stack
                  #   (default 3).
M="-M 2"          # M: Maximum distance (in nucleotides) allowed between stacks
                  #   (default 2).
N="-N 4"          # N: Maximum distance allowed to align secondary reads to
                  #   primary stacks (default: M + 2).
#R="-R"           # R: retain unused reads.
H="-H"            # H: disable calling haplotypes from secondary reads.
p="-p 16"         # p: enable parallel execution with num_threads threads.
r="-r"            # r: enable the Removal algorithm, to drop highly-repetitive
                  #   stacks (and nearby errors) from the algorithm.
d="-d"            # d: enable the Deleveraging algorithm, used for resolving
                  #   over merged tags.
max_locus_stacks="--max_locus_stacks 2"
#k_len="--k_len 31" # --k_len: specify k-mer size for matching between
                    # alleles and loci (automatically calculated by default).
model_type="--model_type bounded" #--model_type: either 'snp' (default), 'bounded', or 'fixed'

#alpha="--alpha 0.05"
bound_low="--bound_low 0"
bound_high="--bound_high 1"
#bc_err_freq="--bc_err_freq 1"

# Launch ustacks for all the individuals
id=1
for file in 04-all_samples/*.fq.gz
do
    echo -e "\n\n##### Treating individual $id: $file\n\n"
    ustacks $t $o $i $m $M $N $R $H $p $r $d $max_locus_stacks $k_len \
        $model_type $alpha $bound_low $bound_high $bc_err_freq -f $file -i $id
    id=$(echo $id + 1 | bc)
done 2>&1 | tee 98-log_files/"$TIMESTAMP"_stacks_1a_ustacks.log

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

