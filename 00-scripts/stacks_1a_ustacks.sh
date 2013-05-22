#!/bin/bash
# Launch ustacks to treat all the samples individually
# If you have a large number of samples and have access to several computers or best computer cluster, you can split us tacks among computers/different nodes 

# Comment out options that you do not wish to use

t="-t fastq"      # t: input file Type. Supported types: fasta, fastq, gzfasta, or gzfastq
#f                # f: input file path
o="-o 05-stacks"  # o: output path to write results
i="-i 1"          # i: SQL ID to insert into the output to identify this sample
m="-m 3"          # m: Minimum depth of coverage required to create a stack (default 3)
M="-M 2"          # M: Maximum distance (in nucleotides) allowed between stacks (default 2)
N="-N 4"          # N: Maximum distance allowed to align secondary reads to primary stacks (default: M + 2)
#R="-R"           # R: retain unused reads
#H="-H"           # H: disable calling haplotypes from secondary reads
p="-p 16"         # p: enable parallel execution with num_threads threads

#-h               # h: display help message                

#Stacks assembly options:
r="-r"            # r: enable the Removal algorithm, to drop highly-repetitive stacks (and nearby errors) from the algorithm
d="-d"            # d: enable the Deleveraging algorithm, used for resolving over merged tags
mls="--max_locus_stacks 3"                 # mls: maximum number of stacks at a single de novo locus (default 3)

#Model options:
mt="--model_type snp"               # either 'snp' (default), 'bounded' or 'fixed' model type

#For the SNP or Bounded SNP model:
alpha="--alpha 0.1"                 # chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001

#For the Bounded SNP model:
#bound_low="--bound_low 0"          # lower bound for epsilon, the error rate, between 0 and 1.0 (default 0)
#bound_high="--bound_high 1"        # upper bound for epsilon, the error rate, between 0 and 1.0 (default 1)

#For the Fixed model:
#bc_err_freq="--bc_err_freq 1"      #specify the barcode error frequency, between 0 and 1.0




# Launch ustacks for all the individuals
# ATTENTION: If ustacks is run on several computers or cluster id=1 need to be modified for each run

id=1
for file in 04-all_samples/*.fq; do echo -e "\n\n##### Treating individual $id: $file\n\n"; ustacks $t $d $r $o $i $m $M $N $p $R $H $mls $mt $alpha $bound_low $bound_high $bc_err_freq -f $file -i $id; id=$(echo $id + 1 | bc); done


