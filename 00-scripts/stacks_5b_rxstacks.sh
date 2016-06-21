#!/bin/bash
# Launch rxstacks 
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# OPTIONS: Comment out options that you do not wish to use
b="-b 1"           #Batch ID to examine when exporting from the catalog.
P="-P 05-stacks"   #path to the Stacks output files.
o="-o 06-stacks_rx"    #output path to write results.
t="-t 1"           #number of threads to run in parallel sections of code.
#v="-v"             #print program version.
#h="-h"             #display this help messsage.

###Filtering options:
lnl_filter="--lnl_filter"     #filter catalog loci based on the mean
                              #log likelihood 
                              #of the catalog locus in the population.
lnl_lim="--lnl_lim -10"       #minimum log likelihood required to keep 
                              #a catalog locus.
#lnl_dist="--lnl_dist"         #print distribution of mean log likelihoods for 
                              #catalog loci.
conf_filter="--conf_filter"   #filter confounded loci.
conf_lim="--conf_lim 0.75"    #between 0.0 and 1.0 (default 0.75), proportion 
                              #of loci in population that must be confounded 
                              #relative to the catalog locus.
prune_haplo="--prune_haplo"   #prune out non-biological haplotypes unlikely to 
                              #occur in the population.
#max_haplo_cnt="--max_haplo_cnt <limit>" #only consider haplotypes for
                                         #pruning if they occur in fewer
                                         #than max_haplo_cnt samples.

###Model options:
model_type="--model_type bounded" #either 'snp' (default), 'bounded', or 
                                  #'fixed'

###For the SNP or Bounded SNP model:
#alpha="--alpha 0.1"  #chi square significance level required to call a 
                     #heterozygote or homozygote, either 0.1 (default), 0.05, 
                     #0.01, or 0.001.

###For the Bounded SNP model:
bound_low="--bound_low 0"      #lower bound for epsilon, the error rate,
                               #between 0 and 1.0 (default 0).
bound_high="--bound_high 1" #upper bound for epsilon, the error rate,
                               #between 0 and 1.0 (default 1).

###Logging Options:
#verbose="--verbose"  #extended logging, including coordinates of all changed 
                     #nucleotides (forces single-threaded execution).

rxstacks $b $P $o $t $v $h $lnl_filter $lnl_lim $lnl_dist $conf_filter $conf_lim $prune_haplo $max_haplo_cnt $model_type $alpha $bound_low $bound_high $verbose | tee 10-log_files/"$TIMESTAMP"_stacks_5b_rxstacks.log

