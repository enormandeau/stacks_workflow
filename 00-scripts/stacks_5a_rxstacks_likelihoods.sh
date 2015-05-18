#!/bin/bash
# Run rxstacks to extract log likelihoods and determine cutoff
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# OPTIONS: Comment out options that you do not wish to use
# rxstacks 1.20
# rxstacks -b batch_id -P path [-o path] [-t threads] [-v] [-h]
b="-b 1"              #  b: Batch ID to examine when exporting from the catalog.
P="-P 05-stacks"      #  P: path to the Stacks output files.
o="-o 05-stacks_rx"   #  o: output path to write results.
t="-t 16"              #  t: number of threads to run in parallel sections of code.
#v="-v"               #  v: print program version.
#h="-h"               #  h: display this help messsage.

#  Filtering options:
#    --lnl_filter: filter catalog loci based on the mean log likelihood of the
#      catalog locus in the population.
#    --lnl_lim <limit>: minimum log likelihood required to keep a catalog locus.
#    --lnl_dist: print distribution of mean log likelihoods for catalog loci.
#    --conf_filter: filter confounded loci.
#    --conf_lim <limit>: between 0.0 and 1.0 (default 0.75), proportion of loci
#      in population that must be confounded relative to the catalog locus.
#    --prune_haplo: prune out non-biological haplotypes unlikely to occur in
#      the population.
#    --max_haplo_cnt <limit>: only consider haplotypes for pruning if they
#      occur in fewer than max_haplo_cnt samples.
#  Model options:
#    --model_type <type>: either 'snp' (default), 'bounded', or 'fixed'
#    For the SNP or Bounded SNP model:
#      --alpha <num>: chi square significance level required to call a
#        heterozygote or homozygote, either 0.1 (default), 0.05, 0.01, or 0.001.
#    For the Bounded SNP model:
#      --bound_low <num>: lower bound for epsilon, the error rate, between 0
#        and 1.0 (default 0).
#      --bound_high <num>: upper bound for epsilon, the error rate, between 0
#        and 1.0 (default 1).
#  Logging Options:
#      --verbose: extended logging, including coordinates of all changed
#        nucleotides (forces single-threaded execution).

# Launch rxstacks
rxstacks $b $P $o $t --lnl_dist --model_type bounded

# Exploring the log likelihoods
cat 05-stacks_rx/batch_1.rxstacks_lnls.tsv | \
    grep -v "^#" | \
    awk '{bucket=(int($2)); lnls[bucket] += 1} END { for (bucket in lnls) print bucket, "\t", lnls[bucket]}' | \
    sort -n > 05-stacks_rx/rxstacks_log_likelihoods.tsv

# Create figure
R -q -e 'source("./00-scripts/utility_scripts/plot_log_likelihoods.r")' 2>&1 > /dev/null

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

