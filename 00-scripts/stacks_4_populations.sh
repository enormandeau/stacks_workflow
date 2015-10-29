#!/bin/bash
# Launch populations
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98-log_files"
INFO_FILES_FOLDER="01-info_files"
POP_MAP="population_map.txt"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"
cp $INFO_FILES_FOLDER/$POP_MAP $LOG_FOLDER/"$TIMESTAMP"_"$POP_MAP"

# IMPORTANT: make sure your read about the available options for 'populations'
# in the STACKS papers

# Comment out the options required to run your analysis

b="-b 1"               # Batch ID to examine when exporting from the catalog
P="-P 05-stacks"       # Path to the Stacks output files.
M="-M 01-info_files/population_map.txt"  # Path to the population map, a
                                         # tab separated file describing
                                         # which individuals belong in
                                         # which population
#s="-s file_for_sql"  # Output a file to import results into an SQL database
#B="-B blacklist_file.txt"   # Specify a file containing Blacklisted markers
                             # to be excluded from the export
#W="-W whitelist_file.txt"   # Specify a file containing Whitelisted markers
                             # to include in the export
#e="-e ENZYME"   # Restriction enzyme, required if generating 'genomic' output
t="-t 16"         # Number of threads to run in parallel sections of code
#v="-v"          # Print program version.
#h="-h"          # Display this help message.

# Merging and Phasing:
#merge_sites="--merge_sites"   # merge loci that were produced from the same
                               # restriction enzyme cutsite
                               # (requires reference-aligned data).
#merge_prune_lim="--merge_prune_lim"   # when merging adjacent loci,
                                      # if at least X% samples posses both loci
                                      # prune the remaining samples out of the analysis.
# Data filtering
r="-r 0.5"          # Minimum percentage of individuals in a population
                    # required to process a locus for that population
p="-p 4"            # Minimum number of populations a locus must be present
                    # in order to process a locus
m="-m 4"            # Specify a minimum stack depth required for individuals
                    # at a locus
f="-f p_value"      # Specify a correction to be applied to Fst values:
min_maf="-a 0.0"         # Specify a minimum minor allele frequency required
                    # before calculating Fst at a locus (0 < a < 0.5)
                    # 'p_value', 'bonferroni_win', or 'bonferroni_gen'
p_value_cutoff="--p_value_cutoff 0.10"  # Required p-value to keep an Fst
                                        # measurement (0.05 by default). Also
                                        # used as base for Bonferroni
                                        # correction
lnl_lim="--lnl_lim -10"   # Filter loci with log likelihood values below this
                         # threshold.
#write_single_snp="--write_single_snp"  # write only the first SNP per locus in
                                        # Genepop and Structure outputs
#write_random_snp="--write_random_snp"  # Restrict data analysis to one random
                                        # SNP per locus.

#Fstats:
#fstats="--fstats"     # Enable SNP and haplotype-based F statistics

# Kernel-smoothing algorithm
#k="-k"                # enable kernel-smoothed Pi, Fis, Fst, Fst', and Phi_st calculations
#window_size="--window_size 150Kb"     # distance over which to average values
                                       #(sigma, default 150Kb)
# Bootstrap resampling
#bootstrap="--bootstrap"               # turn on broostrap resampling
                                       # for all smoothed statistics.
#bootstrap_pifis="--bootstrap_pifis"   # turn on boostrap resampling
                                       # for smoothed SNP-based Pi
                                       # and Fis calculations.
#bootstrap_fst="--bootstrap_fst"       # turn on boostrap resampling
                                       # for smoothed Fst calculations
                                       # based on pairwise population
                                       # comparison of SNPs.
#bootstrap_div="--bootstrap_div"       # turn on boostrap resampling
                                       # for smoothed haplotype diveristy
                                       # and gene diversity calculations
                                       # based on haplotypes.
#bootstrap_phist="--bootstrap_phist"   # turn on boostrap resampling
                                       # for smoothed Phi_st calculations
                                       # based on haplotypes.
#bootstrap_reps="--bootstrap_reps 100" # number of bootstrap resamplings
                                       # to calculate (default 100).
#bootstrap_wl="--bootstrap_wl [path]"  # only bootstrap loci contained
                                       # in this whitelist.

# file output options
#genomic="--genomic"               # Output each nucleotide position
                                   # (fixed or polymorphic) in all population
                                   # members to a file.
#fasta="--fasta"                   # Output full sequence for each allele,
                                   # from each sample locus in FASTA format.
vcf="--vcf"                       # Output results in Variant Call Format (VCF).
vcf_haplotypes="--vcf_haplotypes" # Output haplotypes in Variant Call Format (VCF).
#genepop="--genepop"               # Output results in GenePop format
#structure="--structure"           # Output results in Structure format
#phase="--phase"                   # Output genotypes in PHASE/fastPHASE format.
#fastphase="--fastphase"           # Output genotypes in fastPHASE format
#beagle="--beagle"                 # Output genotypes in Beagle format.
#beagle_phased="--beagle_phased"   # Output haplotypes in Beagle format.
#plink="--plink"                   # Output genotypes in PLINK format.
#phylip="--phylip"                 # Output nucleotides that are fixed-within
                                   # and variant among populations in Phylip
                                   # format for phylogenetic tree construction
#phylip_var="--phylip_var"         # Include variable sites in the phylip output
#hzar="--hzar"                     # Output genotypes in Hybrid Zone Analysis
                                   # using R (HZAR) format.

# Debugging
#verbose="--verbose"            # turn on additional logging.
#log_fst_comp="--log_fst_comp"  # log components of Fst calculations to a file.

# Launch populations

populations $b $P $M $r $m $g $V $B $W $s $e $t $v $h \
    $merge_sites $merge_prune_lim \
    $r $p $m $f $min_maf $p_value_cutoff $lnl_lim $write_single_snp \
    $write_random_snp $fstats $k $window_size $bootstrap $bootstrap_pifis \
    $bootstrap_fst $bootstrap_div $bootstrap_phist $bootstrap_reps \
    $bootstrap_wl $genomic $fasta $vcf $vcf_haplotypes $genepop $structure
    $phase $fastphase $beagle $beagle_phased $plink $phylip $phylip_var $hzar \
    $verbose $log_fst_comp 2>&1 | tee 98-log_files/"$TIMESTAMP"_stacks_4_populations.log

