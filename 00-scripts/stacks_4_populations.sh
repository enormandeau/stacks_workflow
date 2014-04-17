#!/bin/bash
# Launch populations

# Options: Comment out options that you do not wish to use
b="-b 1"            # b: Batch ID to examine when exporting from the catalog
P="-P 05-stacks"    # P: path to the Stacks output files.
M="-M 01-info_files/population_map.txt"  #M: path to the population map, a tab
                    #   separated file describing which individuals belong in
                    #   which population
#s="-s file_for_sql"  #s: output a file to import results into an SQL
                      #   database
#B="-B blacklist_file.txt"   # B: specify a file containing Blacklisted markers
                             #   to be excluded from the export
#W="-W whitelist_file.txt"   # W: specify a file containing Whitelisted markers
                             #   to include in the export
#e="-e ENZYME"      # e: restriction enzyme, required if generating 'genomic'
                    #  output
t="-t 16"           # t: number of threads to run in parallel sections of code

# Data filtering
r="-r 0.2"          # r: minimum percentage of individuals in a population
                    #  required to process a locus for that population
#p="-p 1"           # p: minimum number of populations a locus must be present
                    #   in order to process a locus
m="-m 6"            # m: specify a minimum stack depth required for individuals
                    #   at a locus
#a="-a 0.05"        # a: specify a minimum minor allele frequency required
                    #   before calculating Fst at a locus (0 < a < 0.5)
#f="-f p_value"     # f: specify a correction to be applied to Fst values:
                    #   'p_value', 'bonferroni_win', or 'bonferroni_gen'
#p_value_cutoff="--p_value_cutoff 0.05"   # --p_value_cutoff [num]: required
                    #   p-value to keep an Fst measurement (0.05 by default)
                    #   Also used as base for Bonferroni correction

# Kernel-smoothing algorithm
#k="-k"             # k: enable kernel-smoothed Pi, Fis, and Fst calculations
#window_size="--window_size 150Kb"       # --window_size [num]: distance over
                    #   which to average values (sigma, default 150Kb)


# Bootstrap resampling
#bootstrap="--bootstrap"                 # --bootstrap [exact|approx]: enable
                    #   bootstrap resampling for population statistics
                    #   (reference genome required)
#bootstrap_reps="--bootstrap_reps 100"   #--bootstrap_reps [num]: number of
                    #   bootstrap resamplings to calculate (default 100)

# file output options
#genomic="--genomic"        # --genomic: output each nucleotide position
                            #   (fixed or polymorphic) in all population members
                            #   to a file
vcf="--vcf"                 # --vcf: output results in Variant Call Format (VCF)
#genepop="--genepop"        # --genepop: output results in GenePop format
#structure="--structure"    # --structure: output results in Structure format
#phase="--phase"            # --phase: output genotypes in PHASE/fastPHASE format.
#beagle="--beagle"          # --beagle: output genotypes in Beagle format.
#plink="--plink"            # --plink: output genotypes in PLINK format.

#phylip="--phylip"          # --phylip: output nucleotides that are fixed-within
                            #   and variant among populations in Phylip format
                            #   for phylogenetic tree construction
#phylip_var="--phylip_var"  # --phylip_var: include variable sites in the phylip
                            #   output
#write_single_snp="--write_single_snp"  #--write_single_snp: write only the
                            #   first SNP per locus in Genepop and Structure
                            #   outputs
#log_fst_comp="--log_fst_comp"  # log components of Fst calculations to a file.

# Launch populations
populations $b $P $M $r $m $g $V $B $W $s $e $t $r $p $m $a $f $p_value_cutoff \
    $k $window_size $bootstrap $bootstrap_reps $genomic $vcf $genepop \
    $structure $phase $beagle $plink $phylip $phylip_var \
    $write_single_snp $log_fst_comp

#TODO experiment with the m option
# m: 4 6 8 10

