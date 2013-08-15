#!/usr/bin/python
"""Filter sumstats.tsv files from STACKS

To get detailed list of options, run:
    python 05_filterStacksSNPs.py -h
"""

# Importing modules
from collections import defaultdict
import argparse
import sys

# Defining classes
class SNP(object):
    """Store and manipulate SNP information
    """
    def __init__(self, line):
        l = line.strip().split("\t")
        self.line = l
        self.locus = int(l[1])
        self.position = int(l[4])
        self.population = l[5]
        self.maxFreq = float(l[9])
        self.minFreq = 1 - self.maxFreq
        self.obsHet = float(l[10])
        self.expHet = float(l[12]) 
        self.fis = float(l[17])
        try:
            self.presence = int(l[8])
        except:
            self.presence = 0
        try:
            self.hetero = float(l[10])
        except:
            self.hetero = 0.0
    def __str__(self):
        return "\t".join(self.line) + "\n"

class Locus(object):
    """Locus object containing 1 or more SNPs and with methods to filter them

    self.snps[POSITION][POPULATION] = SNP instance
    """
    def __init__(self, snps):
        self.snps = defaultdict(dict)
        for s in snps:
            self.snps[s.position][s.population] = s
    def __repr__(self):
        to_return = ""
        for pos in self.snps:
            for pop in self.snps[pos]:
                to_return += str(self.snps[pos][pop])
        return to_return

# Defining functions
def sumstats_parser(sumstats_file, pops):
    """Parse sumstats file and yield one locus at a time
    """
    with open(sumstats_file) as f:
        started = False
        previous_locus = -99
        snps = []
        for line in f:
            if not line.strip().startswith("#"):
                snp = SNP(line)
                if started and snp.locus != previous_locus:
                    if len(snps) > 0:
                        yield Locus(snps)
                    snps = []
                    if snp.population in pops:
                        snps.append(snp)
                        previous_locus = snp.locus
                else:
                    started = True
                    if snp.population in pops:
                        snps.append(snp)

# Filter functions
def filter_empty_loci(loci):
    """Remove loci without snps
    """
    for locus in loci:
        if not len(locus.snps) == 0:
            yield locus

def filter_number_individuals(loci, presence, criterion):
    """WARNING: filter individual SNPs, not the locus
    """
    print "Presence wanted:", presence, "Num pops:", criterion
    for locus in loci:
        passing = False
        number_pass = 0
        print "=== New locus ===="
        pos = locus.snps.keys()[0]
        for pop in locus.snps[pos]:
            print "  Population:", pop, "Presence:", locus.snps[pos][pop].presence, "/", presence
            if locus.snps[pos][pop].presence >= presence:
                number_pass += 1
                print "  Locus:", locus.snps[pos][pop].locus, "  Pass:", number_pass
        if number_pass >= criterion:
            passing = True
        if passing:
            print r"  \\\ Passed ///"
            yield locus

def filter_global_maf(loci, maf, criterion):
    """WARNING: filter individual SNPs, not the locus
    """
    for locus in loci:
        yield locus

def filter_locus_maf(loci, maf, criterion):
    """WARNING: filter individual SNPs, not the locus
    """
    for locus in loci:
        yield locus

def filter_heterozygozity(loci, max_hetero, criterion):
    for locus in loci:
        yield locus

def filter_min_fis(loci, min_fis, criterion):
    for locus in loci:
        yield locus

def filter_max_fis(loci, max_fis, criterion):
    for locus in loci:
        yield locus

def filter_snp_number(loci, max_snps):
    for locus in loci:
        yield locus

def write_to_file(loci, output_file):
    with open(output_file, "w") as out_f:
        for locus in loci:
            out_f.write(str(locus))

# Conditional functions
def filter_true(lst): # filter true results
    return [l for l in lst if l]

def at_least_n_true(lst, n):
    return len(filter_true(lst)) >= n

# Main function
if __name__ == "__main__":
    # Parsing user input
    parser = argparse.ArgumentParser(description=
            """Filter SNPs from STACKS sumstats result file""")
    parser.add_argument('-i', '--input_file', type=str, required=True,
            help = 'input sumstats file')
    parser.add_argument('-o', '--output_file', type=str, required=True,
            help = 'output file')
    parser.add_argument('-P', '--population_map_file', type=str, required=True,
            help = 'population_map.txt file from "01-info_files"')
    parser.add_argument('-r', '--remove_pops', action="store_true",
            help = 'population_map.txt file from "01-info_files"')
    parser.add_argument('-k', '--keep_only_pops', action="store_true",
            help = 'population_map.txt file from "01-info_files"')
    parser.add_argument('-l', '--list_of_populations', type=str,
            help = 'list of populations to remove or keep with -r of -k (string, no spaces, separaged by comas, eg: 1,2,3,7)')
    parser.add_argument('-n', '--max_allele_number', type=int, default=2,
            help = 'maximum number of alleles per SNP, (int, 0 to 4, default: 2)')
    parser.add_argument('-p', '--min_presence', type=int, default=0,
            help = 'minimum number of individuals to keep locus (int, 0 to 100, default: 0)')
    parser.add_argument('-H', '--max_hetero', type=float, default=0.5,
            help = 'maximum proportion of heterozygous individuals (float, 0 to 1, default: 0.5)')
    parser.add_argument('-a', '--maf_global', type=float, default=0.05,
            help = 'minimum minor allele frequency, global (float, 0 to 1, default: 0.05)')
    parser.add_argument('-A', '--maf_locus', type=float, default=0.05,
            help = 'minimum minor allele frequency, locus (float, 0 to 1, default: 0.05)')
    parser.add_argument('-f', '--min_fis', type=float, default=-1,
            help = 'minimum Fis value (float, -1 to 1, default: -1)')
    parser.add_argument('-F', '--max_fis', type=float, default=1,
            help = 'maximum Fis value (float, -1 to 1, default: 1)')
    parser.add_argument('-s', '--max_snp_number', type=int, default=999,
            help = 'maximum number of SNPs per locus (int, 0 or more, default: 999)')
    args = parser.parse_args()
    
    # Assert proper values for parameters (do directly with argparse?)
    assert args.max_snp_number > 0, "max_snp_number must be a non-null integer"
    assert args.max_allele_number in [2, 3, 4], "max_allele_number must be 2, 3 or 4"
    assert 0 <= args.min_presence <= 100, "min_presence must be an integer between 0 and 100" 
    assert 0 <= args.max_hetero <= 1, "max_hetero must be a decimal between 0 and 1"
    assert 0 <= args.maf_global <= 0.5, "maf_global must be a decimal between 0 and 0.5"
    assert 0 <= args.maf_locus <= 0.5, "maf_locus must be a decimal between 0 and 0.5"
    assert -1 <= args.min_fis <= 1, "min_fis must be a decimal between -1 and 1"
    assert -1 <= args.max_fis <= 1, "max_fis must be a decimal between -1 and 1"
    assert not (args.remove_pops and args.keep_only_pops), "you must choose either remove_pops or keep_only_pops"

    # Create dictionary of number of individuals per population
    num_ind_per_pop = defaultdict(int)
    with open(args.population_map_file) as f:
        for line in f:
            if line.strip() != "":
                l = line.strip().split()
                pop = l[1]
                num_ind_per_pop[pop] += 1

    # Define which populations should be kept
    pops = num_ind_per_pop.keys()
    if args.remove_pops:
        unwanted = set(args.list_of_populations)
        pops = [p for p in pops if p not in unwanted]
    elif args.keep_only_pops:
        wanted = set(args.list_of_populations)
        pops = [p for p in pops if p in wanted]

    # TODO temporary
    criterion = 2 # len(pops)

    # Parse sumstats file and remove empty loci
    loci = filter_empty_loci(sumstats_parser(args.input_file, pops))

    # Remove loci with too few individuals
    loci = filter_number_individuals(loci, args.min_presence, criterion)

    ## Remove SNPs with too low global MAF
    #loci = filter_global_maf(loci, args.maf_global, criterion)

    ## Remove SNPs with too low locus MAF
    #loci = filter_locus_maf(loci, args.maf_locus, criterion)

    ## Remove loci with too many heterozygous individuals
    #loci = filter_heterozygozity(loci, args.min_presence, criterion)

    ## Remove loci with too low Fis values
    #loci = filter_min_fis(loci, args.min_fis, criterion)

    ## Remove loci with too high Fis values
    #loci = filter_max_fis(loci, args.max_fis, criterion)

    ## Remove loci with too many SNPs
    #loci = filter_snp_number(loci, args.max_snp_number)

    # Write remaining loci to file
    write_to_file(loci, args.output_file)

