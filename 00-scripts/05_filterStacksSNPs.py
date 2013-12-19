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
        self.maxAlleleFreq = float(l[9])
        self.maf = 1 - self.maxAlleleFreq
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
    """Locus object containing 1 or more SNPs
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
        previous_locus = -99
        snps = []
        for line in f:
            if not line.strip().startswith("#"):
                snp = SNP(line)
                if snp.locus != previous_locus:
                    if len(snps) > 0:
                        yield Locus(snps)
                    snps = []
                    if snp.population in pops:
                        snps.append(snp)
                        previous_locus = snp.locus
                else:
                    if snp.population in pops:
                        snps.append(snp)

# Filter functions
def filter_empty_loci(loci):
    """Remove loci without snps
    """
    for locus in loci:
        if len(locus.snps) > 0:
            yield locus

def filter_number_individuals(loci, min_presence, num_ind_per_pop, criterion, use_percent, header):
    """Remove snps that do not have enough individuals for all populations
    """
    with open("06-filters/01_filtered_number_individuals.tsv", "w") as out_f:
        for h in header:
            out_f.write(h)
        for locus in loci:
            pos_to_remove = []
            for pos in locus.snps:
                number_pass = 0
                for pop in locus.snps[pos]:
                    max_num_ind = num_ind_per_pop[pop]
                    presence = locus.snps[pos][pop].presence
                    if use_percent:
                        if 100 * float(presence) / max_num_ind > min_presence:
                            number_pass += 1
                    else:
                        if presence >= min_presence:
                            number_pass += 1
                if number_pass < criterion:
                    pos_to_remove.append(pos)
                    for pop in locus.snps[pos]:
                        out_f.write(str(locus.snps[pos][pop]))
            for p in pos_to_remove:
                locus.snps.pop(p)
            if len(locus.snps) > 0:
                yield locus

def filter_maf(loci, maf_global, maf_population, criterion_global, header):
    """Remove SNPs that do not have high enough global or population-wise MAFs
    """
    with open("06-filters/02_filtered_maf", "w") as out_f:
        for h in header:
            out_f.write(h)
        for locus in loci:
            pos_to_remove = []
            for pos in locus.snps:
                number_pass_global = 0
                number_pass_population = 0
                for pop in locus.snps[pos]:
                    if locus.snps[pos][pop].maf >= maf_global:
                        number_pass_global += 1
                    if locus.snps[pos][pop].maf >= maf_population:
                        number_pass_population += 1
                if number_pass_global < criterion_global and number_pass_population == 0:
                    pos_to_remove.append(pos)
                    for pop in locus.snps[pos]:
                        out_f.write(str(locus.snps[pos][pop]))
            for p in pos_to_remove:
                locus.snps.pop(p)
            if len(locus.snps) > 0:
                yield locus

def filter_heterozygozity(loci, max_hetero, joker, header):
    """Remove all loci where one population in one locus has too many
    heterozygous individuals
    """
    with open("06-filters/03_filtered_heterozygozity", "w") as out_f:
        for h in header:
            out_f.write(h)
        for locus in loci:
            failed = 0
            for pos in locus.snps:
                for pop in locus.snps[pos]:
                    if locus.snps[pos][pop].obsHet > max_hetero:
                        out_f.write(str(locus.snps[pos][pop]))
                        failed += 1
            if failed > joker:
                continue
            else:
                yield locus

def filter_fis(loci, min_fis, max_fis, joker, header):
    """Remove all loci where the Fis value of one population in one locus is
    outside the (min_fis, max_fis) range
    """
    with open("06-filters/04_filtered_fis", "w") as out_f:
        for h in header:
            out_f.write(h)
        for locus in loci:
            failed = 0
            for pos in locus.snps:
                for pop in locus.snps[pos]:
                    if not min_fis <= locus.snps[pos][pop].fis <= max_fis:
                        out_f.write(str(locus.snps[pos][pop]))
                        failed += 1
            if failed > joker:
                continue
            else:
                yield locus

def filter_snp_number(loci, max_snps, header):
    """Remove all loci with too many snps
    """
    with open("06-filters/05_filtered_snp_numbers", "w") as out_f:
        for h in header:
            out_f.write(h)
        for locus in loci:
            if len(locus.snps) > max_snps:
                out_f.write(str(locus))
            else:
                yield locus

def write_to_file(loci, output_file, header):
    num_loci_kept = 0
    with open(output_file, "w") as out_f:
        for h in header:
            out_f.write(h)
        for locus in loci:
            num_loci_kept += 1
            out_f.write(str(locus))
    print num_loci_kept, "loci were retained"

def get_header(input_file):
    header = "" 
    with open(input_file) as f:
        for line in f:
            if line.strip().startswith("#"):
                header = line
            else:
                return header

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
            help = 'use the -l option to list populations that should be removed from the analysis')
    parser.add_argument('-k', '--keep_only_pops', action="store_true",
            help = 'use the -l option to list populations that should be kept from the analysis')
    parser.add_argument('-l', '--list_of_populations', type=str,
            help = 'list of populations to remove or keep with -r of -k (string, no spaces, separaged by comas, eg: 1,2,3,7)')
    parser.add_argument('-p', '--min_presence', type=int, default=0,
            help = 'minimum number of individuals to keep locus (int, default: 0)')
    parser.add_argument('-x', '--min_presence_joker_populations', type=int, default=0,
            help = 'number of populations where it is permitted that the -p threshold does not pass, (int, 0 or more, default: 0')
    parser.add_argument('--use_percent', action="store_true",
            help = 'whether to use percentage (float, 0 to 100, default 0) instead of minimal number of individuals')
    parser.add_argument('-H', '--max_hetero', type=float, default=0.5,
            help = 'maximum proportion of heterozygous individuals (float, 0 to 1, default: 0.5)')
    parser.add_argument('-y', '--max_hetero_joker', type=int, default=0,
            help = 'number of populations where it is permitted that the -H threshold does not pass, (int, 0 or more, default: 0')
    parser.add_argument('-a', '--maf_global', type=float, default=0.05,
            help = 'minimum minor allele frequency that must be respected in all populations to retain locus (float, 0 to 1, default: 0.05)')
    parser.add_argument('-A', '--maf_population', type=float, default=0.1,
            help = 'minimum minor allele frequency that must be found in at least one population to retain locus (float, 0 to 1, default: 0.1)')
    parser.add_argument('-f', '--min_fis', type=float, default=-1,
            help = 'minimum Fis value (float, -1 to 1, default: -1)')
    parser.add_argument('-F', '--max_fis', type=float, default=1,
            help = 'maximum Fis value (float, -1 to 1, default: 1)')
    parser.add_argument('-z', '--fis_joker', type=int, default=0,
            help = 'number of populations where it is permitted that the -f or -F thresholds do not pass, (int, 0 or more, default: 0')
    parser.add_argument('-s', '--max_snp_number', type=int, default=999,
            help = 'maximum number of SNPs per locus (int, 0 or more, default: 999)')
    args = parser.parse_args()

    # Assert proper values for parameters (do directly with argparse?)
    assert args.max_snp_number > 0, "max_snp_number must be a non-null integer"
    assert 0 <= args.min_presence <= 100, "min_presence must be an integer between 0 and 100" 
    assert args.min_presence_joker_populations >= 0, "min_presence_joker_populations must be an integer that is 0 or more"
    assert 0 <= args.max_hetero <= 1, "max_hetero must be a decimal between 0 and 1"
    assert 0 <= args.maf_global <= 1.0, "maf_global must be a decimal between 0 and 1.0"
    assert 0 <= args.maf_population <= 1.0, "maf_population must be a decimal between 0 and 1.0"
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
        unwanted = set(args.list_of_populations.split(","))
        print "Unwanted pops:", list(sorted(unwanted))
        pops = [p for p in pops if p not in unwanted]
    elif args.keep_only_pops:
        wanted = set(args.list_of_populations.split(","))
        print "Wanted pops:", list(sorted(wanted))
        pops = [p for p in pops if p in wanted]

    num_populations = len(pops)

    # Input Filtering Output
    header = get_header(args.input_file)
    loci = filter_empty_loci(sumstats_parser(args.input_file,
        pops))
    loci = filter_number_individuals(loci, args.min_presence,
            num_ind_per_pop,
            num_populations - args.min_presence_joker_populations,
            args.use_percent,
            header)
    loci = filter_maf(loci,
            args.maf_global,
            args.maf_population,
            num_populations,
            header)
    loci = filter_heterozygozity(loci,
            args.max_hetero,
            args.max_hetero_joker,
            header)
    loci = filter_fis(loci,
            args.min_fis,
            args.max_fis,
            args.fis_joker,
            header)
    loci = filter_snp_number(loci,
            args.max_snp_number,
            header)
    write_to_file(loci, args.output_file,
            header)

