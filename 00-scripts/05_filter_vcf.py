#!/usr/bin/env python
"""Filtering VCF files output by STACKS

WARNING! will not work properly with STACKS versions older to v1.30)
"""

# Import
from collections import defaultdict
import argparse
import sys

# Classes
class Sample(object):
    """Genotype informantion of one sample
    """
    def __init__(self, info):
        self.info = info.split(":")
        self.genotype = self.info[0]
        self.depth = int(self.info[1])
        self.ref, self.alt = self.info[2].split(",")

        try:
            self.ref = int(self.ref)
        except:
            self.ref = 0

        try:
            self.alt = int(self.alt)
        except:
            self.alt = 0

        self.likelihood = self.info[3].split(",")[1]

        try:
            self.likelihood = float(self.likelihood)
        except:
            self.likelihood = 0

    def __repr__(self):
        return "\t".join([self.genotype,
            str(self.depth),
            str(self.ref),
            str(self.alt),
            str(self.likelihood)])

class SNP(object):
    def __init__(self, line):
        """Parse SNP information from one SNP line in the VCF file
        """
        self.line = line.split("\t")
        self.chrom = self.line[0]
        self.pos = int(self.line[1])
        self.locus_id = int(self.line[2])
        self.ref = self.line[3]
        self.alt = self.line[4]
        self.qual = self.line[5]
        self.filter = self.line[6]
        self.info = self.line[7]
        self.number_samples, self.global_maf = self.info.split(";")
        self.number_samples = int(self.number_samples.split("=")[1])
        self.global_maf = float(self.global_maf.split("=")[1].split(",")[1])
        self.format = self.line[8]
        self.samples = [Sample(x) for x in self.line[9:]]
        self.flags = Flags()

    def __repr__(self):
        return "\t".join([
                          self.chrom,
                          str(self.pos),
                          str(self.locus_id),
                          self.ref,
                          self.alt,
                          self.qual,
                          self.filter,
                          self.info,
                          self.format,
                          "\t".join([str(x) for x in self.samples])
                         ])
class Flags(object):
    """Flags used for filtering
    """

    # Counters for the number of filtered SNPs by reason
    # These are used by the class method 'report_filters'
    min_allele_coverage_count = 0
    min_presence_count = 0
    maf_global_count = 0
    maf_population_count = 0
    heterozygosity_count = 0
    min_fis_count = 0
    max_fis_count = 0
    max_snp_number_count = 0
    total_snps_count = 0
    total_kept_count = 0
    total_filtered_count = 0

    def __init__(self):
        self.min_allele_coverage = True
        self.min_presence = True
        self.maf_global = True
        self.maf_population = True
        self.heterozygosity = True
        self.min_fis = True
        self.max_fis = True
        self.max_snp_number = True

    def pass_filters(self):
        """Return True if a SNP passes the filters, False otherwise
        """
        Flags.total_snps_count += 1
        passing = (
                   self.min_allele_coverage and
                   self.min_presence and
                   (self.maf_global or self.maf_population) and
                   self.heterozygosity and
                   self.min_fis and
                   self.max_fis and
                   self.max_snp_number
                  )

        if passing:
            Flags.total_kept_count += 1
        else:
            Flags.total_filtered_count += 1

        return passing

    @classmethod
    def report_filters(cls, locus_counter):
        """Print a report about the numbers of SNPs that were filtered by
        each of the filters
        """
        print "--- Filtering results ------------------------------------"
        print "  {}\tGenotypes removed because of {}".format(cls.min_allele_coverage_count, "min_allele_coverage")
        print "  {}\tSNPs failed: {}".format(cls.min_presence_count, "min_presence")
        print "  {}\tSNPs failed: {}".format(cls.maf_global_count, "maf_global")
        print "  {}\tSNPs failed: {}".format(cls.maf_population_count, "maf_population")
        print "  {}\tSNPs failed: {}".format(cls.heterozygosity_count, "heterozygosity")
        #print "  {}\tSNPs failed: {}".format(cls.min_fis_count, "min_fis")
        #print "  {}\tSNPs failed: {}".format(cls.max_fis_count, "max_fis")
        print "  {}\tSNPs failed: {}".format(cls.max_snp_number_count, "max_snp_number")
        print "----------------------------------------------------------"
        print "  {}\tSNPs ({} loci) in input file".format(cls.total_snps_count, locus_counter)
        try:
            print "  {}\tSNPs ({}%) filtered out".format(cls.total_filtered_count, str(100.0 * float(cls.total_filtered_count) / float(cls.total_snps_count))[0:5])
        except:
            print "  {}\tSNPs ({}%) filtered out".format(cls.total_filtered_count, 100)
        print "  {}\tSNPs retained".format(cls.total_kept_count)
        print "----------------------------------------------------------"

    def format_filters(self):
        """Format the flag information for printing in filter file
        """
        return "\t".join([
                          str(self.min_presence),
                          str(self.maf_global),
                          str(self.maf_population),
                          str(self.heterozygosity),
                          str(self.min_fis),
                          str(self.max_fis),
                          str(self.max_snp_number)
                         ]) + "\n"

# Functions
def get_pop_info(line):
    """Return dictionary of population names with a list of their sample
    positions in the lines of the VCF file.
    """

    l = [x.split("_")[0] for x in line.split("\t")[9:]]
    unique_pops = sorted(list(set(l)))
    print "(" + str(len(unique_pops)) + " populations) <<<"
    pop_dict = {}

    for p in unique_pops:
        pop_dict[p] = []

    for i, p in enumerate(l):
        pop_dict[p].append(i)

    return pop_dict

def locus_iterator(input_file):
    """Iterate over the loci, one at a time, and yield the next one as a list
    containing all the SNPs for that locus
    """
    with open(input_file) as in_f:
        locus = []
        first = True
        last_id = -99

        for line in in_f:
            if line.startswith("#"):
                pass
            else:
                snp = SNP(line.strip())
                current_id = snp.locus_id

                if current_id == last_id:
                    locus.append(snp)
                else:
                    if first:
                        first = False
                    else:
                        yield locus
                    locus = [snp]
                    last_id = current_id
        yield locus

def write_filters(locus, handle):
    """Output filters to file
    """
    for snp in locus:
        handle.write("\t".join([
                                str(snp.pos),
                                str(snp.locus_id),
                                snp.flags.format_filters()
                               ]))

def write_locus(locus, handle):
    """Output good SNPs to file
    """
    for snp in locus:
        if snp.flags.pass_filters():
            handle.write("\t".join(snp.line) + "\n")

## Filter functions
def test_min_allele_coverage(locus, pop_info, min_allele_coverage):
    """Test if each genotype is backed by enough reads
    """
    for snp in locus:
        for sample in snp.samples:
            # calculate allele coverage for both alleles
            # change genotypes that do not meet threshold to '0/0'
            if sample.genotype in ["0/1", "1/0"]:
                if sample.alt < min_allele_coverage or sample.ref < min_allele_coverage:
                    sample.genotype = "./."
                    Flags.min_allele_coverage_count += 1

# Min presence
def test_min_presence(locus, pop_info, min_presence, joker, percent):
    """Test if enough populations have enough samples genotyped
    """
    for snp in locus:
        # Get Sample objects for each populations for that SNP
        populations = {}

        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        # Test if SNP meets minimum criterion for all pops
        pops_passing = 0

        for pop in populations:
            passing = 0
            for sample in populations[pop]:
                if sample.genotype != "./.":
                    passing += 1

            # Verifying if enough samples passed
            num_samples = len(populations[pop])
            needed_samples = min_presence
            needed_percent = float(min_presence) / 100.0 * float(num_samples)

            # Using percentage
            if percent and (passing >= needed_percent):
                pops_passing += 1

            # Using number of individuals
            elif passing >= needed_samples:
                pops_passing += 1

        if pops_passing < len(populations) - joker:
            snp.flags.min_presence = False
            Flags.min_presence_count += 1

# Maf global
def test_maf_global(locus, maf_global):
    """Test if the global MAF value is sufficiently high
    """
    for snp in locus:
        if snp.global_maf < maf_global:
            snp.flags.maf_global = False
            Flags.maf_global_count += 1

# Maf population
def test_maf_population(locus, pop_info, maf_population):
    """Test if at least one population has a sufficiently high MAF value
    """
    for snp in locus:
        # Get Sample objects for each populations for that SNP
        populations = {}
        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        # Test if SNP meets minimum criterion for all pops
        pops_passing = 0
        for pop in populations:
            allele_count = defaultdict(int)

            for sample in populations[pop]:
                if sample.genotype != "./.":
                    a1, a2 = sample.genotype.split("/")
                    allele_count[a1] += 1
                    allele_count[a2] += 1
            minimum = min(allele_count.values())
            total = sum(allele_count.values())
            if len(allele_count.values()) > 1:
                maf = float(minimum) / total
            else:
                maf = 0

            if maf >= maf_population:
                pops_passing += 1

        if pops_passing < 1:
            snp.flags.maf_population = False
            Flags.maf_population_count += 1

# Heterozygosity
def test_heterozygosity(locus, pop_info, max_hetero, max_hetero_joker):
    """Test whether too many individuals are heterozygous in at least
    one (or move if the joker value is used) population.
    """
    for snp in locus:
        # Get Sample objects for each populations for that SNP
        populations = {}
        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        # Test if SNP meets minimum criterion for all pops
        pops_passing = 0
        for pop in populations:

            num_hetero = 0
            num_genotypes = 0
            for sample in populations[pop]:
                if sample.genotype != "./.":
                    num_genotypes += 1
                    if sample.genotype in ["0/1", "1/0"]:
                        num_hetero += 1

            proportion = float(num_hetero) / num_genotypes

            if proportion <= max_hetero:
                pops_passing += 1

        if pops_passing < (len(populations) - max_hetero_joker):
            snp.flags.heterozygosity = False
            Flags.heterozygosity_count += 1

# TODO Fis min
def test_fis_min(locus, pop_info, min_fis):
    """NOT IMPLEMENTED YET.

    Test that the minimum Fis value of all the populations is high enough
    """
    pass

# TODO Fis max
def test_fis_max(locus, pop_info, max_fis):
    """NOT IMPLEMENTED YET.

    Test that the maximum Fis value of all the populations is low enough
    """
    pass

# Max number of SNPs (before or after filtering for the others?)
def test_max_snp_number(locus, pop_info, max_snp_number):
    """Test if there are too many SNPs at one locus
    """
    for snp in locus:
        if len(locus) > max_snp_number:
            snp.flags.max_snp_number = False
            Flags.max_snp_number_count += 1

# Main
if __name__ == '__main__':

    # Parsing user input with 'argparse'
    parser = argparse.ArgumentParser(description=
            """Filter SNPs from STACKS VCF result file
            (WARNING! will not work properly with STACKS versions older to v1.30)
            """)

    parser.add_argument("-i", "--input_file", type=str, required=True,
            help = "input VCF file")
    parser.add_argument("-o", "--output_file", type=str, required=True,
            help = "output VCF file")
    parser.add_argument("-c", "--min_allele_coverage", type=int, default=0,
            help = "minimum allele depth to keep a genotype (otherwise modified to '0/0') (int, default: 0)")
    parser.add_argument("-p", "--min_presence", type=int, default=0,
            help = "minimum number of individuals per population to keep locus (int, default: 0)")
    parser.add_argument("-x", "--min_presence_joker_populations", type=int, default=0,
            help = "number of populations where it is permitted that the -p threshold does not pass, (int, 0 or more, default: 0")
    parser.add_argument("--use_percent", action="store_true",
            help = "whether to use percentage (float, 0 to 100, default 0) instead of minimal number of individuals")
    parser.add_argument("-a", "--maf_global", type=float, default=0.0,
            help = "minimum minor allele frequency that must be respected in all populations to retain locus (float, 0 to 1, default: 0)")
    parser.add_argument("-A", "--maf_population", type=float, default=0.0,
            help = "minimum minor allele frequency that must be found in at least one population to retain locus (float, 0 to 1, default: 0)")
    parser.add_argument("-H", "--max_hetero", type=float, default=1,
            help = "maximum proportion of heterozygous individuals (float, 0 to 1, default: 1)")
    parser.add_argument("-y", "--max_hetero_joker", type=int, default=0,
            help = "number of populations where it is permitted that the -H threshold does not pass, (int, 0 or more, default: 0")

    # Not implemented yet
    parser.add_argument("-f", "--min_fis", type=float, default=-1,
            help = "minimum Fis value (float, -1 to 1, default: -1)")
    parser.add_argument("-F", "--max_fis", type=float, default=1,
            help = "maximum Fis value (float, -1 to 1, default: 1)")
    parser.add_argument("-z", "--fis_joker", type=int, default=0,
            help = "number of populations where it is permitted that the -f or -F thresholds do not pass, (int, 0 or more, default: 0")

    parser.add_argument("-s", "--max_snp_number", type=int, default=999,
            help = "maximum number of SNPs per locus (int, 0 or more, default: 999)")
    args = parser.parse_args()

    # Assert proper values for parameters
    assert args.max_snp_number > 0, "max_snp_number must be a non-null integer"
    assert 0 <= args.min_allele_coverage <= 100, "min_allele_coverage must be an integer between 0 and 100"
    assert 0 <= args.min_presence <= 100, "min_presence must be an integer between 0 and 100" 
    assert args.min_presence_joker_populations >= 0, "min_presence_joker_populations must be an integer that is 0 or more"
    assert 0 <= args.max_hetero <= 1, "max_hetero must be a decimal between 0 and 1"
    assert 0 <= args.maf_global <= 1.0, "maf_global must be a decimal between 0 and 1.0"
    assert 0 <= args.maf_population <= 1.0, "maf_population must be a decimal between 0 and 1.0"

    # Not implemented yet
    assert -1 <= args.min_fis <= 1, "min_fis must be a decimal between -1 and 1"
    assert -1 <= args.max_fis <= 1, "max_fis must be a decimal between -1 and 1"

    # Get header from VCF
    header = []
    print ">>> Treating: " + args.input_file, 
    with open(args.input_file) as in_file:
        for line in in_file:
            l = line.strip()
            if l.startswith("##"):
                header.append(line)
            elif l.startswith("#CHROM"):
                header.append(line)
                pop_info = get_pop_info(l)
            else:
                break

    # Open output files handles
    with open(args.output_file, "w") as out_file:
        with open(args.output_file + "_filters.tsv", "w") as filters_file:

            # Printing filtered CSV header
            out_file.writelines(header)

            # Printing filters file header
            filters_file.write("\t".join([
                                          "Position",
                                          "ID",
                                          "MinPresence",
                                          "MafGlobal",
                                          "MafPopulation",
                                          "Hetero",
                                          "FisMin",
                                          "FisMax",
                                          "MaxSnpNumber"
                                         ]) + "\n")

            locus_counter = 1
            report_every = 10000
            max_loci = 9999999999
            #max_loci = 10

            # Iterate over the loci and filter the SNPs
            for locus in locus_iterator(args.input_file):

                if locus_counter >= max_loci:
                    break
                else:
                    locus_counter += 1

                if locus_counter % report_every == 0:
                    print "  Treating locus number: " + str(locus_counter)

                # Run filters
                # The filter functions automatically update SNP flags
                test_min_allele_coverage(
                                         locus,
                                         pop_info,
                                         args.min_allele_coverage
                                        )
                test_min_presence(
                                  locus,
                                  pop_info,
                                  args.min_presence,
                                  args.min_presence_joker_populations,
                                  args.use_percent
                                 )
                
                test_maf_global(locus, args.maf_global)

                test_maf_population(
                                    locus,
                                    pop_info,
                                    args.maf_population
                                   )

                test_heterozygosity(
                                    locus,
                                    pop_info,
                                    args.max_hetero,
                                    args.max_hetero_joker
                                   )

                test_fis_min(locus, pop_info, args.min_fis)

                test_fis_max(locus, pop_info, args.max_fis)

                test_max_snp_number(locus, pop_info, args.max_snp_number)

                # Write output
                write_filters(locus, filters_file)
                write_locus(locus, out_file)

    # Print filtering report
    Flags.report_filters(locus_counter)
