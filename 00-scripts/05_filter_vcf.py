#!/usr/bin/env python3
"""Filtering VCF files output by STACKS

WARNING! will not work properly with STACKS versions older to v1.30 since the
minor allele counts in the VCF file created by STACKS were not reported
properly at that time.
"""

# Modules
import gzip
from collections import defaultdict
import subprocess
import argparse
import pprint
import numpy
import math
import sys
import os

# Classes
class Sample(object):
    """Genotype information of one sample
    """

    def __init__(self, info):
        self.info_full = info
        self.info = info.split(":")
        self.genotype = self.info[0]

        if self.genotype == "./.":
            self.depth = 0
            self.ref = 0
            self.alt = 0

        else:
            self.depth = int(self.info[1])

            if self.info[2] == ".":
                self.ref = self.depth
                self.alt = 0
            else:
                self.ref, self.alt = self.info[2].split(",")

        try:
            self.ref = int(self.ref)
        except:
            self.ref = 0

        try:
            self.alt = int(self.alt)
        except:
            self.alt = 0

        try:
            maximum = max([self.ref, self.alt])
            minimum = min([self.ref, self.alt])
            self.allelic_imbalance = float(maximum) / float(minimum)
        except:
            self.allelic_imbalance = 0.0

        # STACKS 2
        if len(self.info) == 5:
            self.genotype_likelihood = ":".join(self.info[3:5])

        # STACKS 1
        else:
            try:
                self.genotype_likelihood = self.info[3].split(",")[1]
                self.genotype_likelihood = float(self.genotype_likelihood)
            except:
                self.genotype_likelihood = 0.0

    def __repr__(self):
        return self.info_full
        #return ":".join([self.genotype,
        #    str(self.depth),
        #    str(self.ref) + "," + str(self.alt),
        #    str(self.genotype_likelihood)]
        #    )


class SNP(object):
    def __init__(self, line, locus_separator):
        """Parse SNP information from one SNP line in the VCF file
        """
        
        self.line = line.split("\t")
        self.chrom = self.line[0]
        self.pos = int(self.line[1])
        self.locus_full = self.line[2]
        self.locus_id = self.locus_full.split(locus_separator)[0]

        if self.locus_id == ".":
            self.locus_id = self.line[0]

        self.ref = self.line[3]
        self.alt = self.line[4]
        self.qual = self.line[5]
        self.filter = self.line[6]
        self.info = self.line[7]
        self.number_samples = self.info.split(";")[0]
        self.global_maf = self.info.split(";")[1]
        self.number_samples = int(self.number_samples.split("=")[1])

        try:    # Version 1.41-
            self.global_maf = float(self.global_maf.split("=")[1].split(",")[1])
        except: # Version 1.42+
            self.global_maf = float(self.global_maf.split("=")[1])

        self.format = self.line[8]
        self.samples = [Sample(x) for x in self.line[9:]]
        self.flags = Flags()

    def __repr__(self):
        return "\t".join([
                          self.chrom,
                          str(self.pos),
                          self.locus_full,
                          self.ref,
                          self.alt,
                          self.qual,
                          self.filter,
                          self.info,
                          self.format,
                          "\t".join([str(x) for x in self.samples])
                         ])

class Locus(object):
    """A locus is a list of SNP objects and methods to manipulate them
    """

    def __init__(self, snp_list):
        self.snps = snp_list

    def __repr__(self):
        return "There are {} snps in this locus".format(len(self.snps))

class Flags(object):
    """Flags used for filtering each SNP
    """

    # Counters for the number of filtered SNPs by reason
    # These are used by the class method 'report_filters'
    min_allele_coverage_count = 0
    min_depth_count = 0
    max_allele_coverage_count = 0
    min_presence_count = 0
    min_maf_global_count = 0
    min_maf_population_count = 0
    min_mas_count = 0
    max_heterozygosity_count = 0
    min_fis_count = 0
    max_fis_count = 0
    max_allelic_imbalance_count = 0
    max_snp_number_count = 0
    total_snps_count = 0
    total_good_snps_count = 0
    total_filtered_count = 0

    def __init__(self):
        """Initialize flag values for new SNP flags
        """

        self.min_presence = True
        self.maf_global = True
        self.maf_population = True
        self.mas = True
        self.heterozygosity = True
        self.min_fis = True
        self.max_fis = True
        self.max_snp_number = True

    def pass_filters(self, add_count=False):
        """Return True if a SNP passes the filters, False otherwise
        """

        if add_count:
            Flags.total_snps_count += 1

        passing = (
                   self.min_presence and
                   (self.maf_global or self.maf_population) and
                   self.mas and
                   self.heterozygosity and
                   self.min_fis and
                   self.max_fis and
                   self.max_snp_number
                  )

        if passing and add_count:
            Flags.total_good_snps_count += 1
        elif add_count:
            Flags.total_filtered_count += 1

        return passing

    @classmethod
    def report_filters(cls, locus_counter, total_good_loci, report_file):
        """Print report about the numbers of SNPs that didn't pass each of the filters
        """

        report = []
        report.append("===================================================")
        report.append("  {} Genotypes removed  ({}) [{}]".format(pad(cls.min_allele_coverage_count),
                "min_allele_coverage", args.min_allele_coverage))
        report.append("  {} Genotypes removed  ({}) [{}]".format(pad(cls.min_depth_count),
                "min_depth", args.min_depth))
        report.append("  {} Genotypes removed  ({}) [{}]".format(pad(cls.max_allelic_imbalance_count),
                "max_allelic_imbalance", args.max_allelic_imbalance))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.max_allele_coverage_count),
                "max_allele_coverage", args.max_allele_coverage))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.min_presence_count),
                "min_presence", args.min_presence))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.min_maf_global_count),
                "maf_global", args.maf_global))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.min_maf_population_count),
                "maf_population", args.maf_population))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.min_mas_count),
                "mas", args.mas))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.max_heterozygosity_count),
                "heterozygosity", args.max_hetero))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.min_fis_count),
                "min_fis", args.min_fis))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.max_fis_count),
                "max_fis", args.max_fis))
        report.append("  {} SNPs failed        ({}) [{}]".format(pad(cls.max_snp_number_count),
                "max_snp_number", args.max_snp_number))
        report.append("---------------------------------------------------")
        report.append("  {} SNPs ({} loci) in input file".format(pad(cls.total_snps_count), locus_counter))

        try:
            report.append("  {} SNPs ({}%) filtered out".format(pad(cls.total_filtered_count),
                    str(round(100.0 * float(cls.total_filtered_count) / float(cls.total_snps_count), 1))))
        except:
            report.append("  {} SNPs ({}%) filtered out".format(pad(cls.total_filtered_count), 100))

        report.append("  {} SNPs retained".format(pad(cls.total_good_snps_count), total_good_loci))
        report.append("===================================================\n")

        report = "\n".join(report)
        # Print report
        print(report)

        # Write report to file
        report_file.write(report)

    def format_filters(self):
        """Format the flag information for printing in filter file
        """

        return "\t".join([
                          str(self.min_presence),
                          str(self.maf_global),
                          str(self.maf_population),
                          str(self.mas),
                          str(self.heterozygosity),
                          str(self.min_fis),
                          str(self.max_fis),
                          str(self.max_snp_number)
                         ]) + "\n"

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def get_individual_info(line):
    """Return dictionary of positions in the lines of the VCF file and their
    corresponding samples ([population, sample]).
    """

    l = [x.split("_")[0:2] for x in line.split("\t")[9:]]
    individual_dict = {}

    for i in range(len(l)):
        population, sample = l[i]
        individual_dict[i] = [population, sample, 0]

    individual_dict["total"] = 0

    return individual_dict

def get_population_info(line):
    """Return dictionary of population names with a list of their sample
    positions in the lines of the VCF file.
    """

    l = [x.split("_")[0] for x in line.split("\t")[9:]]
    unique_pops = sorted(list(set(l)))
    print("  " + str(len(unique_pops)) + " populations")
    pop_dict = {}

    for p in unique_pops:
        pop_dict[p] = []

    for i, p in enumerate(l):
        pop_dict[p].append(i)

    return pop_dict

def locus_iterator(input_file, locus_separator):
    """Iterate over the loci, one at a time, and yield the next one as a list
    containing all the SNPs for that locus
    """

    with myopen(input_file) as in_f:
        snps = []
        first = True
        last_id = -99

        for line in in_f:
            if line.startswith("#"):
                pass
            else:
                snp = SNP(line.strip(), locus_separator)
                current_id = snp.locus_id

                if current_id == last_id:
                    snps.append(snp)
                else:
                    if first:
                        first = False
                    else:
                        yield Locus(snps)
                    snps = [snp]
                    last_id = current_id

        yield Locus(snps)

## Writing output to files
def write_filters(locus, handle):
    """Output filters to file
    """

    for snp in locus.snps:
        handle.write("\t".join([str(snp.pos), str(snp.locus_id),
                                snp.flags.format_filters()]))

def write_locus(locus, handle):
    """Output good SNPs to file
    """

    for snp in locus.snps:
        if snp.flags.pass_filters(add_count=True):
            handle.write(str(snp) + "\n")

def write_whitelist(locus, handle):
    """Output ID of loci with good SNPs
    """

    global total_good_loci
    passed = False

    for snp in locus.snps:
        passed = False
        if snp.flags.pass_filters():
            passed = True

    if passed:
        handle.write(str(snp.locus_id) + "\n")
        total_good_loci += 1

## General purpose
def median(lst):
    lst = sorted(lst)

    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[int((len(lst)+1)/2)-1]
    else:
            return float(sum(lst[int(len(lst)/2)-1:int(len(lst)/2)+1]))/2.0

def pad(text, char=" ", min_length=6):
    """Pad text with char to min_length
    """

    # Will crash if not possible
    text = str(text)
    char = str(char)

    if len(text) >= min_length:
        return text
    else:
        missing_length = min_length - len(text)
        return missing_length * char + text

## Filter and data gathering functions
# Min allelie coverage
def test_min_allele_coverage(locus, pop_info, min_allele_coverage):
    """Test if each genotype is backed by enough reads
    """

    for snp in locus.snps:
        for sample in snp.samples:
            # calculate allele coverage for both alleles
            # change genotypes that do not meet threshold to './.'
            if sample.genotype in ["0/1", "1/0"]:
                if sample.alt < min_allele_coverage or sample.ref < min_allele_coverage:
                    sample.genotype = "./."
                    Flags.min_allele_coverage_count += 1

            elif sample.genotype in ["0/0", "1/1"]:
                if sample.depth < min_allele_coverage:
                    sample.genotype = "./."
                    Flags.min_allele_coverage_count += 1


# Min depth
def test_min_depth(locus, pop_info, min_depth):
    for snp in locus.snps:
        for sample in snp.samples:
            # calculate global genotype coverage for both alleles
            # change genotypes that do not meet threshold to './.'
            if sample.genotype in ["0/1", "1/0"]:
                if sample.depth < min_depth:
                    sample.genotype = "./."
                    sample.info_full = "./." + sample.info_full[3:]
                    Flags.min_depth_count += 1

            elif sample.genotype in ["0/0", "1/1"]:
                if sample.depth < min_depth:
                    sample.genotype = "./."
                    sample.info_full = "./." + sample.info_full[3:]
                    Flags.min_depth_count += 1

def get_depth_data(graph_dict, locus, pop_info):
    """Get min, median and max depth by pop and globally
    """

    for snp in locus.snps:
        for pop in pop_info:
            # By pop
            samples = [snp.samples[i] for i in pop_info[pop]]
            minimum = min([x.depth for x in samples])
            med = int(median([x.depth for x in samples]))
            maximum = max([x.depth for x in samples])
            graph_dict[pop]["medDepth"][med] += 1
            graph_dict[pop]["maxDepth"][maximum] += 1

            # Globally
            graph_dict["global"]["medDepth"][med] += 1
            graph_dict["global"]["maxDepth"][maximum] += 1

        break

# Max allelic_imbalance
def test_max_allelic_imbalance(locus, pop_info, max_allelic_imbalance):
    """Test if each heterozygote's genotype has a low enough allelic allelic_imbalance
    """

    for snp in locus.snps:
        for sample in snp.samples:
            if sample.genotype in ["0/1", "1/0"]:
                if sample.allelic_imbalance > max_allelic_imbalance:
                    sample.genotype = "./."
                    Flags.max_allelic_imbalance_count += 1

def get_allelic_imbalance_data(graph_dict, locus, pop_info):
    """Get allelic imbalance data by pop and globally
    """

    for snp in locus.snps:
        # By pop
        for pop in pop_info:
            samples = [snp.samples[i] for i in pop_info[pop]]
            imbalance = [x.allelic_imbalance for x in samples if x.allelic_imbalance != 0.0]
            imbalance = [round(math.log(x, 2), 3) for x in imbalance]
            for i in imbalance:
                graph_dict[pop]["allImbalance"][i] += 1
                graph_dict["global"]["allImbalance"][i] += 1

# Min presence
def test_min_presence(locus, pop_info, min_presence, joker, percent):
    """Test if enough populations have enough samples genotyped
    """

    for snp in locus.snps:
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
            if percent:
                if passing >= needed_percent:
                    pops_passing += 1

            # Using number of individuals
            else:
                if passing >= needed_samples:
                    pops_passing += 1

        if pops_passing < len(populations) - joker:
            snp.flags.min_presence = False
            Flags.min_presence_count += 1

def get_presence_data(graph_dict, locus, pop_info):
    """Get presence data by pop and globally
    """

    for snp in locus.snps:
        # Get Sample objects for each populations for that SNP
        populations = {}

        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        for pop in populations:
            num_samples = 0
            passing = 0
            for sample in populations[pop]:
                num_samples += 1
                if sample.genotype != "./.":
                    passing += 1

            proportion_passing = round(float(passing) / num_samples, 4)
            graph_dict[pop]["presence"][proportion_passing] += 1
            graph_dict["global"]["presence"][proportion_passing] += 1

# Maf global
def test_maf_global(locus, maf_global):
    """Test if the global MAF value is sufficiently high
    """

    for snp in locus.snps:
        if snp.global_maf < maf_global:
            snp.flags.maf_global = False
            Flags.min_maf_global_count += 1

def get_maf_global_data(graph_dict, locus, pop_info):
    """Get maf data globally
    """

    for snp in locus.snps:
        maf = round(snp.global_maf, 3)
        graph_dict["global"]["mafGlobal"][maf] += 1

# Mas (minor allele samples)
def test_mas(locus, mas):
    """Test if the global MAF value is sufficiently high
    """

    for snp in locus.snps:
        num_samples_with_minor_allele = len([s for s in snp.samples if "1" in s.genotype])

        if num_samples_with_minor_allele < mas:
            snp.flags.mas = False
            Flags.min_mas_count += 1

# Maf population
def test_maf_population(locus, pop_info, maf_population):
    """Test if at least one population has a sufficiently high MAF value
    """

    for snp in locus.snps:
        # Get Sample objects for each populations for that SNP
        populations = {}
        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        # Test if SNP meets minimum criterion for at least one population
        pops_passing = 0
        for pop in populations:
            allele_count = defaultdict(int)

            for sample in populations[pop]:
                if sample.genotype != "./.":
                    a1, a2 = sample.genotype.split("/")
                    allele_count[a1] += 1
                    allele_count[a2] += 1

            if len(allele_count.values()) > 1:
                minimum = min(allele_count.values())
                total = sum(allele_count.values())
                maf = float(minimum) / total
            else:
                maf = 0

            if maf >= maf_population:
                pops_passing += 1

        if pops_passing < 1:
            snp.flags.maf_population = False
            Flags.min_maf_population_count += 1

def get_maf_population_data(graph_dict, locus, pop_info):
    """Get maf data by pop
    """

    for snp in locus.snps:
        # Get Sample objects for each populations for that SNP
        populations = {}

        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        for pop in populations:
            allele_count = defaultdict(int)

            for sample in populations[pop]:
                if sample.genotype != "./.":
                    a1, a2 = sample.genotype.split("/")
                    allele_count[a1] += 1
                    allele_count[a2] += 1

            if len(allele_count.values()) > 1:
                minimum = min(allele_count.values())
                total = sum(allele_count.values())
                maf = round(float(minimum) / total, 3)
            else:
                maf = 0.0

            graph_dict[pop]["mafPopulation"][maf] += 1
            graph_dict["global"]["mafPopulation"][maf] += 1

# Heterozygosity
def test_heterozygosity(locus, pop_info, max_hetero, max_hetero_joker):
    """Test whether too many individuals are heterozygous in at least
    one (or move if the joker value is used) population.
    """

    for snp in locus.snps:
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

            if num_genotypes > 0:
                proportion = float(num_hetero) / num_genotypes
            else:
                proportion = 0.0

            if proportion <= max_hetero:
                pops_passing += 1

        if pops_passing < (len(populations) - max_hetero_joker):
            snp.flags.heterozygosity = False
            Flags.max_heterozygosity_count += 1

def get_heterozygosity_data(graph_dict, locus, pop_info):
    """Get heterozygosity data by pop and globally
    """

    for snp in locus.snps:
        # Get Sample objects for each populations for that SNP
        populations = {}
        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        # Test if SNP meets minimum criterion for all pops
        for pop in populations:
            num_hetero = 0
            num_genotypes = 0

            for sample in populations[pop]:
                if sample.genotype != "./.":
                    num_genotypes += 1
                    if sample.genotype in ["0/1", "1/0"]:
                        num_hetero += 1

            if num_genotypes > 0:
                proportion = float(num_hetero) / num_genotypes
            else:
                proportion = 0.0

            graph_dict[pop]["heterozygosity"][round(proportion, 3)] += 1
            graph_dict["global"]["heterozygosity"][round(proportion, 3)] += 1

# Fis
def test_fis(locus, pop_info, min_fis, max_fis, fis_joker):
    """Test if each population passes minimum and maximum Fis filter
    """

    for snp in locus.snps:
        # Get Sample objects for each populations for that SNP
        populations = {}

        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        # Test if SNP meets minimum criterion for all pops
        pops_passing_min_fis = 0
        pops_passing_max_fis = 0

        for pop in populations:
            allele_count = defaultdict(int)
            num_samples = 0
            num_hetero = 0

            for sample in populations[pop]:
                if sample.genotype != "./.":
                    num_samples += 1
                    if sample.genotype in ["0/1", "1/0"]:
                        num_hetero += 1

                    a1, a2 = sample.genotype.split("/")
                    allele_count[a1] += 1
                    allele_count[a2] += 1

            if num_samples > 0:
                num_alleles = num_samples * 2.0
                p = float(allele_count["0"]) / num_alleles
                q = float(allele_count["1"]) / num_alleles
                expected_hetero_freq = 2.0 * p * q
                observed_hetero_freq = float(num_hetero) / num_samples

                try:
                    fis = (expected_hetero_freq - observed_hetero_freq) / expected_hetero_freq
                except:
                    fis = 0.0

            else:
                fis = 0.0

            # Modify flag if pop fails min_fis
            if fis >= min_fis:
                pops_passing_min_fis += 1

            # Modify flag if pop fails max_fis
            if fis <= max_fis:
                pops_passing_max_fis += 1

        if pops_passing_min_fis < len(pop_info) - fis_joker:
            snp.flags.min_fis = False
            Flags.min_fis_count += 1

        if pops_passing_max_fis < len(pop_info) - fis_joker:
            snp.flags.max_fis = False
            Flags.max_fis_count += 1

def get_fis_data(graph_dict, locus, pop_info):
    """Get heterozygosity data by pop and globally
    """

    for snp in locus.snps:
        # Get Sample objects for each populations for that SNP
        populations = {}

        for pop in pop_info:
            populations[pop] = [snp.samples[i] for i in pop_info[pop]]

        for pop in populations:
            allele_count = defaultdict(int)
            num_samples = 0
            num_hetero = 0

            for sample in populations[pop]:
                if sample.genotype != "./.":
                    num_samples += 1
                    if sample.genotype in ["0/1", "1/0"]:
                        num_hetero += 1

                    a1, a2 = sample.genotype.split("/")
                    allele_count[a1] += 1
                    allele_count[a2] += 1

            if num_samples > 0:
                num_alleles = num_samples * 2.0
                p = float(allele_count["0"]) / num_alleles
                q = float(allele_count["1"]) / num_alleles
                expected_hetero_freq = 2.0 * p * q
                observed_hetero_freq = float(num_hetero) / num_samples

                try:
                    fis = (expected_hetero_freq - observed_hetero_freq) / expected_hetero_freq
                except:
                    fis = 0.0
            else:
                fis = 0.0
                
            graph_dict[pop]["fis"][round(fis, 3)] += 1
            graph_dict["global"]["fis"][round(fis, 3)] += 1

# Maximum allele coverage
def test_max_allele_coverage(locus, pop_info, max_allele_coverage):
    """Test if median coverage of SNP is above max_allele_coverage
    """

    for snp in locus.snps:
        coverages = []
        for sample in snp.samples:
            # calculate allele coverage for both alleles
            # change genotypes that do not meet threshold to './.'
            if sample.genotype != "./.":
                coverages.append(sample.depth)

        if not coverages or median(coverages) >= max_allele_coverage:
            Flags.max_allele_coverage_count += 1

# Max number of SNPs
def test_max_snp_number(locus, pop_info, max_snp_number):
    """Test if there are too many SNPs at one locus
    """

    for snp in locus.snps:
        if len(locus.snps) > max_snp_number:
            snp.flags.max_snp_number = False
            Flags.max_snp_number_count += 1

def get_numSNP_data(graph_dict, locus, pop_info):
    """Get number of SNPs per locus data
    """

    graph_dict["global"]["numSNP"][len(locus.snps)] += 1

def get_individual_coverage(locus, ind_info):
    """Count the number of missing data per indivual
    """

    for snp in locus.snps:
        ind_info["total"] += 1

        for i in range(len(snp.samples)):
            sample = snp.samples[i]

            if sample.genotype == "./.":
                ind_info[i][2] += 1

        break

# Main
if __name__ == '__main__':

    # Parsing user input with 'argparse'
    parser = argparse.ArgumentParser(description=
            """Filter SNPs from STACKS VCF result file
            (WARNING! Will not work properly with stacks versions older to v1.30)
            """)

    parser.add_argument("-q", "--quiet", action="store_true",
            help = "do not print progress messages")
    parser.add_argument("-i", "--input_file", type=str, required=True,
            help = "input VCF file (can be compressed with gzip, ending in .gz)")
    parser.add_argument("-o", "--output_file", type=str, required=True,
            help = "output VCF file (can be compressed with gzip, ending in .gz) or name of directory if you use -g option")
    parser.add_argument("-g", "--graphs", action="store_true",
            help = "produce parameter distribution graphs instead of filtering")
    parser.add_argument("-c", "--min_allele_coverage", type=int, default=0,
            help = "minimum allele coverage (rare allele for heterozygotes, global coverage for homozygotes) to keep a genotype (or modified to './.') (int, default: 0)")
    parser.add_argument("-m", "--min_depth", type=int, default=0,
            help = "minimum depth (global depth by adding coverage of both alleles) to keep a genotype (or modified to './.') (int, default: 0)")
    parser.add_argument("-I", "--max_allelic_imbalance", type=float, default=1000.0,
            help = "maximum coverage fold change among alleles in heterozygotes (or modified to './.') (float, 0.0 or more, default 1000.0)")
    parser.add_argument("-C", "--max_allele_coverage", type=int, default=10000,
            help = "maximum median allele depth to keep a SNP (int, default: 10000)")
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
    parser.add_argument("-S", "--mas", type=int, default=2,
            help = "minimum number of samples with rare allele to retain snp (int, 1 or more, default: 0)")
    parser.add_argument("-H", "--max_hetero", type=float, default=1.0,
            help = "maximum proportion of heterozygous individuals (float, 0 to 1, default: 1)")
    parser.add_argument("-y", "--max_hetero_joker", type=int, default=0,
            help = "number of populations where it is permitted that the -H threshold does not pass, (int, 0 or more, default: 0")
    parser.add_argument("-f", "--min_fis", type=float, default=-1.0,
            help = "minimum Fis value (float, -1 to 1, default: -1)")
    parser.add_argument("-F", "--max_fis", type=float, default=1.0,
            help = "maximum Fis value (float, -1 to 1, default: 1)")
    parser.add_argument("-z", "--fis_joker", type=int, default=0,
            help = "number of populations where it is permitted that the -f or -F thresholds do not pass, (int, 0 or more, default: 0")
    parser.add_argument("-s", "--max_snp_number", type=int, default=999,
            help = "maximum number of SNPs per locus (int, 1 or more, default: 999)")
    args = parser.parse_args()

    # Assert proper values for parameters
    assert args.max_snp_number > 0, "max_snp_number must be a non-null integer"
    assert 0 <= args.min_allele_coverage <= 100, "min_allele_coverage must be an integer between 0 and 100"
    assert args.max_allelic_imbalance >= 0.0, "max_allelic_imbalance must be a floating point number equal to or greater than 0.0"
    assert args.max_allele_coverage >= 1, "max_allele_coverage must be an integer equal to or greater than 1"
    assert 0 <= args.min_presence <= 100000, "min_presence must be an integer between 0 and 100" 
    assert args.min_presence_joker_populations >= 0, "min_presence_joker_populations must be an integer that is 0 or more"
    assert 0 <= args.max_hetero <= 1, "max_hetero must be a floating point number between 0 and 1"
    assert 0 <= args.maf_global <= 1.0, "maf_global must be a floating point number between 0 and 1.0"
    assert 1 <= args.mas, "mas must be an integer of 1 or more"
    assert 0 <= args.maf_population <= 1.0, "maf_population must be a floating point number between 0 and 1.0"
    assert -1 <= args.min_fis <= 1, "min_fis must be a floating point number between -1 and 1"
    assert -1 <= args.max_fis <= 1, "max_fis must be a floating point number between -1 and 1"
    assert 0 <= args.fis_joker <= 100, "fis_joker must be an integer contained between 0 and 100"
    assert args.max_snp_number >= 1, "max_snp_number must be an integer equal to or greater than 1"

    # Get header from VCF and population information
    header = []
    print("Treating: " + args.input_file)
    stacks_version = None

    with myopen(args.input_file) as in_file:
        for line in in_file:
            l = line.strip()

            if l.startswith("##"):
                header.append(line)

                # Get STACKS version
                if l.startswith("##source="):
                    stacks_version = l.split('"')[1].split(" ")[1]

                    if "v1" in stacks_version:
                        print(f"  STACKS version: {stacks_version}")
                        locus_separator = "_"
                        stacks_version = 1

                    elif "v2" in stacks_version:
                        print(f"  STACKS version: {stacks_version}")
                        locus_separator = ":"
                        stacks_version = 2

            elif l.startswith("#CHROM"):

                if not stacks_version:
                    print("Error: STACKS version not recognized")
                    sys.exit(1)

                header.append(line)
                pop_info = get_population_info(l)
                assert len(pop_info) > 0, "Input file does not contain a header"
                ind_info = get_individual_info(l)
            else:
                break

    # Producing graphs
    if args.graphs:
        print("Collecting data to produce distribution graphs")

        # Initializing dictionary
        graph_dict = {}
        for pop in list(pop_info.keys()) + ["global"]:
            graph_dict[pop] = {}
            graph_dict[pop]["medDepth"] = defaultdict(int)
            graph_dict[pop]["maxDepth"] = defaultdict(int)
            graph_dict[pop]["allImbalance"] = defaultdict(int)
            graph_dict[pop]["genLikelihood"] = defaultdict(int)
            graph_dict[pop]["presence"] = defaultdict(int)
            graph_dict[pop]["mafPopulation"] = defaultdict(int)
            graph_dict[pop]["mafGlobal"] = defaultdict(int)
            graph_dict[pop]["heterozygosity"] = defaultdict(int)
            graph_dict[pop]["fis"] = defaultdict(int)
            graph_dict[pop]["numSNP"] = defaultdict(int)

        # For debugging purposes
        report_every = 100
        locus_counter = 1
        max_loci = 9999999999
        #max_loci = 100

        # Iterate over the loci and filter the snps
        for locus in locus_iterator(args.input_file, locus_separator):

            # For debugging purposes
            if locus_counter >= max_loci:
                break
            else:
                locus_counter += 1

            # Reporting progress
            if not args.quiet:
                if locus_counter % report_every == 0:
                    print("    Treating locus number: " + str(locus_counter))

            # Getting graph data
            get_depth_data(graph_dict, locus, pop_info)
            get_allelic_imbalance_data(graph_dict, locus, pop_info)
            get_presence_data(graph_dict, locus, pop_info)
            get_maf_population_data(graph_dict, locus, pop_info)
            get_maf_global_data(graph_dict, locus, pop_info)
            get_heterozygosity_data(graph_dict, locus, pop_info)
            get_fis_data(graph_dict, locus, pop_info)
            get_numSNP_data(graph_dict, locus, pop_info)
            get_individual_coverage(locus, ind_info)

        # Create folder with graphs
        directory = args.output_file
        subdirectory_1 = os.path.join(directory, "global")
        subdirectory_2 = os.path.join(directory, "populations")
        subdirectory_3 = os.path.join(directory, "missing_data")

        if not os.path.exists(directory):
            os.makedirs(directory)

        if not os.path.exists(subdirectory_1):
            os.makedirs(subdirectory_1)

        if not os.path.exists(subdirectory_2):
            os.makedirs(subdirectory_2)

        if not os.path.exists(subdirectory_3):
            os.makedirs(subdirectory_3)

        # Write missing data info to file and create figure
        missing_file = os.path.join(directory, "missing_data.tsv")
        total = ind_info["total"]
        del ind_info["total"]
        with open(missing_file, "w") as mfile:
            mfile.write("Population\tSample\tProportionMissing\n")
            for i in sorted(ind_info):
                ind_info[i][2] = round(float(ind_info[i][2]) / float(total), 3)
                mfile.write("\t".join([str(x) for x in ind_info[i]]) + "\n")

        # Exporting graph data
        graph_file = os.path.join(directory, "graph_data.tsv")
        graph_data_temp = []
        with open(graph_file, "w") as gfile:
            for pop in graph_dict:
                for param in graph_dict[pop]:
                    for k,v in graph_dict[pop][param].items():
                        graph_data_temp.append("\t".join([param, pop, str(k), str(v)]) + "\n")

            graph_data_temp.sort()

            gfile.write("\t".join(["Parameter", "Population", "Value", "Count"]) + "\n")
            for g in graph_data_temp:
                gfile.write(g)

        # Finished producing graphs, quiting
        with open(".temp_graph_folder", "w") as gf:
            gf.write(directory + "\n")

        subprocess.call("R -q -e 'source(\"00-scripts/utility_scripts/distribution_graphs.r\")' > /dev/null",
                shell=True)

        # Creating figures of missing data per population and individual
        subprocess.call("R -q -e 'source(\"00-scripts/utility_scripts/missing_data_graphs.r\")' > /dev/null",
                shell=True)

        print("Distribution graphs were writen in folder:\n  '{}'".format(directory))
        sys.exit(0)

    # Filtering
    # Open output files handles
    out_file = open (args.output_file, "w")
    filters_file = open(args.output_file + "_filters.tsv", "w")
    whitelist_file = open(args.output_file + "_whitelist.txt", "w")
    report_file = open(args.output_file + "_report.txt", "w")
    report_file.write(" ".join(sys.argv) + "\n\n")

    # Writing header to filtered VCF file
    out_file.writelines(header)

    # Writing header to filters report file
    filters_file.write("\t".join(["Position", "ID", "MaxAlleleCoverage",
                                  "MinPresence", "MafGlobal", "MafPopulation",
                                  "Hetero", "FisMin", "FisMax",
                                  "MaxSnpNumber"]) + "\n")

    # For debugging purposes
    report_every = 100
    locus_counter = 1
    max_loci = 9999999999
    #max_loci = 200

    # Global variable used by Flags.report_filters
    total_good_loci = 0

    # Iterate over the loci and filter the SNPs
    for locus in locus_iterator(args.input_file, locus_separator):

        # For debugging purposes
        if locus_counter >= max_loci:
            break
        else:
            locus_counter += 1

        # Reporting progress
        if not args.quiet:
            if locus_counter % report_every == 0:
                print("  Treating locus number: " + str(locus_counter))

        # Run filters
        # The filter functions automatically update SNP flags
        test_min_allele_coverage(locus, pop_info, args.min_allele_coverage)
        test_min_depth(locus, pop_info, args.min_depth)
        test_max_allelic_imbalance(locus, pop_info, args.max_allelic_imbalance)
        test_max_allele_coverage(locus, pop_info, args.max_allele_coverage)
        test_min_presence(locus, pop_info, args.min_presence,
                          args.min_presence_joker_populations, args.use_percent)
        test_maf_global(locus, args.maf_global)
        test_maf_population(locus, pop_info, args.maf_population)
        test_mas(locus, args.mas)
        test_heterozygosity(locus, pop_info, args.max_hetero, args.max_hetero_joker)
        test_fis(locus, pop_info, args.min_fis, args.max_fis, args.fis_joker)
        test_max_snp_number(locus, pop_info, args.max_snp_number)

        # Write output
        write_filters(locus, filters_file)
        write_locus(locus, out_file)
        write_whitelist(locus, whitelist_file)

    # Print filtering report
    Flags.report_filters(locus_counter, total_good_loci, report_file)

    out_file.close()
    filters_file.close()
    whitelist_file.close()
    report_file.close()
