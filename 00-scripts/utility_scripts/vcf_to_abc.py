#!/usr/bin/env python
"""Create ABC input from filtered haplotype VCF file (from STACKS)

WARNING! will not work properly with STACKS versions older to v1.30 since the
minor allele counts in the VCF file created by STACKS were not reported
properly at that time.
"""

# Import
from collections import defaultdict
from collections import Counter
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
        self.info = info.split(":")
        self.genotype = self.info[0]
        try:
            self.depth = int(self.info[1].split(",")[0])
        except:
            self.genotype = "./."
            self.depth = "0"

    def __repr__(self):
        return self.genotype

class SNP(object):
    def __init__(self, line):
        """Parse SNP information from one SNP line in the VCF file
        """
        
        self.line = line.strip().split("\t")
        self.chrom = self.line[0]
        self.pos = int(self.line[1])
        self.locus_id = self.line[2].split("_")[0]
        self.ref = self.line[3]
        self.alt = self.line[4]
        self.qual = self.line[5]
        self.filter = self.line[6]
        self.info = self.line[7]
        self.global_maf = self.info.split(";")[1]

        try:    # Version 1.41-
            self.global_maf = float(self.global_maf.split("=")[1].split(",")[1])
        except: # Version 1.42+
            self.global_maf = float(self.global_maf.split("=")[1])

        self.format = self.line[8]
        self.samples = [Sample(x) for x in self.line[9:]]

    def __repr__(self):
        return " ".join([
            self.chrom,
            str(self.pos),
            self.locus_id,
            self.ref,
            self.alt,
            self.qual,
            self.filter,
            self.info,
            str(self.global_maf)
            ])

class Locus(object):
    """A locus is a list of SNP objects and methods to manipulate them
    """

    def __init__(self, snp_list):
        self.snps = snp_list

    def __repr__(self):
        #return "There are {} snps in this locus".format(len(self.snps))
        return str(self.snps)

# Functions
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
    print "(" + str(len(unique_pops)) + " populations)"
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
        snps = []
        first = True
        last_id = -99

        for line in in_f:
            if line.startswith("#"):
                continue

            elif "\t\t" in line:
                continue

            elif line.strip():
                snp = SNP(line.strip())
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
            help = "input VCF file")

    parser.add_argument("-p", "--pops", type=str, required=True,
            help = "ID of populations, separated by a comma")

    parser.add_argument("-l", "--locus_length", type=int, default=80,
            help = "length of loci in VCF file (int, default: 80)")

    parser.add_argument("-n", "--effective_size", type=int, default=10000,
            help = "population effective size (int, default: 10000)")

    parser.add_argument("-M", "--min_pop_size", type=int, default=4,
            help = "minimum number of samples in each population (int, default: 4)")

    parser.add_argument("-m", "--mutation_rate", type=float, default=10e-8,
            help = "mutation rate (int, default: 10e-8)")

    # Compile parser
    args = parser.parse_args()

    # Get header from VCF and population information
    header = []
    print "Treating: " + args.input_file, 
    with open(args.input_file) as in_file:
        for line in in_file:
            l = line.strip()

            if l.startswith("##"):
                header.append(line)
            elif l.startswith("#CHROM"):
                header.append(line)
                pop_info = get_population_info(l)
                assert len(pop_info) > 0, "Input file does not contain a header"
                ind_info = get_individual_info(l)
            else:
                break

    # Open output files handles
    output_stub = args.pops.replace(",", "_")
    locus_ms_file = open (output_stub + ".locus.ms", "w")
    bp_file = open (output_stub + ".bpfile", "w")
    spinput_file = open (output_stub + ".spinput.txt", "w")

    # create locus.ms output header
    locus_ms_file.write("../../bin/msnsam tbs 15000 -t 0.08 -r 0.06 100 -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs\n")
    locus_ms_file.write("98128 624518 817721\n")

    # Keep infos for output
    num_sample_per_snp = []
    kept_loci = set()

    # Iterate over the loci and filter the SNPs
    for locus in locus_iterator(args.input_file):
        # Correct genotypes for samples who do not have genotypes for all snps of locus
        
        locus_snps = []
        snp = locus.snps[0]
        alleles = [snp.ref] + snp.alt.split(",")

        # keep only wanted samples for each SNP
        populations = []
        num_samples = 0
        num_samples_per_pop = []

        for pop in args.pops.split(","):

            # also remove samples with missing genotypes
            samples = [snp.samples[i] for i in pop_info[pop] if snp.samples[i].genotype != "./."]
            populations += samples
            num_samples += len(samples)
            num_samples_per_pop.append(len(samples))

        num_snps = len(alleles[0])
        short_genotypes = []

        haplotypes = []
        for sample in populations:
            genotype = sample.genotype
            for allele in genotype.split("/"):
                a = alleles[int(allele)]
                haplotypes.append(a)

        if len(set(haplotypes)) > 1 and min(num_samples_per_pop) >= args.min_pop_size:
            num_sample_per_snp.append(num_samples_per_pop)
            kept_loci.add(snp.locus_id)

            # Encode alleles in 0-1
            encoded = []
            for h in haplotypes:
                encoded.append("")
                for i in xrange(len(h)):
                    position = [x[i] for x in haplotypes]
                    if len(set(position)) > 1:

                        ref_allele = alleles[0][i]

                        if h[i] == ref_allele:
                            encoded[-1] += "0"
                        else:
                            encoded[-1] += "1"

            # Write locus header lines
            locus_ms_file.write("\n//	160	80	80	10.00000	5.00000	1.50000	2.00000	5.00000	5.00000	2.00000\n")
            locus_ms_file.write("segsites: " + str(len(encoded[-1])) + "\n")
            positions = [str(float(x) / 20)[:6] for x in range(1, 20)]
            positions = positions[:len(encoded[-1])]

            locus_ms_file.write("positions: " + " ".join(positions) + "\n")
            locus_ms_file.write("\n".join(encoded) + "\n")

    # Create bpfile
    bp_file.write("#" + " ".join(args.pops.split(",")) + " " + str(args.effective_size) + "\n")
    line1 = []
    line2 = []
    line3 = []
    line4 = []

    for n in num_sample_per_snp:
        line1.append(str(args.locus_length))
        line2.append(str(n[0]*2))
        line3.append(str(n[1]*2))
        line4.append(str(4.0 * args.locus_length * args.effective_size * args.mutation_rate))

    bp_file.write("\t".join(line1) + "\n")
    bp_file.write("\t".join(line2) + "\n")
    bp_file.write("\t".join(line3) + "\n")
    bp_file.write("\t".join(line4) + "\n")
    bp_file.write("\t".join(line4) + "\n")

    # create spinput
    spinput_file.write("\n")
    spinput_file.write(str(len(kept_loci)) + "\n")

    for i in xrange(len(line1)):

        spinput_file.write(str(line2[i]) + "\n")
        spinput_file.write(str(line3[i]) + "\n")
        spinput_file.write(str(line1[i]) + "\n")

    spinput_file.write("1\n")
    spinput_file.write("locus.ms\n")

    # Close file handles
    locus_ms_file.close()
    bp_file.close()
    spinput_file.close()
