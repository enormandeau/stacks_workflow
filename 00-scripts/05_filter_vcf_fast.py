#!/usr/bin/env python3
"""Filtering SNPs in VCF file output by STACKS1 or STACKS2 minimaly

Usage:
    <program> input_vcf min_cov max_missing min_mas output_vcf

Where:
    input_vcf: is the name of the VCF file to filter
    min_cov: minimum allele coverage to keep genotype <int>, eg: 4 or more
    max_missing: maximum proportion of missing data (applies to ALL populations) <float> eg: 0.5 to 1.0
    min_mas: minimum number of samples with rare allele <int> eg: 2 or more
    output_vcf: is the name of the filtered VCF

WARNING:
    The filtering is done purely on a SNP basis. Loci are not taken into account.
"""

# Import
import sys

# Functions
def get_population_info(line):
    """Return dictionary of population names with a list of their sample
    positions in the lines of the VCF file.
    """

    pops = [x.split("_")[0] for x in line.split("\t")[9:]]
    unique_pops = sorted(list(set(pops)))
    print("  " + str(len(unique_pops)) + " populations")
    pop_dict = {}

    for p in unique_pops:
        pop_dict[p] = []

    for i, p in enumerate(pops):
        pop_dict[p].append(i)

    return pop_dict

def correct_genotype(genotype_info, min_cov):
    """Correct genotype to ./. if coverage < min_cov
    """
    infos = genotype_info.split(":")

    if infos[0] == "./.":
        return genotype_info

    cov = int(infos[1])

    if cov >= min_cov:
        return genotype_info

    else:
        return ":".join(["./."] + infos[1:])

# Parse user input
try:
    input_vcf = sys.argv[1]
    min_cov = int(sys.argv[2])
    max_missing = float(sys.argv[3])
    min_mas = int(sys.argv[4])
    output_vcf = sys.argv[5]

except:
    print(__doc__)
    sys.exit(1)

with open(input_vcf) as infile:
    with open(output_vcf, "w") as outfile:
        for line in infile:
            l = line.strip().split("\t")

            # Header info
            if line.startswith("##"):
                outfile.write(line)
                continue

            # Sample names
            if line.startswith("#CHROM"):
                outfile.write(line)
                pop_info = get_population_info(line)
                continue

            # SNP line split into info and genotypes
            infos = l[:9]
            genotypes = l[9:]

            # Correct genotypes with coverage below min_cov
            genotypes = [correct_genotype(x, min_cov) for x in genotypes]

            # Create corrected line
            line = "\t".join(infos + genotypes) + "\n"

            # Remove SNPs with MAS below threshold
            mas = len([1 for x in genotypes if x.split(":")[0] in ["0/1", "1/0", "1/1"]])

            if mas < min_mas:
                continue

            # Remove SNPs with too much missing data in at least one population
            max_missing_failed = False

            for pop in pop_info:
                sample_ids = pop_info[pop]
                num_samples = len(sample_ids)
                samples = [genotypes[i] for i in sample_ids]
                num_missing = len([1 for x in samples if x.split(":")[0] == "./."])
                prop_missing = num_missing / num_samples

                if prop_missing > max_missing:
                    max_missing_failed = True
                    break

            if not max_missing_failed:
                outfile.write(line)
