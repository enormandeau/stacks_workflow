#!/usr/bin/env python3
"""Filtering SNPs in VCF file output by STACKS1 or STACKS2 minimaly

Usage:
    <program> input_vcf min_cov percent_genotypes max_pop_fail min_mas output_vcf

Where:
    input_vcf: is the name of the VCF file to filter (can be compressed with gzip, ending in .gz)
    min_cov: minimum allele coverage to keep genotype <int>, eg: 4 or more
    percent_genotypes: minimum percent of genotype data per population <float> eg: 50, 70, 80, 100
    max_pop_fail: maximum number of populations that can fail percent_genotypes <int> eg: 1, 2, 3
    min_mas: minimum number of samples with rare allele <int> eg: 2 or more
    output_vcf: is the name of the filtered VCF (can be compressed with gzip, ending in .gz)

WARNING:
    The filtering is done purely on a SNP basis. Loci are not taken into account.
"""

# Modules
import gzip
import sys

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def get_population_info(line):
    """Return dictionary of population names with a list of their sample
    positions in the lines of the VCF file.
    """

    pops = [x.split("_")[0] for x in line.split("\t")[9:]]
    unique_pops = sorted(list(set(pops)))
    print("  " + str(len(pops)) + " samples in " + str(len(unique_pops)) + " populations")
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
    percent_genotypes = float(sys.argv[3])
    max_pop_fail = int(sys.argv[4])
    min_mas = int(sys.argv[5])
    output_vcf = sys.argv[6]

except:
    print(__doc__)
    sys.exit(1)

# Validate parameter values
assert min_cov >= 0, "min_cov needs to be zero or a positive integer"
assert percent_genotypes >= 0 and percent_genotypes <= 100.0, "percent_genotypes needs to be a number between zero and 100"
assert max_pop_fail >= 0, "max_pop_fail needs to be a null or positive integer"
assert min_mas >= 1, "min_mas needs to be a non-null positive integer"

# Loop over VCF
counter = 0

with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
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

            # Print progress
            counter += 1
            if not counter % 10000:
                print(f"Treating locus number: {counter}")

            # SNP line split into info and genotypes
            infos = l[:9]
            genotypes = l[9:]

            # Correct genotypes with coverage below min_cov
            genotypes = [correct_genotype(x, min_cov) for x in genotypes]

            # Remove SNPs with MAS below threshold
            mas = len([1 for x in genotypes if x.split(":")[0] in ["0/1", "1/0", "1/1"]])

            # The second part of the test (after the or) is to take into
            # account that we may be filtering a VCF that is a subset of a
            # larger VCF file where a SNP could be 100% homozygote for the rare
            # allele in the samples we kept, even if its MAF was globally less
            # than 0.5 in the original VCF.
            non_null_genotypes = [x for x in genotypes if not x.split(":")[0] in ["./.", "."]]

            if mas < min_mas or mas > len(non_null_genotypes) - min_mas + 1:
                continue

            # Remove SNPs with too much missing data in too many populations
            pops_failed = 0
            max_missing_failed = False

            for pop in pop_info:
                sample_ids = pop_info[pop]
                num_samples = len(sample_ids)
                samples = [genotypes[i] for i in sample_ids]
                num_missing = len([1 for x in samples if x.split(":")[0] == "./."])
                prop_missing = num_missing / num_samples

                if prop_missing > 1 - (percent_genotypes / 100):
                    pops_failed += 1

                    if pops_failed > max_pop_fail:
                        max_missing_failed = True
                        break

            if not max_missing_failed:

                # Create corrected line
                line = "\t".join(infos + genotypes) + "\n"
                outfile.write(line)
