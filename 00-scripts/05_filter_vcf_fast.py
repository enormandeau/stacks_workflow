#!/usr/bin/env python3
"""Filtering SNPs in VCF file output by STACKS1 or STACKS2 minimaly

Usage:
    <program> input_vcf min_cov percent_genotypes max_pop_fail min_mas min_maf output_vcf [group_all]

Where:
    input_vcf: is the name of the VCF file to filter (can be compressed with gzip, ending in .gz)
    min_cov: minimum allele coverage to keep genotype <int>, eg: 4 or more
    percent_genotypes: minimum percent of genotype data per population <float> eg: 50, 70, 80, 100
    max_pop_fail: maximum number of populations that can fail percent_genotypes <int> eg: 1, 2, 3
    min_mas: minimum number of samples with rare allele (MAS) <int> eg: 2 or more
    min_maf: minimum Minor Allele Frequency (MAF) <float> eg: from 0 to 1 inclusively
    output_vcf: is the name of the filtered VCF (can be compressed with gzip, ending in .gz)
    group_all: whether to consider all samples as one population <int> 0, 1 (default: 0)

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

def get_population_info(line, group_all):
    """Return dictionary of population names with a list of their sample
    positions in the lines of the VCF file.
    """

    pops = [x.split("_")[0] for x in line.split("\t")[9:]]

    if group_all:
        for i, p in enumerate(pops):
            pops[i] = "pop1"

    unique_pops = sorted(list(set(pops)))
    print(str(len(pops)) + " samples in " + str(len(unique_pops)) + " populations")
    pop_dict = {}

    for p in unique_pops:
        pop_dict[p] = []

    for i, p in enumerate(pops):
        pop_dict[p].append(i)

    return pop_dict

def correct_genotype(genotype_info, min_cov):
    """Correct genotype to ./. if coverage < min_cov
    """
    global counter_genotypes
    infos = genotype_info.split(":")

    if infos[0] == "./.":
        return genotype_info

    cov = int(infos[1])

    if cov >= min_cov:
        counter_genotypes += 1
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
    min_maf = float(sys.argv[6])
    output_vcf = sys.argv[7]

except:
    print(__doc__)
    sys.exit(1)

# group_all
try:
    group_all = int(sys.argv[8])
except:
    group_all = 0

assert group_all in [0, 1], "Only zero (0) or one (1) are accepted for group_all"

# Validate parameter values
assert min_cov >= 0, "min_cov needs to be zero or a positive integer"
assert percent_genotypes >= 0 and percent_genotypes <= 100.0, "percent_genotypes needs to be a number between zero and 100"
assert max_pop_fail >= 0, "max_pop_fail needs to be a null or positive integer"
assert min_mas >= 1, "min_mas needs to be a non-null positive integer"

# Loop over VCF
counter = 0
counter_missing = 0
counter_mas = 0
counter_maf = 0
counter_genotypes = 0
counter_retained = 0

with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
        print()

        for line in infile:
            l = line.strip().split("\t")

            # Header info
            if line.startswith("##"):
                outfile.write(line)
                continue

            # Sample names
            if line.startswith("#CHROM"):
                outfile.write(line)
                pop_info = get_population_info(line, group_all)
                continue

            # Print progress
            #counter += 1
            #if not counter % 10000:
            #    print(f"Treating SNP number: {counter}")

            # SNP line split into info and genotypes
            infos = l[:9]
            genotypes_raw = l[9:]

            # Correct genotypes_raw with coverage below min_cov
            genotypes_raw = [correct_genotype(x, min_cov) for x in genotypes_raw]
            genotypes = [x.split(":")[0] for x in genotypes_raw]

            # MISSING DATA
            # Remove SNPs with too much missing data in too many populations
            pops_failed = 0
            max_missing_failed = False

            for pop in pop_info:
                sample_ids = pop_info[pop]
                num_samples = len(sample_ids)
                samples = [genotypes[i] for i in sample_ids]
                num_missing = len([x for x in samples if x == "./."])
                prop_missing = num_missing / num_samples

                if prop_missing > 1 - (percent_genotypes / 100):
                    pops_failed += 1

                    if pops_failed > max_pop_fail:
                        max_missing_failed = True
                        counter_missing += 1
                        break

            if max_missing_failed:
                continue

            # MAF
            # Remove SNPs with MAF below min_maf
            num_alleles = 2 * len([x for x in genotypes if not x in ["./.", "."]])
            rare_alleles = len([x for x in genotypes if x in ["0/1", "1/0"]]) \
                    + 2 * len([x for x in genotypes if x == "1/1"])
            maf = rare_alleles / num_alleles

            if maf > 0.5:
                maf = 1 - maf

            if maf < min_maf:
                counter_maf += 1
                continue

            # MAS
            # Remove SNPs with MAS below threshold
            mas = len([1 for x in genotypes if x in ["0/1", "1/0", "1/1"]])

            non_null_genotypes = [x for x in genotypes if not x in ["./.", "."]]

            # We may be filtering a VCF that is a subset of a larger VCF file
            # where a SNP could be 100% homozygote for the rare allele in the
            # samples we kept, even if its MAF was globally less than 0.5 in
            # the original VCF. As a result, we do the adjustment below.
            if mas > len(non_null_genotypes) / 2:
                mas = len([1 for x in genotypes if x in ["0/0", "0/1", "1/0"]])
                if mas < min_mas:
                    print(non_null_genotypes)

            if mas < min_mas:
                counter_mas += 1
                continue

            # Create corrected line
            counter_retained += 1
            line = "\t".join(infos + genotypes_raw) + "\n"
            outfile.write(line)

# Reporting on filtration
print(f"Genotypes with insufficient coverage (min {min_cov}): {counter_genotypes}")
print("Number of SNPs filtered by reason:")
print(f"  > Proportion non-missing data (min {percent_genotypes}): {counter_missing}")
print(f"  > MAS (min {min_mas}): {counter_mas}")
print(f"  > MAF (min {min_maf}): {counter_maf}")
print(f"Number of retained SNPs: {counter_retained}")
print()
