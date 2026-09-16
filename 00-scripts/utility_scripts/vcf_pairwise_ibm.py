#!/usr/bin/env python3
"""Report pairwise ibm distances of samples in a filtered VCF from STACKS

Usage:
    <program> input_vcf output_ibm
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

def distance(geno1, geno2):
    """Return ibm distance between two sets of genotypes

    Using all SNPs, return the number of missing/non-missing cases divided by
    the number of considered SNPs.
    """
    ndiff = 0
    nkept = 0

    for i, g1 in enumerate(geno1):
        g2 = geno2[i]

        # Use only SNPs where at least one sample has missing data
        if (g1 == "./.") or (g2 == "./."):
            nkept += 1

            # Use XOR to test missing
            if (g1 == "./.") ^ (g2 == "./."):
                ndiff += 1

    try:
        return ndiff / nkept
    except:
        return 0.0

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_ibm = sys.argv[2]
except:
    print(__doc__)
    sys.exit()

# Read info from VCF
sample_data = dict()
sample_from_pos = dict()

with myopen(input_vcf, "rt") as infile:
    for line in infile:
        l = line.strip().split()
        data = l[9:]

        if line.startswith("##"):
            continue

        elif line.startswith("#CHROM"):
            for i, sample in enumerate(data):
                sample_data[sample] = []
                sample_from_pos[i] = sample
            continue

        for i, genotype in enumerate(data):
            genotype = genotype.split(":")[0]
            sample_data[sample_from_pos[i]].append(genotype)


# Compute and output ibm
computed = set()
n = len(sample_data)
with open(output_ibm, "wt") as outfile:
    for i, s1 in enumerate(sample_data):

        for s2 in sample_data:

            # Avoid computing twice
            if tuple(sorted([s1, s2])) in computed:
                continue

            computed.add(tuple(sorted([s1, s2])))

            geno1 = sample_data[s1]
            geno2 = sample_data[s2]
            ibm = distance(geno1, geno2)
            outfile.write("\t".join([s1, s2, "{0:0.6f}".format(ibm)]) + "\n")

            if s1 != s2:
                outfile.write("\t".join([s2, s1, "{0:0.6f}".format(ibm)]) + "\n")

        # Report progress
        x = i+1
        print(f"{x}/{len(sample_data)} ({round( 200* (x * (n - x) + x**2 / 2) / n**2, 2)}%): {s1}")

