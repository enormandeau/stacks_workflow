#!/usr/bin/env python3
"""Report pairwise distances of samples in a filtered VCF from STACKS

Usage:
    <program> input_vcf output_distances
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
    """Return distance between two sets of genotypes

    Using only SNPs where neither sample has missing data (./.), return the sum
    of different genotypes divided by the number of considered SNPs.
    """
    ndiff = 0
    nkept = 0

    for i, g1 in enumerate(geno1):
        g2 = geno2[i]
        if "./." not in [g1, g2]:
            nkept += 1

            if g1 != g2:
                ndiff += 1

    return ndiff / nkept

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_distances = sys.argv[2]
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


# Compute and output distances
computed = set()
with open(output_distances, "wt") as outfile:
    for s1 in sample_data:
        print(s1)
        for s2 in sample_data:

            # Avoid computing twice
            if tuple(sorted([s1, s2])) in computed:
                continue

            computed.add(tuple(sorted([s1, s2])))

            geno1 = sample_data[s1]
            geno2 = sample_data[s2]
            dist = distance(geno1, geno2)
            outfile.write("\t".join([s1, s2, "{0:0.6f}".format(dist)]) + "\n")

            if s1 != s2:
                outfile.write("\t".join([s2, s1, "{0:0.6f}".format(dist)]) + "\n")
