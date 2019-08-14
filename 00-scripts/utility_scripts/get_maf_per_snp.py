#!/usr/bin/env python3
"""Compute MAF per SNP from a STACKS2 VCF

Usage:
    <program> input_vcf output_file
"""

# Modules
from collections import Counter
import gzip
import sys

def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read VCF
with open(output_file, "wt") as outfile:
    with myopen(input_vcf) as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            l = line.strip().split()
            info = l[:3]
            data = l[9:]
            genotypes = Counter([x.split(":")[0] for x in data if x.split(":")[0] not in ["./.", "."]])
            num_genotypes = sum(genotypes.values())
            num_rare = genotypes["0/1"] + 2 * genotypes["1/1"]
            maf = num_rare / (2 * num_genotypes)
            outfile.write("\t".join(info + [str(maf)]) + "\n")
