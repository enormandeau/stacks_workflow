#!/usr/bin/env python3
"""Compute Ho, He, and Fis for a whole VCF

Usage:
    <program> input_vcf
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

# Parsing user input
try:
    input_vcf = sys.argv[1]
except:
    print(__doc__)
    sys.exit(1)

# Read VCF and compute values
ho_all_snps = []
he_all_snps = []
fis_all_snps = []

with myopen(input_vcf, "rt") as infile:
    for line in infile:
        if line.startswith("#"):
            continue

        data = [x.split(":")[0] for x in line.strip().split("\t")[9:]]
        data = [x for x in data if not "." in x]
        num_samples = len(data)

        # Compute and update Ho
        het = [x for x in data if x in ["0/1", "1/0"]]
        ho = len(het) / num_samples
        ho_all_snps.append(ho)

        # Compute and update He
        alleles = "".join(data).replace("/", "")
        zeros = alleles.count("0")
        ones = alleles.count("1")
        maf = ones / (zeros + ones)
        he = 2 * maf * (1-maf)
        he_all_snps.append(he)

        # Compute Fis
        try:
            fis = 1 - ho / he
        except:
            fis = 0

        fis_all_snps.append(fis)
        print("\t".join([str(x) for x in [ho, he, fis]]))

#print(ho_all_snps)
#print(he_all_snps)
#print(fis_all_snps)
