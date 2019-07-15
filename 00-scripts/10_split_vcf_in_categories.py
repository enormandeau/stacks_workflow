#!/usr/bin/env python3
"""Split VCF into singleton vs other duplicated categories

Usage:
    <program> input_vcf input_categories
    Input and output VCFs can be compressed with gzip, ending in .gz
"""

# Modules
from collections import defaultdict
from collections import Counter
import gzip
import sys

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

# Parse user input
try:
    input_vcf = sys.argv[1]
    input_categories = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read categories
categories = set()
snps = dict()
loci = defaultdict(list)

with open(input_categories) as infile:
    for line in infile:
        if line.startswith("Scaffold"):
            continue

        scaffold, position, snp, category = line.strip().split()
        categories.add(category)
        locus = snp.split("_")[0]

        snps[(scaffold, position, snp)] = category
        loci[locus].append(category)

# Open output file handles
output_vcfs = dict()

for category in categories:
    if input_vcf.endswith(".gz"):
        compression = ".gz"
    else:
        compression = ""

    output_vcfs[category] = myopen(input_vcf.replace(".vcf", "").replace(".gz", "") + "." + category + ".vcf" + compression, "wt")

# Read and split VCF
with myopen(input_vcf) as infile:
    for line in infile:
        if line.startswith("#"):
            for output in output_vcfs:
                output_vcfs[output].write(line)

            continue

        # Write SNPs in their respective output files
        l = line.strip().split()
        scaffold, position, snp = l[:3]
        category = snps[(scaffold, position, snp)]
        output_vcfs[category].write(line)
