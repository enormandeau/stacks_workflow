#!/usr/bin/env python3
"""Split VCF into singleton vs other duplicated categories

Usage:
    <program> input_vcf input_categories
"""

# Modules
from collections import defaultdict
from collections import Counter
import sys

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

cases = [", ".join(list(set(x))) for x in loci.values() if len(x) > 1]
c = Counter(cases)
for k, v in c.items():
    print(f"{k}: {v}")

# Open output file handles
output_vcfs = dict()

for cat in categories:
    output_vcfs[cat] = open(input_vcf.replace(".vcf", "") + "." + cat + ".vcf", "w")

# Read and split VCF
with open(input_vcf) as infile:
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
