#!/usr/bin/env python3
"""Keep only first SNP for each locus

Usage:
    <program> input_vcf output_vcf
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

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Keep first SNP per locus
seen_ids = set()

with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
        for line in infile:
            l = line.strip().split("\t")

            # Infos lines
            if line.startswith("#"):
                outfile.write(line)
                continue

            # Data lines
            locus_id = l[2].split(":")[0]

            if locus_id in seen_ids:
                continue

            seen_ids.add(locus_id)
            outfile.write(line)
