#!/usr/bin/env python3
"""Remove SNPs from bed ranges from VCF

WARNING:
    - Consumes about 32 GB of RAM per million ranges (avg range size: ~400 bp)

Usage:
    <program> input_vcf input_bed output_vcf
"""

# Modules
from collections import defaultdict
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
    input_bed = sys.argv[2]
    output_vcf = sys.argv[3]
except:
    print(__doc__)
    sys.exit()

# Parse bed file
ranges = defaultdict(set)

counter = 0
with open(input_bed, "rt") as infile:
    for line in infile:
        l = line.strip().split("\t")
        chrom = l[0]
        start = int(l[1])
        stop = int(l[2])

        counter += 1
        if not counter % 100000:
            print(counter, chrom, start, stop)

        # Also exclude preceding and following regions
        ranges[chrom].update(range(start, stop+1))

with myopen(output_vcf, "wt") as outfile:
    with myopen(input_vcf, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                l = line.strip().split("\t")
                chrom = l[0]
                pos = int(l[1])
                if not pos in ranges[chrom]:
                    outfile.write(line)
