#!/usr/bin/env python3
"""Reorder samples in VCF

Usage:
    <program> input_vcf wanted_order output_vcf

Where:
    input_vcf is the name of the VCF file to filter (can be compressed with gzip, ending in .gz)
    wanted_order
    output_vcf is the name of the filtered VCF (can be compressed with gzip, ending in .gz)
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
    wanted_order = sys.argv[2]
    output_vcf = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Get wanted order
order = [x.strip() for x in open(wanted_order).readlines()]

# Filter
with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
        for line in infile:
            l = line.strip().split("\t")

            # Infos lines
            if line.startswith("##"):
                outfile.write(line)
                continue

            # Samples line
            elif line.startswith("#CHROM"):

                samples = l[9:]
                indexes = []

                for o in order:
                    i = samples.index(o)
                    indexes.append(i)

            # Data lines
            infos = l[:9]
            genotypes = l[9:]
            genotypes = [genotypes[i] for i in indexes]
            filtered_line = "\t".join(infos + genotypes) + "\n"
            outfile.write(filtered_line)
