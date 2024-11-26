#!/usr/bin/env python3
"""Fix positions of vcf from split genome

Usage:
    <program> input_vcf chunk_size output_vcf

Input and output VCFs can be compressed with gzip

NOTE: Genome was split with:
    fasta_split_sliding_window.py genome.fasta 1000000 0 genome.split.fasta
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
    chunk_size = int(sys.argv[2])
    output_vcf = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Fix positions
with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
        for line in infile:

            # Skip split contig info
            if line.startswith("##contig"):
                continue

            # Write header info
            elif line.startswith("#"):
                outfile.write(line)
                continue

            l = line.strip().split("\t")

            # Fix one position
            chunk = int(l[0].split("chunk")[1])
            l[0] = l[0].split("_")[0]
            l[1] = str(int(l[1]) + (chunk - 1) * chunk_size)

            outfile.write("\t".join(l) + "\n")
