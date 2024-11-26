#!/usr/bin/env python3
"""Remove ranges around unwanted positions from a bam file piped into this script

Usage:
    samtools view -h <INPUT.bam> | <this_program> bedfile | samtools view -Sb > <OUTPUT.bam>

Example with parallel:
    parallel -j 20 samtools view -h {} \| ../../00-scripts/utility_scripts/bam_exclude_positions.py \
            ../../highcov.bed \| samtools view -Sb \> ../{} ::: *.sorted.bam

NOTE: Don't forget to index with:
    samtools index *.bam
"""

# Modules
from collections import defaultdict
import sys

# Parse user input
try:
    bedfile = sys.argv[1]
except:
    print(__doc__)
    sys.exit(1)

# Load positions into range dict
size = 500
ranges = defaultdict(set)

with open(bedfile, "rt") as infile:
    for line in infile:
        l = line.strip().split("\t")
        chrom = l[0]
        pos = (int(l[1]) + int(l[2])) // 2
        unwanted = pos // size
        ranges[chrom].add(unwanted-1)
        ranges[chrom].add(unwanted)
        ranges[chrom].add(unwanted+1)

# Read bam and remove unwanted regions
for line in sys.stdin:
    if line.startswith("@"):
        print(line.rstrip())

    else:
        l = line.rstrip().split("\t")
        chrom = l[2]
        pos = int(l[3]) // size

        if pos not in ranges[chrom]:
            print(line.rstrip())
