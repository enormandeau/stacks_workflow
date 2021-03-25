#!/usr/bin/env python3

#########################################
#########################################
#########################################
#  WARNING! This script is not finished #
#########################################
#########################################
#########################################

"""Return the median AFD per pairs of populations

Usage:
    <program> input_vcf output_file
"""

# Modules
from collections import defaultdict
import gzip
import sys

# Defining functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

# Parsing user input
try:
    input_vcf = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read VCF
with myopen(input_vcf) as infile:
    for line in infile:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            samples = line.strip().split("\t")[9:]
            pops = enumerate([x.split("_")[0] for x in samples])
            pop_info = defaultdict(list)

            for p in pops:
                pop_info[p[1]].append(p[0])

            print(pops)
            print(pop_info)
            for p in pop_info:
                print(len(pop_info[p]), p)
            break
