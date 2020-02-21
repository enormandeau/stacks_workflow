#!/usr/bin/env python3
"""Rename scaffold names in STACKS VCF for plink

Usage:
    <program> input_vcf output_vcf
"""

# Modules
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
    output_vcf = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read VCF and modify scaffold names
seen_scaffold = set()
scaffold_number = 1000
with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:

        for line in infile:

            if line.startswith("#"):
                outfile.write(line)

            else:
                l = line.strip().split()
                scaffold_name = l[0]

                if scaffold_name in seen_scaffold:
                    l[0] = str(scaffold_number)
                    outfile.write("\t".join(l) + "\n")

                else:
                    seen_scaffold.add(scaffold_name)
                    scaffold_number += 1
                    l[0] = str(scaffold_number)
                    outfile.write("\t".join(l) + "\n")
