#!/usr/bin/env python3
"""Convert STACKS VCF to simplified csv file with genotypes

Usage:
    <program> input_vcf output_csv
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
    output_csv = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Convert
with myopen(input_vcf, "rt") as infile:
    with myopen(output_csv, "wt") as outfile:
        for line in infile:
            if line.startswith("##"):
                continue
            else:
                if line.startswith("#"):
                    line = line[1: ]

                l = line.strip().split("\t")
                to_write = l[: 5]

                data = l[9: ]

                if ":" in line:
                    data = [x.split(":")[0] for x in data]

                to_write += data

                if "," in "".join(to_write):
                    print("ERROR: Comma in data line, problem for csv file")
                    sys.exit(1)

                outfile.write(",".join(to_write) + "\n")
