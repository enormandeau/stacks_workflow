#!/usr/bin/env python3
"""Convert STACKS VCF to simplified csv file with genotypes

Usage:
    <program> input_vcf output_csv [genotype_symbols]

Where:
    genotype_symbols is a string or 3 symbols to be used to represent 0/0, 0/1 or 1/0, and 1/1
    Example1: "012"  (for frequent, homozygote, rare)
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

def correct_symbol(s, s1, s2, s3):
    assert s in ("0/0", "0/1", "1/0", "1/1", "./."), f"Bad genotype encoutered: {s}"

    if s == "0/0":
        return s1
    elif s in ("0/1", "1/0"):
        return s2
    elif s == "1/1":
        return s2
    elif s == "./.":
        return "NA"

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_csv = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

try:

    s1, s2, s3 = list(sys.argv[3])
    print(f"Using the following genotype symbols: {s1}, {s2}, {s3}")
except:
    print(__doc__)
    print()
    print("ERROR: Bad genotype symbols")
    print(s1, s2, s3)
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

                    # Correct genotype symbols
                    if s1:
                        data = [correct_symbol(x, s1, s2, s3) for x in data]

                to_write += data

                if "," in "".join(to_write):
                    print("ERROR: Comma in data line, problem for csv file")
                    sys.exit(1)

                outfile.write(",".join(to_write) + "\n")
