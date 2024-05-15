#!/usr/bin/env python3
"""Fix positions of plink map from split genome

Usage:
    <program> input_plink_map chunk_size output_plink_map
"""

# Modules
import sys

# Parsing user input
try:
    input_plink_map = sys.argv[1]
    chunk_size = int(sys.argv[2])
    output_plink_map = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Fix positions
with open(input_plink_map) as infile:
    with open(output_plink_map, "wt") as outfile:
        for line in infile:

            if line.startswith("#"):
                outfile.write(line)
                continue

            l = line.strip().split("\t")
            chunk = int(l[0].split("chunk")[1])
            l[3] = str(int(l[3]) + (chunk - 1) * chunk_size)
            outfile.write("\t".join(l) + "\n")
