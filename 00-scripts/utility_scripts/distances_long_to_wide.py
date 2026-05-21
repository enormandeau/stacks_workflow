#!/usr/bin/env python3
"""Convert pairwise distance file from long to wide format

Usage:
    <program> input_distances output_distances
"""

# Modules
import sys

# Parsing user input
try:
    input_distances = sys.argv[1]
    output_distances = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Load distances
distances = [x.strip().split("\t") for x in open(input_distances, "rt").readlines()]
samples = sorted(set([x[0] for x in distances]))

distances_dict = {}

for d in distances:
    distances_dict[(d[0], d[1])] = d[2]

# Write wide format to file
with open(output_distances, "wt") as outfile:
    header = ["Name"] + samples
    outfile.write("\t".join(header) + "\n")

    for s1 in samples:
        line = [s1]

        for s2 in samples:
            line.append(distances_dict[(s1, s2)])

        outfile.write("\t".join(line) + "\n")
