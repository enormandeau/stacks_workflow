#!/usr/bin/env python3
"""Prepare data to create figure of missing data proportion as a funciton of reads per sample

Usage:
    <program> input_numreads input_propmissing
"""

# Modules
#import matplotlib.pyplot as plt
#import numpy as np
import sys

# Parse user input
try:
    input_numreads = sys.argv[1]
    input_propmissing = sys.argv[2]
    output_png = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Preamble
sample_dict = {}

# Read numreads data
numreads = [x.strip().split() for x in open(input_numreads).readlines()]
for n in numreads:
    name, count = n
    name = name.split(".")[0]
    count = int(count)
    sample_dict[name] = [count]

# Read propmissing data
propmissing = [x.strip().split() for x in open(input_propmissing).readlines()]
for p in propmissing:
    if p[0] == "Population":
        continue

    name = p[0] + "_" + p[1]
    proportion = float(p[2])

    if name in sample_dict:
        sample_dict[name].append(proportion)

# Create data for figure
data = []

for s in sorted(sample_dict):
    if len(sample_dict[s]) != 2:
        continue

    count, proportion = sample_dict[s]
    print(proportion, count)

## Fixing random state for reproducibility
#x = np.random.rand(N)
#y = np.random.rand(N)
#plt.scatter(x, y, s=area, c=colors, alpha=0.5)
#plt.show()
