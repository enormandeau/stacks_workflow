#!/usr/bin/env python
"""Take a VCF file and create a bitmap representation of the genotypes

Usage:
    <program> input_vcf output_bitmap
"""

# Modules
from PIL import Image
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
    output_bitmap = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read VCF into list of lists
genotypes = []
counter = 0

print("Reading genotypes")
with myopen(input_vcf,'rt') as infile:
    for line in infile:

        if line.startswith("#"):
            continue

        counter += 1
        l = line.strip().split('\t')[9: ]
        l = [x.split(":")[0] for x in l]
        l = [x.count("1") for x in l]

        # Skip lines with too few rare alleles
        if sum([x > 0 for x in l]) < 20:
            continue

        # Skip lines with too many frequent allels
        if sum([x == 0 for x in l]) < len(l) / 2:
            continue

        genotypes.append(l)

        if not counter % 1000:
            print(".", end="", flush=True)
print()

# Creating image
pixels = []
ncol = len(genotypes[0])
nlines = len(genotypes)
im = Image.new('RGB', (ncol, nlines))
print(f"Generating pixels ({ncol} x {nlines})")

for line in genotypes:
    for g in line:
        if g == 0:
            rgb = (10, 10, 10)
        elif g == 1:
            rgb = (150, 150, 255)
        else:
            rgb = (255, 100, 100)

        pixels.append(rgb)

# Create image
print("Writing image")
im.putdata(pixels)
    
# Write image
im.save(output_bitmap)
