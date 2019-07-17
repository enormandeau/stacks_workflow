#!/usr/bin/env python
"""Reduce a STACKS catalog with a whitelist of loci

Usage:
    ./00-scripts/utility_scripts/catalog_filter.py catalog_stub input_whitelist output_suffix
"""

# Modules
import gzip
import sys

# Parsing user input
try:
    catalog_stub = sys.argv[1]
    input_whitelist = sys.argv[2]
    output_suffix = sys.argv[3]
except:
    print __doc__
    sys.exit(1)

# Read whitelist
whitelist = set([x.strip() for x in open(input_whitelist).readlines() if x])

# Filter files
files = ["alleles.tsv.gz", "snps.tsv.gz", "tags.tsv.gz"]

for _file in files:
    print "Treating:", _file
    with gzip.open(catalog_stub + _file) as infile:
        with gzip.open(catalog_stub + output_suffix + "." + _file, "w") as outfile:
            for line in infile:
                l = line.strip().split()
                if line.startswith("#"):
                    outfile.write(line)
                elif l[2] in whitelist:
                    outfile.write(line)
