#!/usr/bin/env python3
"""Create the SampleName and genotype info columns for Rubias

Note: You will still need to add the first 3 metadata columns

Usage:
    <program> input_vcf output_rubias
"""

# Modules
import sys

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_rubias = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read VCF and collect information
print("hello")

# Write rubias file one line (or sample) at a time
