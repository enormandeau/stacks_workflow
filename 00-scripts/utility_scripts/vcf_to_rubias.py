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
samples_info = []
loci_names = ["Sample"]

with open(input_vcf) as infile:
    for line in infile:
        l = line.strip().split("\t")
        locus = "_".join(l[:3])

        if line.startswith("##"):
            continue

        elif line.startswith("#CHROM"):
            for i in range(9, len(l)):
                samples_info.append([l[i]])

        else:
            loci_names.append(locus)
            loci_names.append(locus)
            for i in range(9, len(l)):
                genotype_info = l[i].split(":")[0]
                genotype_info = genotype_info.split("/")
                genotype_info = [x if x != "." else "NA" for x in genotype_info]
                samples_info[i-9] += genotype_info

# Write rubias file one line (or sample) at a time
with open(output_rubias, "w") as outfile:
    outfile.write("\t".join(loci_names) + "\n")

    for sample in samples_info:
        outfile.write("\t".join(sample) + "\n")
