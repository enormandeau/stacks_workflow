#!/usr/bin/env python3
"""Filtering VCF files output by STACKS1 or STACKS2 to remove SNPs

Usage:
    <program> input_vcf filtering_type info_file output_vcf

Where:
    input_vcf is the name of the VCF file to filter
    filtering_type is either "wanted" or "unwanted"
    info_file contains wanted or unwanted SNP ids (3 first columns of the VCF)
    output_vcf is the name of the filtered VCF
"""

# Import
import sys

# Parse user input
try:
    input_vcf = sys.argv[1]
    filtering_type = sys.argv[2]
    info_file = sys.argv[3]
    output_vcf = sys.argv[4]
except:
    print(__doc__)
    sys.exit(1)

if not filtering_type in ["wanted", "unwanted"]:
    print("ERROR: filtering_type must be 'wanted' or 'unwanted'")
    sys.exit(2)

# Get wanted or unwanted SNPs
listed_snps = set([tuple(x.strip().split("\t")) for x in open(info_file).readlines()])

# Filter
with open(input_vcf) as infile:
    with open(output_vcf, "w") as outfile:
        for line in infile:
            l = line.strip().split("\t")

            if line.startswith("#"):
                outfile.write(line)
                continue

            if filtering_type == "wanted" and tuple(l[:3]) in listed_snps:
                outfile.write(line)

            elif filtering_type == "unwanted" and tuple(l[:3]) not in listed_snps:
                outfile.write(line)
