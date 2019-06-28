#!/usr/bin/env python2
"""Use a filtered vcf to filter an haplotypes.vcf file

Usage:
    <program> filtered_stacks_vcf input_haplotypes max_multi_haplotypes min_proportion output_haplotypes

Where:
    filtered_stacks_vcf: STACKS VCF filtered with 05_filter_vcf.py
    input_haplotypes: STACKS batch_1.haplotypes.vcf
    max_multi_haplotypes: Maximum number of samples with 3+ haplotypes
    min_proportion: Minimum proportion of samples with data (1 - prop_missing)
    output_haplotypes: output haplotype.vcf file
"""

# Modules
from collections import defaultdict
import sys

# Parse user input
try:
    filtered_vcf = sys.argv[1]
    input_haplotypes = sys.argv[2]
    max_multi_haplotypes = int(sys.argv[3])
    min_proportion = float(sys.argv[4])
    output_haplotypes = sys.argv[5]
except:
    print __doc__
    sys.exit(1)

# Read filtered VCF to collect wanted loci IDs
wanted_loci = set()
with open(filtered_vcf) as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        else:
            l = line.strip().split("\t")
            locus = l[2]
            if "_" in locus:
                locus = locus.split("_")[0]
                wanted_loci.add(locus)

# Treat input haplotypes
distribution_tri_allelic = defaultdict(int)

with open(output_haplotypes, "w") as outfile:
    with open(input_haplotypes) as infile:
        for line in infile:
            l = line.strip().split("\t")

            # Write header
            if line.startswith("#"):
                outfile.write(line)

            # Treat loci
            else:
                infos = l[:4]
                locus = infos[2]
                
                if locus not in wanted_loci and "consensus" not in line:
                    continue

                else:
                    outfile.write(line.replace("\t\t", "\t./.:0\t").replace("\t\t", "\t./.:0\t"))
