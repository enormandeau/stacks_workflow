#!/usr/bin/env python3
"""Replace missing data in VCF by most common genotype in each group

Groups are defined as samples sharing the same first part of their name before
an underscore. Here are some example groups:

# group pop1
    pop1_01
    pop1_02
    ...
    pop1_30

# group SAU
    SAU_0212
    SAU_0328
    ...
    SAU_1297

Usage:
    <program> input_vcf output_vcf
"""

# Modules
from collections import Counter, defaultdict
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
    output_vcf = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read VCR and impute missing
group_samples = defaultdict(list)
sample_groups = {}

with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
        for line in infile:
            l = line.strip().split("\t")

            if line.startswith("#"):
                outfile.write(line)

                # Create indexes of samples per group and groups for each sample
                if line.startswith("#CHROM"):
                    samples = l[9:]
                    num_group = enumerate([x.split("_")[0] for x in samples])

                    for num, group in num_group:
                        sample_groups[num] = group
                        group_samples[group].append(num)

            else:
                info = l[:9]
                data = l[9:]
                most_common_per_group = {}

                # Get most common genotype per group
                for group in group_samples:
                    sample_data = [data[i] for i in group_samples[group]]

                    genotypes = [x.split(":")[0] for x in sample_data if x.split(":")[0] != "./."]
                    counts = Counter(genotypes)
                    most_common = counts.most_common()[0][0]
                    most_common_per_group[group] = most_common

                new_data = []

                for i, d in enumerate(data):

                    # Replace missing by most common genotype in group
                    if d.startswith("./."):
                        most_common_in_group = most_common_per_group[sample_groups[i]]
                        new_d = ":".join([most_common_in_group] + d.split(":")[1:])
                        new_data.append(new_d)

                    else:
                        new_data.append(d)

                new_line = "\t".join(info + new_data) + "\n"
                outfile.write(new_line)
