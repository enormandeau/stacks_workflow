#!/usr/bin/env python3
"""Replace missing data in VCF by most common genotype globaly

Usage:
    <program> input_vcf output_vcf
"""

# Modules
from collections import Counter
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
with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
        for line in infile:

            if line.startswith("#"):
                outfile.write(line)

            else:
                l = line.strip().split("\t")
                info = l[:9]
                data = l[9:]
                genotypes = [x.split(":")[0] for x in data if x.split(":")[0] != "./."]
                counts = Counter(genotypes)
                most_common = counts.most_common()[0][0]

                new_data = []

                for d in data:

                    if d.startswith("./."):
                        new_d = ":".join([most_common] + d.split(":")[1:])
                        new_data.append(new_d)

                    else:
                        new_data.append(d)

                new_line = "\t".join(info + new_data) + "\n"
                outfile.write(new_line)
