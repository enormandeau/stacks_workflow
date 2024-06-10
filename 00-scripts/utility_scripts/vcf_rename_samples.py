#!/usr/bin/env python3
"""Rename samples in a VCF using a correspondance file

Usage:
    <program> input_vcf correspondance_file output_vcf

Details:
    - Renaming from old name (in first column of tab separated correspondance
      file) to new name. If a sample is not found in the correspondance file,
      its name remains unchanged.
    - The input and output VCF files can be compressed with gzip. In this case,
      their name must end in `.gz`.
"""

# Modules
import gzip
import sys

# Defining functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

# Parsing user input
try:
    input_vcf = sys.argv[1]
    correspondance_file = sys.argv[2]
    output_vcf = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Read correspondance file
corr = dict()
with open(correspondance_file) as infile:
    for line in infile:
        l = line.strip().split("\t")
        corr[l[0]] = l[1]

# Modify VCF
with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:
        for line in infile:
            if line.startswith("#CHROM"):
                l = line.strip().split("\t")
                data = l[:9]
                names = l[9:]

                for i, n in enumerate(names):
                    if n in corr:
                        names[i] = corr[n]
                        
                new_line = "\t".join(data + names) + "\n"
                outfile.write(new_line)

            else:
                outfile.write(line)
