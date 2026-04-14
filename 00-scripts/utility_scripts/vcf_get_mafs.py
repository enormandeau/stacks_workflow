#!/usr/bin/env python3
"""Report MAF values for each population found in VCF

Usage:
    <program> input_vcf output_mafs
"""

# Modules
import gzip
import sys

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def get_maf(genotypes):
    """Return Minor Allele Frequency in genotypes
    """
    g = [x for x in genotypes if x != "./."]
    g = "".join(g).replace("/", "")
    return round(g.count("1") / len(g) ,4)

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_mafs = sys.argv[2]
except:
    print(__doc__)
    sys.exit()

# Read info from VCF
pop_positions = dict()

with open(output_mafs, "wt") as outfile:
    with myopen(input_vcf, "rt") as infile:
        for line in infile:
            l = line.strip().split()
            snp_info = l[:2]
            data = l[9:]

            if line.startswith("##"):
                continue

            elif line.startswith("#CHROM"):
                for i, sample in enumerate(data):
                    pop = sample.split("_")[0]
                    if pop not in pop_positions:
                        pop_positions[pop] = [i]
                    else:
                        pop_positions[pop].append(i)

                outfile.write("Chrom\tPos\t" + "\t".join(pop_positions.keys()) + "\n")
                continue

            # Extract MAFs
            genotypes = [x.split(":")[0] for x in data]
            pop_mafs = dict()

            for pop in pop_positions:
                snp_info.append(get_maf([genotypes[i] for i in pop_positions[pop]]))

            outfile.write("\t".join([str(x) for x in snp_info]) + "\n")
