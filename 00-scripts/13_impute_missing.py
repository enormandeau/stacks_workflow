#!/usr/bin/env python3
"""Impute missing data using admixture membership info

Usage:
    <program> input_vcf input_admixture output_vcf
"""

# Modules
import numpy
import gzip
import sys

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def compute_group_weights(data, memberships):
    """For each admixture group, compute the sum of the weights for each allele
    """
    genotypes = [x.split(":")[0].split("/") for x in data]
    weights = []

    # Sum allele weights per group
    for group in range(len(memberships[0])):
        group_weights = [0.0, 0.0]

        for sample, g in enumerate(genotypes):
            for allele in g:

                if allele in "01":
                    group_weights[int(allele)] += memberships[sample][group]

        # Normalize to a sum of 1
        group_weights = [x / sum(group_weights) for x in group_weights]
        weights.append(group_weights)

    return weights

def impute(sample_membership, group_weights):
    """Use computed allele weights per admixture group and membership to each
    group to find final weight of alleles and randomly draw alleles to form
    the genotypes
    """
    weights = [0.0, 0.0]

    for group, weight in enumerate(group_weights):
        weights[0] += weight[0] * sample_membership[group]
        weights[1] += weight[1] * sample_membership[group]

    weights = [x / sum(weights) for x in weights]

    # Weighted random sampling of two independent allels
    genotype = "/".join(sorted(numpy.random.choice(["0", "1"], 2, p=weights)))

    # Return composed genotype and associated info
    return genotype + ":0:0,0:0:0,0,0"

# Parsing user input
try:
    input_vcf = sys.argv[1]
    input_admixture = sys.argv[2]
    output_vcf = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Read admixture input
memberships = open(input_admixture).read().strip().split("\n")
memberships = [[float(y) for y in x.split(" ")] for x in memberships]

# Process VCF
with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:

        num_genotypes = 0
        num_imputed = 0

        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            l = line.strip().split()
            info = l[:9]
            data = l[9:]
            new_data = []

            group_weights = compute_group_weights(data, memberships)

            for sample, genotype in enumerate(data):
                num_genotypes += 1

                if genotype.startswith("./."):
                    new_data.append(impute(memberships[sample], group_weights))
                    num_imputed += 1

                else:
                    new_data.append(genotype)

            outfile.write("\t".join(info + new_data) + "\n")

percent = 100 * num_imputed / num_genotypes
print(f"Imputed {percent:.2f}% of the genotypes ({num_imputed}/{num_genotypes})")
