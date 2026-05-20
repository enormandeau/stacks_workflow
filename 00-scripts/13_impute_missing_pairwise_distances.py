#!/usr/bin/env python3
"""Impute missing data using sample pairwise distances

Usage:
    <program> input_vcf input_distances number_closest output_vcf
"""

# Modules
from collections import Counter
import numpy
import gzip
import sys

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def impute(closest_genotypes):
    """Impute using most frequent genotype from closely related samples
    """
    genotypes = [x.split(":")[0] for x in closest_genotypes]
    genotypes = [x for x in genotypes if x != "./."]

    if not genotypes:
        genotype = "0/0"
    else:
        # Impute most frequent genotype
        counts = Counter(genotypes)
        genotype = counts.most_common(1)[0][0]

    # Return composed genotype and associated info
    return genotype + ":0:0,0:0:0,0,0"

# Parsing user input
try:
    input_vcf = sys.argv[1]
    input_distances = sys.argv[2]
    number_closest = int(sys.argv[3])
    output_vcf = sys.argv[4]
except:
    print(__doc__)
    sys.exit(1)

# Read admixture input
distances = [x.split("\t") for x in open(input_distances).read().strip().split("\n")]
distances = [[x[0], x[1], float(x[2])] for x in distances]

# Pre-hash liste of closest 11 samples for each sample
samples = set([x[0] for x in distances])

closest_samples = dict()

for s in samples:
    data = [x for x in distances if x[0] == s]

    # Decorate and sort
    data = sorted([[x[2]] + x for x in data])

    # Extract wanted samples
    data = [x[2] for x in data[1: 1 + number_closest]]
    closest_samples[s] = data

# Process VCF
with myopen(input_vcf) as infile:
    with myopen(output_vcf, "wt") as outfile:

        num_genotypes = 0
        num_imputed = 0

        for line in infile:
            l = line.strip().split()

            if line.startswith("#CHROM"):
                outfile.write(line)
                samples = list(enumerate(l[9:]))
                samples_from_id = dict(samples)
                ids_from_sample = dict([x[::-1] for x in samples])
                print(f"Treating {len(samples)} samples")
                continue

            elif line.startswith("#"):
                outfile.write(line)
                continue

            info = l[:9]
            data = l[9:]
            new_data = []

            for i, genotype in enumerate(data):
                sample = samples_from_id[i]
                closest = closest_samples[sample]
                closest_ids = set([ids_from_sample[x] for x in closest])
                closest_genotypes = [x[1] for x in enumerate(data) if x[0] in closest_ids]
                num_genotypes += 1

                if genotype.startswith("./."):
                    new_data.append(impute(closest_genotypes))
                    num_imputed += 1

                else:
                    new_data.append(genotype)

            outfile.write("\t".join(info + new_data) + "\n")

percent = 100 * num_imputed / num_genotypes
print(f"Imputed {percent:.2f}% of the genotypes ({num_imputed}/{num_genotypes})")
