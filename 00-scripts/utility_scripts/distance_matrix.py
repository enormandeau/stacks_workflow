#!/usr/bin/env python
"""Compute a pairwise similarity matrix among samples from a filtered VCF file

Usage:
    ./00-scripts/utility_scripts/similarity_matrix.py input_vcf output_similarity output_missing

Will create an output matrix to use with the R script:
    ./00-scripts/utility_scripts/plot_distance_matrix.R
"""

# Modules
print("Loading modules...")
from collections import defaultdict
from multiprocessing import Process, Value, Array
import numpy as np
import sys

# Functions
def compute_similarity(samples, genotypes):
    num_samples = len(samples)
    similarity = Array('d', [1.0] * (num_samples ** 2))
    missingness = Array('d', [0.4] * (num_samples ** 2))
    sample_count = 0

    def pair_similarity(i, j, similarity, missingness, genotypes, num_samples):
        s1 = genotypes[:,i]
        s2 = genotypes[:,j]
        dist = 0.0
        dist_count = 0
        missing = 0.0
        missing_count = 0
        for pos in xrange(len(s1)):
            g1 = s1[pos].split(":")[0]
            g2 = s2[pos].split(":")[0]
            if g1 != "./." and g2 != "./.":
                dist_count += 1
                missing_count += 1
                if g1 == "1/0":
                    g1 = "0/1"
                if g2 == "1/0":
                    g2 = "0/1"
                if g1 == g2:
                    dist += 1
            elif g1 == "./." or g2 == "./.":
                missing_count += 1
                if g1 != g2:
                    missing += 1

        pos = i * num_samples + j
        rev_pos = j * num_samples + i
        similarity[pos] = float(dist) / float(dist_count)
        similarity[rev_pos] = similarity[pos]
        missingness[pos] = float(missing) / float(missing_count)
        missingness[rev_pos] = missingness[pos]

    for i in xrange(num_samples):
        sample_count += 1
        print("  Treating sample {}/{}: {}".format(sample_count, num_samples, samples[i]))

        for j in xrange(i + 1, num_samples):
            p = Process(target=pair_similarity, args=(i, j, similarity, missingness, genotypes, num_samples))
            p.start()

    similarity = np.reshape(similarity, [num_samples, num_samples])
    missingness = np.reshape(missingness, [num_samples, num_samples])
    return [similarity, missingness]

def compute_missing(samples, genotypes):
    num_samples = len(samples)
    missing = np.ones((num_samples, num_samples))
    return missing

# Main
if __name__ == '__main__':
    try:
        input_vcf = sys.argv[1]
        output_similarity = sys.argv[2]
        output_missing = sys.argv[3]
    except:
        print __doc__
        sys.exit(1)

    # Parse VCF file
    print("Parsing VCF...")
    genotypes = []
    with open(input_vcf) as vfile:
        for line in vfile:
            if not line.startswith("##"):
                if line.startswith("#"):
                    samples = line.strip().split()[9:]
                else:
                    genotypes.append(line.strip().split()[9:])

    num_samples = len(samples)
    genotypes = np.array(genotypes)

    # Calculate similarity and missingness genotypes
    # Distance is the proportion of genotypes that differ
    print("Computing similarity and missingness...")
    similarity, missingness = compute_similarity(samples, genotypes)

    # Writing similarity matrix to file
    with open(output_similarity + ".csv", "w") as outf:
        outf.write("\t".join(["Sample"] + samples) + "\n")
        for i in xrange(num_samples):
            data = [str(x) for x in similarity[i,]]
            outf.write("\t".join([samples[i]] + data) + "\n")

    # Writing missingness matrix to file
    with open(output_missing + ".csv", "w") as outf:
        outf.write("\t".join(["Sample"] + samples) + "\n")
        for i in xrange(num_samples):
            data = [str(x) for x in missingness[i,]]
            outf.write("\t".join([samples[i]] + data) + "\n")

    # Plotting
    # Similarity
    from matplotlib import pyplot as plt
    heatmap = plt.pcolor(similarity, cmap=plt.cm.get_cmap('Blues'), vmax=1)
    plt.savefig(output_similarity + ".png", dpi=300, bbox_inches='tight')
    ##plt.show()

    # Similarity
    from matplotlib import pyplot as plt
    heatmap = plt.pcolor(missingness, cmap=plt.cm.get_cmap('Blues'))
    plt.savefig(output_missing + ".png", dpi=300, bbox_inches='tight')
    #plt.show()
