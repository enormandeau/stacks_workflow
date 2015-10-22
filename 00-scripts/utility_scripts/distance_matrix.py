#!/usr/bin/env python
"""Compute a pairwise similarity matrix among samples from a filtered VCF file

Usage:
    ./00-scripts/utility_scripts/similarity_matrix.py input_vcf output_matrix

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
    sample_count = 0

    def pair_similarity(i, j, similarity, genotypes, num_samples):
        s1 = genotypes[:,i]
        s2 = genotypes[:,j]
        dist = 0
        count = 0
        for pos in xrange(len(s1)):
            g1 = s1[pos].split(":")[0]
            g2 = s2[pos].split(":")[0]
            if g1 != "./." and g2 != "./.":
                count += 1
                if g1 == "1/0":
                    g1 = "0/1"
                if g2 == "1/0":
                    g2 = "0/1"
                if g1 == g2:
                    dist += 1

        pos = i * num_samples + j
        rev_pos = j * num_samples + i
        similarity[pos] = float(dist) / float(count)
        similarity[rev_pos] = similarity[pos]

    for i in xrange(num_samples):
        sample_count += 1
        print("  Treating sample {}/{}: {}".format(sample_count, num_samples, samples[i]))

        for j in xrange(i + 1, num_samples):
            p = Process(target=pair_similarity, args=(i, j, similarity, genotypes, num_samples))
            p.start()

    return np.reshape(similarity, [num_samples, num_samples])

def compute_missing(samples, genotypes):
    num_samples = len(samples)
    missing = np.ones((num_samples, num_samples))
    return missing

# Main
if __name__ == '__main__':
    try:
        input_vcf = sys.argv[1]
        output_matrix = sys.argv[2]
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

    # Distance is the proportion of genotypes that differ
    print("Computing similarity...")
    similarity = compute_similarity(samples, genotypes)
    print similarity.shape

    # Missing is the proportion of genotypes found in only one sample
    print("Computing missing...")
    missing = compute_missing(samples, genotypes)

    # Writing matrix to file
    with open(output_matrix, "w") as outf:
        outf.write("\t".join(["Sample"] + samples) + "\n")
        for i in xrange(num_samples):
            print "Sample:", i
            #outf.write("\t".join([samples[i]] + [str(x) for x in similarity[i:]]) + "\n")

    # Plotting
    from matplotlib import pyplot as plt
    heatmap = plt.pcolor(similarity, cmap=plt.cm.get_cmap('Blues'), vmax=1)
    plt.savefig(output_matrix + ".png", dpi=300, bbox_inches='tight')
    #plt.show()

