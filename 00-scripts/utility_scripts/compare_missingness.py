#!/usr/bin/env python
"""Find populations with similar information missingness across loci

Usage:
    ./00-scripts/utility_scripts/compare_missingness.py sumstats_file output_file min_num_pop
"""

# Modules
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
import sys

# Classes
class Population(object):
    def __init__(self, pop_num, pop_name, num_samples):
        self.pop_num = pop_num
        self.pop_name = pop_name
        self.num_samples = num_samples

    def __repr__(self):
        return "\t".join([self.pop_num, self.pop_name, str(self.num_samples)])

# Functions
def get_pop_counts(handle):
    pop_info = dict()
    with open(handle) as infile:
        for line in infile:
            if "Batch ID" in line or not line.startswith("#"):
                break
            else:
                l = line.split()
                pop_num = l[1]
                samples = l[2].split(",")
                pop_name = samples[0].split("_")[0]
                num_samples = len(samples)
                p = Population(pop_num, pop_name, num_samples)
                pop_info[pop_num] = p

    return pop_info

def get_missingness(handle, pop_info):
    missingness = defaultdict(dict)
    with open(handle) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            else:
                l = line.split("\t")
                snp_id = l[3]
                pop_num = str(int(l[5]) + 1)
                num_genotyped = int(l[8])
                missingness[snp_id][pop_num] = float(num_genotyped) / pop_info[pop_num].num_samples

    return missingness

def get_num_pop(missingness):
    pops = set()
    for snp in missingness:
        for pop in missingness[snp]:
            pops.add(pop)

    print pops
    return pops

def compute_similarity(missingness, pop_info):
    n = len(pop_info)
    similarity = np.diag([0.1] * n)
    pop_numbers = sorted([int(x) for x in pop_info])
    for pop1 in range(len(pop_numbers)):
        for pop2 in range(pop1 + 1, len(pop_numbers)):
            n_loci = 0
            sim = 0.0
            p1 = str(pop_numbers[pop1])
            p2 = str(pop_numbers[pop2])
            for snp in missingness:
                try:
                    m1 = missingness[snp][p1]
                except:
                    continue
                try:
                    m2 = missingness[snp][p2]
                except:
                    continue

                n_loci += 1
                sim += abs(m1 - m2)

            total_sim = sim / n_loci
            similarity[pop1, pop2] = total_sim
            similarity[pop2, pop1] = total_sim

    return similarity

# Main
if __name__ == '__main__':
    # Parse user input
    try:
        sumstats_file = sys.argv[1]
        output_file = sys.argv[2]
        min_num_pop = int(sys.argv[3])
    except:
        print __doc__
        sys.exit(1)

    # Get population informations
    pop_info = get_pop_counts(sumstats_file)

    # Extract missingness by population
    print("Extracting missingness info by population")
    missingness = get_missingness(sumstats_file, pop_info)

    # Get number of populations
    populations = get_num_pop(missingness)
    print("There were {} populations".format(len(populations)))

    ## Remove absent populations from pop_info
    #to_remove = []
    #for pop in pop_info:
    #    if pop not in populations:
    #        to_remove.append(pop)

    #for pop in to_remove:
    #    pop_info.pop(pop)

    # Remove SNPs with less than <min_num_pop> populations
    to_remove = []
    for snp in missingness:
        if len(missingness[snp]) <  min_num_pop:
            to_remove.append(snp)

    print("  {}% of the SNPs had less than {} populations present".format(
        round(100.0 * float(len(to_remove)) / float(len(missingness)), 2),
        min_num_pop))

    for snp in to_remove:
        missingness.pop(snp)

    print("  Keeping {} SNPs".format(len(missingness)))

    # Compute pairwise similarity among populations
    print("Computing similarity...")
    similarity = compute_similarity(missingness, pop_info)

    # Print population names
    pop_numbers = sorted([int(x) for x in pop_info])
    population_names = []
    for p in pop_numbers:
        name = pop_info[str(p)].pop_name
        print p, name
        population_names.append(name)

    # Write similarity to file
    with open(output_file, "w") as outfile:
        outfile.write("\t".join(population_names) + "\n")
        for i in range(similarity.shape[0]):
            outfile.write("\t".join([str(x) for x in list(similarity[i,])]) + "\n")

    ## Plotting
    #heatmap = plt.pcolor(similarity, cmap=plt.cm.get_cmap('Blues'), vmax=1)
    #plt.savefig(output_file + ".png", dpi=300, bbox_inches='tight')
    ##plt.show()

