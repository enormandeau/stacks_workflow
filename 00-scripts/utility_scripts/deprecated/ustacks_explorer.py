#!/usr/bin/env python
"""Get stats from ustacks results

Usage:
    ./01_scripts/01_ustacks_explorer.py input_folder output_file
"""

# Modules
from collections import defaultdict
import gzip
import sys
import os

# Functions
def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

# Parse user input
try:
    input_folder = sys.argv[1]
    output_file = sys.argv[2]
except:
    print __doc__
    sys.exit(1)

# Identify samples to treat
folders = sorted(list(set( os.listdir(input_folder))))

with open(output_file, "w") as outfile:
    outfile.write("Run\tm\tM\tN\tSample\tNumLoci\tMedCov\tMeanCov\tNumSNP\tNumHetSNP\n")

    for folder in folders:
        samples = sorted(list(set([x.split(".")[0] for x in os.listdir(os.path.join(input_folder, folder))])))

        # For each sample, get stats
        for sample in samples:
            stats = [folder]
            params = list(folder.split("_")[1])
            stats += params
            stats.append(sample)
            loci = set()
            snps = set()
            heterozygote = 0
            loci_dict = defaultdict(float)

            ## number of loci and coverage per locus
            with gzip.open(os.path.join(input_folder, folder, sample + ".alleles.tsv.gz")) as infile:

                for line in infile:
                    if line.startswith("#"):
                        continue

                    l = line.strip().split()
                    locus = l[2]
                    coverage = float(l[5])
                    loci.add(locus)
                    loci_dict[locus] += coverage

            stats.append(len(loci))
            stats.append(median(loci_dict.values()))
            stats.append(float(sum(loci_dict.values())) / len(loci_dict.values()))

            ## number of SNPs and number of heterozygote SNPs
            with gzip.open(os.path.join(input_folder, folder, sample + ".snps.tsv.gz")) as infile:

                for line in infile:
                    if line.startswith("#"):
                        continue

                    l = line.strip().split()
                    snp = l[2]
                    het = l[4]
                    snps.add(snp)
                    if het == "E":
                        heterozygote += 1

            stats.append(len(snps))
            stats.append(heterozygote)
            outfile.write("\t".join([str(x) for x in stats]) + "\n")
            print "\t".join([str(x) for x in stats])
