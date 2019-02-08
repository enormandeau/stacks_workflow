#!/usr/bin/env python2
"""Use a filtered vcf to filter an haplotypes.tsv file

Usage:
    <program> filtered_stacks_vcf input_haplotypes max_multi_haplotypes min_proportion output_haplotypes

Where:
    filtered_stacks_vcf: STACKS VCF filtered with 05_filter_vcf.py
    input_haplotypes: STACKS batch_1.haplotypes.tsv
    max_multi_haplotypes: Maximum number of samples with 3+ haplotypes
    min_proportion: Minimum proportion of samples with data (1 - prop_missing)
    output_haplotypes: output haplotype.tsv file
"""

# Modules
from collections import defaultdict
import sys

# Parse user input
try:
    filtered_stacks_vcf = sys.argv[1]
    input_haplotypes = sys.argv[2]
    max_multi_haplotypes = int(sys.argv[3])
    min_proportion = float(sys.argv[4])
    output_haplotypes = sys.argv[5]
except:
    print __doc__
    sys.exit(1)

# Read filtered VCF to collect wanted loci IDs
wanted_loci = set()
with open(filtered_stacks_vcf) as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        else:
            l = line.strip().split("\t")
            locus = l[2]
            if "_" in locus:
                locus = locus.split("_")[0]
                wanted_loci.add(locus)

# Treat input haplotypes
distribution_tri_allelic = defaultdict(int)

with open(output_haplotypes, "w") as outfile:
    with open(input_haplotypes) as infile:
        for line in infile:
            l = line.strip().split("\t")

            # Write header
            if line.startswith("Catalog"):
                outfile.write(line)

            # Treat loci
            else:
                infos = l[:2]
                locus = infos[0]
                
                if locus not in wanted_loci and "consensus" not in line:
                    continue

                haplotypes = l[2:]

                num_missing = len([x for x in haplotypes if x == "-"])
                num_tri_allelic = len([x for x in haplotypes if x.count("/") > 1])
                distribution_tri_allelic[num_tri_allelic] += 1
                multi_haplotypes = [x for x in haplotypes if x.count("/") > 1]

                # If too many, go to next
                if len(multi_haplotypes) > max_multi_haplotypes:
                    continue

                # Correct haplotypes containing Ns
                if "N" in "".join(haplotypes):
                    bad_haplotypes = [x for x in haplotypes if "N" in x]
                    bad_positions = set()
                    for bad in bad_haplotypes:
                        for b in bad.split("/"):
                            for i in range(len(b)):
                                if b[i] == "N":
                                    bad_positions.add(i)

                    bad_positions = sorted(list(bad_positions), reverse=True)

                    corrected_haplotypes = []
                    for h in haplotypes:

                        corrected_alleles = []
                        for allele in h.split("/"):
                            if allele == "-":
                                corrected_alleles.append(allele)
                                continue

                            for p in bad_positions:
                                allele = allele[:p] + allele[p+1:]

                            corrected_alleles.append(allele)

                        corrected_haplotypes.append("/".join(corrected_alleles))

                    haplotypes = corrected_haplotypes

                # Else, correct the genotypes to null ("-")
                else:
                    for i in xrange(len(haplotypes)):
                        if haplotypes[i].count("/") > 1:
                            haplotypes[i] = "-"

                    prop_missing = float(num_missing) / len(haplotypes)
                    prop_with_data = 1 - prop_missing
                    if prop_with_data >= min_proportion:
                        outfile.write("\t".join(infos + haplotypes) + "\n")
# Report distribution of tri-allelic samples
print("Distribution of tri-allelic samples")
print("---")
print("NumTri\tFrequency")
for n in sorted(distribution_tri_allelic):
    print str(n) + "\t" + str(distribution_tri_allelic[n])
print("---")
print("Note to future self:")
print("  Add R script to produce figure if needed")
