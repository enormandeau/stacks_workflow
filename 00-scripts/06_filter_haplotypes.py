#!/usr/bin/env python
"""Use a filtered vcf to filter an haplotypes.tsv file

Usage:
    ./00-scripts/06_filter_haplotypes.py filtered_vcf input_haplotypes max_multi_haplotypes min_proportion output_haplotypes
"""

# Modules
import sys

# Parse user input
try:
    filtered_vcf = sys.argv[1]
    input_haplotypes = sys.argv[2]
    max_multi_haplotypes = int(sys.argv[3])
    min_proportion = float(sys.argv[4])
    output_haplotypes = sys.argv[5]
except:
    print __doc__
    sys.exit(1)

# Read filtered VCF to collect wanted loci IDs
wanted_loci = set()
with open(filtered_vcf) as infile:
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
                
                if locus not in wanted_loci:
                    continue

                haplotypes = l[2:]
                multi_haplotypes = [x for x in haplotypes if x.count("/") > 1]

                # If too many, go to next
                if len(multi_haplotypes) > max_multi_haplotypes:
                    continue

                # Else, correct the genotypes to null ("-")
                else:
                    for i in xrange(len(haplotypes)):
                        if haplotypes[i].count("/") > 1:
                            haplotypes[i] = "-"

                    num_missing = len([x for x in haplotypes if x == "-"])
                    prop_missing = float(num_missing) / len(haplotypes)
                    prop_with_data = 1 - prop_missing
                    if prop_with_data >= min_proportion:
                        outfile.write("\t".join(infos + haplotypes) + "\n")
