#!/usr/bin/env python3
"""Get allele frequencies for each population and each SNP

Usage:
    <program> input_vcf output_file exponent num_snps

Where:
    input_vcf: input vcf (can be gzipped compressed)
    exponent: x**e, where e is above zero (preferably above 0.5 and below or at 2.0)
    num_snps: number of SNPs for the least differentiated group
    output_vcf: output vcf (can be gzipped compressed)

"""

# Modules
from collections import defaultdict
from math import sqrt
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
    exponent = float(sys.argv[2])
    num_snps = int(sys.argv[3])
    output_file = sys.argv[4]
except:
    print(__doc__)
    sys.exit(1)

assert 0 < exponent, "exponent: x**e, where e is above zero (preferably above 0.5 and below or at 2.0)"
assert 1 < num_snps, "num_snps must be an integer greater than 1"

# Read VCF
pop_ids = defaultdict(list)
pairwise_afds = defaultdict(list)
best_afds = defaultdict(list)
best_snps = set()

with myopen(input_vcf, "rt") as infile:
#with myopen(output_file, "wt") as outfile:
    for line in infile:

        # Skip comment lines
        if line.startswith("##"):
            continue

        l = line.strip().split("\t")

        # Gather sample names and population info
        if line.startswith("#CHROM"):
            l[0] = l[0].replace("#", "")
            samples = l[9: ]

            for s in list(enumerate(l))[9: ]:
                ID = s[0]
                pop = s[1].split("_")[0]

                pop_ids[pop].append(ID)

            pops = sorted(pop_ids)
            snp = l[:3]
            infos = snp + pops + ["AFD"]

            for i, p1 in enumerate(pops):
                for p2 in pops[i+1:]:
                    infos.append(p1 + "-" + p2)

            #outfile.write("\t".join(infos) + "\n")

        # Extract allele frequencies per SNP
        else:
            snp = l[:3]
            afs = []
            afs_dict = {}
            for p in pops:

                alleles = "".join([l[i].split(":")[0].replace("/", "")
                    for i in pop_ids[p] if not l[i].startswith("./.")])

                af = alleles.count("1") / len(alleles)
                afs.append(af)
                afs_dict[p] = af

            # Add allele frequency range
            afd = max(afs) - min(afs)
            infos = snp + [f"{a:.4f}" for a in afs] + [f"{afd:.4f}"]

            # Get allele differences (do second loop)
            for i, p1 in enumerate(pops):
                for p2 in pops[i+1:]:
                    diff = abs(afs_dict[p1] - afs_dict[p2])
                    infos.append(f"{diff:.4f}")
                    pair = p1 + "-" + p2
                    pairwise_afds[pair].append((diff**exponent, snp))

            #outfile.write("\t".join(infos) + "\n")

# Select best SNPs
# This process is iterative. You choose number of SNPs for the group with the
# smallest AFD and all the other pairs will get a number of SNPs that leads to
# the same sum of AFD (or sum of AFD**2) for all pairwise comparisions. If you
# are aiming for a specific number of SNPs, adjust the number for the least
# differentiated group up or down until you get the right amount of SNPs.

afd_sums = []
for p in pairwise_afds:
    # Squared: 214 SNPs for a total of 500
    # Not squared: 127 SNPs for a total of 500
    best_afds[p] = sorted(pairwise_afds[p], reverse=True)[:num_snps]

    afd_sums.append((sum([x[0] for x in best_afds[p]]), p))

min_afd = sorted(afd_sums)[0][0]

# Shorten SNP list to have same AFD sum as min_afd
print("Number of snps kept per pairwise comparison:")
for p in best_afds:
    cumulative_sum = 0.0

    num_kept = 0
    for afd in best_afds[p]:
        diff, snp = afd
        cumulative_sum += diff

        if cumulative_sum <= min_afd:
            num_kept += 1
            best_snps.add(tuple(snp))
    print(f"  {p}, {num_kept}")

print(f"--\nKept a total of {len(best_snps)} SNPs")

with myopen(output_file + ".best_snps.tsv", "wt") as outfile:
    for snp in sorted(best_snps):
        outfile.write("\t".join(snp) + "\n")
