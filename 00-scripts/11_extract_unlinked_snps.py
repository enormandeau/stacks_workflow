#!/usr/bin/env python3
"""Keep only one SNP per linked group within a locus.

If there are multiple SNPs in a linked group, the one on the left is kept.

Usage:
    <program> input_vcf diff_threshold output_vcf

Where:
    input_vcf: name of a VCF from STACKS 1.4x (not vcftools...)

    diff_threshold: minimum difference between 0.0 and 1.0 to keep a second
    SNP. The recommended value is 0.5. A good value is anywhere between 0.2 and
    0.5. A value of 0.2 is more permissive and will retain a few false
    positives, ie: SNPs with low MAFs that have exactly the same information
    but where some one or a few samples were mis-genotyped. A value of 0.5 will
    get you 99.9% of what is different without false positives. Values above
    0.5 will lose you true positives.

    output_vcf: name of the output filtered VCF
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

def write_snp_line(l, outfile):
    outfile.write("\t".join(l) + "\n")

def encode_genotype(g):
    """Change genotypes into 0, 1, and 2 from 0/0, 0/1 or 1/0, and 1/1
    """
    code = {
            "./.": -99,
            "0/0": 0,
            "0/1": 1,
            "1/0": 1,
            "1/1": 2
            }

    return code[g]

def invert_genotypes(s):
    return [-1 * (g - 2) if g != -99 else g for g in s]
    
def diff_pair(s1, s2):
    """Compute difference [0.0 - 1.0] of genotypes for two SNPs

    Use only samples without missin data at both SNPs and where at least one of
    the two genotypes contains the rare variant.
    """

    genotypes = [[], []]

    for i in range(len(s1)):
        g1 = s1[i]
        g2 = s2[i]

        # Find differences between two original SNPs
        if (g1 != 0 or g2 != 0) and not (g1 == -99 or g2 == -99):
            genotypes[0].append(g1)
            genotypes[1].append(g2)

    differences = 0

    for i in range(len(genotypes[0])):
        if genotypes[0][i] != genotypes[1][i]:
            differences += 1

    try:
        return (differences / len(genotypes[0]))
    except:
        # Special cases "LOW-MAF-AND-MISSING-DATA"
        #
        # either:
        # 1. One SNP has an MAF of 0.0
        # 2. Both SNPs have a low MAF and all the samples in the SNP without the
        #    rare allele have missing data
        #
        # In both cases, return 0.0 to keep only the first SNP

        print("Found special case LOW-MAF-AND-MISSING-DATA")
        print(genotypes)
        return 0.0

def difference(s1, s2):
    """Return minimum expected difference
    """

    s1 = [encode_genotype(x.split(":")[0]) for x in s1[9:]]
    s2 = [encode_genotype(x.split(":")[0]) for x in s2[9:]]
    s2inversed = invert_genotypes(s2)

    minimum = min(diff_pair(s1, s2), diff_pair(s1, s2inversed))
    return minimum

def keep_all_different(snps, diff_threshold, outfile):
    """Keep only SNPs with different information)

    - always keep first SNP (write to file directly)
    - go through all remaining and remove unwanted
    - While 2+ SNPs remaining, repeat
    """

    # If zero SNP, do nothing (for recursion)
    if len(snps) == 0:
        pass

    # Write lone SNPs without further work
    elif len(snps) == 1:
        write_snp_line(snps[0], outfile)

    # If 2+ SNPs, write first, prune others, and recurse on remaining
    else:
        first_snp = snps[0]
        write_snp_line(first_snp, outfile)
        other_snps = snps[1:]

        different_snps = [snp for snp in other_snps if difference(first_snp, snp) > diff_threshold]

        keep_all_different(different_snps, diff_threshold, outfile)

# Main
if __name__ == '__main__':
    # Parse user input
    try:
        input_vcf = sys.argv[1]
        diff_threshold = float(sys.argv[2])
        output_vcf = sys.argv[3]
    except:
        print(__doc__)
        sys.exit(1)

    # Open file handles
    infile = myopen(input_vcf)
    outfile = myopen(output_vcf, "wt")

    # Iterate through SNPs
    snps = []

    for line in infile:

        # Treat header lines
        if line.startswith("#"):

            # Get STACKS version
            if line.startswith("##source="):
                stacks_version = line.split('"')[1].split(" ")[1]

                if "v1" in stacks_version:
                    print(f"# STACKS version: {stacks_version}")
                    locus_separator = "_"
                    stacks_version = 1

                elif "v2" in stacks_version:
                    print(f"# STACKS version: {stacks_version}")
                    locus_separator = ":"
                    stacks_version = 2

            elif line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                populations = set([x.split("_")[0] for x in samples])
                print(f"# VCF contains {len(samples)} samples in {len(populations)} populations")

            # Write header to output file
            outfile.write(line)
            continue

        # Parse SNP data
        l = line.strip().split()
        locus = l[2].split(locus_separator)[0]
        if locus == ".":
            locus = l[0]

        # First SNP of VCF
        if snps == []:
            snps.append(l)
            previous_locus = locus

        # SNP in same locus, add to list
        elif locus == previous_locus:
            snps.append(l)

        else:
            # Treat locus
            keep_all_different(snps, diff_threshold, outfile)

            # Re-initialize snps and previous_locus
            snps = [l]
            previous_locus = locus

    ## Treat last locus
    keep_all_different(snps, diff_threshold, outfile)

    # Close file handles
    infile.close()
    outfile.close()
