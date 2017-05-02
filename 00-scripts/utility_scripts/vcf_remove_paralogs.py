#!/usr/bin/env python
"""Remove or clean loci with paralogs

Usage:
    ./vcf_remove_paralogs.py input_vcf paralog_info threshold output_vcf

Where:
    input_vcf is the VCF file to treat
    paralog_info contains information about paralogous loci (see format below)
    threshold is the maximum number of individuals to keep the locus by removing
        only the bad genotypes. Above this threshold, the entire locus is removed
    output_vcf is the output VCF file

Format of paralog_info file:
    NOTE: Header lines MUST start with a '#' sign

#LOCUS  POP_ID  INDIVIDUALS  HAPLOTYPES
6      JCF     JCF_14_Sf    GCA/GTA/TTA
7      ABE     ABE_14_Sf    AAATGTG/ACATATG/ACATGTG/GAATGTG
7      ABE     ABE_16_Sf    ACAGATG/ACATGTG/GCATATG/GCATGTG
7      ABE     ABE_5_Sf     AAATATG/ACATATG/ACATGTG/GAATGTG
7      ABE     ABE_7_Sf     ACAGATG/ACATATG/ACATGTG
7      AMA     AMA_30_Sf    ACAGATG/ACATATG/ACATGTG
7      ARB     ARB_5_Sf     ACATATG/ACATGTG/GAATGTG
7      BRO     BRO_10_Sf    AACGGTG/AACTGTG/ACAGGTG/ACATGTG
7      BRU     BRU_10_Sf    AAAGATG/AACGATG/AACTATG/AACTGTG/ACAGATG/ACATATG/ACATGTG
"""

# Modules
from collections import defaultdict
from collections import Counter
import sys

# Classes

# Functions

# Parsing user input
try:
    input_vcf = sys.argv[1]
    paralog_info = sys.argv[2]
    threshold = int(sys.argv[3])
    output_vcf = sys.argv[4]
except:
    print __doc__
    sys.exit(1)

# Parsing paralog_info file
loci = defaultdict(set)

with open(paralog_info) as hapfile:
    for line in hapfile:
        if not line.startswith("#"):
            locus, pop, sample, haplotype = line.strip().split("\t")
            loci[locus].add(sample)
            
num_paralog_per_locus = []
for i in loci.values():
    num_paralog_per_locus.append(len(i))

paralog_frequency = sorted(Counter(num_paralog_per_locus).items())
with open(paralog_info + ".paralog_frequency", "w") as freqfile:
    for i in paralog_frequency:
        freqfile.write(str(i[0]) + "\t" + str(i[1]) + "\n")

# Treating the VCF
vcffile = open(input_vcf)
outputfile = open(output_vcf, "w")

sample_pos = dict()

for line in vcffile:
    # Treat header lines
    if line.startswith("#"):
        outputfile.write(line)

        # Get sample positions in VCF
        if line.startswith("#CHROM"):
            infos = line.strip().split()[9:]
            for index, sample in enumerate(infos):
                sample_pos[index] = sample

    # Treate loci infos
    else:
        l = line.strip().split()
        locus = l[2]
                    
        if "_" in locus:
            locus = locus.split("_")[0]

        # Locus contains paralogs
        if locus in loci:
            num_bad = len(loci[locus])

            # If num paralogous samples > threshold, flush
            if num_bad > threshold:
                continue

            # If num paralogous samples <= threshold, correct genotypes
            else:
                meta = l[:9]
                infos = l[9:]
                good_infos = []
                for index, genotype in enumerate(infos):
                    sample = sample_pos[index]
                    if sample in loci[locus]:
                        genotype = "./." + genotype[3:]

                    good_infos.append(genotype)

                line = "\t".join(meta + good_infos) + "\n"
                outputfile.write(line)

        # Locus did not contain paralogs
        else:
            outputfile.write(line)

# Close files
vcffile.close()
outputfile.close()
