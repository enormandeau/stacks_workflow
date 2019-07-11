#!/usr/bin/env python3
"""Get informations about overmerged / CNV loci

Usage:
    <program> input_vcf output_file
    Input and output VCFs can be compressed with gzip, ending in .gz
"""

# Modules
import statistics
import gzip
import sys

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def duplicated_likelihood(avg_ratio, total_coverage_heterozygotes, med_coverage_heterozygotes, fis):
    return -1

# Read VCF and compute allelic imbalance for each SNP
with myopen(input_vcf) as infile:
    with myopen(output_file, "wt") as outfile:

        # Write header
        outfile.write("Scaffold\tPosition\tID\tMedRatio\tAvgRatio\tMedCovHet\tTotCovHet\tMedCovHom\tNumHet\tPropHomFreq\tPropHet\tPropHomRare\tNumRare\tFis\n")

        # Iterate over loci and SNPs
        for line in infile:
            l = line.strip().split("\t")

            if line.startswith("#"):
                continue

            # Get locus information
            scaffold, position, locus_id = l[:3]
            locus = locus_id.split("_")[0]
            num_samples = len([x for x in l[9:] if x.split(":")[0] != "./."])

            # Get only coverage data per sample for heterozygote samples
            data_homozygotes_freq = [(int(i[0]), int(i[1])) for i in [x.split(":")[2].split(",")
                for x in l[9:] if x.split(":")[0] == "0/0" and "," in x.split(":")[2]]
                    ]

            data_heterozygotes = [(int(i[0]), int(i[1])) for i in [x.split(":")[2].split(",")
                    for x in l[9:] if x.split(":")[0] in ["0/1", "1/0"] and "," in x.split(":")[2]]
                    ]

            data_homozygotes_rare = [(int(i[0]), int(i[1])) for i in [x.split(":")[2].split(",")
                    for x in l[9:] if x.split(":")[0] == "1/1" and "," in x.split(":")[2]]
                    ]

            # allele ratio of snp
            allele_ratios = [x[1] / sum(x) for x in data_heterozygotes]
            try:
                med_ratio = statistics.median(allele_ratios)
                avg_ratio = statistics.mean(allele_ratios)
            except:
                med_ratio = 0.0
                avg_ratio = 0.0

            # median and stdev coverage
            coverages_heterozygotes = [sum(x) for x in data_heterozygotes]
            total_coverage_heterozygotes = sum(coverages_heterozygotes)
            try:
                med_coverage_heterozygotes = statistics.median(coverages_heterozygotes)
            except:
                med_coverage_heterozygotes = 0.0

            coverages_homozygotes = [sum(x) for x in data_homozygotes_freq + data_homozygotes_rare]
            coverages_homozygotes_rares = [sum(x) for x in data_homozygotes_rare]

            try:
                med_coverage_homozygotes = statistics.median(coverages_homozygotes)
            except:
                med_coverage_homozygotes = 0.0

            ## Skip high coverage SNPs
            #coverages_total = coverages_heterozygotes + coverages_homozygotes
            #med_coverage_total = statistics.median(coverages_total)

            # proportion heterozygotes (only these with a genotype)
            num_heterozygotes = len(coverages_heterozygotes)
            num_rare = len(coverages_homozygotes_rares) + len(coverages_heterozygotes)

            prop_heterozygotes = len(coverages_heterozygotes) / num_samples

            # proportion homozygotes for the frequent and rare allele
            prop_homozygotes_freq = len(data_homozygotes_freq) / num_samples
            prop_homozygotes_rare = len(data_homozygotes_rare) / num_samples

            # Compute Fis
            Hobs = prop_heterozygotes
            p = prop_homozygotes_freq + prop_heterozygotes / 2
            q = prop_homozygotes_rare + prop_heterozygotes / 2
            Hexp = 2 * p * q

            fis = 1 - Hobs / Hexp

            # Compute likelihood of being a single / duplicated SNP
            likelihood = duplicated_likelihood(
                    avg_ratio,
                    total_coverage_heterozygotes,
                    med_coverage_heterozygotes,
                    fis
                    )

            # Get all infos
            snp_info = [
                    scaffold,
                    position,
                    locus_id,
                    med_ratio,
                    avg_ratio,
                    med_coverage_heterozygotes,
                    total_coverage_heterozygotes,
                    med_coverage_homozygotes,
                    num_heterozygotes,
                    prop_homozygotes_freq,
                    prop_heterozygotes,
                    prop_homozygotes_rare,
                    num_rare,
                    fis
                    ]

            info = [str(x) for x in snp_info]
            outfile.write("\t".join(info) + "\n")
