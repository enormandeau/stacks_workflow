#!/usr/bin/env python3
"""Get informations about overmerged / CNV loci

Usage:
    <program> input_vcf output_file
"""

# Modules
import statistics
import sys

# Parse user input
try:
    input_vcf = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read VCF and compute allelic imbalance for each SNP
with open(input_vcf) as infile:

        loci = []
        previous_locus = None

        # Iterate over loci and SNPs
        for line in infile:
            l = line.strip().split("\t")

            if line.startswith("#"):
                continue

            # Get locus information
            scaffold, position, locus_id = l[:3]
            locus = locus_id.split("_")[0]

            # Get only coverage data per sample for heterozygote samples
            data_homozygotes_freq = [(int(i[0]), int(i[1])) for i in [x.split(":")[2].split(",")
                    for x in l[9:] if x.split(":")[0] in ["0/0"]]
                    ]

            data_heterozygotes = [(int(i[0]), int(i[1])) for i in [x.split(":")[2].split(",")
                    for x in l[9:] if x.split(":")[0] in ["0/1", "1/0"]]
                    ]

            data_homozygotes_rare = [(int(i[0]), int(i[1])) for i in [x.split(":")[2].split(",")
                    for x in l[9:] if x.split(":")[0] in ["1/1"]]
                    ]

            if len(data_heterozygotes) < 2:
                continue

            # allele ratio of snp
            allele_ratios = [x[1] / sum(x) for x in data_heterozygotes]
            med_ratio = statistics.median(allele_ratios)

            # median and stdev coverage
            coverages_heterozygotes = [sum(x) for x in data_heterozygotes]
            med_coverage_heterozygotes = statistics.median(coverages_heterozygotes)
            std_coverage_heterozygotes = statistics.stdev(coverages_heterozygotes)

            coverages_homozygotes = [sum(x) for x in data_homozygotes_freq + data_homozygotes_rare]
            med_coverage_homozygotes = statistics.median(coverages_homozygotes)
            std_coverage_homozygotes = statistics.stdev(coverages_homozygotes)

            # Skip high coverage SNPs
            coverages_total = coverages_heterozygotes + coverages_homozygotes
            med_coverage_total = statistics.median(coverages_total)

            if med_coverage_total > 80:
                continue

            # proportion heterozygotes (only these with a genotype)
            prop_heterozygotes = len(coverages_heterozygotes) / (
                    len(coverages_heterozygotes) + len(coverages_homozygotes))

            # proportion homozygotes for the frequent allele
            prop_homozygotes_freq = len(data_homozygotes_freq) / (
                    len(data_homozygotes_freq) + len(data_homozygotes_rare) + len(data_heterozygotes))

            # proportion homozygotes for the rare allele
            prop_homozygotes_rare = len(data_homozygotes_rare) / (
                    len(data_homozygotes_freq) + len(data_homozygotes_rare) + len(data_heterozygotes))

            # Get all infos
            snp_info = [
                    scaffold,
                    position,
                    locus_id,
                    med_ratio,
                    med_coverage_heterozygotes,
                    std_coverage_heterozygotes,
                    med_coverage_homozygotes,
                    std_coverage_homozygotes,
                    prop_heterozygotes,
                    prop_homozygotes_freq,
                    prop_homozygotes_rare
                    ]

            if locus != previous_locus:
                loci.append([])

            loci[-1].append(snp_info)
            previous_locus = locus

with open(output_file, "w") as outfile:

    # Write header
    outfile.write("#Scaffold\tPosition\tID\tMedRatio\tMedCovHet\tStdCovHet\tMedCovHom\tStdCovHom\tPropHet\tPropHomFreq\tPropHomRare\tNumSnpsLocus\tMedLocusRatioLocus\tMedLocusCovLocus\n")

    for locus in sorted(loci):
        num_snps = len(locus)

        med_allele_ratio = statistics.median([x[3] for x in locus])

        med_allele_cov = statistics.median([x[4] for x in locus])

        for snp in locus:
            snp += [num_snps, med_allele_ratio, med_allele_cov]

            info = [str(x) for x in snp]
            outfile.write("\t".join(info) + "\n")
