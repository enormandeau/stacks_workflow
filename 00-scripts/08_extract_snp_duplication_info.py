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

def allele_depths(sample):
    """Read the expected third (AD) field; missing depths are not zero reads.

    This retains the script's existing GT:DP:AD field-order assumption.
    """
    fields = sample.split(":")
    if len(fields) < 3 or fields[2] in ("", ".", ".,."):
        return None
    values = fields[2].split(",")
    if len(values) != 2:
        raise ValueError("Expected two biallelic AD values")
    if "." in values:
        return None
    depths = tuple(int(x) for x in values)
    if min(depths) < 0:
        raise ValueError("Allele depths must be non-negative")
    return depths

def median_or_na(values):
    return statistics.median(values) if values else None

def proportion_or_na(count, total):
    return count / total if total else None

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
            called = {g: [] for g in ("0/0", "0/1", "1/0", "1/1")}
            for sample in l[9:]:
                gt = sample.split(":")[0]
                if gt in called:
                    called[gt].append(sample)
            num_samples = sum(len(samples) for samples in called.values())

            def depths_for(genotypes):
                depths = [allele_depths(sample) for g in genotypes for sample in called[g]]
                return [x for x in depths if x is not None]

            data_homozygotes_freq = depths_for(("0/0",))
            data_heterozygotes = depths_for(("0/1", "1/0"))
            data_homozygotes_rare = depths_for(("1/1",))

            # allele ratio of snp
            allele_ratios = [x[1] / sum(x) for x in data_heterozygotes if sum(x) > 0]
            med_ratio = median_or_na(allele_ratios)
            avg_ratio = statistics.mean(allele_ratios) if allele_ratios else None

            # median and stdev coverage
            coverages_heterozygotes = [sum(x) for x in data_heterozygotes]
            total_coverage_heterozygotes = sum(coverages_heterozygotes) if coverages_heterozygotes else None
            med_coverage_heterozygotes = median_or_na(coverages_heterozygotes)

            coverages_homozygotes = [sum(x) for x in data_homozygotes_freq + data_homozygotes_rare]
            coverages_homozygotes_rares = [sum(x) for x in data_homozygotes_rare]

            med_coverage_homozygotes = median_or_na(coverages_homozygotes)

            ## Skip high coverage SNPs
            #coverages_total = coverages_heterozygotes + coverages_homozygotes
            #med_coverage_total = statistics.median(coverages_total)

            # proportion heterozygotes (only these with a genotype)
            # Genotype counts must not depend on AD availability.
            num_heterozygotes = len(called["0/1"]) + len(called["1/0"])
            num_rare = len(called["1/1"]) + num_heterozygotes

            prop_heterozygotes = proportion_or_na(num_heterozygotes, num_samples)

            # proportion homozygotes for the frequent and rare allele
            prop_homozygotes_freq = proportion_or_na(len(called["0/0"]), num_samples)
            prop_homozygotes_rare = proportion_or_na(len(called["1/1"]), num_samples)

            # Compute Fis
            fis = None
            if num_samples:
                Hobs = prop_heterozygotes
                p = prop_homozygotes_freq + prop_heterozygotes / 2
                q = prop_homozygotes_rare + prop_heterozygotes / 2
                Hexp = 2 * p * q
                if Hexp > 0:
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

            info = ["NA" if x is None else str(x) for x in snp_info]
            outfile.write("\t".join(info) + "\n")
