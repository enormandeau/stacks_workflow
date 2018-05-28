#!/usr/bin/env python3
"""Summarize reproducibility of RADSeq using duplicates

Usage:
    <program> input_vcf separator_character output_filename
"""

# Module
from collections import defaultdict
import sys

# Parse user input
try:
    input_vcf = sys.argv[1]
    separator_character = sys.argv[2]
    output_filename = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Find sample names and columns of pairs
with open(input_vcf) as infile:
    for line in infile:
        if line.startswith("#CHROM"):
            l = line.strip().split("\t")
            sample_names = l[9:]
            break

sample_name_stubs = sorted(list(set([x.split(separator_character)[0] for x in sample_names])))
pairs = dict()

for stub in sample_name_stubs:
    sample1, sample2 = sorted([x for x in sample_names if x.startswith(stub)])
    column1 = sample_names.index(sample1)
    column2 = sample_names.index(sample2)
    pairs[stub] = [sample1, column1, sample2, column2]

# Read input_vcf and collect informations about pairs
geno_dictionary = {
        "./.": -9,
        "0/0": 0,
        "0/1": 1,
        "1/0": 1,
        "1/1": 2}

# duplicates_dict[pair][reproducibility][count]
duplicates_dict = defaultdict(lambda: defaultdict(int))
duplicates_per_coverage_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
duplicate_cases = ["both_missing",
        "one_missing",
        "both_homo_frequent",
        "both_hetero",
        "both_homo_rare",
        "hetero_homo",
        "diff_homo"]

with open(input_vcf) as infile:
    with open(output_filename, "w") as outfile:
        header = "Contig\tcontig_pos\tlocus_id\tsnp_pos\tallele1\tallele2\tnum_samples\tmaf\tsample1\tgeno1\tcov1\ta1_cov1\ta2_cov1\tsample2\tgeno2\tcov2\ta1_cov2\ta2_cov2\n"
        outfile.write(header)

        for line in infile:

            # Ignore VCF header
            if line.startswith("#"):
                continue

            # Gather SNP infos
            infos = line.strip().split("\t")[:8]
            contig, contig_pos, locus_id, allele1, allele2, _, _, num_and_maf = infos
            snp_pos = locus_id.split("_")[1]
            num_and_maf = num_and_maf.split(";")
            number_samples = num_and_maf[0]
            allele_freq = num_and_maf[1]
            number_samples = number_samples.split("=")[1]
            allele_freq = allele_freq.split("=")[1]

            data = line.strip().split("\t")[9:]

            # Gather infos for each sample/duplicate pair
            for stub in pairs:
                sample1, column1, sample2, column2 = pairs[stub]

                data_column1 = data[column1].split(":")
                geno1 = data_column1[0]
                cov1 = data_column1[1]
                allele_cov1 = data_column1[2]
                geno1 = geno_dictionary[geno1]
                cov1_a1, cov1_a2 = allele_cov1.split(",")

                data_column2 = data[column2].split(":")
                geno2 = data_column2[0]
                cov2 = data_column2[1]
                allele_cov2 = data_column2[2]
                geno2 = geno_dictionary[geno2]
                cov2_a1, cov2_a2 = allele_cov2.split(",")

                min_cov = str(min([int(cov1), int(cov2)]))

                if geno1 == geno2 == -9:
                    duplicates_dict[stub]["both_missing"] += 1
                    duplicates_per_coverage_dict[stub]["both_missing"][min_cov] += 1

                elif geno1 == -9 or geno2 == -9:
                    duplicates_dict[stub]["one_missing"] += 1
                    duplicates_per_coverage_dict[stub]["one_missing"][min_cov] += 1

                elif geno1 == geno2 == 0:
                    duplicates_dict[stub]["both_homo_frequent"] += 1
                    duplicates_per_coverage_dict[stub]["both_homo_frequent"][min_cov] += 1

                elif geno1 == geno2 == 1:
                    duplicates_dict[stub]["both_hetero"] += 1
                    duplicates_per_coverage_dict[stub]["both_hetero"][min_cov] += 1

                elif geno1 == geno2 == 2:
                    duplicates_dict[stub]["both_homo_rare"] += 1
                    duplicates_per_coverage_dict[stub]["both_homo_rare"][min_cov] += 1

                elif geno1 == 1 and geno2 in [0, 2]:
                    duplicates_dict[stub]["hetero_homo"] += 1
                    duplicates_per_coverage_dict[stub]["hetero_homo"][min_cov] += 1

                elif geno1 in [0, 2] and geno2 == 1:
                    duplicates_dict[stub]["hetero_homo"] += 1
                    duplicates_per_coverage_dict[stub]["hetero_homo"][min_cov] += 1

                elif (geno1 == 0 and geno2 == 2) or (geno1 == 2 and geno2 == 0):
                    duplicates_dict[stub]["diff_homo"] += 1
                    duplicates_per_coverage_dict[stub]["diff_homo"][min_cov] += 1

                else:
                    print("You missed a case, dumass!")
                    print(geno1, geno2)
                    sys.exit(1)

                # Report all infos
                report = "\t".join([contig, contig_pos, locus_id, snp_pos,
                        allele1, allele2, number_samples, allele_freq,
                        sample1, str(geno1), cov1, cov1_a1, cov1_a2,
                        sample2, str(geno2), cov2, cov2_a1, cov2_a2]) + "\n"

                outfile.write(report)

# Summary
with open(output_filename + "_summary.tsv", "w") as outfile:
    for sample in pairs:
        for case in duplicate_cases:
            outfile.write(
                    "\t".join([sample, case,
                        str(sum(duplicates_per_coverage_dict[sample][case].values()))]) +
                    "\n")

# By coverage
with open(output_filename + "_summary_by_coverage.tsv", "w") as outfile:
    for sample in pairs:
        for case in duplicate_cases:
            for cov in range(4, 101):
                outfile.write(
                        "\t".join([sample, case, str(cov),
                            str(duplicates_per_coverage_dict[sample][case][str(cov)])]) +
                        "\n")
