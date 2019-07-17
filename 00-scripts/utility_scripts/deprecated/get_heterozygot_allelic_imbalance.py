#!/usr/bin/env python3
"""Get allelic imbalance for heterozygote samples for each SNP
in order to find paralogs

Usage:
    <program> input_vcf output_file
"""

# Modules
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
    with open(output_file, "w") as outfile:

        # Write header
        outfile.write("Locus\tPosition\tID\tNumSamples\tMedCov\tMedRatio\n")

        # Iterate over SNPs
        for line in infile:
            l = line.strip().split("\t")

            if line.startswith("#"):
                continue

            # Get locus information
            locus, position, locus_id = l[:3]

            # Get only coverage data per sample for heterozygote samples
            data = [(int(i[0]), int(i[1])) for i in [x.split(":")[2].split(",")
                    for x in l[9:] if x.split(":")[0] in ["0/1", "1/0"]]
                    ]

            coverages = sorted([sum(x) for x in data])

            num_samples = len(data)

            if num_samples == 0:
                continue

            ratios = sorted([x[1] / sum(x) for x in data])
            cov = coverages[int(len(coverages)/2)]
            
            # Mean
            #data = sum(ratios) / len(ratios)

            # Median
            med_ratio = ratios[int(len(ratios)/2)]

            output = "\t".join([
                locus, position, locus_id, str(num_samples), str(cov), str(med_ratio)
                ]) + "\n"

            outfile.write(output)
