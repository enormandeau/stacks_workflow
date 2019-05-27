#!/usr/bin/env python
"""Keep only the first SNP per locus in vcf file

Usage:
    ./extract_first_snp.py input_vcf output_vcf
"""

# Modules
import sys

# Functions

# Main
if __name__ == '__main__':
    # Parse user input
    try:
        input_vcf = sys.argv[1]
        output_vcf = sys.argv[2]
    except:
        print(__doc__)
        sys.exit(1)

    # Open file handles
    infile = open(input_vcf)
    outfile = open(output_vcf, "w")

    # Extract header and SNPs
    extracted_snps = set()
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
        else:
            snp_id = line.strip().split()[2]
            if "_" in snp_id:
                snp_id = snp_id.split("_")[0]

            if snp_id not in extracted_snps:
                outfile.write(line)
                extracted_snps.add(snp_id)

    # Close file handles
    infile.close()
    outfile.close()
