#!/usr/bin/env python2
"""Keep only the SNP with the highest maf in each locus in a STACKS 1.x vcf file

Usage:
    ./extract_snp_with_max_maf.py input_vcf output_vcf
"""

# Modules
import sys

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
    extracted_snps = dict()

    line_counter = 0
    for line in infile:
        line_counter += 1

        # Header
        if line.startswith("#"):
            outfile.write(line)

        # Loci
        else:
            snp_id = line.strip().split()[2]
            maf = float(line.strip().split()[7].split(";")[1].split("=")[1])

            # Correct for snp id format
            if "_" in snp_id:
                snp_id = snp_id.split("_")[0]

            # Decide if we keep the snp
            if snp_id not in extracted_snps:
                extracted_snps[snp_id] = (line_counter, maf, line)
                #print "new", snp_id, maf
                #print maf

            else:
                if maf > extracted_snps[snp_id][1]:
                    #print "better", snp_id, ">>", maf
                    extracted_snps[snp_id] = (line_counter, maf, line)

    snps = sorted(extracted_snps.values())
    for snp in snps:
        line = snp[2]
        outfile.write(line)

    # Close file handles
    infile.close()
    outfile.close()
