#!/bin/bash

# Based on 07_populations_low_ram.sh by pop VCF files, filter catalog
for i in 09-vcfs_by_pop/*.vcf
do
    grep -v "^#" "$i" |
    cut -f 3 |
    cut -d "_" -f 1 |
    uniq
done | sort -nu > 09-vcfs_by_pop/wanted_loci_for_catalog.txt

# Filter all VCF files

# Get list of wanted loci

# Filter Catalog

# Backup full catalog

# Put filtered catalog in place

# Report

echo "---"
echo "Catalog filtered: $n loci kept"
echo "You can now re-run sstacks and populations"
