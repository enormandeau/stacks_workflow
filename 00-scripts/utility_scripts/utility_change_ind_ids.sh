#!/bin/bash

# alleles col 2
echo "  >>> treating alleles"
id=1
for file in $(ls -1 04-stacks/*.alleles.tsv | grep -v catalog); do echo $file; awk -F "\t" -v id="${id}" '{ $2=id; print }' $file | perl -pe 's/ /\t/g' > $file".corrected"; id=$(echo $id + 1 | bc); done

# matches col 4
echo "  >>> treating matches"
id=1
for file in $(ls -1 04-stacks/*.matches.tsv | grep -v catalog); do echo $file; awk -F "\t" -v id="${id}" '{ $4=id; print }' $file | perl -pe 's/ /\t/g' > $file".corrected"; id=$(echo $id + 1 | bc); done

# snps col 2
echo "  >>> treating snps"
id=1
for file in $(ls -1 04-stacks/*.snps.tsv | grep -v catalog); do echo $file; awk -F "\t" -v id="${id}" '{ $2=id; print }' $file | perl -pe 's/ /\t/g' > $file".corrected"; id=$(echo $id + 1 | bc); done

# tags col 2
echo "  >>> treating tags"
id=1
for file in $(ls -1 04-stacks/*.tags.tsv | grep -v catalog); do echo $file; awk -F "\t" -v id="${id}" '{ $2=id; print }' $file | perl -pe 's/ /\t/g' > $file".corrected"; id=$(echo $id + 1 | bc); done

