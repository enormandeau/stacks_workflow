# stacks_workflow revamp
Features and documentation updates.

## STACKS v1.30 (1.31 upcoming soon, mostly minor bugfixes)
* Make sure it works for STACKS v1.30
* Give sstacks all samples at once

## Introduction
- Other papers to read? (01-stacks.rst)
- Other options to use with ./configure for stacks? (02-Step_2.rst)

## Data preparation
- FastQC
  - FastQC on raw data
  - FastQC on cleaned data
  . Installing FastQC on MacOS
* Parallelize cutadapt

## What to check for at this point
. Number of sequences per individual
. Number of good individuals per population
* Estimation of M parameter with Dan Ilut's code (can we distribute it?)

## Using a reference genome
* Revive the stacks_1b_pstacks.sh script
* Test genome assisted STACKS
  - Use mapping options that work with Illumina AND Ion Proton

## Filters
# TODO list for 05_filter_vcf.py
* Add filters:
  * Fis (calculate from VCF)
  * Output a Whitelist (to rerun populations)
  - Allele imbalance (correct genotypes?)
  . Genotype likelihood (>10)

- Look for differences of maf, Fis, Het among SNPs of a same locus (?)
- Remove individuals with more than 2 haplotypes in haplotype file

## Post STACKS
. Preparing files for population genetics analyses
. Add Laura's software suggestions and code examples
. Format conversion (PGD Spider et. al.)
. Discuss with STACKS group for other ideas

## Documentation
. Add sanity checkpoints with plots (ask Laura)
. Add know-how sections where needed

## Ion Proton
* Discuss impact of that technology
  - Indels
  - Do we lose a high proportion of reads?
  - How to solve this potential problem

