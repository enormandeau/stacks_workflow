# stacks_workflow revamp

Implement new features and update documentation.

## STACKS v1.30
* Make sure it works
  - Especially sstacks (give all samples at once)
  - Confirm number of SNPs

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

## Scripts and parameters
* Test the effect of --lnl-lim (-20, -10, -5)
. Discuss choosing and testing values

## Filters
* Implement vcf filtering (use pyVCF)
. General purpose filtering script (with Thierry)

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

