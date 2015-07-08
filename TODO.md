# Stacks Workflow development
Features and documentation updates

## STACKS v1.32
* Copy and timestamp `population_map.csv` and `sample_information.csv` when used
* Modify `03_rename_samples_complex.sh` to create links for samples with only
  one source to save lots of space in some cases and remove need for 'simple'
  and 'complex' scripts.
* Give sstacks all samples at once (faster since catalogue is read only once)
* Add v1.32 tag

## Reference genome
* Revive the stacks_1b_pstacks.sh script
* Test genome assisted STACKS
  - Use mapping options that work with Illumina AND Ion Proton

## Introduction
- Other papers to read? (01-stacks.rst)

## Data preparation
- FastQC
  - FastQC on raw data
  - FastQC on cleaned data
  . Installing FastQC on MacOS
* Parallelize cutadapt

## What to check for at this point
. Number of sequences per individual
. Number of good individuals per population
. (?) Estimation of M parameter with Dan Ilut's code (can we distribute it?)

## Filters
# TODO list for 05_filter_vcf.py
* Add figures (eg: distribution of parameters by pop)
* Add filters:
  * Fis (calculate from VCF)
  * Output a Whitelist (to rerun populations)
  - Allele imbalance (correct genotypes?)
  - Genotype likelihood threshold

- (?) Look for differences of maf, Fis, Het among SNPs of a same locus
- Remove individuals with more than 2 haplotypes in haplotype file

## Post STACKS
* Preparing files for population genetics analyses
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

