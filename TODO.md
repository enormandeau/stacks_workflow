# Stacks Workflow development
Features and documentation updates

## Filters
* Update documentation
  - How to create distribution graphs and what graphs are produced where
  - How to combine two distribution graphs folders
  - Required python libraires (numpy and PIL --> how to install anaconda and PIL)

## Maybe
. (?) Look for differences of maf, Fis, Het among SNPs of a same locus
. (?) Remove individuals with more than 2 haplotypes in haplotype file
. Remove SNPs (haplotypes?) with more than 2 alleles (not sure how to do this)

## Reference genome
* Revive `stacks_1b_pstacks.sh`
* Test genome assisted STACKS
  - Use mapping options that work well with Illumina AND Ion Proton

## Paralelization
- Paralelize long steps that can use only one CPU or are ineficient
. Give sstacks all samples at once (faster since catalogue is read only once)

## Data preparation
* FastQC
  - FastQC on raw data (use small subset)
  - FastQC on cleaned data (use small subset)
  . Installing FastQC on MacOS

## Post STACKS
* Preparing files for population genetics analyses
. Add Laura's software suggestions and code examples
. Format conversion (PGD Spider et. al.)
. Discuss with STACKS group for other ideas

## Ion Proton
* Discuss impact of that technology
  - Indels
  - Do we lose a high proportion of reads?
  - How to solve this potential problem

