# Stacks Workflow development
Features and documentation updates

## Filters
* Distribution graphs
  - Improve graphs
  - Update documentation
. (?) Look for differences of maf, Fis, Het among SNPs of a same locus
. (?) Remove individuals with more than 2 haplotypes in haplotype file
. Remove SNPs with more than 2 alleles (not sure how to do this)

## Reference genome
* Revive `stacks_1b_pstacks.sh`
* Test genome assisted STACKS
  - Use mapping options that work well with Illumina AND Ion Proton

## Paralelization
- Paralelize long steps that can use only one CPU
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

