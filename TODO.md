# Stacks Workflow development
Features and documentation updates

## STACKS v1.34
- Read changelog of 1.33 and 1.34
- Validate all the options are the same as in `stacks_workflow`
- Modify as needed
- Create v1.34 tag and push

## Filters (`05_filter_vcf.py`)
* Add figures (eg: distribution of parameters by pop, MAF, He, Fis, GL, nbSNPs)
  - Create and improve graphs
  - Update documentation
. (?) Look for differences of maf, Fis, Het among SNPs of a same locus
. (?) Remove individuals with more than 2 haplotypes in haplotype file
. Remove SNPs with more than 2 alleles (not sure how to do this)

## Reference genome
* Revive `stacks_1b_pstacks.sh`
* Test genome assisted STACKS
  - Use mapping options that work with Illumina AND Ion Proton

## Paralelization
- Paralelize long steps that can use only one CPU
  * Cutadapt
. Give sstacks all samples at once (faster since catalogue is read only once)

## Data preparation
- FastQC
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

