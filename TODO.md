# Stacks Workflow development
Features and documentation updates

## Filters
# TODO list for 05_filter_vcf.py
* Add figures (eg: distribution of parameters by pop, MAF, He, Fis, GL, nbSNPs)
* Add filters:
  * Fis (calculate from VCF)
  * Output a Whitelist (to rerun populations)
  - Maximum allele coverage
  - Allele imbalance (correct genotypes?)
  - Genotype likelihood threshold
  - SNPs with more than 2 alleles

- (?) Look for differences of maf, Fis, Het among SNPs of a same locus
- (?) Remove individuals with more than 2 haplotypes in haplotype file

## Paralelization
- Paralelize long steps that can use only one CPU
. Give sstacks all samples at once (faster since catalogue is read only once)

## Reference genome
* Revive the stacks_1b_pstacks.sh script
* Test genome assisted STACKS
  - Use mapping options that work with Illumina AND Ion Proton

## Data preparation
- FastQC
  - FastQC on raw data
  - FastQC on cleaned data
  . Installing FastQC on MacOS
* Parallelize cutadapt

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

