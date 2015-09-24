# Stacks Workflow development
Features and documentation updates

## Filters
* Fix number of retained loci (currently wrong in both filter and graphics mode)
* Remove individuals with more than 2 haplotypes in haplotype file (write utility script)
* Update documentation
  - How to create distribution graphs, what graphs are produced, and where they are
  - How to combine two distribution graphs folders
  - Required python libraires (numpy and PIL --> how to install anaconda and PIL)

## Reference genome
* Test genome assisted STACKS
  - Use mapping options that work well with Illumina AND Ion Proton

## Maybe
* Revive `stacks_1b_pstacks.sh`

## Paralelization
- Paralelize long steps that can use only one CPU or are ineficient
  - (done) cutadapt
  - (done) process radtags
  - ustacks ?
  - sstacks ?
. Give sstacks all samples at once (faster since catalogue read only once)

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

