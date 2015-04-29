# stacks_workflow revamp

Implement new features and update documentation.

## Data preparation
- FastQC on raw data
- Cutadapt
  - Choosing the right primers and cut site
- FastQC on cleaned data
- Process radtags
- Combining samples
- What to check for at this point
  - Number of sequences per individual
  - Number of good individuals per population
  - Estimation of M parameter with Dan Ilut's code
- Testing different trim lengths to keep more reads

## Scripts and parameters
- Change defaults to currently used values
- Discuss choosing and testing values
- Timestamp the logs
- Copy scripts as they were when run and timestamp

## Filters
- General purpose filtering script (with Thierry)

## Post STACKS
- Preparing files for population genetics analyses
- Add Laura's software suggestions and code examples
- Format conversion (PGD Spider et. al.)
- Discuss with STACKS group for other ideas

## Documentation
- STACKS version 1.29
- Linux and MacOS versions
- Add sanity checkpoints with plots
- Add know-how sections where needed

## Ion Proton
- Impact of that technology
  - Indels
  - Do we lose a high proportion of reads?
  - How to solve this potential problem

