# stacks_workflow revamp

Implement new features and update documentation.

## Introduction
* Other papers to read? (01-stacks.rst)
* Other options to use with ./configure for stacks? (02-Step_2.rst)

## Data preparation
. Installing FastQC on MacOS
- FastQC on raw data
- FastQC on cleaned data

* Combining samples (simple and complex)
  * Create and add `complex` script

. What to check for at this point
  . Number of sequences per individual
  . Number of good individuals per population
  . Estimation of M parameter with Dan Ilut's code

### Stop here for now ###

## Timestamp the logs
* Add timestamp to log files
* Add script as run, also with timestamp

## Scripts and parameters
* Test the effect of --lnl-lim (-20, -10, -5)
- Discuss choosing and testing values
- Copy scripts as they were when run and timestamp

## Filters
- General purpose filtering script (with Thierry)

## Post STACKS
- Preparing files for population genetics analyses
- Add Laura's software suggestions and code examples
- Format conversion (PGD Spider et. al.)
- Discuss with STACKS group for other ideas

## Documentation
- Add sanity checkpoints with plots
- Add know-how sections where needed

## Ion Proton
- Impact of that technology
  - Indels
  - Do we lose a high proportion of reads?
  - How to solve this potential problem

