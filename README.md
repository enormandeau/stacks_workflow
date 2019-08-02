# Stacks Workflow

RADseq workflow using [STACKS](http://creskolab.uoregon.edu/stacks/)

Developed by [Eric Normandeau](https://github.com/enormandeau) in
[Louis Bernatchez](http://www.bio.ulaval.ca/louisbernatchez/presentation.htm)'s
laboratory.

**Warning!**: this software is provided "as is", without warranty of any kind,
express or implied, including but not limited to the warranties of
merchantability, fitness for a particular purpose and noninfringement. In no
event shall the authors or copyright holders be liable for any claim, damages
or other liability, whether in an action of contract, tort or otherwise,
arising from, out of or in connection with the software or the use or other
dealings in the software.

stack_workflow was developed with the needs of our research group in mind. We
make no claim about its usefulness to other groups or in other contexts, but we
know it has been and continues to be useful to other groups.

**NOTE**: Works with STACKS2 when using a reference genome (2019-06-17)
**NOTE**: Implementing STACKS2 *without* a reference genome (2019-07-29)

## Licence

stacks_workflow is licensed under the gpl3 license. See the LICENCE file
provided with stacks_workflow for more details.

## Downloading

[click here to
download](https://github.com/enormandeau/stacks_workflow/archive/master.zip)
the current version of the stacks workflow. alternatively, you can clone this
repository with:

```bash
git clone https://github.com/enormandeau/stacks_workflow
```

## About STACKS

The [stacks analysis pipeline](http://creskolab.uoregon.edu/stacks/) is used
for restriction-site associated DNA sequencing (RADseq) studies, with and
without a reference genome.

Before starting to use STACKS, you should read the official STACKS papers
found at the bottom of the official STACKS page (see link above).

# Stacks workflow tutorial

The goal of this workflow is to simplify the use of the STACKS pipeline and
make the analyses more reproducible. Special care is taken to have a
standardized SNP filtering procedure to produce quality SNPs at the end.

## Overview of the steps

1. Install stacks_workflow and STACKS (or STACKS2 if you have a reference genome) with its dependencies
1. Download your raw data files (illumina lanes or ion proton chips)
1. Clean the reads and assess their quality
1. Extract individual data with process_radtags
1. Rename the sample files
1. Align reads to a reference genome (optional)
1. Stacks pipeline
1. Filtering the results

## Where to find this manual

To read an HTML version of this file online, go to the
[stacks_workflow project page](https://github.com/enormandeau/stacks_workflow)
on GitHub.

## Installing stacks_workflow

### Download and install the most recent version of this workflow

It is recommended to download the most recent version of stacks_workflow **for
each new project**, as opposed to re-using the same directory for multiple
projects and naming the outputs differently.  This is **central to
`stack_workflow`'s philosophy** of reproducibility.

**One stacks_workflow folder should contain only one analysis**

Deviate from this at your own risk ;)

#### From the terminal

```bash
cd ~/desktop
wget https://github.com/enormandeau/stacks_workflow/archive/master.zip
unzip master.zip
```

#### Using git

If `git` is installed on your computer, you can run the following command
instead to get a complete `stacks_workflow` git repository.

```bash
git clone https://github.com/enormandeau/stacks_workflow
```

For the rest of the project, use the extracted or cloned folder as your working
directory. **All the commands in this manual are launched from that directory.**

### Download and install [STACKS](http://creskolab.uoregon.edu/stacks/)

**NOTE**: `stacks_workflow` works with STACKS2 when using a reference genome (2019-06-17)
**NOTE**: Implementing STACKS2 *without* a reference genome (2019-07-29)

Follow the instructions on the website to install and test the installation with:

```bash
populations
which populations
```

This will output the help of the populations program and where it is located
on your computer. You will also be able to confirm the version number of
your STACKS installation.

#### Installing Cutadapt

There are different ways you can install Cutadapt. If you have `pip` (a Python
package manager) installed, you can use the following command:

```bash
sudo pip install --user --upgrade cutadapt
```

Otherwise, [visit their website to download it and install
it](https://pypi.python.org/pypi/cutadapt/)

## Prepare your raw data files

### Downloading your data

Download your raw Illumina or Ion Proton data files from your sequencing
service provider.

### Copy your raw files

Put a copy of (or a link to) your raw data files in the `02-raw` folder of
stacks_workflow.

All file names **must** end with **.fastq.gz** for the following scripts to
work.

### Preparing the `lane_info.txt` file

This file will contain the names of the raw data files and is used by
stacks_workflow later. From the stacks_workflow folder, run:

```bash
./00-scripts/00_prepare_lane_info.sh
```

### Running Cutadapt

We trim our data using Cutadapt **in single-end mode** with the following
command:

```bash
./00-scripts/01_cutadapt.sh numCPUs
```

Where `numCPUs` is the number of CPUs you wish to use in parallel. If you do
not put a number, the script will use only one CPU.

### Scan the cutadapt logs

The cutadapt log files can be found in the `10-log_files` folder. Scan them
to confirm that cutadapt has done an appropriate job.

There may be differences in adapters and filter parameters to use with data
produced by Illumina and Ion Proton sequencers.

## Extract individual data with `process_radtags`

### Prepare a file named `sample_information.csv`

Use the same format found in the `example_sample_information.csv` file in
the `01-info_files` folder.

Save this file in the `01-info_files` folder and name it exactly
`sample_information.csv`.

The `sample_information.csv` file will be used to extract the samples and
rename the extracted sample files automatically.

The first column **MUST** contain the **EXACT** name of the data file for the
lane of each sample.

**Notes**:

- The **columns are separated by tabulations** (even if the extension is .csv)
- The second column contains the barcode sequence of each sample.
- The third column contains the population name of each sample.
- The fourth column contains the name of the sample (do not include the
population name or abbreviation in the sample name).
- Neither the population name nor the sample name should contain underscores `_`
- The fifth column contains a number identifying the populations.
- The sixth column contains the plate well identifier.

Columns three and four are treated as text, so they can contain either text or
numbers. Other columns can be present after the fifth one and will be ignored.
However, it is crucial that the six first columns respect the format in the
example file exactly. Be especially careful not to include errors in this file,
for example by mixing lower and capital letters in population or sample names
(e.g.: Pop01 and pop01), since these will be treated as two different
populations.

### Launch process_radtags

#### One restriction enzyme

```bash
./00-scripts/02_process_radtags.sh <trimLength> <enzyme>
```

Where:

 - **trimLength** = length to trim all the sequences. This should be the length
   of the Illumina reads minus the length of the longest tag or MID.
 - **enzyme** = name of enzyme (run `process_radtags`, without options, for a
   list of the supported enzymes)

#### Two restriction enzymes

```bash
./00-scripts/02_process_radtags_2_enzymes.sh <trimLength> <enzyme1> <enzyme2>
```

Where:

 - **trimLength** = length to trim all the sequences. This should be the length
   of the Illumina reads minus the length of the longest tag or MID.
 - **enzyme1** = name of the first enzyme (run `process_radtags`, without
   options, for a list of the supported enzymes)
 - **enzyme2** = name of the second enzyme (run `process_radtags`, without
   options, for a list of the supported enzymes)

#### Two restriction enzymes in parallel over multiple CPUs

```bash
./00-scripts/02_process_radtags_2_enzymes_parallel.sh <trimLength> <enzyme1> <enzyme2> <numCPUs>
```

Where:

 - **trimLength** = length to trim all the sequences. This should be the length
   of the Illumina reads minus the length of the longest tag or MID.
 - **enzyme1** = name of the first enzyme (run `process_radtags`, without
   options, for a list of the supported enzymes)
 - **enzyme2** = name of the second enzyme (run `process_radtags`, without
   options, for a list of the supported enzymes)
 - **numCPUs** = number of CPUs to use

#### Testing different trim lengths

If you are using Ion Proton data, the effect of the trimLength parameter used
above on the number of usable SNPs you recover at the end may not be trivial.
As a rule of thumb, a trimmed length of 80bp should produce good results. We
suggest you run tests with a smaller group of samples to determine what length
to trim to. For highly species with high genetic variability, short loci will
be more likely to contain SNPs and long loci to contain more than one SNP,
which is not always informative. Thus, trimming to shorter lengths may be more
interesting for highly variant species or when coverage is limiting.

### Rename samples

#### Rename and copy the samples

We provide a script to rename the extracted samples and move them into the
`04-all_samples` folder. The script behaves differently for samples that are
present only once in the `01-info_files/sample_information.csv` and for those
that are present more than once. If a sample is present only once, a link is
created, using no additional disk space. If it is present more than once,
all the copies are concatenated into one file, doubling the amount of disk
space taken by this sample (all the individual files PLUS the combined one).

```bash
./00-scripts/03_rename_samples.sh
```

#### Deleting samples with too few reads

If after splitting your samples you notice that some have too few reads, you
can remove these from the `04-all_samples` folder. The threshold for the
minimum number of reads will depend on your project, including on the number of
expected cut sites generated by your library preparation protocol. Keep samples
with low coverages if you are not sure what threshold to use at this point. We
will filter the VCF for this later and will then have better information.

#### Align reads to a reference genome (optional)

#### Install `bwa <http://bio-bwa.sourceforge.net>`

#### Download reference genome to the `08-genome`

#### Index the reference genome

```bash
bwa index -p 08-genome/genome ./08-genome/<genome reference>
```

### Align samples

Different bwa alignment scripts are available in 00-scripts.

**NOTE**: The two scripts for single-end reads (ie: not the one with `PE` in
its name) have options that are specific for IonProton data. To align Illumina
data, remove the `-O 0,0` and `-E 2,2` options.

```bash
00-scripts/bwa_mem_align_reads.sh
00-scripts/bwa_mem_align_reads_by_n_samples.sh
00-scripts/bwa_mem_align_reads_PE.sh
```

## STACKS pipeline

### Prepare population info file

```bash
./00-scripts/04_prepare_population_map.sh
```

## Edit script parameters

You will need to go through the scripts named `stacks_*` in the `00-scripts
folder` and edit the options to suit your needs. Depending on your project
(eg: *de novo* vs reference), you will not use all the scripts.

**Warning!** This step is very important. Choosing appropriate parameters for
your study is crucial in order to generate meaningful and optimal results. Read
the STACKS documentation on their website to learn more about the different
options.

## Run the STACKS1 programs

### Without a reference genome

```bash
./00-scripts/stacks_1a_ustacks.sh
./00-scripts/stacks_2_cstacks.sh
./00-scripts/stacks_3_sstacks.sh
./00-scripts/stacks_4_populations.sh
./00-scripts/stacks_5a_rxstacks_likelihoods.sh
```

Visualize the distribution of log likelihoods in
`rxstacks_log_likelihoods.png` and choose a cutoff to use in the next script
(`./00-scripts/stacks_5b_rxstacks.sh`). Then launch:

```bash
./00-scripts/stacks_5b_rxstacks.sh
./00-scripts/stacks_6_cstacks_rx.sh
./00-scripts/stacks_7_sstacks_rx.sh
./00-scripts/stacks_8_populations_rx.sh
```

### With a reference genome

```bash
./00-scripts/bwa_mem_align_reads.sh
./00-scripts/stacks_1a_pstacks.sh
./00-scripts/stacks_2_cstacks.sh
./00-scripts/stacks_3_sstacks.sh
./00-scripts/stacks_4_populations.sh
```

**NOTE**: To align Illumina data, remove the `-O 0,0` and `-E 2,2` options from
the bwa script.

## Run the STACKS2 programs

### Without a reference genome

  - ustacks (params: )
  - cstacks (params: )
  - sstacks (params: )
  - tsv2bam (params: )
  - gstacks (params: )
  - populations (params: )

```bash
./00-scripts/stacks2_ustacks.sh
./00-scripts/stacks2_cstacks.sh
./00-scripts/stacks2_sstacks.sh
./00-scripts/stacks2_tsv2bam.sh
./00-scripts/stacks2_gstacks.sh
./00-scripts/stacks2_populations.sh
```

### With a reference genome

```bash
./00-scripts/stacks2_gstacks.sh
./00-scripts/stacks2_populations.sh
```
## Filtering the results

**NOTE**: All the filtering scripts that take a VCF for input or output can read
and write compressed VCF files. The files must be compressed with gzip and end
with the `.gz` extension. This is how the Python scripts recognize them. As a
result, it is recommended to compress your original VCF files from populations
with gzip as well as any further steps in order to save disk space.

### 1. Filter STACKS or STACKS2 VCF minimally and create graphs

#### Fast filter

This new filter script (2019-07-08) is recommended over the older, slower one.

- Less parameters
- Faster (5-10X depending on dataset)
- Uses only needed parameters
- Recommended for all analyses and much faster for big datasets

Here is the help from this script:

```bash
# Filtering SNPs in VCF file output by STACKS1 or STACKS2 minimaly
#
# Usage:
#     <program> input_vcf min_cov percent_genotypes max_pop_fail min_mas output_vcf
#
# Where:
#     input_vcf: is the name of the VCF file to filter
#     min_cov: minimum allele coverage to keep genotype <int>, eg: 4 or more
#     percent_genotypes: minimum percent of genotype data per population <float> eg: 50, 70, 80, 100
#     max_pop_fail: maximum number of populations that can fail percent_genotypes <int> eg: 1, 2, 3
#     min_mas: minimum number of samples with rare allele <int> eg: 2 or more
#     output_vcf: is the name of the filtered VCF
#
# WARNING:
#     The filtering is done purely on a SNP basis. Loci are not taken into account.

# Filtering (STACKS1)
./00-scripts/05_filter_vcf_fast.py 05-stacks/batch_1.vcf 4 70 0 2 filtered_m4_p70_S2

# Filtering (STACKS2)
./00-scripts/05_filter_vcf_fast.py 05-stacks/populations.snps.vcf 4 70 0 2 filtered_m4_p70_S2

# Graphs
./00-scripts/05_filter_vcf -i filtered_m4_p70_S2 -o graphs_filtered_m4_p70_S2 -g
```

**Note**: The last option filters on the **MAS**, which is akin to the MAF and
MAC. It keeps only SNPs where the rare allele has been found in *at least* a
certain number of samples. For example: `2` means that at least two samples
have the rare alleles. For RADseq data, the MAS is better than the MAF and MAC,
which are boosted by genotyping errors where one heterozygote sample is
genotyped as a rare-allele homozygote.

#### Slow filter

- More parameters but they are not needed with this new filtering procedure
- Slower (5-10X depending on dataset)
- Keeping only for backward compatibility

```bash
# Filtering (STACKS1)
./00-scripts/05_filter_vcf -i 05-stacks/batch_1.vcf -m 4 -p 70 --use_percent -S 2 -o filtered_m4_p70_S2

# Filtering (STACKS2)
./00-scripts/05_filter_vcf -i 05-stacks/populations.snps.vcf -m 4 -p 70 --use_percent -S 2 -o filtered_m4_p70_S2

# Graphs
./00-scripts/05_filter_vcf -i filtered_m4_p70_S2 -o graphs_filtered_m4_p70_S2 -g
```

**Note**: The `-S` option filters on the **MAS**, which is akin to the MAF and MAC.
It keeps only SNPs where the rare allele has been found in *at least* a certain
number of samples. For example: `-S 2` means that at least two samples have the
rare alleles. For Radseq data, the MAS is better than the MAF and MAC, which are
often boosted by genotyping errors where one heterozygote sample is genotyped as
a rare-allele homozygote.

### 2. Identify bad samples in lightly filtered VCF

#### 2.1. Too much missing data

  - Use data from `missing_data.png` and `missing_data.txt` from the graph step just above
  - Decide on a threshold and create a file with unwanted samples (one sample name per line)
  - Remove these bad samples from original populations VCF with `06_filter_samples_with_list.py`
    **BEFORE** you proceed to the next steps. Samples with a lot of missing data will create
    strange relatedness patterns.
  - Filter original populations VCF again with `05_filter_vcf_fast.py`

#### 2.2. Relatedness

  - Run `vcftools --relatedness --vcf <INPUT_VCF> --out samples` (use `--gzvcf` for comressed VCF files)
    to identify samples with potential errors / problems
  - Plot graph with `./00-scripts/utility_scripts/plot_relatedness_graphs.R samples.relatedness 0.5`
  - Decide on a threshold and create a file with unwanted samples (one sample name per line)

#### 2.3. Heterozygosity

  - Use `vcftools --het --vcf <INPUT_VCF> --out samples` (use `--gzvcf` for comressed VCF files)
  - Plot heterozygosity graph (see steps below)
  - Decide on a threshold and create a file with unwanted samples (one sample name per line)
  - Format data with:

```bash
awk '{print $5,$1,$1}' samples.het | cut -d "_" -f 1,2 > samples.het.data
```
  - Plot graph with `./00-scripts/utility_scripts/plot_heterozygozity.R samples.het.data`
  - Decide on a threshold and create a file with unwanted samples (one sample name per line)
  - Extract samples below that threshold with:

```bash
awk '$1 < -0.4 {print $2}' samples.het.data > samples.het.ids
```

#### 2.4. Remove bad samples

  - Create list of all unwanted samples from subsections  2.2, and 2.3 (one sample name per line)
  - Filter original populations VCF with `06_filter_samples_with_list.py`
  - This will create an unfiltered VCF where the bad samples are removed

### 3. If needed, make bigger groups of samples
#### TODO Permit the use of a population map with fast filter script to avoid this
  - If your dataset contains many small populations, regroup samples into fewer and bigger
    groups to avoid strict and overly stochastic filtering
  - STACKS1: Make a copy of `05-stacks/batch_1.vcf` or `06-stacks_rx/batch_1.vcf`
  - STACKS1: Make a copy of `05-stacks/populations.snps.vcf`
  - Modify sample names (`POP1_sample` -> `Group1_POP1-sample`. Note that the underscore `_` becomes a dash `-`

### 4. Filter new VCF

**NOTE**: You can launch the `05_filter_vcf_fast.py` without options to see documentation.

```bash
./00-scripts/05_filter_vcf_fast.py batch_1_grouped.vcf 4 70 0 2 filtered_bad_samples_removed_m4_p70_x0_S2
```

### 5. Explore SNP duplication using the following scripts

```bash
./00-scripts/08_extract_snp_duplication_info.py
./00-scripts/09_classify_snps.R
./00-scripts/10_split_vcf_in_categories.py

```

  - The following criteria are used by in `09_classify_snps.R`. Modify these in the script to fit your data.
    - Low Confidence: Extreme allele ratios (< 0.4 and > 0.6) with least one rare homozygote
    - Duplicated: Fis < -0.1
    - Duplicated: Fis + MedRatio / 3 < 0.11
    - Diverged: Fis < -0.6
    - Low Confidence: Fis > 0.6
    - High Coverage: MedCovHom > 40 or MedCovHet > 40
    - Minor Allele Sample (MAS): NumRare <= 2

### 6. Keep all unlinked SNPs
  - Using the singleton SNPs, keep only unlinked SNPs using:

```bash
00-scripts/11_extract_unlinked_snps.py
```

### 7. Onwards!

You should now have a very clean SNP dataset for your project. Analyze singletons,
duplicated, diverged, and high coverage SNPs separately.

  - Run population genomics analyses
  - Publish a paper!

## What to write in the Methods section of your paper

Here is a summary of things to write in the Methods section of your paper:

### Sample preparation

- Data prepared and SNPs VCF generated using [stacks_workflow](https://github.com/enormandeau/stack_workflow) <version> (eg: 2.4)
- Analyses done with STACKS <version> (eg: 1.48 or 2.4)
- Raw data cleaned with Cutadapt <version> (eg: 2.1)
- Samples extracted with `process_radtags`(part of STACKS)

### With a reference genome

- Cleaned and demultiplexed reads aligned to genome with bwa <version> (eg: 0.7.17-r1188)
  and samtools <version> (eg: 1.8)

- STACKS1 pipeline
  - pstacks (params: )
  - cstacks (params: )
  - sstacks (params: )
  - populations (params: )

- STACKS2 pipeline
  - gstacks (params: )
  - populations (params: )

### Withoug a reference genome (denovo)

- STACKS1 pipeline
  - ustacks (params: )
  - cstacks (params: )
  - sstacks (params: )
  - populations (params: )
  - rxstacks (params: )
  - cstacks (params: )
  - sstacks (params: )
  - populations (params: )

- STACKS2 pipeline
  - ustacks (params: )
  - cstacks (params: )
  - sstacks (params: )
  - tsv2bam (params: )
  - gstacks (params: )
  - populations (params: )

### Filtering

- STACKS VCF filtered a first time with `05_filter_vcf_fast.py` (params: 4 60 2 3)
- Create graphs to find samples with high missing data `05_filter_vcf.py` (params: -g)
- Decide missing data threshold and remove these samples with `06_filter_samples_with_list.py`
- Look for sample relatedness and heterozygosity problems in new VCF with vcftools
- Remove them with `06_filter_samples_with_list.py`
- If needed, regroups
- Filter this VCF with `05_filter_vcf_fast.py` (params: 4 60 0 3)
- Classify SNPs into singleton, duplicated, diverged, high coverage, low confidence, MAS with
  - `./00-scripts/08_extract_snp_duplication_info.py`
  - `./00-scripts/09_classify_snps.R`
  - `./00-scripts/10_split_vcf_in_categories.py`
- Keep only SNPS that are unlinked within loci with `11_extract_unlinked_snps.py`

## TODO

- Look for shared patterns of missing data caused by the sequencing
  - `plink --vcf <INPUT_VCF> --cluster missing --out <OUTPUT_VCF> --mds-plot 4 --allow-extra-chr`
  - Create figure using strata file to color samples

- Missing data imputation
  - Look for structure
  - Impute within each differentiable group
  - Assign genotypes randomly using observed genotype frequencies

### Running into problems

1. Consider joining the [STACKS Google
   group](https://groups.google.com/forum/#!forum/stacks-users)
1. [Biostar](https://www.biostars.org) is a useful bioinformatics forum.
1. [Stack Overflow (no link with STACKS)](https://stackoverflow.com/) is an essential programming forum.
