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

**NOTE**: Working on providing STACKS2 support (2019-05-29)

## Downloading

[click here to
download](https://github.com/enormandeau/stacks_workflow/archive/master.zip)
the current version of the stacks workflow. alternatively, you can clone this
repository with:

```
git clone https://github.com/enormandeau/stacks_workflow
```

## Licence

The stacks_workflow is licensed under the gpl3 license. See the licence file
for more details.

# Stacks workflow tutorial

The goal of this workflow is to simplify the use of the STACKS pipeline and to
create a folder architecture and code structure to help people analysing RADseq
data in a reproductible manner. It was developed with the needs of our research
group in mind. We make no claim about its usefulness to other groups or in other
contexts, but we still believe it may be of use to some.

## About STACKS

The [stacks analysis pipeline](http://creskolab.uoregon.edu/stacks/) is used
for restriction-site associated DNA sequencing (RADseq) studies, with and
without a reference genome.

Before starting to use STACKS, you should read the official STACKS papers
found at the bottom of the official STACKS page (see just above).

## Overview of the steps

1. Install stacks_workflow and STACKS with its dependencies
1. Download your raw data files (illumina lanes or ion proton chips)
1. Clean the reads and assess their quality
1. Extract individual data with process_radtags
1. Rename the sample files
1. Align reads to a reference genome (optional)
1. Stacks pipeline
1. Filtering the results

## Where to find this manual

To read a html version of this file online, go to the
[stacks_workflow project page](https://github.com/enormandeau/stacks_workflow)
on GitHub.

## Installing stacks_workflow

### Download and install the most recent version of this workflow

It is recommended to do this **for each new project**, as opposed to re-using
the same directory for multiple projects and naming the outputs differently.
This is **central to `stack_workflow`'s phylosophy** of reproducibility.

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

Follow the instructions on the website.

```bash
cstacks
which cstacks
```

This will output the help of the cstacks program and where it is located
on your computer. You will also be able to confirm the version number of
your STACKS installation.

#### Installing Cutadapt

There are different ways you can install Cutadapt. If you have `pip` (a
Python package installer) installed, you can use the following command:

```bash
sudo pip install --user --upgrade cutadapt
```

Otherwise, [visit their website to download it and install
it](https://pypi.python.org/pypi/cutadapt/)

#### Installing FastQC

We use FastQC to assess read quality before and after filtering with Cutadapt.
To install FastQC, visit [this
page](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
and download the version that works with your system (Linux or MacOS). Then
launch these commands:

##### Installing on Linux

```bash
# Installing
unzip fastqc_v0.11.3_source.zip
cd FastQC
chmod 755 fastqc

# Look at the options
./fastqc --help

# Run FastQC
./fastqc
```

## Prepare your raw datafiles

### Downloading your data

Download your raw Illumina or Ion Proton data files from your sequencing
service provider.

### Copy your raw files

Put a copy of (or a link to) your raw data files in the `02-raw` folder of
stacks_workflow.

All file names **must** end with **.fastq.gz** for the following scripts to
work.

### Preparing the `lane_info.txt` file

This file will contains the names of the raw data files and is used by
stacks_workflow later.  From the stacks_workflow folder, run:

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
not put a number, the script will use only one CPU (for backward
compatibility).

### Scan the cutadapt logs

The cutadapt log files can be found in the `10-log_files` folder. Scan them
to confirm that cutadapt has done an appropriate job.

There may be difference in adapters and filter parameters to use with data
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

**Notes:**

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
for example mixing lower and capital letters in population or sample names
(e.g.: Pop01 and pop01), since these will be treated as two different
populations.

### Launch process_radtags with:

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

#### Testing different trim lengths

If you are using Ion Proton data, the effect of the trimLength parameter used
above on the number of usable SNPs you recover at the end may not be trivial.
As a rule of thumb, a trimmed length of 80bp should produce good results. We
suggest you run tests with a smaller group of samples to determine what length
to trim to. For highly variant species, short loci will be more likely to
contain SNPs and long loci to contain more than one SNP, which is not always
informative. Thus, trimming to shorter lengths may be more interesting for
highly variant species or when coverage is limiting.

### Rename samples

#### Rename and copy the samples

We provide a scripts to rename the extracted samples and move them into the
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

If after splitting your samples you notice that some of have too few reads, you
can remove these from the `04-all_samples` folder. The threshold for the
minimum number of reads will depend on your project, including on the number of
expected cut sites generated by your library preparation protocol. Keep samples
with low coverages if you are not sure at this point what threshold to use. We
will filter the output later for this and will then have better informations.

#### Align reads to a reference genome (optional)

#### Install `bwa <http://bio-bwa.sourceforge.net>`_

#### Download reference genome to the `01-info_files`

#### Index the reference genome

```bash
bwa index -p genome -a bwtsw ./01-info_files/<genome reference>
mv genome.* 01-info_files
```

### Align samples

```bash
./00-scripts/bwa_mem_align_reads.sh
```

## STACKS pipeline

### Prepare population info file

```bash
./00-scripts/04_prepare_population_map.sh
```

## Edit script parameters

You will need to go through the scripts named `stacks_*` in the `00-scripts
folder` and edit the options to suite your needs. Depending on your project
(eg: *de novo* vs reference), you will not use all the scripts.

**Warning!** This step is very important. Choosing appropriate parameters for
your study is crucial in order to generate meaninful and optimal results. Read
the STACKS documentation on their website to learn more about the different
options.

## Run the STACKS programs

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

**Warning!** The documentation and scripts used with a reference genome have
not been updated in a long time. We believe they should not be used at the
moment.


```bash
./00-scripts/bwa_mem_align_reads.sh
./00-scripts/stacks_1a_pstacks.sh
./00-scripts/stacks_2_cstacks.sh
./00-scripts/stacks_3_sstacks.sh
./00-scripts/stacks_4_populations.sh
```

## Filtering the results

1. Filter STACKS VCF a minimum and create graphs:
```
./00-scripts/05_filter_vcf -i 05-stacks/batch_1.fcf -m 4 -p 70 --use_percent -S 2 -o filtered_m4_p70_S2
./00-scripts/05_filter_vcf -i filtered_m4_p70_S2 -o graphs_filtered_m4_p70_S2 -g
```

2. Identify bad samples
  - Identify samples with too much missing data from figure and `missing_data.txt`
  - Run `vcftools --relatedness` and identify potential errors / problems
  - Remove the files from these samples from `05-stacks` or `06-stacks_rx` (put them in a subfolder)
  - Re-run `population`

3. If needed, make bigger groups of samples
  - If many small populations, group samples into fewer groups to avoid strict and stochastic filtering
    - Copy `05-stacks/batch_1.vcf` (or `06-stacks_rx/batch_1.vcf`)
    - Modify names of samples (`POP1_sample` -> `Group1_POP1-sample`)

4. Filter new VCF
```
./00-scripts/05_filter_vcf -i batch_1_grouped.vcf -m 4 -p 70 --use_percent -S 2 -o filtered_bad_samples_removed_m4_p70_S2
```

5. Explore SNP duplication using the following scripts
```
./00-scripts/08_extract_snp_duplication_info.py
./00-scripts/09_classify_snps.R
./00-scripts/10_split_vcf_in_categories.py
```

  - The following criteria are used by in `09_classify_snps.R`. Modify these for your data.
    - Low Confidence: Extreme allele ratios (< 0.4 and > 0.6) with least one rare homozygote
    - Duplicated: Fis < -0.1
    - Duplicated: Fis + MedRatio / 3 < 0.11
    - Diverged: Fis < -0.6
    - Low Confidence: Fis > 0.49
    - High Coverage: MedCovHom > 40 or MedCovHet > 40
    - Minor Allele Sample: NumRare <= 2

6. Keep all unlinked SNPs
  - Using the singleton SNPs, keep only unlinked SNPs using
```
00-scripts/utility_scripts/extract_unlinked_snps.py
```

7. Onwards!

You should now have a very clean SNP dataset for your project. Analyze singletons,
duplicated, diverged, and high coverage SNPs separately.

  - Run population genomics analyses
  - Publish a paper!

### Running into problems

1. Consider joining the [Stacks Google
   group](https://groups.google.com/forum/#!forum/stacks-users)
1. [Biostar](https://www.biostars.org) is a useful bioinformatics forum.
1. [SEQanswers](http://seqanswers.com) is another useful forum for all things
   related to Next Generation Sequencing.
1. Of course, you can always ask [Google](https://www.google.com) for help.
