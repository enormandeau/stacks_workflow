# STACKS Workflow2 

RADseq workflow using [STACKS2](http://creskolab.uoregon.edu/stacks/)

Developed by [Eric Normandeau](https://github.com/enormandeau) in
[Louis Bernatchez](http://www.bio.ulaval.ca/louisbernatchez/presentation.htm)'s
laboratory.

**NOTE!**: stacks_workflow no longer supports STACKS1. For the latest version
of stack_workflow that does, find `v2.5.2_last_version_supporting_STACKS1` in
the releases on the stack_workflow GitHub page.

**Warning!**: this software is provided "as is", without warranty of any kind,
express or implied, including but not limited to the warranties of
merchantability, fitness for a particular purpose and noninfringement. In no
event shall the authors or copyright holders be liable for any claim, damages
or other liability, whether in an action of contract, tort or otherwise,
arising from, out of or in connection with the software or the use or other
dealings in the software.

stack_workflow was developed with the needs of our research group in mind. We
make no claim about its usefulness to other groups or in other contexts, but it
has been and continues to be useful to other groups.

## Licence

stacks_workflow is licensed under the gpl3 license. See the LICENCE file
provided with stacks_workflow for more details.

## About STACKS

The [stacks analysis pipeline](http://creskolab.uoregon.edu/stacks/) is used
for restriction-site associated DNA sequencing (RADseq) studies, with and
without a reference genome.

Before starting to use STACKS, you should read the official STACKS papers
found at the bottom of the official STACKS page (see link above) as well
as the official STACKS documentation found at the link above.

# STACKS workflow tutorial

The goal of this workflow is to simplify the use of the STACKS pipeline and
make the analyses more reproducible. One of the major contributions is the
standardized SNP filtering procedures used to produce high quality SNP
datasets.

## Overview of the steps

1. Install stacks_workflow and STACKS2
1. Download your raw data files (illumina lanes or ion proton chips)
1. Clean the reads and assess their quality
1. Extract individual data with process_radtags
1. Rename the sample files
1. Align reads to a reference genome (optional)
1. Use the STACKS pipeline
1. Filter the results
1. Impute missing data if needed

## Installing stacks_workflow

### Download and install the most recent version of this workflow

It is recommended to download the most recent version of stacks_workflow **for
each new project**, as opposed to re-using the same directory for multiple
projects and naming the outputs differently. This is **central to
`stack_workflow`'s philosophy** of reproducibility.

**One stacks_workflow folder should contain only one analysis**

Deviate from this at your own risk ;)

[click here to
download](https://github.com/enormandeau/stacks_workflow/archive/master.zip)

#### Clone stacks_workflow with git

If `git` is installed on your computer, you can run the following command
instead to get a complete stacks_workflow git repository.

```bash
git clone https://github.com/enormandeau/stacks_workflow
```

For the rest of the project, use the extracted stacks_workflow archive or
cloned folder as your working directory. It may be a good idea to rename it by
adding information about your project, the date or any other useful
information.

**All the commands in this manual are launched from that directory.**

### Download and install [STACKS](http://creskolab.uoregon.edu/stacks/)

Follow the instructions on the STACKS website to install and test the
installation with:

```bash
populations --version
which populations
```

This will output the version of the populations program (part of STACKS) and
where it is located on your computer.

### Dependencies

- STACKS (latest 2.x version is recommended)
- Linux or MacOS
- gnu parallel
- cutadapt (install with conda or pip)
- Python3 and packages (best managed with conda):
  - matplotlib
  - numpy
  - pandas
  - PIL, xlrd, and xlutil (optional: for inter-chip normalization)
- admixture (optional: for missing data imputation)
- plink (optional: for missing data imputation and exploration)
- R
  - adegenet (optional: for admixture plots)
- Imagemagick (optional: to join admixture plots)

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

The cutadapt log files can be found in the `10-log_files` folder. Scan them or
look at the file sizes to confirm that cutadapt has done an appropriate job.

There may be differences in adapters and filter parameters to use with data
produced by Illumina and Ion Proton sequencers.

## Extract individual data with `process_radtags`

### Prepare a file named `sample_information.csv`

Use the same format found in the `example_sample_information.csv` file located
in the `01-info_files` folder.

Save this file in the `01-info_files` folder and name it exactly
`sample_information.csv`. This file will be used to extract the samples and
rename the extracted sample files automatically.

The first column **MUST** contain the **EXACT** name of the data file for the
lane/chip of each sample.

**Notes**:

- The **columns are separated by tabulations** (even if the extension is .csv)
- The second column contains the barcode sequence of each sample.
- The third column contains the population name of each sample.
- The fourth column contains the name of the sample (do not include the
  population name or abbreviation in the sample name).
- Neither the population name nor the sample name should contain underscores `_`
- The fifth column contains a number or string identifying the populations. you
  can use the same as in the the the third column.
- The sixth column contains the plate well identifier.

Columns three, four, and five are treated as text, so they can contain either
text or numbers. Other columns can be present after the fifth one and will be
ignored. However, it is crucial that the six first columns respect the format
in the example file exactly. Be especially careful not to include errors in
this file, for example by mixing lower and capital letters in population or
sample names (e.g.: Pop01 and pop01), since these will be treated as two
different populations.

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
As a rule of thumb, a trimmed length of 80bp should produce good results in
most projects. We suggest you run tests with a smaller group of samples to
determine what length to trim to. For highly species with high genetic
variability, short loci will be more likely to contain SNPs and long loci to
contain more than one SNP, which is not always informative. Thus, trimming to
shorter lengths may be more interesting for highly variant species or when
coverage is limiting. On the other hand, trimming to keep longer sequences
(for example 120pb) can be more interesting if the read coverage is very
good and the genetic varability is low.

### Rename samples

#### Rename and copy the samples

We provide a script to rename the extracted samples and move them into the
`04-all_samples` folder. The script behaves differently for samples that are
present only once in the `01-info_files/sample_information.csv` and for those
that are present more than once. If a sample is present only once, a link is
created, using no additional disk space. If it is present more than once, it
means that this sample has been sequenced on multiple lanes/chips and all the
copies are concatenated into one file, doubling the amount of disk space taken
by this sample (all the individual files PLUS the combined one).

```bash
./00-scripts/03_rename_samples.sh
```

#### Assessing the quality of your reads

After this step, you will want to run FastQC on the read sequences found in
`04-all_samples`. A nice way of visualizing them is to use `multiqc` to create
a unique report for all the reads. Pay special attention to the duplication
level. You probably want to have high duplication in the 10-50X range, but if a
high proportion of your data is in the 100X+ range, then maybe your library
suffers from lower complexity than is ideal. This is up to you to judge given
what you know of your species (genome), enzyme(s) used, and sequencing
coverage.

#### Deleting samples with too few reads

If after splitting your samples you notice that some have too few reads, you
can remove these from the `04-all_samples` folder. The threshold for the
minimum number of reads will depend on your project, including the number of
expected cut sites generated by your library preparation protocol and the
number of reads per sample. Keep samples with low coverages if you are not sure
what threshold to use at this point. We will filter the VCF for this later and
will then have better information then.

#### Align reads to a reference genome (optional)

##### Install `bwa <http://bio-bwa.sourceforge.net>`

##### Download reference genome to the `08-genome`

Decompress the genome if needed and make a copy of it named `genome.fasta`.

##### Index the reference genome

```bash
bwa index ./08-genome/genome.fasta
```

### Align samples

Different bwa alignment scripts are available in 00-scripts.

**IMPORTANT NOTE**: The two scripts for single-end reads (ie: not the one with
`PE` in its name) have options that are specific for IonProton data. To align
Illumina data, remove the `-O 0,0` and `-E 2,2` options.

```bash
00-scripts/bwa_mem_align_reads.sh
00-scripts/bwa_mem_align_reads_by_n_samples.sh
00-scripts/bwa_mem_align_reads_PE.sh
```

**NOTE**: To align Illumina data, remove the `-O 0,0` and `-E 2,2` options from
the bwa script.

## STACKS pipeline

### Prepare population info file

```bash
./00-scripts/04_prepare_population_map.sh
```

## Edit script parameters

You will need to go through the scripts named `stacks2_*` in the
`00-scripts folder` and edit the options to suit your needs. Depending on your
project (eg: *de novo* vs reference), you will not use all the scripts.

**Warning!** This step is very important. Choosing appropriate parameters for
your study is crucial in order to generate meaningful and optimal results. Read
the STACKS documentation on their website to learn more about the different
options.

## Run the STACKS2 programs

### Without a reference genome

  - ustacks
  - cstacks
  - sstacks
  - tsv2bam
  - gstacks
  - populations

```bash
./00-scripts/stacks2_ustacks.sh
./00-scripts/stacks2_cstacks.sh
./00-scripts/stacks2_sstacks.sh
./00-scripts/stacks2_tsv2bam.sh
./00-scripts/stacks2_gstacks.sh
./00-scripts/stacks2_populations.sh
```

### With a reference genome

After the reads are aligned with bwa, run:

```bash
./00-scripts/stacks2_gstacks.sh
./00-scripts/stacks2_populations.sh
```
## Filtering the results

**NOTE**: All the filtering scripts that take a VCF for input or output can
read and write compressed VCF files. The files must be compressed with gzip and
end with the `.gz` extension.  This is how the Python scripts recognize them.
As a result, it is recommended to compress your original VCF files from
populations with gzip as well as any further steps in order to save disk space,
especially for big projects.

### 1. Filter the VCF minimally and create graphs

#### Fast filter

This new filter script (2019-07-08) is recommended instead of the older, slower
one.

Reasons to use the faster filter script:

- Less parameters
- Uses only needed parameters
- Faster (5-10X depending on dataset)
- Recommended for all analyses and much faster for big datasets

Here is the documentation from this script:

```bash
# Filtering SNPs in VCF file output by STACKS1 or STACKS2 minimaly
#
# Usage:
#     <program> input_vcf min_cov percent_genotypes max_pop_fail min_mas output_vcf
#
# Where:
#     input_vcf: is the name of the VCF file to filter
#     min_cov: minimum allele coverage to keep genotype <int>, eg: 4 or more
#     percent_genotypes: minimum percent of genotype data per population <float> eg: 50, 70, 80, 90, 100
#     max_pop_fail: maximum number of populations that can fail percent_genotypes <int> eg: 1, 2, 3
#     min_mas: minimum number of samples with rare allele <int> eg: 2 or more
#     output_vcf: is the name of the filtered VCF
#
# WARNING:
#     The filtering is done purely on a SNP basis. Loci are not taken into account.

# Filtering (STACKS1)
./00-scripts/05_filter_vcf_fast.py 05-stacks/batch_1.vcf 4 70 0 2 filtered_m4_p70_x0_S2.vcf

# Filtering (STACKS2)
./00-scripts/05_filter_vcf_fast.py 05-stacks/populations.snps.vcf 4 70 0 2 filtered_m4_p70_x0_S2.vcf

# Graphs
./00-scripts/05_filter_vcf.py -i filtered_m4_p70_x0_S2 -o graphs_filtered_m4_p70_x0_S2 -g
```

**Note**: The last option filters on the **MAS**, which is akin to the MAF and
MAC. It keeps only SNPs where the rare allele has been found in *at least* a
certain number of samples. For example: `2` means that at least two samples
have the rare alleles. For RADseq data, the MAS is better than the MAF and MAC,
which are artificially boosted by genotyping errors where one heterozygote
sample is falsely genotyped as a rare-allele homozygote. Given the nature of RADseq,
these errors are quite frequent.

#### Slow filter

- More parameters but they are not needed with this new filtering procedure.
  They are a relic of an "early era" in the exploration of quality filtering.
- Slower (5-10X depending on dataset)
- Keeping only for backward compatibility and to generate descriptive graphs

```bash
# Filtering (STACKS1)
./00-scripts/05_filter_vcf.py -i 05-stacks/batch_1.vcf -m 4 -p 70 --use_percent -S 2 -o filtered_m4_p70_x0_S2

# Filtering (STACKS2)
./00-scripts/05_filter_vcf.py -i 05-stacks/populations.snps.vcf -m 4 -p 70 --use_percent -S 2 -o filtered_m4_p70_x0_S2

# Graphs
./00-scripts/05_filter_vcf.py -i filtered_m4_p70_x0_S2 -o graphs_filtered_m4_p70_x0_S2 -g
```

**Note**: The `-S` option filters on the **MAS**, which is akin to the MAF and
MAC. It keeps only SNPs where the rare allele has been found in *at least* a
certain number of samples. For example: `2` means that at least two samples
have the rare alleles. For RADseq data, the MAS is better than the MAF and MAC,
which are artificially boosted by genotyping errors, where one heterozygote
sample is genotyped as a rare-allele homozygote. Given the nature of RADseq,
these errors are quite frequent.

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
  - If your dataset contains many small populations, regroup samples into fewer
    and bigger groups to avoid strict and overly stochastic filtering
  - Make a copy of `05-stacks/populations.snps.vcf`
  - Modify sample names (eg: `POP1_sample` -> `Group1_POP1-sample`). Note that the
    underscore `_` becomes a dash `-`.
  - Use bcftools to do that:
    - `bcftools reheader -s names.txt input.vcf > renamed.vcf`
    - The `names.txt` file contains current sample names in the first column
and desired sample names in a second column.
    - The columns are separated by a tabulation.

### 4. Filter new VCF

**NOTE**: You can launch the `05_filter_vcf_fast.py` without options to see documentation.

```bash
./00-scripts/05_filter_vcf_fast.py populations.snps.grouped.vcf 4 70 0 2 filtered_bad_samples_removed_m4_p70_x0_S2
```

### 5. Explore SNP duplication using the following scripts

```bash
./00-scripts/08_extract_snp_duplication_info.py
./00-scripts/09_classify_snps.R
./00-scripts/10_split_vcf_in_categories.py

```

  - The following criteria are used by in `09_classify_snps.R`. Modify these in
      the script to fit your data.
    - Low Confidence: Extreme allele ratios (< 0.4 and > 0.6) with least one
    rare homozygote
    - Duplicated: Fis < -0.1
    - Duplicated: Fis + MedRatio / 3 < 0.11
    - Diverged: Fis < -0.6
    - Low Confidence: Fis > 0.6
    - High Coverage: MedCovHom > 40 or MedCovHet > 40
    - Minor Allele Sample (MAS): NumRare <= 2

### 6. Keep all unlinked SNPs

It is often thought that SNPs appearing within the same STACKS locus are 100%
linked because they are really close. However, this is often not the case.
Frequently, you will find SNPs that are not linked within the same locus. In
order to filter and keep as much genetic information as possible, while
avoiding close by SNPs with high Linkage Disequilibrium, you can keep all the
SNPs that we refer to as unlinked in all the loci.

The procedure is as follows:
  - Keep the first SNP and remove all the other ones appear linked to it
  - If you have SNPs remaining, repeat

Two SNPs are linked when sample genotypes are highly correlated for these two
SNPs. Since RADseq data has 1) missing data and 2) mostly SNPs with low MAF
values, we need to be careful when comparing sample genotypes between two SNPs.
As a result, when comparing two SNPs, we only use samples that have no missing
data in both SNPs and who possess the rare allele in at least one of the SNPs.

Using the singleton SNPs, keep only unlinked SNPs using one the following
scripts. The input parameters are described by the scripts themselves.

```bash
# Denovo
00-scripts/11_extract_unlinked_snps.py

# Reference
00-scripts/11_extract_unlinked_snps_genome.py
```

### 7. Missing data imputation

Impute missing data in a VCF using Admixture ancestry relationships

/!\ **WARNING** /!\  

Whatever the method of choice, missing data imputation cannot impute **CORRECT
GENOTYPES**, only **GENOTYPES THAT MINIMIZE BIASES** in a dataset. You should
use imputation **ONLY** when you really need it. For example when some piece
of software will not accept missing data in its input VCF.

#### Limitations of the ancestry-based missing data imputation

- Light: Admixture is slow with big datasets. You can thin down your SNP dataset
  if this becomes problematic (see admixture manual).
- Light: Using all the SNPs versus using only neutral SNPs with admixture can
  change the ancestry estimation of samples. For example, the CV could vary
  differently as a function of K.
- Light: Even using cross-validation in Admixture (CV values), the best K value is
  chosen by the user and so the groups and ancestry will vary. This will have
  an impact on the imputation but the approach should be fairly robust around K
  values that make biological sense.
- **Moderate**: Admixture requires that the individuals be unrelated. Some level
  of half-sibs or full sibs is probably OK, but watch out for datasets with a lot
  of related samples. You can use the relatedness part of the filtration steps
  listed above to check that.
- **Moderate**: Identity by missing data, where patterns of similarity among
  samples is the result of non-random missing data within groups of samples, is
  problematic for admixture. You need to assert that this pattern is not present
  in your dataset (using plink) or remove the loci succeptible to this from
  your VCF before using vcf_impute. See details in the procedure below.
- **IMPORTANT**: Admixture is a poor choice for samples with a continuous
  genetic gradient, a pattern of isolation by distance or a dataset with a lot
  of populations with very low or unequal sample numebrs.  Using a k-nearest
  neighbors approach may be better in this case.
- **VERY IMPORTANT**: Large genomic features, such as big inversions, can create
  strong groupings in admixture but that group structure would only apply to
  local parts of the genome, or even none at all for complex cases. If you feel
  like different parts of your genomes could lead to a very different
  hierarchical population structure, using a k-nearest approach may be better.

#### Advantages of the ancestry-based missing data imputation

- Major: Avoid using overfitted models that depend on information from other
  loci to impute genotypes in the current locus. It is our belief that, in most
  RADseq studies, apparent correlation among loci exists because of stochastic
  rather than biological reasons. For that reason, using information from loci
  that seem correlated is not a good choice to infer missing genotypes. This
  is because the genotypes at these other pseudo-correlated loci have a low
  probability to be informative for the imputation of the missing genotype.

#### Running the imputation

1. Format contig/scaffold names

In order to use admixture, contig/scaffold names (refered to as chromosomes in
admixture) must be integers. We use the following script to correct this. Make
sure the output vcf is EXACTLY named `input_renamed.vcf`. The input VCF can be
compressed with gzip.

```bash
./00-scripts/12_rename_vcf_scaffolds_for_plink.py <input.vcf> input_renamed.vcf
```

And check for patterns of identify by missing and potentialy filter the VCF to
remove the SNPs responsible of any such pattern (not covered in this document).

```bash
./00-scripts/utility_scripts/plink_cluster_missing.sh
./00-scripts/utility_scripts/plink_cluster_missing_figure.R input_renamed.mds
```

2. Use plink to create bed file

```bash
plink --vcf input_renamed.vcf --make-bed --out input_renamed --allow-extra-chr
```

3. Use admixture and find a good K value

```bash
# Run admixture
# Adjust the `seq 10` value for your dataset. This number is the highest number
# of groups (K) that will be tried with admixture
seq 10 | parallel admixture input_renamed.bed {} -j4 --cv -C 0.1 \> 11-admixture/input_renamed.{}.log
mv *.P *.Q 11-admixture/

# Explore CV values and choose an appropriate K value
grep -h CV 11-admixture/*.log | sort -V  # May not work on MacOs or BSD descendents because of the -V option
grep -h CV 11-admixture/*.log | cut -d " " -f 4,3 | awk '{print $2,$1}' | sort -n

# Look at (crude) graphs of group memberships to assist in choosing the K value
# (Thanks to Nicolas Leroux for the original plot R script!)
# The .png files will be found in the 11-admixture folder
# Requires the adegenet package
parallel ./00-scripts/utility_scripts/plot_admixture.R ::: 11-admixture/*.Q

# If you have imagemagick installed, you can combine all the graphs in one to
# help choose the best K value
convert $(ls -1 11-admixture/input_renamed.*.png | sort -V) -trim -border 0x4 -gravity center -append all_admixture_figures.png
```

4. Impute missing genotypes using sample related groups

```bash
# Replace K by the desired number of clusters
./00-scripts/13_impute_missing.py input_vcf 11-admixture/input_renamed.K.Q output_vcf
```
### 8. Onwards!

You should now have a very clean SNP dataset for your project. Analyze only
singletons or analyse the different categories of SNPs separately.

  - Run population genomics analyses
  - Publish a paper!

## For the Methods section of your paper

Here is a summary of informations that should go in the Methods section of your paper.

### Sample preparation

- Data prepared, SNPs VCF generated and filtered using:
  - [STACKS](http://catchenlab.life.illinois.edu/stacks/) <version> (eg: 1.48 or 2.54)
  - [stacks_workflow](https://github.com/enormandeau/stack_workflow) <version> (eg: 2.5.4)
- Raw data cleaned with Cutadapt <version> (eg: 2.1)
- Samples extracted with `process_radtags` (part of STACKS)

### Denovo

- **STACKS2 pipeline (Denovo)**
  - ustacks (ex. params: -m 4, -M 3, -N 5)
  - cstacks (ex. params: -n 3)
  - sstacks (ex. params: na)
  - tsv2bam (ex. params: na)
  - gstacks (ex. params: --max-clipped 0.1)
  - populations (ex. params: -p 2, -r 0.6, --fasta-loci, --vcf)

### Reference genome

- Cleaned and demultiplexed reads aligned to genome with:
  - bwa <version> (eg: 0.7.17-r1188)
  - samtools <version> (eg: 1.8)

- **STACKS2 pipeline (Reference genome)**
  - gstacks (ex. params: --max-clipped 0.1)
  - populations (ex. params: -p 2, -r 0.6, --ordered-export, --fasta-loci, --vcf)

### Filtering

- STACKS VCF filtered a first time with `05_filter_vcf_fast.py` (ex. params: 4 60 2 3)
- Create graphs to find samples with high missing data `05_filter_vcf.py` (ex. params: -g)
- Decide missing data threshold and remove these samples with `06_filter_samples_with_list.py`
- Look for sample relatedness and heterozygosity problems in new VCF with vcftools
- Remove them with `06_filter_samples_with_list.py`
- If needed, regroup populations into larger groups to prevent spurious filtering
- Filter this new VCF with `05_filter_vcf_fast.py` (ex. params: 4 60 0 3)
- Classify SNPs into singleton, duplicated, diverged, high coverage, low confidence, MAS with
  - `./00-scripts/08_extract_snp_duplication_info.py`
  - `./00-scripts/09_classify_snps.R`
  - `./00-scripts/10_split_vcf_in_categories.py`
- Keep only SNPS that are unlinked within loci with `11_extract_unlinked_snps.py`
- Impute missing data with admixture

### Running into problems

1. Consider joining the [STACKS Google
   group](https://groups.google.com/forum/#!forum/stacks-users)
1. [Biostar](https://www.biostars.org) is a useful bioinformatics forum.
1. [Stack Overflow (no link with STACKS)](https://stackoverflow.com/) is an essential programming forum.
