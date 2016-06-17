# Stacks Workflow

An integrated workflow around the [STACKS software
pipeline](http://creskolab.uoregon.edu/stacks/) to help streamline RAD and GBS
analyses

stacks workflow is developped in [louis bernatchez'
lab](http://www.bio.ulaval.ca/louisbernatchez/presentation.htm).

**warning** 

the software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. in no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the
software.

## downloading

[click here to
download](https://github.com/enormandeau/stacks_workflow/archive/master.zip)
the current version of the stacks workflow. alternatively, you can clone this
repository with:

```
git clone https://github.com/enormandeau/stacks_workflow
```

## licence

the stacks_workflow is licensed under the gpl3 license. see the licence file
for more details.

# stacks workflow tutorial

the goal of this workflow is to simplify the use of the stacks pipeline and to
create a folder architecture and code structure to help people analysing gbs or
rad data, especially those with limited unix/linux experience, to jump on the
analysis wagon faster. 

it was developed with the needs of our research group in mind as well as with
an emphasis on non-model species studies. we make no claim about its usefulness
to other groups or in other contexts, but we still believe it may be of use to
some.

## about stacks
 
the [stacks analysis pipeline](http://creskolab.uoregon.edu/stacks/) is used
for snp discovery in genotyping by sequencing (gbs) and restriction-site
associated dna sequencing (rad) studies, with and without a reference genome.

upon starting to use stacks, it is highly suggested to read the two official
stacks papers:

[catchen, j. m., amores, a., hohenlohe, p. a., cresko, w. a., postlethwait, j.
h., & de koning, d. j. (2011). stacks: building and genotyping loci de novo
from short-read sequences. g3, 1(3), 171–182.
doi:10.1534/g3.111.000240](http://www.g3journal.org/content/1/3/171.full)

[catchen, j. m., hohenlohe, p. a., bassham, s., amores, a., & cresko, w. a.
(2013). stacks: an analysis tool set for population genomics. molecular
ecology, 22(11), 3124–3140.
doi:10.1111/mec.12354](http://onlinelibrary.wiley.com/doi/10.1111/mec.12354/abstract)

also very useful paper to read before attempting to run stacks:

[mastretta yanes a, arrigo n, alvarez n et al. (2014) rad sequencing,
genotyping error estimation and de novo assembly optimization for population
genetic inference. molecular ecology resources,
n/a–n/a.](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12291/abstract;jsessionid=a32722e1462a2a2714ee53a6fd4c7194.f04t04)
 
## overview of the steps

1. install stacks_workflow and stacks with its dependencies
1. download your raw data files (illumina lanes or ion proton chips)
1. clean the reads and assess their quality
1. extract individual data with process_radtags
1. rename the sample files
1. align reads to a reference genome (optional)
1. stacks pipeline
1. filtering the results
1. further analyses

## where to find this manual

to read a html version of this file, simply go to the github home page
for [stacks_workflow](https://github.com/enormandeau/stacks_workflow)

## installing stacks_workflow

### download and install the most recent version of this workflow

#### from the terminal

```bash
cd ~/desktop
wget https://github.com/enormandeau/stacks_workflow/archive/master.zip
unzip master.zip
```

#### using git

if `git` is installed on your computer, you can run the following command
instead to get a complete `stacks_workflow` git repository.

```bash
git clone https://github.com/enormandeau/stacks_workflow
```

for the rest of the project, use the extracted or cloned folder as your working
directory. **all the commands in this manual are launched from that
directory.**

### set up for mac osx

if you have a macos computer, you can prepare your computer to run gbs
analyses. thierry gosselin prepared a companion tutorial about gbs for cloud
computing and you should read the section about [setting up your macos
computer](http://gbs-cloud-tutorial.readthedocs.org/en/latest/03_computer_setup.html#mac-osx>)

### download and install [stacks](http://creskolab.uoregon.edu/stacks/)

#### installing google's sparsehash

google sparsehash is a hash table implementation with a small memory footprint.
stacks can (and should) be compiled using it. this will make stacks use much
less memory.

```bash
wget http://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz

tar -xvf sparsehash-2.0.2.tar.gz
cd sparsehash-2.0.2
./configure
make  # add '-j n' to use n cpus during the compilation
sudo make install
```

#### Installing stacks

```bash
# Modify the version number as needed
wget http://creskolab.uoregon.edu/stacks/source/stacks-1.29.tar.gz

tar -xvf stacks-1.29.tar.gz
cd stacks-1.29

# Install the binaries in /usr/local/bin
./configure --enable-sparsehash
make  # Add '-j n' to use n CPUs during the compilation
sudo make install

# Remove the temporary install folders and archives
cd ..
sudo rm -R stacks-1.29 stacks-1.29.tar.gz sparsehash-2.0.2 sparsehash-2.0.2.tar.gz
```
 
#### Test the STACKS installation

```bash
cstacks
```

This will output the help of the cstacks program. You will also be able to
confirm the version number of your STACKS installation.

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

##### Installing on MacOS

You're on your own, cowboy ;)

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

### Removing adapters

We use Cutadapt in order to remove the adapters present in the raw data. This
step is done before we extract the reads for each individual.

In order to clean your reads properly, you must determine which adapters are
present in your data and whether a cut site can also be found with it.

Before running cutadapt, we must create a fasta file containing the adapters
that are present in our data. This file must be found in the `01-info_files/`
folder and be names `adapters.fasta`. The file named `adapters.fasta`, found in
the `01-info_files/` folder should be used as a template.

### Running Cutadapt

Now that the `01-info_files/adapters.fasta` file contains all the adapters we
wish to remove from our data, we can launch Cutadapt **in single-end mode**
with the following command:

```bash
./00-scripts/01_cutadapt.sh numCPUs
```

Where `numCPUs` is the number of CPUs you wish to use in parallel. If you do
not put a number, the script will use only one CPU (for backward
compatibility).

### Scan the cutadapt logs

The cutadapt log files can be found in the `98-log_files` folder. Scan them
to confirm that cutadapt has done an appropriate job.

There may be difference in adapters and filter parameters to use with data
produced by Illumina and Ion Proton sequencers.

### Assess sequence quality with FastQC

**TODO:**

- Before and after Cutadapt
- Use the script
  - Command line interface
  - Run only on first 1,000,000 reads or less
  - Gathers important figures together

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
- The fifth column contains a number identifying the populations. 

Columns three and four are treated as text, so they can contain either text or
numbers. Other columns can be present after the fifth one and will be ignored.
However, it is crucial that the five first columns respect the format in the
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
We suggest you run tests with a smaller group of samples to determine what
length to trim to. For highly variant species, short loci will be more likely
to contain SNPs and long loci to contain more than one SNP, which is not always
informative.  Thus, trimming to shorter lengths may be more interesting for
highly variant species.

### Rename samples

#### Rename and copy the samples

We provide a scripts to rename the extracted samples and move them into the
`04-all_samples` folder. The script behaves differently for samples that are
present only once in the `01-info_files/sample_information.csv` and for those
that are present more than once. If a sample is present only once, a link is
created, thus using no additional disk space. If it is present more than once,
all the copies are concatenated into one file, doubling the amount of disk
space taken by this sample (all the individual files PLUS the combined one).

```bash
./00-scripts/03_rename_samples.sh
```

#### Deleting samples with too few reads

If after splitting your samples you notice that some of have too few reads, you
can remove these from the `04-all_samples` folder and modify your
`01-info_files/sample_information.csv` file accordingly. The threshold for the
minimum number of reads will depend on your project, including on the number of
expected cut sites generated by your library preparation protocol.

#### Align reads to a reference genome (optional)

#### Install `bwa <http://bio-bwa.sourceforge.net>`_

#### Download reference genome to the `01-info_files`

#### Index the reference genome

```bash
bwa index -p genome -a bwtsw ./01-info_files/<genome reference>
mv genome.* 01-info_files
```

###Align samples

```bash
./00-scripts/bwa_commands.sh
```

## STACKS pipeline

### Prepare population info file

```bash
./00-scripts/04_prepare_population_map.sh
```

## Edit script parameters

You will need to go through all the scripts named `stacks_*` in the `00-scripts
folder` and edit the options to suite your needs.

**Warning!** This step is most important. Choosing appropriate parameters for
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
./00-scripts/stacks_1b_ustacks.sh
./00-scripts/stacks_2_cstacks.sh
./00-scripts/stacks_3_sstacks.sh
./00-scripts/stacks_4_populations.sh
./00-scripts/stacks_5a_rxstacks_likelihoods.sh
```

Visualize the distribution of log likelihoods in `rxstacks_log_likelihoods.png`
and choose a cutoff to use in the next script
(`./00-scripts/stacks_5b_rxstacks.sh`). Then launch:

```bash
./00-scripts/stacks_5b_rxstacks.sh
./00-scripts/stacks_6_cstacks_rx.sh
./00-scripts/stacks_7_sstacks_rx.sh
./00-scripts/stacks_8_populations_rx.sh
```

## Filtering the results

`Stacks_workflow` includes a script to filter the VCF file output by STACKS.
To print the short documentation of the filtering script, launch the script
without options.

```bash
./00-scripts/05_filter_vcf.py
```

For the long documentation, use the -h option.

```bash
./00-scripts/05_filter_vcf.py -h
```

The following example shows how to use the script with some of the options.
These parameter values are only for demonstration purpose. Choose your
own values carefully.

```bash
./00-scripts/05_filter_vcf.py \  
   -i 05-stacks/batch_1.sumstats.tsv \  
   -o filtered.tsv \  
   -c 1 -m 7 -l 10 -I 8 -p 70 --use_percent \
   -a 0.01 -A 0.05 -H 0.5 -f -0.3 -F 0.3 -s 10
```

## Conclusion

### Running into problems

1. Consider joining the [Stacks Google
   group](https://groups.google.com/forum/#!forum/stacks-users)
1. To help you solve computer related problems and start focusing back on
   biology, you should consider getting this book: [Practical Computing for
   Biologists](http://practicalcomputing.org) (go for the paper edition). 
1. [Biostar](https://www.biostars.org) is a useful bioinformatics forum.
1. [SEQanswers](http://seqanswers.com) is another useful forum for all things
   related to Next Generation Sequencing.
1. Of course, you can always ask [Google](https://www.google.com) for help.

