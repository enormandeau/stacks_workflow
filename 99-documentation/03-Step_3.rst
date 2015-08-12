Prepare your raw datafiles
**************************

Downloading your data 
=====================

Download your raw Illumina or Ion Proton data files from your sequencing
service provider.

Copy your raw files
===================

Put a copy of (or a link to) your raw data files in the ``02-raw`` folder of
stacks_workflow.

.. Note::

 All file names **must** end with **.fastq.gz** for the following scripts to
 work.

Preparing the ``lane_info.txt`` file
====================================

This file will contains the names of the raw data files and is used by
stacks_workflow later.  From the stacks_workflow folder, run:

.. code-block:: bash

 ./00-scripts/00_prepare_lane_info.sh

Removing adapters
=================

We use Cutadapt in order to remove the adapters present in the raw data. This
step is done before we extract the reads for each individual.

Creating a file with adapters
-----------------------------

In order to clean your reads properly, you must determine which adapters are
present in your data and whether a cut site can also be found with it.

Before running cutadapt, we must create a fasta file containing the adapters
that are present in our data. This file must be found in the ``01-info_files/``
folder and be names ``adapters.fasta``. The file named
``example_adapters.fasta``, found in the ``01-info_files/`` folder should be
used as a template.

The best way to know which adapters and potential restriction enzyme cut sites
are present in your data is probably to ask the people who created your
sequencing libraries.

In order to identify or confirm the potential adapters found in your data
directly, you can use the following command (you will need to have installed
``jellyfish``):

.. code-block:: bash

 cat $(ls -1 02-raw/*.fastq.gz | head -1) | \
    gunzip -c | grep -E "^[ACTG]+$" | head -10000 | cut -c 20- | \
    awk 'length > 30 {print ">seq_" NR "\n" $_}' > temp.fasta

 rm kmer.count.*
 for k in $(seq 12 30)
 do
    jellyfish count -m $k -s 10000 -t 8 --text -o kmer.count."$k" temp.fasta
    perl -pe 's/^.+\}[^ACTG]+//' kmer.count."$k" | awk '{print $2,$1}' | \
    sort -nr | head -20 > kmer.count.sorted."$k"; done

Then look at the results with the following command. From the number of
identical kmers and how they vary from one kmer length to the next, you should
be able to spot your adapters and potentially a version containing a cut site
just before it if you used two restriction enzymes in your library preparation.

.. code-block:: bash

 for i in kmer.count.sorted.*; do echo $i; cat $i; done | less

 # Delete all these files
 rm kmer.count*

Running Cutadapt
----------------

Now that the ``01-info_files/adapters.fasta`` file contains all the adapters we
wish to remove from our data, we can launch Cutadapt **in single-end mode**
with the following command:

.. code-block:: bash

 ./00-scripts/01_cutadapt.sh numCPUs

Where `numCPUs` is the number of CPUs you wish to use in parallel. If you do
not put a number, the script will use only one CPU (for backward
compatibility).

Scan the cutadapt logs
----------------------

The cutadapt log files can be found in the ``98-log_files`` folder. Scan them
to confirm that cutadapt has done an appropriate job.

.. Note::

 Ther may be difference in adapters and filter parameters to use with data
 produced by Illumina and Ion Proton sequencers.

Assess sequence quality with FastQC
===================================

**TODO:**

- Before and after Cutadapt
- Use the script
  - Command line interface
  - Run only on first 1,000,000 reads or less
  - Gathers important figures together

