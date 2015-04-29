Prepare your raw datafiles
**************************

Download Illumina lanes or Ion Proton chips
===========================================

Download the raw data files from your sequencing service provider.

Put your raw files in ``02-raw`` folder
=======================================

.. Note::

 All file names **must** end with **.fastq.gz** for the following scripts to
 work.

Prepare ``lane_info.txt`` file
==============================

From the stacks_workflow folder, run:

.. code-block:: bash

 ./00-scripts/01_prepare_lane_info.sh

Use Cutadapt to remove adapters
===============================
TODO

- Illumina vs. Ion Proton
- Create file with adapters (script to confirm)
- Run ``Cutadapt`` script
- Scan the logs to assess the damage

Assess sequence quality with FastQC
===================================
TODO

- Before and after Cutadapt
- Use the script
  - Command line interface
  - Run only on first 1,000,000 reads or less
  - Gathers important figures together

