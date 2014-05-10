Extract individual data with ``process_radtags``
************************************************

Prepare a file named ``sample_information.csv``
===============================================

Using the same format found in the ``example_sample_information.csv`` file in
the ``01-info_files`` folder. 

Save this file in the ``01-info_files`` folder.

The ``sample_information.csv`` file will be used to extract the samples and
rename the extracted sample files automatically. 

The first column **MUST** contain the **EXACT** name of the data file for the
lane of each sample. 

The second column contains the barcode sequence of each sample. 

The third column contains the population name of each sample. 

The fourth column contains the name of the sample (do not include the
population name or abbreviation in the sample name). 

The fifth column contains a number identifying the populations. 

Columns three and four are treated as text, so they can contain either text or
numbers. Other columns can be present after the fifth one and will be ignored.
However, it is crucial that the five first columns respect the format in the
example file exactly. Be especially careful not to include errors in this file,
for example mixing lower and capital letters in population or sample names
(e.g.: Pop01 and pop01), since these will be treated as two different
populations.

Launch process_radtags with:
============================

.. code-block:: bash

 ./00-scripts/02_process_radtags.sh <trimLength> <enzyme>

Where:  

 - **trimLength** = length to trim all the sequences. This should be the length
   of the Illumina reads minus the length of the longest tag or MID.  
 - **enzyme** = name of enzyme (run ``process_radtags``, without options, for a
   list of the supported enzymes)

