Stacks workflow tutorial
========================

This workflow aims at making the use of the STACKS pipeline easier and more structured so that people tasked with analysing GBS or RAD data and possessing limited UNIX/Linux experience can jump on the analysis wagon faster. 

It was developed with the needs of our research group in mind with an emphasis on non-model species studies. We make no claim about its usefulness to other groups or in other contexts, but we still believe it may be of use to some.

This workflow has been tested with version 1.19 and earlier versions of STACKS under Linux (Ubuntu 12.04 to 13.10) and MacOSX (10.9 Mavericks).


.. Note::

 **About STACKS**
 
 The `STACKS analysis pipeline <http://creskolab.uoregon.edu/stacks/>`_ is the de facto tool for SNP discovery in Genotyping-by-Sequencing (GBS) and Restriction-site Associated DNA sequencing (RAD) studies, with and without a reference genome. Here is diagram of Stacks pipeline:
 
 .. image:: stacks_diagram.png
    :width: 500pt

 **Upon starting to use STACKS, it is highly suggested to read the two official STACKS papers:**
 
 `Catchen, J. M., Amores, A., Hohenlohe, P. A., Cresko, W. A., Postlethwait, J. H., & De Koning, D. J. (2011). Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. G3, 1(3), 171–182. doi:10.1534/g3.111.000240 <http://www.g3journal.org/content/1/3/171.full>`_
 
 `Catchen, J. M., Hohenlohe, P. A., Bassham, S., Amores, A., & Cresko, W. A. (2013). Stacks: an analysis tool set for population genomics. Molecular Ecology, 22(11), 3124–3140. doi:10.1111/mec.12354 <http://onlinelibrary.wiley.com/doi/10.1111/mec.12354/abstract>`_
 



**Overview of the steps**

Step 1: Read the papers!

Step 2: Install and prepare Stacks_workflow  

Step 3: Download your raw datafiles (Illumina lanes)

Step 4: Extract individual data with process_radtags

Step 5: Rename samples

Step 6: Align reads to a reference genome (optional)

Step 7: STACKS pipeline

Step 8: Filtering the results


**The tutorial can be found in 3 places:**

 - `Github <https://github.com/enormandeau/stacks_workflow>`_
 - `Read the docs <>`_ ** will come soon **

**Licence:**

The Stacks_workflow is licensed under the GPL3 license. See the LICENCE file for more details.


Made by Eric Normandeau with the collaboration of Thierry Gosselin, Jeremy Gaudin, Laura Benestan and Charles Perrier ***(C'Est ok avec toi??)***

