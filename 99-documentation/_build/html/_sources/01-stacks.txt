Stacks workflow tutorial
************************

This workflow aims at making the use of the STACKS pipeline easier and more
structured so that people tasked with analysing GBS or RAD data and possessing
limited UNIX/Linux experience can jump on the analysis wagon faster. 

It was developed with the needs of our research group in mind with as well as
with an emphasis on non-model species studies. We make no claim about its
usefulness to other groups or in other contexts, but we still believe it may be
of use to some.

This workflow has been tested with version 1.24 and earlier versions of STACKS
under Linux (Ubuntu 12.04 to 13.10) and MacOSX (10.9 Mavericks and 10.10
Yosemite).

.. warning::

 The workflow will not work with versions of STACKS older than 1.24 since tags
 of different lengths can now be split together and intermediary files are kept
 in compressed (.gz) format.

About STACKS
============
 
 The `STACKS analysis pipeline <http://creskolab.uoregon.edu/stacks/>`_ is the
 de facto tool for SNP discovery in Genotyping-by-Sequencing (GBS) and
 Restriction-site Associated DNA sequencing (RAD) studies, with and without a
 reference genome. Here is diagram of the STACKS pipeline:
 
 .. image:: stacks_diagram.png

 Upon starting to use STACKS, it is highly suggested to read the two official
 STACKS papers:
 
 `Catchen, J. M., Amores, A., Hohenlohe, P. A., Cresko, W. A., Postlethwait, J.
 H., & De Koning, D. J. (2011). Stacks: Building and Genotyping Loci De Novo
 From Short-Read Sequences. G3, 1(3), 171–182. doi:10.1534/g3.111.000240
 <http://www.g3journal.org/content/1/3/171.full>`_
 
 `Catchen, J. M., Hohenlohe, P. A., Bassham, S., Amores, A., & Cresko, W. A.
 (2013). Stacks: an analysis tool set for population genomics. Molecular
 Ecology, 22(11), 3124–3140. doi:10.1111/mec.12354
 <http://onlinelibrary.wiley.com/doi/10.1111/mec.12354/abstract>`_
 
 Also very useful paper to read before attempting to run Stacks:
 
 `Mastretta Yanes A, Arrigo N, Alvarez N et al. (2014) RAD sequencing,
 genotyping error estimation and de novo assembly optimization for population 
 genetic inference. Molecular Ecology Resources, n/a–n/a.
 <http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12291/abstract;jsessionid=A32722E1462A2A2714EE53A6FD4C7194.f04t04>`_
 
Overview of the steps
=====================

 #. Install and prepare Stacks_workflow  
 #. Download your raw datafiles (Illumina lanes)
 #. Extract individual data with process_radtags
 #. Rename samples
 #. Align reads to a reference genome (optional)
 #. STACKS pipeline
 #. Filtering the results

Where to find this tutorial
===========================

 - In the pipeline folder itself. Go to `Github
   <https://github.com/enormandeau/stacks_workflow>`_, download the whole
   repository and open the file named ``MANUAL.html`` in your web browser.

Licence
=======

The Stacks_workflow is licensed under the GPL3 license. See the LICENCE file
for more details.

