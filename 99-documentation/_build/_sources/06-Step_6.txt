Align reads to a reference genome (optional)
********************************************

Install `bwa <http://bio-bwa.sourceforge.net>`_
===============================================

Download reference genome to the ``01-info_files``
==================================================

Index the reference genome
==========================

.. code-block:: bash

 bwa index -p genome -a bwtsw ./01-info_files/<genome reference>

copy files:

.. code-block:: bash

 cp genome.* 01-info_files

Align samples
=============

.. code-block:: bash

 for i in $(ls -1 04-all_samples/*.fq)
 do
     name=$(basename $i)
     bwa aln -n 5 -k 3 -t 2 ./01-info_files/genome $i | \
     bwa samse -r "@RG\tID:'$name'\tSM:'$name'\tPL:Illumina" \
         ./01-info_files/genome - $i ./04ln-all_samples/$name.sam; \
 done

