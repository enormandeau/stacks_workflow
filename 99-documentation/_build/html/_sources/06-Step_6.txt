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

Where `<genome reference>` is the name of the genome reference file.

Move the index files:

.. code-block:: bash

 mv genome.* 01-info_files

Align samples
=============

.. code-block:: bash

 ./00-scripts/bwa_commands.sh

