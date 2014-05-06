Download your raw datafiles (Illumina lanes)
============================================

Put your raw files in the ``02-raw`` folder of the stacks_workflow folder
-------------------------------------------------------------------------

.. Note::

 All file names MUST end with **.fastq.gz**


Prepare the ``lane_info.txt`` file automatically
------------------------------------------------

From the stacks_workflow folder, run:

.. code-block:: bash

 ./00-scripts/01_prepare_lane_info.sh

