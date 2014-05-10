STACKS pipeline
***************

Prepare population info file
============================

.. code-block:: bash

 ./00-scripts/04_prepare_population_map.sh

Edit script parameters
======================

You will need to go through the scripts named ``stacks_*`` in the ``00-scripts
folder`` and edit the options.

.. warning::

 Choosing appropriate parameters for your study is crucial in order to generate
 meaninful and optimal results.

Run the STACKS programs
=======================

If you do not have access to a reference genome, launch the following commands:

.. code-block:: bash

 ./00-scripts/stacks_1a_ustacks.sh
 ./00-scripts/stacks_2_cstacks.sh
 ./00-scripts/stacks_3_sstacks.sh
 ./00-scripts/stacks_4_populations.sh

Or, if you have a reference genome:

.. code-block:: bash

 ./00-scripts/stacks_1b_pstacks.sh
 ./00-scripts/stacks_2_cstacks.sh
 ./00-scripts/stacks_3_sstacks.sh
 ./00-scripts/stacks_4_populations.sh

