STACKS pipeline
***************

Prepare population info file
============================

.. code-block:: bash

 ./00-scripts/04_prepare_population_map.sh

Edit script parameters
======================

You will need to go through all the scripts named ``stacks_*`` in the
``00-scripts folder`` and edit the options to suite your needs.

.. warning::

 This step is crucial. Choosing appropriate parameters for your study is
 crucial in order to generate meaninful and optimal results. Read the STACKS
 documentation on their website to learn more about the different options.

Run the STACKS programs
=======================

Without a reference genome
--------------------------

.. code-block:: bash

 ./00-scripts/stacks_1a_ustacks.sh
 ./00-scripts/stacks_2_cstacks.sh
 ./00-scripts/stacks_3_sstacks.sh
 ./00-scripts/stacks_4_populations.sh
 ./00-scripts/stacks_5a_rxstacks_likelihoods.sh

Visualize the distribution of log likelihoods in
``rxstacks_log_likelihoods.png`` and choose a cutoff to use in the next script
(``./00-scripts/stacks_5b_rxstacks.sh``). Then launch:

.. code-block:: bash

 ./00-scripts/stacks_5b_rxstacks.sh
 ./00-scripts/stacks_6_cstacks_rx.sh
 ./00-scripts/stacks_7_sstacks_rx.sh
 ./00-scripts/stacks_8_populations_rx.sh

With a reference genome
-----------------------

.. warning::

 The documentation and scripts used with a reference genome have not been
 updated in a long time. We believe they should not be used at the moment.


.. code-block:: bash

 ./00-scripts/stacks_1b_ustacks.sh
 ./00-scripts/stacks_2_cstacks.sh
 ./00-scripts/stacks_3_sstacks.sh
 ./00-scripts/stacks_4_populations.sh
 ./00-scripts/stacks_5a_rxstacks_likelihoods.sh

Visualize the distribution of log likelihoods in
``rxstacks_log_likelihoods.png`` and choose a cutoff to use in the next script
(``./00-scripts/stacks_5b_rxstacks.sh``). Then launch:

.. code-block:: bash

 ./00-scripts/stacks_5b_rxstacks.sh
 ./00-scripts/stacks_6_cstacks_rx.sh
 ./00-scripts/stacks_7_sstacks_rx.sh
 ./00-scripts/stacks_8_populations_rx.sh

