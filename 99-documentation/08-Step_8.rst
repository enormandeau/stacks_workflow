Filtering the results
*********************

``Stacks_workflow`` includes a script to filter your SNPs. To print the
documentation of the filtering script, launch it without options:

.. code-block:: bash
 ./00-scripts/05_filterStacksSNPs.py --help

The following example shows how to use the script with some of the options.
These options are only for demonstration purpose. Choose your threshold values
carefully.

.. code-block:: bash

 ./00-scripts/05_filterStacksSNPs.py \  
    -i 05-stacks/batch_1.sumstats.tsv \  
    -o filtered.tsv \  
    -P 01-info_files/population_map.txt \  
    -p 2 -x 1 -H 0.7 -a 0.05 -A 0 -f -0.3 -F 0.8 -s 10


We also recommend you look at your data with Stacks Web-based interface.
`The installation instruction are also in the cloud tutorial <http://gbs-cloud-tutorial.readthedocs.org/en/latest/09_stacks_web_interface.html>`_
