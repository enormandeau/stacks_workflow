Filtering the results
*********************

``Stacks_workflow`` includes a script to filter the VCF file output by STACKS.
To print the short documentation of the filtering script, launch it without
options.

.. code-block:: bash
 ./00-scripts/05_filter_vcf.py

For the long documentation, use the -h option.

.. code-block:: bash
 ./00-scripts/05_filter_vcf.py -h

The following example shows how to use the script with some of the options.
These options are only for demonstration purpose. Choose your threshold values
carefully.

.. code-block:: bash

 ./00-scripts/05_filterStacksSNPs.py \  
    -i 05-stacks/batch_1.sumstats.tsv \  
    -o filtered.tsv \  
    -p 2 -x 1 -H 0.7 -a 0.05 -A 0 -s 10

