Filtering the results
*********************

``Stacks_workflow`` includes a script to filter the VCF file output by STACKS.
To print the short documentation of the filtering script, launch the script
without options.

.. code-block:: bash

 ./00-scripts/05_filter_vcf.py

For the long documentation, use the -h option.

.. code-block:: bash

 ./00-scripts/05_filter_vcf.py -h

The following example shows how to use the script with some of the options.
These options are only for demonstration purpose. Choose your threshold values
carefully.

.. code-block:: bash

 ./00-scripts/05_filter_vcf.py \  
    -i 05-stacks/batch_1.sumstats.tsv \  
    -o filtered.tsv \  
    -c 2 -l 10 -I 8 -C 50 -p 70 --use_percent \
    -a 0.01 -A 0.05 -H 0.5 -f -0.2 -F 0.1 -s 10

