Rename samples
**************

Rename and copy the samples
===========================

.. code-block:: bash

 ./00-scripts/03_rename_samples.sh

Join samples that should go together
====================================

If you need to combine some samples, for example if some samples were run
multiple times on different lanes, you will need to combine their different
files. This step is done manually.

 - Go to ``04-all_samples`` and join the .fq files that should go together with
   the ``cat`` command
 - Remove partial .fq files that have been joined
 - Remove individuals with too few sequences if needed (optional)

