Installing Stacks_workflow
**************************

Download and install the most recent version of this workflow
=============================================================

From the terminal
-----------------

.. code-block:: bash

 cd ~/Desktop
 wget https://github.com/enormandeau/stacks_workflow/archive/master.zip
 unzip master.zip

Using git
---------

If ``git`` is installed on your computer, you can run the following command
instead to get a complete ``git`` stacks_workflow repository.

.. code-block:: bash

 git clone https://github.com/enormandeau/stacks_workflow

This will make it easier to update to the most recent version of the workflow,
which can be accomplished by running the following command from the
``stacks_workflow`` directory:

.. code-block:: bash

 git pull

For the rest of the project, use the extracted or cloned folder as your working
directory. **All the commands in this manual are launched from that
directory.**

Set up for Mac OSX
=======================

If you have a MacOS computer, make sure you prepare your computer to run you
GBS analyses. Thierry Gosselin prepared a companion tutorial about GBS for
cloud computing and you should read the section about setting up your MacOS
computer:
<http://gbs-cloud-tutorial.readthedocs.org/en/latest/03_computer_setup.html#mac-osx>`_

Download and install `STACKS <http://creskolab.uoregon.edu/stacks/>`_
=====================================================================

Installing Google's SparseHash
------------------------------

Google SparseHash is a hash table implementation with a small memory footprint.
STACKS can (and should) be compiled using it. This will make STACKS use much
less memory.

.. code-block:: bash

 wget http://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz

 tar -xvf sparsehash-2.0.2.tar.gz
 cd sparsehash-2.0.2
 ./configure
 make  # Add '-j n' to use n CPUs during the compilation
 sudo make install

Installing STACKS
-----------------

.. code-block:: bash

 # Modify the version number as needed
 wget http://creskolab.uoregon.edu/stacks/source/stacks-1.29.tar.gz

 tar -xvf stacks-1.29.tar.gz
 cd stacks-1.29
 
 # Install the binaries in /usr/local/bin
 ./configure --enable-sparsehash
 make  # Add '-j n' to use n CPUs during the compilation
 sudo make install
 
 # Remove the temporary install folders and archives
 cd ..
 sudo rm -R stacks-1.29 stacks-1.29.tar.gz sparsehash-2.0.2 sparsehash-2.0.2.tar.gz
 
Test the STACKS installation
----------------------------
 
.. code-block:: bash

 cstacks

This will output the help of the cstacks program. You will also be able to
confirm the version number of your STACKS installation.

Installing Cutadapt
-------------------

There are different ways you can install Cutadapt. If you have ``pip`` (a
Python package installer) installed, you can use the following command:

.. code-block:: bash

 sudo pip install --user --upgrade cutadapt

Otherwise, visit their website to download it and install it:
<https://pypi.python.org/pypi/cutadapt/>`_

Installing jellyfish
--------------------

Jellyfish is a kmer counter. We will use it to identify the presence of
adapters in our raw sequences in order to remove them with Cutadapt. To
install Jellyfish, launch the following commands:

.. code-block:: bash

 # Getting the latest version
 wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.0/jellyfish-2.2.0.tar.gz

 # Installing
 tar xvfz jellyfish-2.2.0.tar.gz
 cd jellyfish-2.2.0
 ./configure
 make  # Add '-j n' to use n CPUs during the compilation
 sudo make install

 # Cleanup
 cd ..
 rm -r jellyfish-2.2.0.tar.gz jellyfish-2.2.0


