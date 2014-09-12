Install Stacks_workflow
***********************

Download and install the most recent version of this workflow
=============================================================

From the terminal, run:

.. code-block:: bash

 cd ~/Desktop
 wget https://github.com/enormandeau/stacks_workflow/archive/master.zip
 unzip master.zip

If ``git`` is installed on your computer (Mac user have it pre-installed on
Mavericks), you can run the following command instead to get a full ``git`` repository.

.. code-block:: bash

 git clone https://github.com/enormandeau/stacks_workflow

This will make it easier to update to the most recent version of the workflow, which can be accomplished by running the following command from the ``stacks_workflow`` directory:

.. code-block:: bash

 git pull

For the rest of the project, use the extracted or cloned folder as your working
directory. All the commands in this manual are launched from that directory.

Download and install `STACKS <http://creskolab.uoregon.edu/stacks/>`_
=====================================================================

.. code-block:: bash

 # You can modify the version number as needed
 wget http://creskolab.uoregon.edu/stacks/source/stacks-1.20.tar.gz

 tar -xvf stacks-1.20.tar.gz
 cd stacks-1.20
 
 # Install the binaries in /usr/local/bin
 ./configure
 make
 sudo make install
 
 # Remove the folder and gz file
 cd ..
 sudo rm -R stacks-1.20 stacks-1.20.tar.gz
 
Test the installation
---------------------
 
.. code-block:: bash

 cstacks

This will output the help of the cstacks program. You will also be able to
confirm the version number of your STACKS installation.

