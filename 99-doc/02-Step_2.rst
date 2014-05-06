Install and prepare Stacks_workflow
===================================

Download and install the most recent version of this workflow
-------------------------------------------------------------

From Terminal, run:

.. code-block:: bash

 cd ~/Desktop
 wget https://github.com/enormandeau/stacks_workflow/archive/master.zip
 unzip master.zip


.. Note::

 Mac user have **git** pre-installed with Mavericks, run the following command:

.. code-block:: bash

 git clone https://github.com/enormandeau/stacks_workflow


Use the extracted or cloned folder as your working directory for the rest of the project. All the commands in this manual are launched from that directory.

`2. Download and install STACKS <http://creskolab.uoregon.edu/stacks/>`_
------------------------------------------------------------------------

.. code-block:: bash

 wget http://creskolab.uoregon.edu/stacks/source/stacks-1.19.tar.gz # you have to use version 1.14 instead? Just change 1.19 to 1.14.
 tar -xvf stacks-1.19.tar.gz
 cd stacks-1.19
 ./configure
 make
 sudo make install # this will install the binaries in /usr/local/bin
 cd ..
 sudo rm -R stacks-1.19 stacks-1.19.tar.gz # to remove the folder and gz file
 
 
Test the installation with this command:
 
.. code-block:: bash

 populations
