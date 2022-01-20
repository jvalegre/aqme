.. _basic:

===============
Basic Structure
===============

This documents provides steps for running the ``aqme`` program.

.. contents::

Command line
------------

The ``aqme`` program can be run from the command line by providing the necessary arguments. In the following example,
an input file is provided and the compute option in invoked.

.. code:: bash

  $ python -m aqme --input [file] --compute


YAML file
---------

The ``aqme`` program can also be run by using a yaml file which contains the necessary arguments. In the following example,
a yaml file named ``params.yaml`` is created where an input file name is provided and the compute option in invoked.

.. code-block:: yaml

  # INPUT FILE
  input : 'file.smi' # input files

  # COMPUTE OPTION TRUE
  compute : True

To run ``aqme`` program with the yaml file, the following command is used

.. code:: bash

  $ python -m aqme --varfile params.yaml


Additional parameters which are needed to be changed can be added in the yaml file. A complete YAML file is
present it the defaults section for reference.

.. note::  If no specific arguments are provided, the default setting would be used. The default arguments can be found in the default section.
