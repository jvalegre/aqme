Solubility descriptors
======================

Along the steps of this example workflow we will show how to: 

i)   RDKit conformer sampling
ii)  xTB porperty calculations to determine molecular and atomic
     properties
iii) Generate a GNN model to predict solubility

Specifically in this example we will calculate descriptors of a set of molecules
represented by smiles with the aim of training a Graph Neural Network (GNN) 
capable of predicting their solubility.

The set of molecules that we will use are in a file named 'solubility.csv' 
that has the following contents:

.. highlight:: none

.. literalinclude:: ../../chemfiles/end_to_end_3_inp.csv

.. highlight:: default

.. note::

   A jupyter notebook containing all the steps shown in this example can be found 
   in the AQME repository in `Github  <https://github.com/jvalegre/aqme>`__ or in 
   `Figshare <https://figshare.com/articles/dataset/AQME_paper_examples/20043665/11>`__


.. contents:: Steps
   :local:
   :depth: 2


Step 1: Run CSEARCH (RDKit) on the CSV
--------------------------------------

.. code:: shell

   python -m aqme --csearch --program rdkit --input solubility.csv --ewin_csearch 1

Step 2 : Run xTB calculations using QDESCP
------------------------------------------

.. code:: shell

   python -m aqme --qdescp --files "CSEARCH/*.sdf" --program xtb


Step 3 : Create the CSV file with descriptors for the GNN model
---------------------------------------------------------------

This step can be run in bash, python or manually using a spreadsheet editor.
The aim is to generate a file named 'solubility_xtb.csv' containing an extra 
column with the filepaths to the QDESCP generated .json files.

Here we include the python commands to do the step.

.. code:: python

   import pandas as pd 

   data =  pd.read_csv(file)
   code_to_filepath = 'QDESCP/boltz/{}_rdkit_boltz.json'.format
   data['xtbjson'] = data['code_name'].apply(code_to_filepath)
   data.to_csv('solubility_xtb.csv',index=False)


Step 4: Run the gnn.py to get results
-------------------------------------

This step requires other files that do not use aqme itself. The other files 
are available at 
`Figshare <https://figshare.com/articles/dataset/AQME_paper_examples/20043665/11>`__ .
Specifically the files gnn.py and gnn_functions.py require to be on the same 
directory as the file 'solubility.csv' and these two files depend on packages 
that aqme does not depend on. So before continuing please make sure that you 
have all the packages required installed as well as the specified files.

We can execute directly the gnn.py script to obtain the results

.. code:: shell

    python gnn.py

For more details on the contents of the gnn.py script please go to 
:doc:`python example <../../examples_python/end_to_end/example_3>` ) steps 5a, 
5b and 5c. 

