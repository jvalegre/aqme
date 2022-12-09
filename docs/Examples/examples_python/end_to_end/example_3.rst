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

::

   code_name,Molecule,ESOL predicted log solubility in mols per litre,Minimum Degree,Molecular Weight,Number of H-Bond Donors,Number of Rings,Number of Rotatable Bonds,Polar Surface Area,measured log solubility in mols per litre,smiles
   mol_1,Amigdalin,-0.974,1,457.432,7,3,7,202.32,-0.77,OCC3OC(OCC2OC(OC(C#N)c1ccccc1)C(O)C(O)C2O)C(O)C(O)C3O 
   mol_2,Fenfuram,-2.885,1,201.225,1,2,2,42.24,-3.3,Cc1occc1C(=O)Nc2ccccc2
   mol_3,citral,-2.579,1,152.237,0,0,4,17.07,-2.06,CC(C)=CCCC(C)=CC(=O)
   mol_4,Picene,-6.618,2,278.354,0,5,0,0,-7.87,c1ccc2c(c1)ccc3c2ccc4c5ccccc5ccc43
   mol_5,Thiophene,-2.232,2,84.143,0,1,0,0,-1.33,c1ccsc1
   mol_6,benzothiazole,-2.733,2,135.191,0,2,0,12.89,-1.5,c2ccc1scnc1c2 
   mol_7,"2,2,4,6,6'-PCB",-6.545,1,326.437,0,2,1,0,-7.32,Clc1cc(Cl)c(c(Cl)c1)c2c(Cl)cccc2Cl
   mol_8,Estradiol,-4.138,1,272.388,2,4,0,40.46,-5.03,CC12CCC3C(CCc4cc(O)ccc34)C2CCC1O
   mol_9,Dieldrin,-4.533,1,380.913,0,5,0,12.53,-6.29,ClC4=C(Cl)C5(Cl)C3C1CC(C2OC12)C3C4(Cl)C5(Cl)Cl
   mol_10,Rotenone,-5.246,1,394.423,0,5,3,63.22,-4.42,COc5cc4OCC3Oc2c1CC(Oc1ccc2C(=O)C3c4cc5OC)C(C)=C 
   mol_11,2-pyrrolidone,0.243,1,85.106,1,1,0,29.1,1.07,O=C1CCCN1
   mol_12,2-Chloronapthalene,-4.063,1,162.619,0,2,0,0,-4.14,Clc1ccc2ccccc2c1
   mol_13,1-Pentene ,-2.01,1,70.135,0,0,2,0,-2.68,CCCC=C

.. note::

   A jupyter notebook containing all the steps shown in this example can be found 
   in the aqme repository in `Github  <https://github.com/jvalegre/aqme>`__ or in 
   `Figshare <https://figshare.com/articles/dataset/AQME_paper_examples/20043665/11>`__


.. contents:: Steps
   :local:
   :depth: 2


Step 1: Import AQME and other python modules, and the required CSV
------------------------------------------------------------------

.. code:: python

   import glob
   from aqme.csearch import csearch
   from aqme.qdescp import qdescp
   
   file = 'solubility.csv'

Step 2: Run CSEARCH (RDKit) on the CSV
--------------------------------------

.. code:: python

   csearch(program='rdkit',
           input=file,
           ewin_csearch=1)

Step 3 : Run xTB calculations using QDESCP
------------------------------------------

.. code:: python

   sdf_rdkit_files = glob.glob(f'CSEARCH/*.sdf')
   qdescp(files=sdf_rdkit_files,
          boltz=True,
          program='xtb')


Step 4 : Create the CSV file with descriptors for the GNN model
---------------------------------------------------------------

.. code:: python

   data =  pd.read_csv(file)
   code_to_filepath = 'QDESCP/boltz/{}_rdkit_boltz.json'.format
   data['xtbjson'] = data['code_name'].apply(code_to_filepath)
   data.to_csv('solubility_xtb.csv',index=False)

Step 5: Run the gnn.py to get results
-------------------------------------

This step requires other files that do not use aqme itself. The other files 
are available at 
`Figshare <https://figshare.com/articles/dataset/AQME_paper_examples/20043665/11>`__ .
Specifically the files gnn.py and gnn_functions.py require to be on the same 
directory as the file 'solubility.csv' and these two files depend on packages 
that aqme does not depend on. So before continuing please make sure that you 
have all the packages required to continue as well as the specified files.

we can execute directly the gnn.py script to obtain the results

.. code:: shell

    python gnn.py

But the following will show the three main steps included in the gnn.py file. 
We will start by importing the necessary modules.

.. code:: python 

   import pandas as pd
   import numpy as np
   import matplotlib.pyplot as plt
   import seaborn as sns
   from gnn_functions import *
   from sklearn.metrics import r2_score
   import sklearn.metrics as metrics
   import tensorflow as tf

Step 5a: Load the solubility CSV file and split the data into training, validation and test sets
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. code:: python

    sol = pd.read_csv('solubility_xtb.csv')
    valid, test, train = np.split(sol[['smiles','xtbjson']].sample(frac=1., random_state=41), [50, 100])

Step 5b: Set up the GNN dataset
+++++++++++++++++++++++++++++++

.. code:: python

    train_dataset, valid_dataset, test_dataset = gnn_data(valid, test, train, sol)
    inputs, outputs = next(train_dataset.as_numpy_iterator())


Step 5c: Set up the GNN model
+++++++++++++++++++++++++++++

.. code:: python

    model = gnn_model()
    model.compile(loss='mae', optimizer=tf.keras.optimizers.Adam(1E-3))
    model.fit(train_dataset, validation_data=valid_dataset, epochs=500)

Step 5d: Predict solubities of external test set using the GNN model
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. code:: python

    # Predict solubility of the external test set
    test_predictions = model.predict(test_dataset)
    test_db_values = sol.set_index('smiles').reindex(test.smiles)['measured log solubility in mols per litre'].values
    
    # Plot the results
    fig = plt.subplots(figsize=(3,3))
    
    ax1 = sns.scatterplot(test_db_values,test_predictions.flatten(),s=30,marker='o',color='b',alpha=0.5)
    ax1.set_xlabel(r'Measured',fontsize=10)
    ax1.set_ylabel(r'Predicted',fontsize=10)
    ax1.grid(linestyle='--', linewidth=1)
    
    mae = metrics.mean_absolute_error(test_db_values,test_predictions.flatten())
    r2 = metrics.r2_score(test_db_values,test_predictions.flatten())
    
    plt.annotate(f"$R^2$ = {round(r2,1)} \nMAE = {round(mae,1)} ", xy=(-1.5, -5.9), fontsize=10)
    plt.savefig('solubility-gnn.jpg',dpi=400)
    plt.show()

