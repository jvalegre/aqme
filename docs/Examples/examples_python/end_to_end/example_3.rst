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

This step requires other files that do not use AQME itself. The other files 
are available at 
`Figshare <https://figshare.com/articles/dataset/AQME_paper_examples/20043665/11>`__ .
Specifically the files gnn.py and gnn_functions.py require to be on the same 
directory as the file 'solubility.csv' and these two files depend on packages 
that AQME does not depend on. So before continuing please make sure that you 
have all the packages required installed as well as the specified files.

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

