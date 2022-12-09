.. |nocov_bimol_chemdraw| image:: ../../images/nocov_bimol_chem.png
   :width: 400

.. |nocov_bimol_3D| image:: ../../images/Quinine-3D-balls.png
   :width: 400

Noncovalent bimolecular complex generated from SMILES
=====================================================

In the following example we will generate conformers of the 
isopentane-water complex using CREST

+-----------------------------------------------+
| .. centered:: **SMILES**                      |
+-----------------------------------------------+
| .. centered:: CCC(C)C.O                       |
+--------------------------+--------------------+
|  |nocov_bimol_chemdraw|  |  |nocov_bimol_3D|  |
+--------------------------+--------------------+

We start by importing the packages

.. code:: python

    from pathlib import Path
    from aqme.csearch import csearch

Compared with the :doc:`organic molecule example <organic_molecule>` this time 
we are going to use the default location for the output file. As a consequence
we can now proceed to run the conformer search from the smiles string:

.. code:: python

    smiles = 'CCC(C)C.O'
    csearch(smi=smiles,
            name='isopent-water-complex',
            program='crest',
            crest_keywords='--nci',     # indicate that it is a non-covalent complex
            cregen=True,                # Include CREGEN post-analysis
            cregen_keywords='--ewin 3') # energy window for CREGEN == 3.0 kcal/mol

