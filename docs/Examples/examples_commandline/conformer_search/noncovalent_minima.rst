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

Compared with the :doc:`organic molecule example <organic_molecule>` this time 
we are going to use the default location for the output file. As a consequence
we can now proceed to run the conformer search from the smiles string:

.. code:: shell

   python -m aqme --csearch --smi "CCC(C)C.O" --name "isopent-water-complex" --program crest --crest_keywords "--nci" --cregen --cregen_keywords "--ewin 3"


