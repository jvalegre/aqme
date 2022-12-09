.. |nocov_ts_chemdraw| image:: ../../images/nocov_ts_chem.png
   :width: 400

.. |nocov_ts_3D| image:: ../../images/Quinine-3D-balls.png
   :width: 400

.. |mapping| image:: ../../images/nocov_TS_map.png
   :width: 400

TS involving a trimolecular complex generated from SMILES
=========================================================

In the following example we will generate conformers of the SN2 transition state 
involving Water, Cloride ion and 2-Fluoro-2-methylpropane. 

+-----------------------------------------------+
| .. centered:: **SMILES**                      |
+-----------------------------------------------+
| .. centered:: O.FC(C)(C)C.[Cl-]               |
+--------------------------+--------------------+
|  |nocov_ts_chemdraw|     |  |nocov_ts_3D|     |
+--------------------------+--------------------+

In the following example we need a mapped SMILES to set up the constraints. We 
can either use rdkit for it in python, in a jupyter notebook (see 
:doc:`python example <../../examples_python/conformer_search/noncovalent_TS>` )
or we can write the SMILES by ourselves (not recommended). 

The mapped smiles that we will be using is the following one: 

|mapping|

:: 

   [Cl-:7].[F:2][C:3]([C:4]([H:10])([H:11])[H:12])([C:5]([H:13])([H:14])[H:15])[C:6]([H:16])([H:17])[H:18].[O:1]([H:8])[H:9]


Now that we have the mapping, we can easily proceed to set up the constraints.
In this case we want the C-F and the Cl-C bond distances to be constrained and 
equal to 1.8 angstroms and we want the angle Cl-C-F to be of 180º. 

As a consequence the distance constraints will be :code:`"[[2,3,1.8],[3,7,1.8]]"`
and the dihedral constrains :code:`"[[2,3,7,180]]"`

Finally we proceed to the conformer generation using CREST

.. code:: shell

   python -m aqme --name "TS-example" --program crest --cregen --crest_keywords "--nci" --constraints_dist "[[2,3,1.8],[3,7,1.8]]" --constraints_angle "[[2,3,7,180]]" --smi "[Cl-:7].[F:2][C:3]([C:4]([H:10])([H:11])[H:12])([C:5]([H:13])([H:14])[H:15])[C:6]([H:16])([H:17])[H:18].[O:1]([H:8])[H:9]"




