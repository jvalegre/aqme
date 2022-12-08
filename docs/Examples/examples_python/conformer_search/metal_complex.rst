.. |metal_comp_chemdraw| image:: ../../images/metal_comp_chemdraw.png
   :width: 300

.. |metal_comp_3D| image:: ../../images/metal_comp_3D.png
   :width: 400

Metal complex from SMILES input
===============================

In the following example we will: 

1) Generate various conformers of a Pd complex using RDKit using a template 
   for square planar complexes.
2) Generate Gaussian input files using an ECP for each of the conformers.

+-----------------------------------------------+
| .. centered:: **SMILES**                      |
+-----------------------------------------------+
| .. centered:: I[Pd]([PH3+])(F)Cl              |
+--------------------------+--------------------+
|  |metal_comp_chemdraw|   |  |metal_comp_3D|   |
+--------------------------+--------------------+

As with the other examples, we start by importing the necesary packages

.. code:: python

    from pathlib import Path
    from aqme.csearch import csearch

Next we specify the output directory where the sdf with conformations will be 
generated

.. code:: python

    w_dir_main = Path.cwd()
    sdf_path = w_dir_main/'Pd_sdf_files'

Finally we proceed to the conformational search using the smiles string of the 
molecule. 

.. code:: python

    smiles = 'I[Pd]([PH3+])(F)Cl'
    csearch(destination=sdf_path,
            smi=smiles,
            name='Pd_complex',
            program='rdkit',
            metal_atoms=['Pd',],         # Symbol of transition metal atoms included
            metal_oxi=[2,],              # Oxidation number per metal_atom
            mult=1,                      # multiplicity   
            complex_type='squareplanar') # Template geometry to use
