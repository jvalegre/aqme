.. |metal_comp_chemdraw| image:: ../../images/metal_comp_chemdraw.png
   :width: 300

.. |metal_comp_3D| image:: ../../images/Quinine-3D-balls.png
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



Step 2: Writing Gaussian input files with the SDF obtained from CSEARCH
.......................................................................

.. code:: python

    # set SDF filenames and directory where the new com files will be created
    com_path = sdf_path.joinpath(f'com_files')
    sdf_rdkit_files = glob.glob(f'{sdf_path}/*.sdf')
    
    # run QPREP input files generator, with:
    # 1) PATH to create the new SDF files (destination=com_path)
    # 2) Files to convert (files=sdf_rdkit_files)
    # 3) QM program for the input (program='gaussian')
    # 4) Keyword line for the Gaussian inputs (qm_input='B3LYP/genecp opt freq')
    # 5) Basis set to use in the atoms included in genECP (bs_gen='def2svp')
    # 6) Basis set to use in the rest of the atoms (bs='6-31G*')
    # 7) Atoms to include as genECP (gen_atoms=['Pd'])
    # 8) Memory to use in the calculations (mem='24GB')
    # 9) Processors to use in the calcs (nprocs=8)
    qprep(destination=com_path,files=sdf_rdkit_files,program='gaussian',qm_input='B3LYP/genecp opt freq',
            bs_gen='def2svp',bs_nogen='6-31G*',gen_atoms=['Pd'],mem='24GB',nprocs=8)
