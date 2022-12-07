Metal complex from SMILES input
-------------------------------

In the following example we will: 

1) Generate various conformers of a Pd complex using RDKit using a template 
   for square planar complexes.
2) Generate Gaussian input files using an ECP for each of the conformers.

Step 1: CSEARCH conformational sampling (creates SDF files)
...........................................................

.. code:: python

    import os, glob
    from pathlib import Path                                                                                                                                                          
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    
    # set working directory and SMILES string
    w_dir_main = Path(os.getcwd())
    sdf_path = w_dir_main.joinpath('Pd_sdf_files')
    smi_metal = 'I[Pd]([PH3+])(F)Cl'
    
    # run CSEARCH conformational sampling, specifying:
    # 1) PATH to create the new SDF files (destination=sdf_path)
    # 2) Simple RDKit sampling (program='rdkit')
    # 3) SMILES string (smi=smi_metal)
    # 4) Name for the output SDF files (name='Pd_complex')
    # 5) The metal is Pd (metal=['Pd'])
    # 6) Oxidation number +2 (metal_oxi=[2])
    # 7) The complex is squareplanar (complex_type='squareplanar')
    csearch(destination=sdf_path,program='rdkit',smi=smi_metal,name='Pd_complex',
            metal_atoms=['Pd'],metal_oxi=[2],mult=1,complex_type='squareplanar')

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
     

ToRemove
........

Bonus 1: If you want to use the same functions using a YAML file that stores all the variables
                                                                                              

.. code:: python

    # to load the variables from a YAML file, use the varfile option
    csearch(varfile='FILENAME.yaml')
    
    # for each option, specify it in the YAML file as follows:
    # program='rdkit' --> program: 'rdkit'
    # name='Pd_complex' --> name: 'Pd_complex'
    # etc

Bonus 2: If you want to use the same functions through command lines
                                                                    

.. code:: python

    csearch(w_dir_main=w_dir_main,destination=sdf_path,program='rdkit',smi=smi_metal,name='Pd_complex',
            metal=['Pd'],metal_oxi=[2],mult=1,complex_type='squareplanar')
    
    # for each option, specify it in the command line as follows:
    # program='rdkit' --> --program 'rdkit'
    # name='Pd_complex' --> --name Pd_complex
    # etc
    # for example: python -m aqme --program rdkit --smi I[Pd](Cl)([PH3+])[N+]1=CC=CC=C1 --name Pd_complex
