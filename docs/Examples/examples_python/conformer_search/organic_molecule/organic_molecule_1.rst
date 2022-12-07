Using RDKit and ANI
-------------------

In the following example we will: 

1) Generate various conformers of the quinine molecule using RDKit
2) Minimize the generated conformers using ANI
3) Generate Gaussian input files for each of the conformers

Step 1: CSEARCH conformational sampling (creates SDF files)
...........................................................

.. code:: python

    import os, glob
    from pathlib import Path                                                                                                                                                          
    from aqme.csearch import csearch
    from aqme.cmin import cmin
    from aqme.qprep import qprep
    from aqme.cmin import cmin
    
    # set working directory and SMILES string
    w_dir_main = Path(os.getcwd())
    sdf_rdkit_path = w_dir_main.joinpath('quinine_rdkit')
    smi = 'COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O'
    
    # run CSEARCH conformational sampling, specifying:
    # 1) PATH to create the new SDF files (destination=sdf_rdkit_path)
    # 2) RDKit sampling (program='rdkit')
    # 3) SMILES string (smi=smi)
    # 4) Name for the output SDF files (name='quinine')
    csearch(destination=sdf_rdkit_path,
            smi=smi,name='quinine',program='rdkit')

Step 2: Doing CMIN with the SDF obtained from CSEARCH
.....................................................

.. code:: python

    sdf_cmin_path = w_dir_main.joinpath('quinine_ani')
    sdf_rdkit_files = glob.glob(f'{sdf_rdkit_path}/*.sdf')
    
    # run CMIN refiner, specifying:
    # 1) PATH to create the new SDF files (destination=sdf_cmin_path)
    # 2) SDF files from RDKit (files=sdf_rdkit_files)
    # 3) ANI re-optimization (program='ani')
    cmin(destination=sdf_cmin_path,files=sdf_rdkit_files,program='ani')

Step 3: Writing Gaussian input files with the SDF obtained from CMIN
....................................................................

.. code:: python

    import os, glob
    from pathlib import Path
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    from aqme.cmin import cmin
    
    # set SDF filenames and directory where the new com files will be created
    com_path = w_dir_main.joinpath(f'ani_com_files')
    sdf_cmin_files = glob.glob(f'{sdf_cmin_path}/quinine_rdkit_ani.sdf')
    
    # run QPREP input files generator, with:
    # 1) PATH to create the new SDF files (destination=com_path)
    # 2) Files to convert (files=sdf_cmin_files)
    # 3) QM program for the input (program='gaussian')
    # 4) Keyword line for the Gaussian inputs (qm_input='wb97xd/6-31+G* opt freq')
    # 5) Memory to use in the calculations (mem='24GB')
    # 6) Processors to use in the calcs (nprocs=8)
    qprep(destination=com_path,files=sdf_cmin_files,program='gaussian',
            qm_input='wb97xd/6-31+G* opt freq',mem='24GB',nprocs=8)
