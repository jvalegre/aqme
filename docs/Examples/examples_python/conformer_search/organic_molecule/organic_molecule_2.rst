Using Fullmonte
---------------

In the following example we will: 

1) Generate various conformers of the quinine molecule using Fullmonte
2) Generate Gaussian input files for each of the conformers

Step 1: CSEARCH conformational sampling (creates SDF files)
...........................................................

.. code:: python

    import os, glob
    from pathlib import Path                                                                                                                                                          
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    
    # set working directory and SMILES string
    w_dir_main = Path(os.getcwd())
    sdf_path = w_dir_main.joinpath('quinine')
    smi = 'COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O'
    
    # run CSEARCH conformational sampling, specifying:
    # 1) PATH to create the new SDF files (destination=sdf_path)
    # 2) SMILES string (smi=smi)
    # 3) Name for the output SDF files (name='quinine')
    # 4) Fullmonte sampling (program='fullmonte')
    csearch(destination=sdf_path,
            smi=smi,name='quinine',program='fullmonte')

Step 2: Writing Gaussian input files with the SDF obtained from CSEARCH
.......................................................................

.. code:: python

    # set SDF filenames and directory where the new com files will be created
    com_path = sdf_path.joinpath(f'com_files')
    sdf_fullmonte_files = glob.glob(f'{sdf_path}/*.sdf')
    
    # run QPREP input files generator, with:
    # 1) PATH to create the new SDF files (destination=com_path)
    # 2) Files to convert (files=sdf_fullmonte_files)
    # 3) QM program for the input (program='gaussian')
    # 4) Keyword line for the Gaussian inputs (qm_input='wb97xd/6-31+G* opt freq')
    # 5) Memory to use in the calculations (mem='24GB')
    # 6) Processors to use in the calcs (nprocs=8)
    qprep(destination=com_path,files=sdf_fullmonte_files,program='gaussian',
            qm_input='wb97xd/6-31+G* opt freq',mem='24GB',nprocs=8)

ToRemove
........

Bonus 1: If you want to use the same functions using a YAML file that stores all the variables
                                                                                              

.. code:: python

    # to load the variables from a YAML file, use the varfile option
    csearch(varfile='FILENAME.yaml')
    
    # for each option, specify it in the YAML file as follows:
    # program='fullmonte' --> program: 'fullmonte'
    # name='quinine' --> name: 'quinine'
    # etc

Bonus 2: If you want to use the same functions through command lines
                                                                    

.. code:: python

    csearch(destination=sdf_path,
            smi=smi,name='quinine',program='fullmonte')
    
    # for each option, specify it in the command line as follows:
    # program='fullmonte' --> --program 'fullmonte'
    # name='quinine' --> --name quinine
    # etc
    # for example: python -m aqme --program fullmonte --smi COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O --name quinine

