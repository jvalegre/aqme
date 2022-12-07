
Using a csv as input
--------------------

In the following example we will: 

1) Generate various conformers of the quinine molecule using RDKit starting 
   from a .csv file.
2) Minimize the generated conformers using ANI
3) Generate Gaussian input files for each of the conformers

Lets asume we have in our working directory the file 'ML_test.csv' with the 
following contents: 

.. ::

   name,SMILES,
   quinine,COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O,

It is important that it contains the columns 'name' and 'SMILES'. The values 
under the column of name will be used in the filenames of the generated .sdf 
files containing the conformers. 

Step 1: CSEARCH conformational sampling (creates SDF files)
...........................................................

.. code:: python

    import os, glob
    from pathlib import Path
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    
    # Ideal for ML or big data projects, only need to replace smi and name with CSV input
    csv_input = 'ML_test.csv'
    sdf_folder = 'ML_test'
    w_dir_main = Path(os.getcwd())
    sdf_path = w_dir_main.joinpath(sdf_folder)
    
    # create conformers for all the entries in the CSV
    csearch(destination=sdf_path,program='rdkit',input=csv_input)

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
    # 4) Keyword line for the Gaussian inputs (qm_input='wb97xd/6-31+G* opt freq')
    # 5) Memory to use in the calculations (mem='24GB')
    # 6) Processors to use in the calcs (nprocs=8)
    qprep(destination=com_path,files=sdf_rdkit_files,program='gaussian',
            qm_input='wb97xd/6-31+G* opt freq',mem='24GB',nprocs=8)
    
To Remove
.........

Bonus 1: using ORCA instead of Gaussian in QPREP
                                                

.. code:: python

    # Only need to change the qm_input and program options.
    # Multiple lines are allowed. For example, this is the input file of a TS calculation:
    ORCA_input = 'BP86 def2-SVP def2/J\n'
    ORCA_input += '%geom\n'
    ORCA_input += 'Calc_Hess true\n'
    ORCA_input += 'Recalc_Hess 5\n'
    ORCA_input += 'end'
    
    qprep(destination=com_path,files=sdf_rdkit_files,program='orca',
            qm_input=ORCA_input,mem='4GB',nprocs=8)

Bonus 2: using a CSV with many SMILES
                                     

.. code:: python

    import os, glob
    from pathlib import Path
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    
    # Ideal for ML or big data projects, only need to replace smi and name with CSV input
    csv_input = 'ML_test.csv'
    sdf_folder = 'ML_test'
    w_dir_main = Path(os.getcwd())
    sdf_path = w_dir_main.joinpath(sdf_folder)
    
    # create conformers for all the entries in the CSV
    csearch(destination=sdf_path,program='rdkit',input=csv_input)
    
    # set SDF filenames and directory where the new com files will be created
    com_path = sdf_path.joinpath(f'com_files')
    sdf_rdkit_files = glob.glob(f'{sdf_path}/*.sdf')
    
    # create COM files
    qprep(destination=com_path,files=sdf_rdkit_files,program='gaussian',
            qm_input='wb97xd/6-31+G* opt freq',mem='24GB',nprocs=8)

Bonus 3: If you want to use the same functions using a YAML file that stores all the variables
                                                                                              

.. code:: python

    # to load the variables from a YAML file, use the varfile option
    csearch(varfile='FILENAME.yaml')
    
    # for each option, specify it in the YAML file as follows:
    # program='rdkit' --> program: 'rdkit'
    # name='quinine' --> name: 'quinine'
    # etc

Bonus 4: If you want to use the same functions through command lines
                                                                    

.. code:: python

    csearch(destination=sdf_path,smi=smi,name='quinine',program='rdkit')
    
    # for each option, specify it in the command line as follows:
    # program='rdkit' --> --program 'rdkit'
    # name='quinine' --> --name quinine
    # etc
    # for example: python -m aqme --program rdkit --smi COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O --name quinine

