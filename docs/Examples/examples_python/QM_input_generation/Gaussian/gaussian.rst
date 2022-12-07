Generate Gaussian Inputs
========================

Enforce Charge and Multiplicity
-------------------------------

.. code:: python

    # calculate reduction potential from sdf files
    import os
    from aqme.qprep import qprep
    
    # folder with input sdf files and their names (*.sdf to include all the sdf files in the folder)
    sdf_files = os.getcwd()+'/sdf_files/*.sdf'
    
    
    # run the QPREP module, with:
    # 1) Names of the files to get atoms and coordinates (files=sdf_files)
    # 2) Set charge for the input files (charge=-1)
    # 3) Set multiplicity for the input files (mult=2)
    # 4) Keyword line(s) used in the inputs (qm_input='wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)')
    # 5) Suffix to add to the file names (suffix='wb97xd-reduced')
    # 6) Program for the input file format (program='gaussian')
    # 7) Memory to use in the calculations (mem='16GB')
    # 8) Processors to use in the calcs (nprocs=8)
    print(f'o  Creating input files with suffix "wb97xd-reduced" \n')
    qprep(files=sdf_files, charge=-1, mult=2,
          qm_input='wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)', suffix='wb97xd-reduced',
          program='gaussian', mem='16GB', nprocs=8)

Gaussian Inputs with genecp
---------------------------

.. code:: python

    # do a single point with genECP included
    
    import os
    from aqme.qprep import qprep
    
    # folder with input json files and their names (*.json to include all the json files in the folder)
    json_files = os.getcwd()+'/json_files/*.json'
    
    # run the QPREP module, with:
    # 1) Names of the files to get atoms and coordinates (files=json_files)
    # 2) Basis set to use in the atoms included in genECP (bs_gen='def2svp')
    # 3) Basis set to use in the rest of the atoms (bs='6-31G*')
    # 4) Atoms to include as genECP (gen_atoms=['C'])
    # 5) Keyword line(s) used in the inputs (qm_input='wb97xd/genecp scrf=(smd,solvent=acetonitrile)')
    # 6) Suffix to add to the file names (suffix='wb97xd-genecp')
    # 7) Program for the input file format (program='gaussian')
    # 8) Memory to use in the calculations (mem='16GB')
    # 9) Processors to use in the calcs (nprocs=8)
    print(f'o  Creating input files with suffix "wb97xd-genecp"\n')
    qprep(files=json_files,
                bs_gen='def2svp', bs_nogen='6-31G*', gen_atoms=['C'],
                qm_input='wb97xd/genecp scrf=(smd,solvent=acetonitrile)', suffix='wb97xd-genecp',
                program='gaussian', mem='16GB', nprocs=8)

Gaussian Inputs with nbo
------------------------

.. code:: python

    # do an NBO calc with a final line after coords (requiring Wyberg bond orders)
    
    import os
    from aqme.qprep import qprep
    
    # folder with input json files and their names (*.json to include all the json files in the folder)
    json_files = os.getcwd()+'/json_files/*.json'
    
    # run the QPREP module, with:
    # 1) Names of the files to get atoms and coordinates (files=json_files)
    # 2) Final line after the coordinates section (qm_end='$nbo bndidx $end')
    # 3) Keyword line(s) used in the inputs (qm_input='wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)')
    # 4) Suffix to add to the file names (suffix='wb97xd-nbo')
    # 5) Program for the input file format (program='gaussian')
    # 6) Memory to use in the calculations (mem='16GB')
    # 7) Processors to use in the calcs (nprocs=8)
    print(f'o  Creating input files with suffix "wb97xd-nbo" \n')
    qprep(files=json_files, qm_end='$nbo bndidx $end',
                qm_input='pop=(nbo6read,savenbos) wb97xd/def2svp', suffix='wb97xd-nbo',
                program='gaussian', mem='16GB', nprocs=8)