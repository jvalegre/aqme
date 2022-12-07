Enforce Charge and Multiplicity
===============================

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

