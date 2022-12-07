Gaussian Inputs with nbo
===========================

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

