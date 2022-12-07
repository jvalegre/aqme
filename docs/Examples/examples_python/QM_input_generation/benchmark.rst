Preparation of QM input files
=============================

QPREP input files preparation of a benchmarking containing 3 levels of theory for Gaussian and 1 level of theory for ORCA input files
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

.. code:: python

    import os
    from aqme.qprep import qprep
    
    # folder with input log files and their names (*.log to include all the log files in the folder)
    log_files = os.getcwd()+'/log_files/*.log'
    
    # specify a list of lists with level of theory, suffix and program used to generate input files
    # 1) Three levels of theory for Gaussian calculations
    lot_suffix_program = [['wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)','wb97xd','gaussian']]
    lot_suffix_program.append(['m062x/def2qzvpp emp=gd3 scrf=(smd,solvent=acetonitrile)','m062x','gaussian'])
    lot_suffix_program.append(['b3lyp/6-31G*','b3lyp','gaussian'])
    
    # 2) A DLPNO example for ORCA calculations
    ORCA_SP = 'Extrapolate(2/3,cc) def2/J cc-pVTZ/C DLPNO-CCSD(T) NormalPNO TightSCF RIJCOSX GridX7\n'
    ORCA_SP += '%cpcm\n'
    ORCA_SP += 'smd true\n'
    ORCA_SP += 'SMDsolvent \"CH2Cl2\"\n'
    ORCA_SP += 'end\n'
    ORCA_SP += '%method\n'
    ORCA_SP += 'Grid 3\n'
    ORCA_SP += 'FinalGrid 5\n'
    ORCA_SP += 'end\n'
    ORCA_SP += '%scf maxiter 500\n'
    ORCA_SP += 'end\n'
    ORCA_SP += '% mdci\n'
    ORCA_SP += 'Density None\n'
    ORCA_SP += 'end\n'
    ORCA_SP += '% elprop\n'
    ORCA_SP += 'Dipole False\n'
    ORCA_SP += 'end'
    
    lot_suffix_program.append([ORCA_SP,'DLPNO','orca'])
    
    # run the QPREP module, with:
    # 1) Names of the files to get atoms and coordinates (files=log_files)
    # 2) Keyword line(s) used in the inputs (qm_input=level[0])
    # 3) Suffix to add to the file names (suffix=level[1])
    # 4) Program for the input file format (program=level[2])
    # 5) Memory to use in the calculations (mem='4GB')
    # 6) Processors to use in the calcs (nprocs=2)
    for level in lot_suffix_program:
        print(f'o  Creating input files with suffix "{level[1]}" \n')
        qprep(files=log_files, 
              qm_input=level[0], suffix=level[1], program=level[2], mem='4GB', nprocs=2)


