.. |Strychnine| image:: ../../images/strychnine_chemdraw.png
   :width: 400


Strychnine
==========

Along the steps of this example workflow we will show how to: 

i)   Generate different conformers of the molecule using `csearch`
ii)  Generate the inputs for the QM geometry optimization
iii) Fix error terminations and imaginary frequencies of the output files
iv)  Calculation and analyze the NMR chemical shifts for the conformers
     generated.
v)   Use GoodVibes to calculate the Boltzmann distributions using Gibbs free
     energies at 298.15 K

Specifically, in this example we will calculate the NMR chemical shifts of the strychnine
starting from the smiles representation of said molecule that we can see below. 

+---------------------------------------------------------------------------------------+
|                         .. centered:: **SMILES**                                      |
+---------------------------------------------------------------------------------------+
| .. centered:: C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75  |
+---------------------------------------------------------------------------------------+
|                      .. centered::  |Strychnine|                                      |
+---------------------------------------------------------------------------------------+

.. note::

   A jupyter notebook containing all the steps shown in this example can be found 
   in the aqme repository in `Github  <https://github.com/jvalegre/aqme>`__ or in 
   `Figshare <https://figshare.com/articles/dataset/AQME_paper_examples/20043665/11>`__

.. contents:: Steps
   :local:


Step 1: Importing AQME and other python modules
-----------------------------------------------

.. code:: python

    import os, subprocess,shutil, glob
    from pathlib import Path 
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    from aqme.qcorr import qcorr
    from aqme.qdescp import qdescp

Step 2: CSEARCH conformational sampling
---------------------------------------

.. code:: python

    name = 'Strychnine'
    smi = 'C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75'
    program = 'rdkit'
    
    # folder where the SDF files are generated
    sdf_path = f'{os.getcwd()}/{name}_sdf-files' 
    
    csearch(destination=sdf_path,
            program=program,
            smi=smi,
            name=name)

Step 3: Creating Gaussian input files for optimization and frequency with QPREP
-------------------------------------------------------------------------------

.. code:: python

    program = 'gaussian'
    qm_input = 'B3LYP/6-31+G(d,p) opt freq'
    mem='24GB'
    nprocs=12
    
    # SDF files from Step 2
    sdf_rdkit_files = f'{sdf_path}/*.sdf' 

    # folder where the COM files are generated
    com_path = f'{os.getcwd()}/{name}_com-files' 
    
    qprep(destination=com_path,
          files=sdf_rdkit_files,
          program=program,
          qm_input=qm_input,
          mem=mem,
          nprocs=nprocs)

Step 4: Running Gaussian inputs for optimization and frequency calcs externally
-------------------------------------------------------------------------------

Now that we have generated our gaussian input files (in the com_path location 
of Step 3) we need to run the gaussian calculations. If you do not know how to 
run the Gaussian calculations in your HPC please refer to your HPC manager. 

As an example, for a single calculation in Gaussian 16 through the terminal we 
would run the following command on a Linux-based system: 

.. code:: shell

    g16 myfile.com



Step 5: QCORR analysis including isomerization filter
-----------------------------------------------------

.. code:: python

    log_files=f'{com_path}/*.log' # LOG files from Step 4
    
    qcorr(files=log_files,
          freq_conv='opt=(calcfc,maxstep=5)',
          isom_type='com',
          isom_inputs=com_path,
          nprocs=24,
          mem='96GB')

Step 6: Resubmission of unsuccessful calculations (if any) with suggestions from AQME
-------------------------------------------------------------------------------------

Now we need to run the generated COM files (in fixed_inp_folder) with Gaussian 
like we did in Step 4

Step 7: Creating Gaussian input files for NMR calcs with QPREP
--------------------------------------------------------------

.. code:: python

    program = 'gaussian'
    qm_input = 'B3LYP/6-311+G(2d,p) scrf=(solvent=chloroform,smd) nmr=giao'
    mem='24GB'
    nprocs=12
    
    # folder where the successful LOG files are stored during the QCORR cycles 
    # (Steps 5 and 6)

    success_folder = com_path+'/success' 
    log_files = f'{success_folder}/*.log'

    # folder to store the new COM inputs for single point NMR calcs
    sp_path = f'{os.getcwd()}/{name}_sp-files'
    
    qprep(w_dir_main=success_folder,
          destination=sp_path,
          files=log_files,
          program=program,
          qm_input=qm_input,
          mem=mem,
          nprocs=nprocs,
          suffix='SP')

Step 8: Running Gaussian NMR calcs
----------------------------------

Now we need to run the generated COM files (in sp_path) with Gaussian 
like we did in Step 4

Step 9: Obtaining Boltzmann weighted NMR shifts with QDESCP
-----------------------------------------------------------

.. code:: python

    # Create JSON files with QCORR to store the information from the resulting LOG files
    log_files=f'{sp_path}/*.log'
    qcorr(files=log_files)
    
    # Analyze the JSON files to calculate the Boltzmann averaged shielding tensors

    ## folder where the JSON files were just created with QCORR
    json_folder = sp_path+'/success/SP_calcs/json_files'
    json_files=f'{json_folder}/*.json'

    ## folder to store the results from QDESCP
    nmr_path = f'{os.getcwd()}/{name}_nmr-files' 
    
    qdescp(program='nmr',
           boltz=True,
           files=json_files,
           destination=nmr_path,
           nmr_slope=[-1.0537, -1.0784],
           nmr_intercept=[181.7815,31.8723], 
           nmr_experim='Experimental_NMR_shifts.csv')

Step 10: Calculating conformer populations with GoodVibes
---------------------------------------------------------

.. code:: python

    log_files = glob.glob(f'{success_folder}/*.log')
    
    w_dir_main  = Path(os.getcwd())
    GV_folder = w_dir_main.joinpath('Strychine_GoodVibes-analysis')
    GV_folder.mkdir(exist_ok=True, parents=True)
    
    for file in log_files:
    	shutil.copy(file, GV_folder)
    
    # run GoodVibes
    os.chdir(GV_folder)
    subprocess.run(['python', '-m', 'goodvibes', '--xyz','-c','1', '*.log','--boltz'])
    os.chdir(w_dir_main)


