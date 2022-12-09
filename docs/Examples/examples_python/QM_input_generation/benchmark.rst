.. |QPREP_scheme| image:: ../../images/QPREP_scheme.png
   :width: 500

Benchmarking methods 
====================

In this example we provide the steps needed to use aqme to generate the input 
for different levels of theory and softwares which is a common occurrence in 
method benchmarking.

.. centered:: |QPREP_scheme|

For this example we are going to asume that we have a folder named 'log_files' 
with several .log files corresponding to gaussian output files of optimization 
and frequency calculations. 

First we start by importing the modules

.. code:: python

    from pathlib import Path
    from aqme.qprep import qprep


Next we list all the files that contain the geometries 

.. code:: python

    log_files = [str(filepath for filepath in Path('log_files').glob('/*.log')]

We now are going to create a list with the different calculations that we want 
to use. 

.. code:: python 

   # We are going to use namedtuples to make the code easier to understand
   from collections import namedtuple 

   Calculation = namedtuple('Calculation','qm_input suffix program')

   inp_1 = Calculation(qm_input='wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)',
                       suffix='wb97xd',
                       program='gaussian')
   inp_2 = Calculation(qm_input='m062x/def2qzvpp emp=gd3 scrf=(smd,solvent=acetonitrile)',
                       suffix='m062x',
                       program='gaussian')
   inp_3 = Calculation(qm_input='b3lyp/6-31G*',
                       suffix='b3lyp',
                       program='gaussian')

   orca_qm_inp = r'''
   Extrapolate(2/3,cc) def2/J cc-pVTZ/C DLPNO-CCSD(T) NormalPNO TightSCF RIJCOSX GridX7
   %cpcm
   smd true
   SMDsolvent "CH2Cl2"
   end
   %method
   Grid 3
   FinalGrid 5
   end
   %scf maxiter 500
   end
   % mdci
   Density None
   end
   % elprop
   Dipole False
   end'''.lstrip()

   inp_4 = Calculation(qm_input=orca_qm_inp,
                       suffix='DLPNO',
                       program='orca')

Finally we proceed to generate all the input files: 

.. code:: 
   
    calculations = [inp_1,inp_2,inp_3,inp_4]

   for qm_inp,suf,prog in calculations: 
      print(f'o  Creating input files with suffix "{suf}" \n')
      qprep(files=log_files, 
            qm_input=qm_inp, 
            suffix=suf, 
            program=prog, 
            mem='4GB', 
            nprocs=2)

