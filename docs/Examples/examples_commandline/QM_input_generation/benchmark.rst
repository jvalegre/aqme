.. |QPREP_scheme| image:: ../../images/QPREP_scheme.png
   :width: 500

Benchmarking methods 
====================

In this example we provide the steps needed to use AQME to generate the input 
for different levels of theory and softwares which is a common occurrence in 
method benchmarking.

.. centered:: |QPREP_scheme|

For this example we are going to asume that we have a folder named 'log_files' 
with the .log files shown in the scheme corresponding to gaussian output files 
of optimization and frequency calculations. Our aim is to generate the input 
files on the right hand side of the scheme to run SP calculations. 

First we are going to generate the SP inputs for the wb97xd functional using 
gaussian: 

.. code:: shell

    python -m aqme --qprep --mem 4GB --nprocs 2 --program gaussian --suffix wb97xd --files "log_files/*.log" --qm_input "wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)"

Then we proceed to generate the SP inputs for the m062x functional:

.. code:: shell

    python -m aqme --qprep --mem 4GB --nprocs 2 --program gaussian --suffix m062x --files "log_files/*.log" --qm_input "m062x/def2qzvpp emp=gd3 scrf=(smd,solvent=acetonitrile)"

Then we move to the b3lyp functional. 

.. code:: shell

    python -m aqme --qprep --mem 4GB --nprocs 2 --program gaussian --suffix b3lyp --files "log_files/*.log" --qm_input "b3lyp/6-31G*"

And finally we proceed to generate the SP inputs for orca. 

.. code:: shell

   python -m aqme --qprep --mem 4GB --nprocs 2 --program orca --suffix DLPNO --files "log_files/*.log" --qm_input "Extrapolate(2/3,cc) def2/J cc-pVTZ/C DLPNO-CCSD(T) NormalPNO TightSCF RIJCOSX GridX7
   %cpcm
   smd true
   SMDsolvent \"CH2Cl2\"
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
   end"

.. note:: 

   See how to generate the four different inputs we only need to 
   change the :code:`program`, :code:`qm_input` and the :code:`suffix` but most 
   of the command line can be re-used.

