.. |QCORR_scheme| image:: ../images/QCORR_scheme.png
   :width: 600

=======================
Correction of QM Output
=======================

The QCORR module focuses on the analysis and correction of the output files of 
QM calculations. Here we refer to correction as: 

*  Generate new inputs from calculations that have an error termination. 
*  Generate new inputs for minima containing a small imaginary frequency
*  Ensure that all provided files have the same level of theory, grid size, 
   program, version, etc.

The following scheme shows how QCORR works and how it sorts the calculations.

|QCORR_scheme|


Analyzing the output files
--------------------------

We start by importing the modules. 

.. code:: python

    from pathlib import Path
    from aqme.qcorr import qcorr

Then we list the files that we want to analyze. In this case we are going to 
analyze Gaussian16 output files. We are going to assume that we have our 
files in the folder 'calculations'

.. code:: python 

    folder = Path('calculations')
    files=[str(filepath) for filepath in folder.glob('*.log')]

Finally we run the analysis of the files.

.. code:: python

    qcorr(files=files,
          freq_conv='opt=(calcfc,maxstep=5)')

Here we specify the `freq_conv` which will attempt to fix calculations whose 
optimization ended normally but whose frequency calculation did not.

Optionally we may include the extension of the initial input files and their 
location: 

.. code:: python

    qcorr(files=files,
          freq_conv='opt=(calcfc,maxstep=5)',
          isom_type='com',         # Extension of the initial input files
          isom_inputs=folder)      # Folder with the initial input files

If we want to check the .json files that we have generated in a separatedly we 
can use the full_check function which will only check if all calculations were 
done using the same theory level, grid size, program and version, solvation ...

.. code:: python

    from aqme.qcorr import full_check

    files = [str(filepath) for filepath in (folder/'json_files').glob('*.json')]
    full_check(destination_fullcheck=folder,files=files)

If we instead wanted to skip the checks and generate the .json files containing 
information about our calculations we can use the `fullcheck` keyword.

.. code:: python

    qcorr(files=files,
          fullcheck=False)
