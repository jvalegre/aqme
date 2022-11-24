.. aqme-banner-start

.. |aqme_banner| image:: ./Logos/AQME_logo.jpg

|aqme_banner|

.. aqme-banner-end

.. badges-start

.. |CircleCI| image:: https://img.shields.io/circleci/build/github/jvalegre/aqme?label=Circle%20CI&logo=circleci
   :target: https://app.circleci.com/pipelines/github/jvalegre/aqme

.. |Codecov| image:: https://img.shields.io/codecov/c/github/jvalegre/aqme?label=Codecov&logo=codecov
   :target: https://codecov.io/gh/jvalegre/aqme

.. |CodeFactor| image:: https://img.shields.io/codefactor/grade/github/jvalegre/aqme/master?label=Codefactor%20grade&logo=Codefactor
   :target: https://www.codefactor.io/repository/github/jvalegre/aqme/overview/master

.. |Codacy| image:: https://img.shields.io/codacy/grade/3a4cc7c7705e46129c7ea0fca58af846?label=Codacy%20grade&logo=Codacy
   :target: https://www.codacy.com/gh/jvalegre/aqme/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jvalegre/aqme&amp;utm_campaign=Badge_Grade

.. |lgtm| image:: https://img.shields.io/lgtm/grade/python/github/jvalegre/aqme?label=LGTM%20grade&logo=lgtm 
   :target: https://lgtm.com/projects/g/jvalegre/aqme/context:python

|CircleCI|
|Codecov|
|CodeFactor|
|Codacy|
|lgtm|

.. badges-end

.. checkboxes-start

.. |check| raw:: html

    <input checked=""  type="checkbox">

.. |check_| raw:: html

    <input checked=""  disabled="" type="checkbox">

.. |uncheck| raw:: html

    <input type="checkbox">

.. |uncheck_| raw:: html

    <input disabled="" type="checkbox">

.. checkboxes-end

================================================
Automated Quantum Mechanical Environments (AQME)
================================================


What is AQME?
-------------

.. introduction-start

The code is an ensemble of automated QM workflows that can be run through jupyter notebooks, command lines and yaml files. Some of the most popular workflows include:  
   |uncheck| RDKit- and CREST-based conformer generator leading to 
   ready-to-submit QM input files starting from individual files or SMILES 
   databases  

   |uncheck| Post-processing of QM output files to fix convergence errors, 
   extra imaginary frequencies, unfinished jobs, duplicates and error 
   terminations, as well as to detect spin contamination, isomerization issues, 
   and more optimization problems  

   |uncheck| Analysis of homogeneity of QM calculations (same level of theory, 
   grid size, program and version, solvation models, etc)  

   |uncheck| Generation of xTB, DFT and RDKit descriptors in json and csv files 
   that are ready to use in machine-learning models or used to predict NMR spectra  

   |uncheck| More other useful workflows  

Don't miss out the latest hands-on tutorials from our 
`YouTube channel <https://www.youtube.com/channel/UCHRqI8N61bYxWV9BjbUI4Xw>`_  

.. introduction-end

.. installation-start

Installation
------------

Check our `AQME installation in 2 mins <https://youtu.be/VeaBzqIZHbo>`_ video 
for a quick installation guide. In a nutshell, AQME and its dependencies are 
installed as follows:

1. Using conda-forge (fastest, one-command install) 

.. code-block:: shell 
   
   conda install -c conda-forge aqme

2. Using the source code (Latest version, recommended)

.. code-block:: shell

   python -m pip install .

3. Using the Python Package Index (pip)

.. code-block:: shell 

   python -m pip install aqme

Installation of the extra requirements
++++++++++++++++++++++++++++++++++++++

If the installation was carried out using pip: 

.. code-block:: shell

   conda install -c conda-forge rdkit openbabel

If the `cmin` module with torchani will be used (torch-related dependencies): 

.. code-block:: shell 

   pip install torch torchvision torchani

.. warning:: *Known incompatibilities:*
   
   -  RDKit cannot be installed through `pip install rdkit` in Windows when 
      Anaconda prompts are used

.. installation-end 

.. requirements-start

Requirements
------------

Python and Python libraries
+++++++++++++++++++++++++++

*  Python >= 3.7
*  pandas
*  Numpy
*  PyYAML
*  progress
*  ase (Atomic Simulation Environment)
*  cclib (Computational Chemistry Library)
*  cffi
*  matplotlib 
*  seaborn
*  goodvibes
*  (opt) torch, torchvision and torchani

Other requirements
++++++++++++++++++

*  RDKit
*  Openbabel
*  XTB
*  CREST

.. requirements-end

.. workflows-start

Example Workflows
-----------------

The inputs to run pre-defined AQME end-to-end workflows are available in the 
"/Example_workflows/End-to-end_Workflows" folder. Choose the workflow and run the inputs.

Automated protocols for individual modules and tasks are provided in the 
/Example_workflows/ folder inside subfolders with the corresponding module names.

.. workflows-end

.. tests-start

Running the tests
-----------------

Requires the pytest library. 

.. code-block:: shell

   cd path/to/aqme/source/code
   cd tests
   pytest --v

.. tests-end

.. features-modules-start

Features and modules
--------------------

csearch
+++++++

Module on charge of conformational sampling starting from multiple input types (SMILES, csv, sdf, xyz, etc). Options:

RDKit-based conformational sampling
...................................

Faster sampling, suitable especially for unimolecular systems. Options:  

   |uncheck| RDKit standard sampling  
   
   |uncheck| Systematic Unbounded Multiple Minimum search (SUMM)  
   
   |uncheck| FullMonte sampling  

CREST-based conformational sampling
...................................

Slower sampling, suitable for all types of systems (including noncovalent 
complexes and constrained systems such as transition states)

cmin
++++

Module used to refine conformers generated in CSEARCH through new geometry 
optimizations. Options:  

   |uncheck| xTB (GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF)  

   |uncheck| ANI (ANI-1x, ANI-1ccx, ANI-2x)  

qprep
+++++

Generator of input files for QM calculations. Options:  

   |uncheck| Gaussian  

   |uncheck| ORCA  

   |uncheck| pySCF (loading parameters in jupyter notebook)  


qcorr
+++++

cclib-based analyzer of output files from multiple QM programs. This module:  

   |uncheck| Separates normally terminated files with no errors, extra imaginary 
   frequencies, duplicates, isomerization to other systems and spin contamination  

   |uncheck| Automatically generates new com files to "fix" the different issues 
   of the calculations with strategies that are optimal for each type of issue 
   (Gaussian and ORCA)  

   |uncheck| Checks that all the calculations are homogeneous (i.e. using the 
   same level of theory, same grid size, same program and version, 
   solvation model, etc)  

qdescp
++++++

Descriptor generator from multiple input types such as SMILES, log files, xyz, etc. Descriptors generated with:  

   |uncheck| RDKit descriptors (i.e. number of polar H, number of aromatic rings, etc)  

   |uncheck| xTB (i.e. atomic charges, molecular dipole, solvation energy, etc)  

   |uncheck| QM programs (i.e. descriptors from cclib) 

.. features-modules-end

Quickstart
----------

.. quickstart-start

Using AQME in Jupyter Notebooks
+++++++++++++++++++++++++++++++

There are multiple ready-to-use workflows presented as jupyter notebooks in the 
'Example workflows' folder. Some examples are: 

  * CSEARCH_CMIN_conformer_generation folder --> CSEARCH/CMIN conformational 
    sampling from SMILES and creation of QM input files  

  * QCORR_processing_QM_outputs --> QCORR analysis of Gaussian output files, 
    generation of JSON files with all the information and creation of new QM input 
    files  

  * QPREP_generating_input_files --> QPREP preparation of input files for 
    Gaussian, ORCA and PySCF from LOG/OUT, SDF and JSON files

Using AQME through the command line
+++++++++++++++++++++++++++++++++++

csearch examples
................

Conformer generation with one SMILES and name: 

.. code-block:: shell

   python -m aqme --csearch --program rdkit --smi CCC --name proprane

Conformer generation with multiple SMILES and names:

.. code-block:: shell 

   python -m aqme --csearch --program rdkit --input FILENAME.csv

.. note:: 
   
   The csv file must contain the list of SMILES in a column called "SMILES" and 
   the corresponding names in a column called "code_name" 
   (see Example_workflows for more information)

Conformer generation using a YAML file containing constraints:

.. code-block:: shell

   python -m aqme --varfile FILENAME.yaml


The YAML file must contain the following parameters 


::

   input : 'smi.csv' #name of input
   output_name : 'csearch' #name for output
   csearch : True #activate CSEARCH
   program : 'rdkit' #program used in CSEARCH


qcorr example
.............

analysis of Gaussian output files and json file generation:  

.. code-block:: shell

   python -m aqme --qcorr --program gaussian --freq_conv "opt=(calcfc,maxstep=5)" --files=*.log


qprep examples
..............

Input file generation from SDF files (coming from CSEARCH for example):  

.. code-block:: shell

   python -m aqme --qprep --program gaussian --qm_input "M062x def2tzvp opt freq" --files *.sdf


Input file generation from last geometry of output files (log or out files):  

.. code-block:: shell

   python -m aqme --qprep --program gaussian--qm_input "M062x def2tzvp opt freq" --files *.log --suffix M062X


Input file generation from json files:  

.. code-block:: shell

   python -m aqme --qprep --program orca --qm_input "BP86 def2-SVP def2/J" --files *.json --suffix BP86

.. quickstart-end

Extended documentation
----------------------

** ReadTheDocs page in process **

.. developers-start

Developers and help desk
------------------------

List of main developers and contact emails:  

|uncheck| Shree Sowndarya S. V. [
`ORCID <https://orcid.org/0000-0002-4568-5854>`__ , 
`Github <https://github.com/shreesowndarya>`__ , 
`email <svss@colostate.edu>`__]
main developer of the CSEARCH, CMIN, QDESCP and VISMOL modules. 

|uncheck| Juan V. Alegre-Requena [
`ORCID <https://orcid.org/0000-0002-0769-7168>`__ , 
`Github <https://github.com/jvalegre>`__ , 
`email <jvalegre@unizar.es>`__ ]
main developer of the QCORR and QPREP modules.   

|uncheck| Turki Alturaifi [
`webpage <https://www.chem.pitt.edu/person/turki-alturaifi>`__ ,
`Github <https://github.com/turkiAlturaifi>`__ , 
`email <turki0@rams.colostate.edu>`__] 
worked in benchmarking the parameters for RDKit-based conformer generation. 

|uncheck| Raúl Pérez-Soto [
`ORCID <https://orcid.org/0000-0002-6237-2155>`__ ,
`Github <https://github.com/rperezsoto>`__ ,
`email <rperezsoto.research@gmail.com>`__ ] 
worked in refactoring the code.

|uncheck| Robert S. Paton [
`ORCID <https://orcid.org/0000-0002-0104-4166>`__ ,
`Github <https://github.com/bobbypaton>`__ , 
`email <robert.paton@colostate.edu>`__]
research group supervisor and code advisor.

For suggestions and improvements of the code (greatly appreciated!), please 
reach out through the issues and pull requests options of Github.

.. developers-end

License
-------

.. license-start 

AQME is freely available under an `MIT License <https://opensource.org/licenses/MIT>`_  

.. license-end

Reference
---------

.. reference-start

AQME v1.3, Alegre-Requena, J. V.; Sowndarya, S.; Pérez-Soto, R.; Alturaifi, T. M.; 
Paton, R. S., 2022. https://github.com/jvalegre/aqme  

.. reference-end
