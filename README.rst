.. aqme-banner-start

.. |aqme_banner| image:: ./Logos/AQME_logo.jpg

|aqme_banner|

.. aqme-banner-end

.. badges-start

.. |CircleCI| image:: https://img.shields.io/circleci/build/github/jvalegre/aqme?label=Circle%20CI&logo=circleci
   :target: https://app.circleci.com/pipelines/github/jvalegre/aqme

.. |Codecov| image:: https://img.shields.io/codecov/c/github/jvalegre/aqme?label=Codecov&logo=codecov
   :target: https://codecov.io/gh/jvalegre/aqme

.. |Downloads| image:: https://img.shields.io/conda/dn/conda-forge/aqme?label=Downloads&logo=Anaconda
   :target: https://anaconda.org/conda-forge/aqme

.. |ReadtheDocs| image:: https://img.shields.io/readthedocs/aqme?label=Read%20the%20Docs&logo=readthedocs
   :target: https://aqme.readthedocs.io
   :alt: Documentation Status

|CircleCI|
|Codecov|
|Downloads|
|ReadtheDocs|

.. badges-end

.. checkboxes-start

.. |check| raw:: html

    <input checked=""  type="checkbox">

.. |check_| raw:: html

    <input checked=""  disabled="" type="checkbox">

.. *  raw:: html

    <input type="checkbox">

.. |uncheck_| raw:: html

    <input disabled="" type="checkbox">

.. checkboxes-end

================================================
Automated Quantum Mechanical Environments (AQME)
================================================

.. contents::
   :local:

What is AQME?
-------------

.. introduction-start

The code is an ensemble of automated QM workflows that can be run through 
jupyter notebooks, command lines and yaml files. Some of the most popular 
workflows include:

   *  RDKit- and CREST-based conformer generator leading to 
      ready-to-submit QM input files starting from individual files or SMILES 
      databases  
   *  Post-processing of QM output files to fix convergence errors, 
      extra imaginary frequencies, unfinished jobs, duplicates and error 
      terminations, as well as to detect spin contamination, isomerization issues, 
      and more optimization problems  
   *  Analysis of homogeneity of QM calculations (same level of theory, 
      grid size, program and version, solvation models, etc)  
   *  Generation of xTB, DFT and RDKit descriptors in json and csv files 
      that are ready to use in machine-learning models or used to predict 
      NMR spectra  
   *  More other useful workflows  

Don't miss out the latest hands-on tutorials from our 
`YouTube channel <https://www.youtube.com/channel/UCHRqI8N61bYxWV9BjbUI4Xw>`_  

.. introduction-end

.. installation-start

Installation
------------

Check our `AQME installation in 2 mins <https://youtu.be/VeaBzqIZHbo>`_ video 
for a quick installation guide. In a nutshell, AQME and its dependencies are 
installed as follows:

1. Install using conda-forge  

.. code-block:: shell 
   
   conda install -c conda-forge aqme

2. Install using pip:  

.. code-block:: shell

   pip install aqme

3. Update to the latest version:  

.. code-block:: shell

   pip install aqme --upgrade

Installation of the extra requirements
++++++++++++++++++++++++++++++++++++++

Extra requirements if xTB or CREST are used (compatible with MacOS and Linux only):  

.. code-block:: shell 

   conda install -y -c conda-forge xtb

.. code-block:: shell 

   conda install -y -c conda-forge crest

Extra requirements if `CMIN` is used with ANI models:  

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

*These libraries (except opt) are installed during the initial conda-forge installation.*

Other requirements
++++++++++++++++++

*  RDKit
*  Openbabel
*  xTB
*  CREST

*RDKit and Openbabel are installed during the initial conda-forge installation.*

.. requirements-end

.. workflows-start

Example Workflows
-----------------

The inputs to run pre-defined AQME end-to-end workflows are available in the 
"/Example_workflows/End-to-end_Workflows" folder. Choose the workflow and run the inputs.

Automated protocols for individual modules and tasks are provided in the 
"/Example_workflows" folder inside subfolders with the corresponding module names.

.. workflows-end

.. tests-start

Running the tests
-----------------

Requires the pytest library. 

.. code-block:: shell

   cd path/to/aqme/source/code
   cd tests
   pytest -v

.. tests-end

.. features-modules-start

Features and modules
--------------------

CSEARCH
+++++++

Module on charge of conformational sampling starting from multiple input types
(SMILES, csv, sdf, xyz, etc). Options:

RDKit-based conformational sampling
...................................

Faster sampling, suitable especially for unimolecular systems. Options:  

   *  RDKit standard sampling  
   *  Systematic Unbounded Multiple Minimum search (SUMM)  
   *  FullMonte sampling  

CREST-based conformational sampling
...................................

Slower sampling, suitable for all types of systems (including noncovalent 
complexes and constrained systems such as transition states)

CMIN
++++

Module used to refine conformers generated in CSEARCH through new geometry 
optimizations. Options:  

   *  xTB (GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF, etc.)  
   *  ANI (ANI-1x, ANI-1ccx, ANI-2x, etc.)  

QPREP
+++++

Generator of input files for QM calculations. Options:  

   *  Gaussian  
   *  ORCA  
   *  pySCF (loading parameters in jupyter notebook)  


QCORR
+++++

cclib-based analyzer of output files from multiple QM programs. This module:  

   *  Separates normally terminated files with no errors, extra imaginary 
      frequencies, duplicates, isomerization to other systems and spin contamination  
   *  Automatically generates new com files to "fix" the different issues 
      of the calculations with strategies that are optimal for each type of issue 
      (Gaussian and ORCA)  
   *  Checks that all the calculations are homogeneous (i.e. using the 
      same level of theory, same grid size, same program and version, 
      solvation model, etc)  

QDESCP
++++++

Descriptor generator from multiple input types such as SMILES, log files, xyz, etc. Descriptors generated with:  

   *  RDKit descriptors (i.e. number of polar H, number of aromatic rings, etc)  
   *  xTB (i.e. atomic charges, molecular dipole, solvation energy, etc)  
   *  QM programs (i.e. descriptors from cclib) 

.. features-modules-end

Quickstart
----------

.. quickstart-start

Using AQME in Jupyter Notebooks
+++++++++++++++++++++++++++++++

There are multiple ready-to-use workflows presented as jupyter notebooks in the 
in the aqme repository in 
`Example_Workflows  <https://github.com/jvalegre/aqme/Example_workflows>`__ 
folder. Some examples are: 

  * CSEARCH_CMIN_conformer_generation folder --> CSEARCH/CMIN conformational 
    sampling from SMILES and creation of QM input files  
  * QCORR_processing_QM_outputs --> QCORR analysis of Gaussian output files, 
    generation of JSON files with all the information and creation of new QM input 
    files  
  * QPREP_generating_input_files --> QPREP preparation of input files for 
    Gaussian, ORCA and PySCF from LOG/OUT, SDF and JSON files

.. note::
   
   For a more examples please see the 'Examples/Examples Python' section 
   in our `ReadtheDocs <https://aqme.readthedocs.io>`__ page. 

Using AQME through the command line
+++++++++++++++++++++++++++++++++++

CSEARCH examples
................

Conformer generation with one SMILES and name using RDKit or CREST (use rdkit or crest in --program): 

.. code-block:: shell

   python -m aqme --csearch --program rdkit --smi "CCC" --name proprane

Conformer generation with multiple SMILES and names (i.e. from a database in CSV format):

.. code-block:: shell

   python -m aqme --csearch --program rdkit --input FILENAME.csv

.. note:: 
   
   The csv file must contain the list of SMILES in a column called "SMILES" and 
   the corresponding names in a column called "code_name" 
   (see Example_workflows for more information)

CMIN examples
................

Geometry optimization with xTB or ANI (use xtb or ani in --program; use sdf, xyz, com/gjf or pdb in --files):

.. code-block:: shell

   python -m aqme --cmin --program xtb --files "*.sdf"

QPREP examples
..............

Input file generation from SDF, JSON and LOG/OUT files (replace "\*.sdf" for the corresponding format):

.. code-block:: shell

   python -m aqme --qprep --program gaussian --qm_input "M062x def2tzvp opt freq" --files "*.sdf"

QCORR examples
..............

Analysis of Gaussian output files and JSON file generation:  

.. code-block:: shell

   python -m aqme --qcorr --program gaussian --freq_conv "opt=(calcfc,maxstep=5)" --files "*.log"

.. quickstart-end

Extended documentation
----------------------

More detailed examples, an API reference and the extended list of currently 
avaliable parameters can be found at 
`https://aqme.readthedocs.io <https://aqme.readthedocs.io>`__ 

Developers and help desk
------------------------

.. developers-start 

List of main developers and contact emails:  

*  Juan V. Alegre-Requena [
   `ORCID <https://orcid.org/0000-0002-0769-7168>`__ , 
   `Github <https://github.com/jvalegre>`__ , 
   `email <jv.alegre@csic.es>`__ ]
   main developer of the CSEARCH, CMIN, QCORR, QPREP and QDESCP modules.  
*  Shree Sowndarya S. V. [
   `ORCID <https://orcid.org/0000-0002-4568-5854>`__ , 
   `Github <https://github.com/shreesowndarya>`__ , 
   `email <svss@colostate.edu>`__]
   main developer of the CSEARCH, CMIN, QDESCP and VISMOL modules. 
*  Raúl Pérez-Soto [
   `ORCID <https://orcid.org/0000-0002-6237-2155>`__ ,
   `Github <https://github.com/rperezsoto>`__ ,
   `email <rperezsoto.research@gmail.com>`__ ] 
   worked in refactoring the code and creating the documentation.
*  Turki Alturaifi [
   `webpage <https://www.chem.pitt.edu/person/turki-alturaifi>`__ ,
   `Github <https://github.com/turkiAlturaifi>`__ , 
   `email <tma53@pitt.edu>`__] 
   worked in benchmarking the parameters for RDKit-based conformer generation. 
*  Robert S. Paton [
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

If you use any of the AQME modules, please include this citation:  
AQME v1.4, Alegre-Requena, J. V.; Sowndarya, S.; Alturaifi, T.; Pérez-Soto, R.; 
Paton, R. ChemRxiv 2022, 
DOI: `10.26434/chemrxiv-2022-dnc48 <https://doi.org/10.26434/chemrxiv-2022-dnc48>`__ .  
  
Additionally, please include the corresponding references for the following programs:  
  * If you used CSEARCH with RDKit methods: `RDKit <https://www.rdkit.org/>`__ 
  * If you used CSEARCH with CREST methods: `CREST <https://crest-lab.github.io/crest-docs/>`__ 
  * If you used CMIN with xTB: `xTB <https://xtb-docs.readthedocs.io/en/latest/contents.html>`__ 
  * If you used CMIN with ANI: `ANI <https://github.com/isayev/ASE_ANI>`__ 
  * If you used QCORR: `cclib <https://cclib.github.io/>`__ 
  * If you used QDESCP with xTB: `xTB <https://xtb-docs.readthedocs.io/en/latest/contents.html>`__ 

.. reference-end
