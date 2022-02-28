![](Logos/AQME_logo.jpg)
[![Build Status](https://img.shields.io/travis/com/jvalegre/aqme?label=Linux%20CI&logo=Travis)](https://travis-ci.com/github/jvalegre/aqme)
[![Tests](https://img.shields.io/static/v1?label=Tests&message=104&color=green&logo=Travis)](https://travis-ci.com/github/jvalegre/aqme)
[![Codecov](https://img.shields.io/codecov/c/github/jvalegre/aqme?label=Codecov&logo=codecov)](https://codecov.io/gh/jvalegre/aqme)

[![CodeFactor](https://img.shields.io/codefactor/grade/github/jvalegre/aqme?label=Codefactor%20grade&logo=codefactor)](https://www.codefactor.io/repository/github/jvalegre/aqme/overview/master)
[![Codacy](https://img.shields.io/codacy/grade/3a4cc7c7705e46129c7ea0fca58af846?label=Codacy%20grade&logo=Codacy)](https://www.codacy.com/gh/jvalegre/aqme/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jvalegre/aqme&amp;utm_campaign=Badge_Grade)
[![lgtm](https://img.shields.io/lgtm/grade/python/github/jvalegre/aqme?label=LGTM%20grade&logo=lgtm)](https://lgtm.com/projects/g/jvalegre/aqme/context:python)
[![CodeClimate](https://img.shields.io/codeclimate/maintainability-percentage/jvalegre/aqme?label=Code%20climate%20maintainability&logo=code%20climate)](https://codeclimate.com/github/jvalegre/aqme)

![Size](https://img.shields.io/github/languages/code-size/jvalegre/aqme)
[![PyPI](https://img.shields.io/pypi/v/aqme?color=blue&label=PyPI&logo=pypi)](https://pypi.org/project/aqme)

___
# <p align="center">Automated Quantum Mechanical Environments (AQME)</p>
## <p align="center">-- Table of contents --</p>
### <p align="center">[What is AQME?](https://github.com/jvalegre/aqme#what-is-aqme) &nbsp; &nbsp; [Installation](https://github.com/jvalegre/aqme#installation) &nbsp; &nbsp; [Requirements](https://github.com/jvalegre/aqme#requirements)</p>
### <p align="center">[Features and modules](https://github.com/jvalegre/aqme#features-and-modules) &nbsp; &nbsp; [Quickstart](https://github.com/jvalegre/aqme#quickstart) &nbsp; &nbsp; [Extended documentation](https://github.com/jvalegre/aqme#extended-documentation-installation-use-examples-etc)</p>
### <p align="center">[Developers and help desk](https://github.com/jvalegre/aqme#developers-and-help-desk) &nbsp; &nbsp; [License](https://github.com/jvalegre/aqme#license) &nbsp; &nbsp; [Reference](https://github.com/jvalegre/aqme#reference)</p>
___
## What is AQME?
The code is an ensemble of automated QM workflows that can be run through jupyter notebooks, command lines and yaml files. Some of the most popular workflows include:
  - [ ] RDKit- and CREST-based conformer generator leading to ready-to-submit QM input files starting from individual files or SMILES databases
  - [ ] Post-processing of QM output files to fix convergence errors, extra imaginary frequencies, unfinished jobs, duplicates and error terminations, as well as to detect spin contamination, isomerization issues, and more optimization problems
  - [ ] Analysis of homogeneity of QM calculations (same level of theory, grid size, program and version, solvation models, etc)
  - [ ] Generation of xTB, DFT and RDKit descriptors in json and csv files that are ready to use in machine-learning models or used to predict NMR spectra
  - [ ] More other useful workflows!

## Installation
AQME and its dependencies are installed as follows:
1. Installing RDKit and Open Babel through conda. For shortcuts:
  * `conda install -c rdkit rdkit` or `conda install -c conda-forge rdkit` (compatible with Python >=3.8)
  * `conda install -c conda-forge openbabel`
2. Install AQME and its dependencies:
  * `pip install aqme` or `python setup.py install`

## Requirements
* Python 3
* RDKit and Open Babel
* Dependencies installed with setup.py (automatic install)

## Features and modules
### CSEARCH
Module on charge of conformational sampling starting from multiple input types (SMILES, csv, sdf, xyz, etc). Options:
#### RDKit-based conformational sampling
Faster sampling, suitable especially for unimolecular systems. Options:
  - [ ] RDKit standard sampling
  - [ ] Systematic Unbounded Multiple Minimum search (SUMM)
  - [ ] FullMonte sampling
#### CREST-based conformational sampling
Slower sampling, suitable for all types of systems (including noncovalent complexes and constrained systems such as transition states)

### CMIN
Module used to refine conformers generated in CSEARCH through new geometry optimizations. Options:
  - [ ] xTB (GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF)
  - [ ] ANI (ANI-1x, ANI-1ccx, ANI-2x)

### QPREP
Generator of input files for QM calculations. Options:
  - [ ] Gaussian
  - [ ] ORCA
  - [ ] pySCF (loading parameters in jupyter notebook)

### QCORR
cclib-based analyzer of output files from multiple QM programs. This module:
  - [ ] Separates normally terminated files with no errors, extra imaginary frequencies, duplicates, isomerization to other systems and spin contamination
  - [ ] Automatically generates new com files to "fix" the different issues of the calculations with strategies that are optimal for each type of issue (Gaussian and ORCA)
  - [ ] Checks that all the calculations are homogeneous (i.e. using the same level of theory, same grid size, same program and version, solvation model, etc)

### QDESC
Descriptor generator from multiple input types such as SMILES, log files, xyz, etc. Descriptors generated with:
  - [ ] RDKit descriptors (i.e. number of polar H, number of aromatic rings, etc)
  - [ ] xTB (i.e. atomic charges, molecular dipole, solvation energy, etc)
  - [ ] QM programs (i.e. descriptors from cclib)

## Quickstart
### Using AQME in Jupyter Notebooks
There are multiple ready-to-use workflows presented as jupyter notebooks in the 'Example_workflows.zip' file. Some examples are:
  * CSEARCH_CMIN_workflows.ipynb (CSEARCH/CMIN conformational sampling from SMILES and creation of QM input files):
    1. RDKit-based automated generation of quinine conformers, followed by geometry reoptimization with xTB and creation of Gaussian input files for QM optimization
    2. RDKit-based conformer generation of a Pd complex, followed by creation of Gaussian input files containing a genECP section to account for Pd atoms
    3. CREST-based conformer sampling of a noncovalent isopentane--water complex, followed by creation of Gaussian input files
    4. CREST-based conformer generation using a transition state (TS) with 3 separated molecules, followed by generation of ORCA input files of TS

  * QCORR_workflows.ipynb (QCORR analysis of Gaussian output files and creation of QM input files):
    1. Analyze optimized QM ground and transition states and create json files of normally terminated files with no errors, extra imaginary frequencies, duplicates, etc.
    2. Use json files to generate single-point energy corrections with multiple levels of theory, charge and multiplicity through for loops:
      3. For Gaussian files, genECP and NBO files containing final extra lines are provided.
      4. For ORCA input files, a DLPNO-CCSD(T) example with multiple % sections is provided.
      5. For pySCF, a calculation is set with DeepMind 21 (DM21).

### Using AQME through command lines in terminals
AQME can also be run through command lines. Some examples are:
  * CSEARCH/CMIN for conformer generation with one SMILES and name:
    ```
    python -m aqme --smi CCCC --name butane --CSEARCH rdkit
    ```
  * CSEARCH/CMIN for conformer generation with multiple SMILES and names:
    ```
    python -m aqme --input smi.csv --CSEARCH rdkit
    ```
    ** The csv file must contain the list of SMILES in a column called "SMILES" and the corresponding names in a column called "code_names" (see Example_workflows for more information)\
    ** Include `--CMIN xtb` or `--CMIN ani` to use the CMIN geometry refiner
    
  * CSEARCH/CMIN for conformer generation using a YAML file containing constrains:
    ```
    python -m aqme --varfile varfile.yaml
    ```
    ** The YAML file must contain the following parameters:
      ```
      input : 'smi.csv' #name of input
      output_name : 'csearch' #name for output
      verbose : True

      CSEARCH : 'rdkit' #csearch option
      qprep : 'gaussian' #program for QM input file generation

      nprocs : 8 #number of processors
      mem : '24GB' #memory
      qm_input : 'B3LYP/6-31G**' #keywords line for the QM inputs
      suffix : 'rdkit' #suffix for the QM inputs
      ```
  * QCORR analysis of Gaussian output files and json file generation:
    ```
    python -m aqme --qcorr --files=*.log
    ```
  * Input file generation from last geometry of output files (log or out files):
    ```
    python -m aqme --qprep --qm_input "M062x def2tzvp opt freq" --files *.log --suffix m062x
    ```
  * Input file generation from json files:
    ```
    python -m aqme --qprep --qm_input "M062x def2tzvp opt freq" --files *.json --suffix m062x
    ```

## Extended documentation (installation, use, examples, etc)
** ReadTheDocs page in process **  
- [ ] QPREP arguments:  
  **files : mol object, str or list of str, default=None**  
		Optionally, starts from: 1) a mol object(s), 2) Gaussian or ORCA output file(s) or 3) json file(s) to prepare the input QM file(s). Options: mol object, FILENAME.log, FILENAME.out or FILENAME.json. Also, lists can be parsed (i.e. [FILE1.log, FILE2.log] or *.FORMAT for glob.glob('*.FORMAT') such as *.json).  
  **atom_types : list of str, default=[]**  
    (If input_file is None) List containing the atoms of the system  
  **cartesians : list of str, default=[]**  
		(If input_file is None) Cartesian coordinates used for further processing  
  **w_dir_main : str, default=os.getcwd()**  
    Working directory  
  **destination : str, default=None,**  
    Directory to create the input file(s)  
  **molecule : str, default="",**  
    Name of the input file without the extension (i.e. NAME for NAME.com)  
  **varfile : str, default=None**  
    Option to parse the variables using a yaml file (specify the filename)  
  **program : str, default='gaussian'**  
    Program required to create the new input files. Current options: 'gaussian', 'orca'  
  **qm_input : str, default=''**  
    Keywords line for new input files (i.e. 'B3LYP/6-31G opt freq')  
  **qm_end : str, default=''**  
    Final line(s) in the new input files  
  **charge : int, default=0**  
		Charge of the calculations used in the following input files  
  **mult : int, default=1**  
		Multiplicity of the calculations used in the following input files  
	**suffix : str, default=''**  
    Suffix for the new input files (i.e. FILENAME_SUFFIX.com for FILENAME.log)  
  **chk : bool, default=False**  
    Include the chk input line in new input files for Gaussian calculations  
  **mem : str, default='4GB'**  
    Memory for the QM calculations (i) Gaussian: total memory; (ii) ORCA: memory per processor  
  **nprocs : int, default=2**  
    Number of processors used in the QM calculations  
  **gen_atoms : list of str, default=[]**  
    Atoms included in the gen(ECP) basis set (i.e. ['I','Pd'])  
  **bs_gen : str, default=''**  
    Basis set used for gen(ECP) atoms	(i.e. 'def2svp')  
  **bs : str, default=''**  
    Basis set used for non gen(ECP) atoms in gen(ECP) calculations (i.e. '6-31G\*')  

- [ ] QCORR arguments:  
  **files : list of str, default=''**  
    Filenames of QM output files to analyze. If \*.log (or other strings that are not lists such as \*.out) are specified, the program will look for all the log files in the working directory through glob.glob(\*.log)  
  **w_dir_main : str, default=os.getcwd()**  
    Working directory  
  **fullcheck : bool, default=True**  
    Perform an analysis to detect whether the calculations were done homogeneously (i.e. same level of theory, solvent, grid size, etc)  
  **varfile : str, default=None**  
    Option to parse the variables using a yaml file (specify the filename)  
  **ifreq_cutoff : float, default=0.0**  
    Cut off for to consider whether a frequency is imaginary (absolute of the specified value is used)  
  **amplitude_ifreq : float, default=0.2**  
    Amplitude used to displace the imaginary frequencies to fix  
  **freq_conv : str, default=None**  
    If a string is defined, it will remove calculations that converged during optimization but did not convergence in the subsequent frequency calculation. Options: opt keyword as string (i.e. 'opt=(calcfc,maxstep=5)'). If readfc is specified in the string, the chk option must be included as well.  
  **s2_threshold : float, default=10.0**  
  	Cut off for spin contamination during analysis in % of the expected value (i.e. multiplicity 3 has an the expected <S\*\*2> of 2.0, if s2_threshold = 10, the <S\*\*2> value is allowed to be 2.0 +- 0.2). Set s2_threshold = 0 to deactivate this option.  
  **dup_threshold : float, default=0.0001**  
    Energy (in hartree) used as the energy difference in E, H and G to detect duplicates  
  **isom : str, default=None**  
    Check for isomerization from the initial input file to the resulting output files. It requires the extension of the initial input files (i.e. isom='com') and the folder of the input files must be added in the isom_inputs option  
  **isom_inputs : str, default=os.getcwd()**  
    Folder containing the initial input files to check for isomerization  
  **vdwfrac : float, default=0.50**  
    Fraction of the summed VDW radii that constitutes a bond between two atoms in the isomerization filter  
  **covfrac : float, default=1.10**  
    Fraction of the summed covalent radii that constitutes a bond between two atoms in the isomerization filter  
  
  **--- Options related to file generation to fix issues found by QCORR ---**  
  New input files are generated through the QPREP module and, therefore, all QPREP arguments can be used when calling QCORR and will overwrite default options. For example, if the user specifies qm_input='wb97xd/def2svp', all the new input files generated to fix issues will contain this keywords line. See examples in the 'Example_workflows' folder for more information.

## Developers and help desk
List of main developers and contact emails:
  - [ ] [Shree Sowndarya S. V.](https://orcid.org/0000-0002-4568-5854), main developer of the CSEARCH and CMIN modules. Contact: [svss@colostate.edu](mailto:svss@colostate.edu)
  - [ ] [Juan V. Alegre-Requena](https://orcid.org/0000-0002-0769-7168), main developer of the QCORR and QPREP modules. Contact: [juanvi89@hotmail.com](mailto:juanvi89@hotmail.com)
  - [ ] [Turki Alturaifi](https://www.chem.pitt.edu/person/turki-alturaifi), worked in benchmarking the parameters for RDKit-based conformer generation. Contact: [turki0@rams.colostate.edu](mailto:turki0@rams.colostate.edu)
  - [ ] [Raúl Pérez-Soto](https://orcid.org/0000-0002-6237-2155), worked in refactoring the code. Contact: [rperez@iciq.es](mailto:rperez@iciq.es)
  - [ ] [Robert S. Paton](https://orcid.org/0000-0002-0104-4166), research group supervisor and code advisor. Contact: [robert.paton@colostate.edu](mailto:robert.paton@colostate.edu)

For suggestions and improvements of the code (greatly appreciated!), please reach out through the issues and pull requests options of Github.

## License
AQME is freely available under an [MIT](https://opensource.org/licenses/MIT) License

## Reference
AQME v1.0, Alegre-Requena, J. V.; Sowndarya, S.; Pérez-Soto, R.; Alturaifi, T. M.; Paton, R. S., 2021. https://github.com/jvalegre/aqme
