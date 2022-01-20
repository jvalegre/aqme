![](Logos/AQME_logo.jpg)
[![Build Status](https://img.shields.io/travis/com/jvalegre/aqme?label=Linux%20CI&logo=Travis)](https://travis-ci.com/github/jvalegre/aqme)
[![Tests](https://img.shields.io/static/v1?label=Tests&message=104&color=green&logo=Travis)](https://travis-ci.com/github/jvalegre/aqme)
[![Codecov](https://img.shields.io/codecov/c/github/jvalegre/aqme?label=Codecov&logo=codecov)](https://codecov.io/gh/jvalegre/aqme)
[![CodeFactor](https://img.shields.io/codefactor/grade/github/jvalegre/aqme?label=Codefactor%20grade&logo=codefactor)](https://www.codefactor.io/repository/github/jvalegre/aqme/overview/master)
[![Codacy](https://img.shields.io/codacy/grade/047e9c6001a84713a82e180669e14c98?label=Codacy%20grade&logo=codacy)](https://www.codacy.com/manual/jvalegre/aqme?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jvalegre/aqme&amp;utm_campaign=Badge_Grade)
[![lgtm](https://img.shields.io/lgtm/grade/python/github/jvalegre/aqme?label=LGTM%20grade&logo=lgtm)](https://lgtm.com/projects/g/jvalegre/aqme/context:python)

![Size](https://img.shields.io/github/repo-size/jvalegre/aqme?label=Size)
![PyPI](https://img.shields.io/pypi/v/aqme?label=PyPI&logo=pypi)
![PyPI Downloads](https://img.shields.io/pypi/dm/aqme?label=PyPI%20downloads&logo=pypi&style=social)

# Automated Quantum Mechanical Environments (AQME)

The code is an ensemble of automated QM workflows that can be run through jupyter notebooks, command lines and yaml files. Some of the most popular workflows include:
* RDKit- and CREST-based conformer generator leading to ready-to-submit QM input files starting from individual files or SMILES databases
* post-processing of QM output files to fix convergence errors, extra imaginary frequencies, unfinished jobs, duplicates and error terminations, as well as to detect spin contamination, isomerization issues, and more optimization problems
* analysis of homogeneity of QM calculations (same level of theory, grid size, program and version, solvation models, etc)
* generation of xTB, DFT and RDKit descriptors in json and csv files that are ready to use in machine-learning models or used to predict NMR spectra
* more other useful workflows!

## Installation
1) Installing RDKit and openbabel through conda:
`conda install -c rdkit rdkit` or `conda install -c conda-forge rdkit` (compatible with Python >=3.8)
`conda install -c conda-forge openbabel`

2) Install AQME and its dependencies
`pip install aqme`
or
`python setup.py install`

## Requirements
* Python 3

## Features and modules
### CSEARCH
Module on charge of conformational sampling starting from multiple input types (SMILES, csv, sdf, xyz, etc). Options:
#### RDKit-based conformational sampling
Faster sampling, suitable especially for unimolecular systems.
* RDKit standard sampling
* Systematic Unbounded Multiple Minimum search (SUMM)
* FullMonte sampling
#### CREST-based conformational sampling
Slower sampling, suitable for all types of systems (including noncovalent complexes and constrained systems such as transition states)

### CMIN
Module used to refine conformers generated in CSEARCH through new geometry optimizations. Options:
* xTB (GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF)
* ANI (ANI-1x, ANI-1ccx, ANI-2x)

### QPREP
Generator of input files for QM calculations. Options:
* Gaussian
* ORCA
* pySCF (loading parameters in jupyter notebook)

### QCORR
cclib-based analyzer of output files from multiple QM programs. This module:
* Separates Normal terminated files with no errors, extra imaginary frequencies, duplicates, isomerization to other systems and spin contamination
* Automatically generates new com files to "fix" the different issues of the calculations with strategies that are optimal for each type of issue (Gaussian and ORCA)
* Checks that all the calculations are homogeneous (i.e. using the same level of theory, same grid size, same program and version, solvation model, etc)

### QDESC
Descriptor generator from multiple input types such as SMILES, log files, xyz, etc. Descriptors generated with:
* RDKit
* xTB
* QM programs (cclib-compatible)

## Quickstart
### Using AQME in Jupyter Notebooks
There are multiple ready-to-use workflows presented as jupyter notebooks in the 'Example_workflows' folder. Some examples are:
* QCORR_workflows.ipynb (QCORR analysis of Gaussian output files):
1) Analyze optimized QM ground and transition states and create json files of normally terminated files with no errors, extra imaginary frequencies, duplicates, etc. 
2) Use json files to generate single-point energy corrections with multiple levels of theory, charge and multiplicity through for loops:
2a) For Gaussian files, genECP and NBO files containing final extra lines are provided. 
2b) For ORCA input files, a DLPNO-CCSD(T) example with multiple % sections is provided. 
2c) For pySCF, a calculation is set with DeepMind 21 (DM21).

### Using AQME through command lines in terminals
AQME can also be run through command lines. Some examples are:
* QCORR analysis of Gaussian output files and json file generation:
1) cd into a folder with output files (in this case, LOG files but other formats such as OUT are compatible as well)
2) Run: `python -m aqme --qcorr --qm_files=*.log`
* Input file generation from json files:
1) cd into a folder with json files
2) Run: `python -m aqme --json2input --qm_input "M062x def2tzvp opt freq" --json_files *.json --suffix m062x`

## Extended documentation (installation, use, examples, etc)
** In process **
XXX LINK READTHEDOCS WEBPAGE XXX

## Reference
AQME v1.0, Alegre-Requena, J. V.; Sowndarya, S.; Pérez-Soto, R.; Alturaifi, T. M.; Paton, R. S., 2021. https://github.com/jvalegre/aqme
