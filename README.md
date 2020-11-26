![](Logos/DBGEN%20logo.tif)
[![Build Status](https://img.shields.io/travis/com/jvalegre/pyconfort?label=Linux%20CI&logo=Travis)](https://travis-ci.com/github/jvalegre/pyCONFORT)
[![Tests](https://img.shields.io/static/v1?label=Tests&message=104&color=green&logo=Travis)](https://travis-ci.com/github/jvalegre/pyCONFORT)
[![Codecov](https://img.shields.io/codecov/c/github/jvalegre/pyconfort?label=Codecov&logo=codecov)](https://codecov.io/gh/jvalegre/pyCONFORT)
[![CodeFactor](https://img.shields.io/codefactor/grade/github/jvalegre/pyconfort?label=Codefactor%20grade&logo=codefactor)](https://www.codefactor.io/repository/github/jvalegre/pyconfort/overview/master)
[![Codacy](https://img.shields.io/codacy/grade/047e9c6001a84713a82e180669e14c98?label=Codacy%20grade&logo=codacy)](https://www.codacy.com/manual/jvalegre/pyCONFORT?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jvalegre/pyCONFORT&amp;utm_campaign=Badge_Grade)
[![lgtm](https://img.shields.io/lgtm/grade/python/github/jvalegre/pyCONFORT?label=LGTM%20grade&logo=lgtm)](https://lgtm.com/projects/g/jvalegre/pyCONFORT/context:python)

![Size](https://img.shields.io/github/repo-size/jvalegre/pyconfort?label=Size)
![PyPI](https://img.shields.io/pypi/v/pyconfort?label=PyPI&logo=pypi)
![PyPI Downloads](https://img.shields.io/pypi/dm/pyconfort?label=PyPI%20downloads&logo=pypi&style=social)

# pyCONFORT

Conformer generator followed by generation of com files for Gaussian starting from 1D (SMILES), 2D (CDX) and 3D (sdf,xyz,etc) inputs.

XXX COPY DESCRIPTION FROM INTRO OF PAPER XXX


## Features
XXX MODIFY BASED ON WHAT WE SAY IN THE PAPER XXX
Allows for creation of different input files by varying the following:
1. Level of Theory
2. Basis Set (accounts for transition metals by invoking genecp)
2. Basis Set genecp atoms (only one type can be specified as of now)
3. Solvation
4. Full optimization vs Single Point Calculations
5. Automatic submission of QM input files once they are created if run on summit.
6. Creation of files based on the number of conformers we want to consider (lowest vs energy gap vs all)

Analysis the .log file features (all files are either in the folders /gaussian or /sp):
1. organizes the log files in each level of theory/ basis set combination into folders namely Finished, Imaginary_frequencies, Failed_Error, Failed_Unfinished
2. If not Normal termination, the creates a new folder with all .com to re-run.
3. If all .log files are moved to Finished, the can create new NMR input files.

Analysis for Boltzmann averaging and combining files
1. Respective files for each molecule are grabbed and outputs for each molecule are written to a .csv files
2. All the .csv for each molecule are grabbed and all thermodynamic data are written to three different .csv (all data, average data, comparison of lowest vs avg G)

## Requirements
* Python 3.6, or 3.7 (true?)

## Conda and PyPI (`pip`)
- Basic Instructions go here

## Extended documentation (installation, use, examples, etc)
XXX LINK READTHEDOCS WEBPAGE XXX

## Acknowledgements
go here

## Reference
XXX DOI FOR PAPER XXX
