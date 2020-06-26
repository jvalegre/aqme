![](Logos/DBGEN%20logo.tif)
[![Build Status](https://img.shields.io/travis/com/jvalegre/pyconfort?label=Linux%20CI&logo=Travis)](https://travis-ci.com/github/jvalegre/pyCONFORT)
[![Tests](https://img.shields.io/static/v1?label=Tests&message=96&color=green&logo=Travis)](https://travis-ci.com/github/jvalegre/pyCONFORT)
[![Codecov](https://img.shields.io/codecov/c/github/jvalegre/pyconfort?label=Codecov&logo=codecov)](https://codecov.io/gh/jvalegre/pyCONFORT)
[![CodeFactor](https://img.shields.io/codefactor/grade/github/jvalegre/pyconfort?label=Codefactor%20grade&logo=codefactor)](https://www.codefactor.io/repository/github/jvalegre/pyconfort/overview/master)
[![Codacy](https://img.shields.io/codacy/grade/047e9c6001a84713a82e180669e14c98?label=Codacy%20grade&logo=codacy)](https://www.codacy.com/manual/jvalegre/pyCONFORT?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jvalegre/pyCONFORT&amp;utm_campaign=Badge_Grade)
[![lgtm](https://img.shields.io/lgtm/grade/python/github/jvalegre/pyCONFORT?label=lgtm%20grade&logo=lgtm)](https://lgtm.com/projects/g/jvalegre/pyCONFORT/context:python)

![Size](https://img.shields.io/github/repo-size/jvalegre/pyconfort?label=Size)
![PyPI](https://img.shields.io/pypi/v/pyconfort?label=PyPI&logo=pypi)
![PyPI Downloads](https://img.shields.io/pypi/dm/pyconfort?label=PyPI%20downloads&logo=pypi&style=social)

# DBGEN
## Description
Conformer generator followed by generation of .com files for Gaussian starting from smiles, sdf, csv or cdx files.
The program allows for two round of optimizations if error terminated in the first round.
The thermodynamics of these conformers for each molecule are compiled and returned as .csv files

As of now, all files are created and analyzed by the program (user only required to run the files)

Allows for creation of different input files by varying the following:
1. Level of Theory
2. Basis Set (accounts for transition metals by invoking genecp)
2. Basis Set genecp atoms (only one type can be specified as of now)
3. Solvation
4. Full optimization vs Single Point Calculations(all single points do an NMR calculation).
5. Automatic submission of .com files once they are created if run on summit.
6. Creation of files based on the number of conformers we want to consider (lowest vs energy gap vs all)

Analysis the .log file features (all files are either in the folders /gaussian or /sp):
1. organizes the log files in each level of theory/ basis set combination into folders namely Finished, Imaginary_frequencies, Failed_Error, Failed_Unfinished
2. If not Normal termination, the creates a new folder with all .com to re-run.
3. If all .log files are moved to Finished, the can create new NMR input files.

Analysis for Boltzmann averaging and combining files
1. Respective files for each molecule are grabbed and outputs for each molecule are written to a .csv files
2. All the .csv for each molecule are grabbed and all thermodynamic data are written to three different .csv (all data, average data, comparison of lowest vs avg G)

## Limitations
Currently, these are the limitations of pyCONFORT:
1. Salts or complexes with more than one molecule wont work. (RDkit doesn't know how to handle multiple molecules, need to figure out this!)
2. Transition states don't work (need to figure how to generalise templates)

## Using octahedral, trigonal bipyramidal, square-planar and square-based pyramidal structures
RDKit uses iodine to generate:
Octahedral I+ complex (coordination number 6), XXX I+ complex (coordination number 5), Tetrahedral I+ complex (coordination number 4)
\*\*Otherwise, RDKit gives error of coordination if you use the central metals directly.

Then, the code replaces the central I for the metal that you choose with the oxidation state and multiplicity that you want. After this replacement, it runs xTB and then you get a geometry that looks like what you would expect for that complex (i.e. you start from a I+ tetrahedral complex from RDKIT, you replace I+ for Pd and run xTB, then xTB gives you a square-planar geometry with Pd). This "cleaned up" geometry is then used to create the final COM input file.

## Possible methods of invoking the script
These are the two pssible methods to run the scripts:
1. (PREFERRED) python -m DBGEN --varfile params.py (it reads all the variables from a file of .py format)
2. python -m DBGEN --compute --input FILENAME.smi args (command line arguments)

## Examples
### (1) File with SMILES
python -m DBGEN --compute --input FILENAME.smi
(where FILENAME.smi has the format /SMILES NAME/:

CCCCC pentane

CCCCCC hexane

### (2) SDF file with 3D molecules
python -m DBGEN --compute --input FILENAME.sdf
(where FILENAME.sdf contains all the molecules to use)

### (3) Multiple SMILES or SDF files
python -m DBGEN --compute --input \*.smi

python -m DBGEN --compute --input \*.sdf

### (4) Multiple SDF files with paramaters adjusted for a certain DFT level
python -m DBGEN --compute --input \*.sdf

*** First, make sure that (1) you have the params.py file in the folder you are running the script and (2) you edit the params.py with the level of theory and type of calculation that you want

## To Do list
Work on progress:
1. Add ENSO conformer generation
2. Automate the work flow including the job running on the cluster.
3. Check how runtime scales with number of atoms and rotatable bonds. Provide some examples.
4. Make the program work with multiple molecules in the same calc (i.e. noncovalent complexes)

## Installation

    (1) Install the python modules below (they are widely used modules, you can use "pip install" or "conda install")

    (2) Download DBGEN folder (there is a DBGEN subfolder inside)

    (3a) If you don't use DBGEN as a module through your PYTHONPATH, you can run the program from the DBGEN main folder

    (3b) You can run DBGEN from other folders if you add the location of the DBGEN directory to the $PYTHONPATH environment variable

## Requirements
(1) Python 3

(2) Python modules:

    (a) General:

        NumPy
        periodictable
        pandas
        openbabel
        pyyaml

     (b) If you use the compute option (conformer generation):

        RDKit

        (b.1) If you use xTB optimizations (mandatory for metal complexes):

            xTB (only if xTB is used for conformer generation)

        (b.2) If you use AN1 optimizations:

            ase
            ase.optimize
            torch
            torchani
            argparse

     (c) If you use the analyze option (post-processing of output files):

     (c) If you use the analyze option (post-processing of output files):
