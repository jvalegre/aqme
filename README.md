![](Logos/AQME_logo.jpg)
[![CircleCI](https://img.shields.io/circleci/build/github/jvalegre/aqme?label=Circle%20CI&logo=circleci)](https://app.circleci.com/pipelines/github/jvalegre/aqme)
[![Codecov](https://img.shields.io/codecov/c/github/jvalegre/aqme?label=Codecov&logo=codecov)](https://codecov.io/gh/jvalegre/aqme)
[![Downloads](https://img.shields.io/conda/dn/conda-forge/aqme?label=Downloads&logo=Anaconda)](https://anaconda.org/conda-forge/aqme)
[![Read the Docs](https://img.shields.io/readthedocs/aqme?label=Read%20the%20Docs&logo=readthedocs)](https://aqme.readthedocs.io/)

___
# <p align="center">Automated Quantum Mechanical Environments (AQME)</p>
## <p align="center">-- Table of contents --</p>

### <p align="center">[What is AQME?](https://github.com/jvalegre/aqme#what-is-aqme) &nbsp; &nbsp; [Installation](https://github.com/jvalegre/aqme#installation) &nbsp; &nbsp; [Requirements](https://github.com/jvalegre/aqme#requirements)</p>
### <p align="center">[Example workflows and running tests](https://github.com/jvalegre/aqme#example-workflows-and-running-tests-1)</p>
### <p align="center">[Features and modules](https://github.com/jvalegre/aqme#features-and-modules) &nbsp; &nbsp; [Quickstart](https://github.com/jvalegre/aqme#quickstart) &nbsp; &nbsp; [Extended documentation](https://github.com/jvalegre/aqme#extended-documentation-installation-use-examples-etc)</p>
### <p align="center">[Developers and help desk](https://github.com/jvalegre/aqme#developers-and-help-desk) &nbsp; &nbsp; [License](https://github.com/jvalegre/aqme#license) &nbsp; &nbsp; [Reference](https://github.com/jvalegre/aqme#reference)</p>
___
## What is AQME?  
The code is an ensemble of automated QM workflows that can be run through jupyter notebooks, command lines and yaml files. Some of the most popular workflows include:  
  - [ ] RDKit- and CREST-based conformer generator leading to ready-to-submit QM input files starting from individual files or SMILES databases  
  - [ ] Post-processing of QM output files to fix convergence errors, extra imaginary frequencies, unfinished jobs, duplicates and error terminations, as well as to detect spin contamination, isomerization issues, and more optimization problems  
  - [ ] Analysis of homogeneity of QM calculations (same level of theory, grid size, program and version, solvation models, etc)  
  - [ ] Generation of xTB, DFT and RDKit descriptors in json and csv files that are ready to use in machine-learning models or used to predict NMR spectra  
  - [ ] Other useful workflows  

Don't miss out the latest hands-on tutorials from our [YouTube channel](https://www.youtube.com/channel/UCHRqI8N61bYxWV9BjbUI4Xw)!  

## Installation
Check our [AQME installation in 2 mins](https://youtu.be/VeaBzqIZHbo) video for a quick installation guide. In a nutshell, AQME and its dependencies are installed as follows:  
1. Using conda-forge: `conda install -c conda-forge aqme`  
2. Using pip: `pip install aqme`  
  
Extra requirements if xTB or CREST are used (MacOS and Linux only):  
  * xTB: `conda install -y -c conda-forge xtb`  
  * CREST: `conda install -y -c conda-forge crest`  

Extra requirements if CMIN is used with ANI models:  
  * torch-related modules: `pip install torch torchvision torchani`  
  
Known incompatibilities:  
  * RDKit cannot be installed through `pip install rdkit` in Windows when Anaconda prompts are used     

## Requirements
* Python 3  
* Any of the AQME installation options as detailed in the installation section
* Torch-related modules if CMIN is used (shown in the installation section)

## Example workflows and running tests
- The inputs to run pre-defined AQME end-to-end workflows are available in the "/Example_workflows/End-to-end_Workflows" folder. Choose the workflow and run the inputs.
- Automated protocols for individual modules and tasks are provided in the /Example_workflows/ folder inside subfolders with the corresponding module names.
- To run the tests, run pytest in a terminal as follows `pytest --v` from the main AQME folder or `pytest --v PATH` using the PATH where AQME is installed.

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
There are multiple ready-to-use workflows presented as jupyter notebooks in the 'Example workflows' folder. Some examples are:  
  * CSEARCH_CMIN_conformer_generation folder --> CSEARCH/CMIN conformational sampling from SMILES and creation of QM input files  

  * QCORR_processing_QM_outputs --> QCORR analysis of Gaussian output files, generation of JSON files with all the information and creation of new QM input files  

  * QPREP_generating_input_files --> QPREP preparation of input files for Gaussian, ORCA and PySCF from LOG/OUT, SDF and JSON files  

### Using AQME through command lines in terminals
AQME can also be run through command lines. Some examples are:  
  * CSEARCH for conformer generation with one SMILES and name using RDKit or CREST (use rdkit or crest in --program):  
    ```
    python -m aqme --csearch --program rdkit --smi "CCC" --name proprane
    ```  
  * CSEARCH for conformer generation with multiple SMILES and names (i.e. from a database in CSV format):  
    ```
    python -m aqme --csearch --program rdkit --input FILENAME.csv
    ```  
    ** The csv file must contain the list of SMILES in a column called "SMILES" and the corresponding names in a column called "code_name" (see Example_workflows for more information)  

  * CMIN geometry optimization with xTB or ANI (use xtb or ani in --program; use sdf, xyz, com/gjf or pdb in --files):  
    ```
    python -m aqme --cmin --program xtb --files "*.sdf"
    ```  

  * QPREP input file generation from SDF, JSON and LOG/OUT files (replace "*.sdf" for the corresponding format):  
    ```
    python -m aqme --qprep --program gaussian --qm_input "M062x def2tzvp opt freq" --files "*.sdf"
    ```  

  * QCORR analysis of Gaussian output files and json file generation:  
    ```
    python -m aqme --qcorr --program gaussian --freq_conv "opt=(calcfc,maxstep=5)" --files "*.log"
    ```  

## Extended documentation (installation, use, examples, etc)
** ReadTheDocs page in process **  
- [ ] CSEARCH arguments:  
    **input : str, default=''**  
        (If smi is None) Optionally, file containing the SMILES strings and names of the molecules. Current file extensions: .smi, .sdf, .cdx, .csv, .com, .gjf, .mol, .mol2, .xyz, .txt, .yaml, .yml, .rtf  
        For .csv files (i.e. FILENAME.csv), two columns are required, 'code_name' with the names and 'SMILES' for the SMILES string  
    **program : str, default=None**  
        Program required in the conformational sampling. Current options: 'rdkit', 'summ', 'fullmonte', 'crest'  
    **smi : str, default=None**  
        Optionally, define a SMILES string as input  
    **name : str, default=None**  
    (If smi is defined) optionally, define a name for the system  
    **w_dir_main : str, default=os.getcwd()**  
        Working directory  
    **varfile : str, default=None**  
        Option to parse the variables using a yaml file (specify the filename)  
    **max_workers : int, default=4**  
        Number of simultaneous RDKit jobs run with multiprocessing (WARNING! More than 12 simultaneous jobs might collapse your computer!)  
    **charge : int, default=None**  
        Charge of the calculations used in the following input files. If charge isn't defined, it automatically reads the charge of the SMILES string  
    **mult : int, default=None**  
        Multiplicity of the calculations used in the following input files. If mult isn't defined, it automatically reads the multiplicity of the mol object created with the SMILES string. Be careful with the automated calculation of mult from mol objects when using metals!  
    **prefix : str, default=''**  
        Prefix added to all the names  
    **suffix : str, default=''**  
        Suffix added to all the names  
    **stacksize : str, default='1G'**  
        Controls the stack size used (especially relevant for xTB/CREST calculations of large systems, where high stack sizes are needed)  

    *-- Options for RDKit-based methods (RDKit, SUMM and Fullmonte), organic and organometallic molecules --*  
    **sample : int, default='auto'**  
        Number of conformers used initially in the RDKit sampling. If this option isn't specified, AQME automatically calculates (previously benchmarked) an approximate number based on number of rotatable bonds, XH (i.e. OH) groups, saturated cycles, etc (see the auto_sampling() function in csearch.py for more information)  
    **auto_sample : int, default=20**  
        Base multiplicator number used in the sample option  
    **ff : str, default='MMFF'**  
        Force field used in RDKit optimizations and energy calculations. Current options: MMFF and UFF (if MMFF fails, AQME tries to use UFF automatically)  
    **ewin_csearch : float, default=5.0**  
        Energy window in kcal/mol to discard conformers (i.e. if a conformer is more than the E window compared to the most stable conformer)  
    **initial_energy_threshold : float, default=0.0001**  
        Energy difference in kcal/mol between unique conformers for the first filter of only E  
    **energy_threshold : float, default=0.25**  
        Energy difference in kcal/mol between unique conformers for the second filter of E + RMS  
    **rms_threshold : float, default=0.25**  
        RMS difference between unique conformers for the second filter of E + RMS  
    **opt_steps_rdkit : int, default=1000**  
        Max cycles used in RDKit optimizations  
    **heavyonly : bool, default=True**  
        Only consider heavy atoms during RMS calculations for filtering (in the Chem.rdMolAlign.GetBestRMS() RDKit function)  
    **max_matches_rmsd : int, default=1000**  
        Max matches during RMS calculations for filtering (maxMatches option in the Chem.rdMolAlign.GetBestRMS() RDKit function)  
    **max_mol_wt : int, default=0**  
        Discard systems with molecular weights higher than this parameter (in g/mol). If 0 is set, this filter is off  
    **max_torsions : int, default=0**  
        Discard systems with more than this many torsions (relevant to avoid molecules with many rotatable bonds). If 0 is set, this filter is off  
    **seed : int, default=62609**  
        Random seed used during RDKit embedding (in the Chem.rdDistGeom.EmbedMultipleConfs() RDKit function)  

    *-- Options for RDKit-based methods (RDKit, SUMM and Fullmonte), organometallic molecules only --*  
    **metal_atoms : list of str, default=[]**  
        Specify metal atom(s) of the system as [ATOM_TYPE]. Multiple metals can be used simultaneously (i.e. ['Pd','Ir']).  This option is important to calculate the charge of metal complexes based on SMILES strings. Requires the use of metal_oxi.  
    **metal_oxi : list of int, default=[]**  
        Specify metal oxidation state as [NUMBER]. Multiple metals can be used simultaneously (i.e. [2,3]).  
    **complex_type : str, default=''**  
        Forces the metal complexes to adopt a predefined geometry. This option is especially relevant when RDKit predicts wrong complex geometries or gives a mixture of geometries. Current options: squareplanar, squarepyramidal, linear, trigonalplanar  

    *-- Options for the SUMM method only --*  
    **degree : float, default=120.0**  
        Interval of degrees to rotate dihedral angles during SUMM sampling (i.e. 120.0 would create 3 conformers for each dihedral, at 0, 120 and 240 degrees)  

    *-- Options for the Fullmonte method only --*  
    **ewin_fullmonte : float, default=5.0**  
        Energy window in kcal/mol to discard conformers (i.e. if a conformer is more than the E window compared to the most stable conformer)  
    **ewin_sample_fullmonte : float, default=2.0**  
        Energy window in kcal/mol to use conformers during the Fullmonte sampling (i.e. conformers inside the E window compared to the most stable conformer are considered as unique in each step of the sampling)  
    **nsteps_fullmonte : int, default=100**  
        Number of steps (or conformer batches) to carry during the Fullmonte sampling  
    **nrot_fullmonte : int, default=3**  
        Number of dihedrals to rotate simultaneously (picked at random) during each step of the Fullmonte sampling  
    **ang_fullmonte : float, default=30**  
        Available angle interval to use in the Fullmonte sampling. For example, if the angle is 120.0, the program chooses randomly between 120 and 240 degrees (picked at random) during each step of the sampling  

    *-- Options for CREST --*  
    **nprocs : int, default=2**  
        Number of processors used in CREST optimizations  
    **constraints_atoms : list, default=[]**  
        Specify constrained atoms as [AT1,AT2,AT3]. An example of multiple constraints (atoms 1, 2 and 5 are frozen: [1,2,5]  
    **constraints_dist : list of lists, default=[]**  
        Specify distance constraints as [AT1,AT2,DIST]. An example of multiple constraints (atoms 1 and 2 with distance 1.8 Å, and atoms 4 and 5 with distance 2.0 Å): [[1,2,1.8],[4,5,2.0]]  
    **constraints_angle : list of lists, default=[]**  
        Specify angle constraints as [AT1,AT2,AT3,ANGLE]. An example of multiple constraints (atoms 1, 2 and 3 with an angle of 180 degrees, and atoms 4, 5 and 6 with an angle of 120): [[1,2,3,180],[4,5,6,120]]  
    **constraints_dihedral : list of lists, default=[]**  
        Specify dihedral constraints as [AT1,AT2,AT3,AT4,DIHEDRAL]. An example of multiple constraints (atoms 1, 2, 3 and 4 with a dihedral angle of 180 degrees, and atoms 4, 5, 6 and 7 with a dihedral angle of 120): [[1,2,3,4,180],[4,5,6,7,120]]  
    **crest_force : float, default=0.5**  
        Force constant for constraints in the .xcontrol.sample file for CREST jobs  
    **crest_keywords : str, default=None**  
        Define additional keywords to use in CREST that are not included in --chrg, --uhf, -T and -cinp. For example: '--alpb ch2cl2 --nci --cbonds 0.5'  
    **cregen : bool, default=False**  
        If True, perform a CREGEN analysis after CREST (filtering options below)  
    **cregen_keywords : str, default=None**  
        Additional keywords for CREGEN (i.e. cregen_keywords='--ethr 0.02')  
    **xtb_keywords : str, default=None**  
        Define additional keywords to use in the xTB pre-optimization that are not included in -c, --uhf, -P and --input. For example: '--alpb ch2cl2 --gfn 1'   

- [ ] CMIN arguments:  
    **files : str or list of str, default=None**  
        Input files. Formats accepted: XYZ, SDF, GJF, COM and PDB. Also, lists can be used (i.e. [FILE1.sdf, FILE2.sdf] or \*.FORMAT such as \*.sdf).  
    **program : str, default=None**  
        Program required in the conformational refining. Current options: 'xtb', 'ani'  
    **w_dir_main : str, default=os.getcwd()**  
        Working directory  
    **destination : str, default=None,**  
        Directory to create the output file(s)  
    **varfile : str, default=None**  
        Option to parse the variables using a yaml file (specify the filename)  
    **nprocs : int, default=2**  
        Number of processors used in the xTB optimizations  
    **charge : int, default=None**  
        Charge of the calculations used in the xTB calculations. If charge isn't defined, it automatically reads the charge from the input SDF files (if the files come from CSEARCH, which adds the property "Real charge") or calculates it from the generated mol object  
    **mult : int, default=None**  
        Multiplicity of the calculations used in the xTB calculations. If charge isn't defined, it automatically reads the charge from the input SDF files (if the files come from CSEARCH, which adds the property "Mult") or calculates it from the generated mol object. Be careful with the automated calculation of mult from mol objects when using metals!  
    **metal_atoms : list of str, default=[]**  
        Specify metal atom(s) of the system as [ATOM_TYPE]. Multiple metals can be used simultaneously (i.e. ['Pd','Ir']). This option is useful to calculate molecular charges automatically (i.e. from metal databases). Requires the use of metal_oxi.  
    **metal_oxi : list of int, default=[]**  
        Specify metal oxidation state as [NUMBER]. Multiple metals can be used simultaneously (i.e. [2,3]).  
    **ewin_cmin : float, default=5.0**  
        Energy window in kcal/mol to discard conformers (i.e. if a conformer is more than the E window compared to the most stable conformer)  
    **initial_energy_threshold : float, default=0.0001**  
        Energy difference in kcal/mol between unique conformers for the first filter of only E  
    **energy_threshold : float, default=0.25**  
        Energy difference in kcal/mol between unique conformers for the second filter of E + RMS  
    **rms_threshold : float, default=0.25**  
        RMS difference between unique conformers for the second filter of E + RMS  
    **stacksize : str, default='1G'**  
        Controls the stack size used (especially relevant for xTB/CREST calculations of large systems, where high stack sizes are needed)  

    *-- Options for xTB --*  
    **xtb_keywords : str, default=None**  
        Define additional keywords to use in xTB that are not included in -c, --uhf, -P and --input. For example: '--alpb ch2cl2 --gfn 1'  
    **constraints_atoms : list, default=[]**  
        Specify constrained atoms as [AT1,AT2,AT3]. An example of multiple constraints (atoms 1, 2 and 5 are frozen: [1,2,5]  
    **constraints_dist : list of lists, default=[]**  
        Specify distance constraints as [AT1,AT2,DIST]. An example of multiple constraints (atoms 1 and 2 with distance 1.8 Å, and atoms 4 and 5 with distance 2.0 Å): [[1,2,1.8],[4,5,2.0]]  
    **constraints_angle : list of lists, default=[]**  
        Specify angle constraints as [AT1,AT2,AT3,ANGLE]. An example of multiple constraints (atoms 1, 2 and 3 with an angle of 180 degrees, and atoms 4, 5 and 6 with an angle of 120): [[1,2,3,180],[4,5,6,120]]  
    **constraints_dihedral : list of lists, default=[]**  
        Specify dihedral constraints as [AT1,AT2,AT3,AT4,DIHEDRAL]. An example of multiple constraints (atoms 1, 2, 3 and 4 with a dihedral angle of 180 degrees, and atoms 4, 5, 6 and 7 with a dihedral angle of 120): [[1,2,3,4,180],[4,5,6,7,120]]  

    *-- Options for ANI --*  
    **opt_steps : int, default=1000**  
        Maximum number of steps used in the ase.optimize.BFGS optimizer.  
    **opt_fmax : float, default=0.05**  
        Maximum force value to determine convergence in the ase.optimize.BFGS optimizer.  
    **ani_method : str, default='ANI2x'**  
        ANI model used in the ase.optimize.BFGS optimizer.  

- [ ] QPREP arguments:  
    **files : mol object, str or list of str, default=None**  
        Files used to prepare input QM file(s). Formats accepted: mol object(s), Gaussian or ORCA LOG/OUT output files, JSON, XYZ, SDF, PDB. Also, lists can be used (i.e. [FILE1.log, FILE2.log] or \*.FORMAT such as \*.json).  
    **program : str, default=None**  
        Program required to create the new input files. Current options: 'gaussian', 'orca'  
    **atom_types : list of str, default=[]**  
        (If files is None) List containing the atoms of the system  
    **cartesians : list of str, default=[]**  
        (If files is None) Cartesian coordinates used for further processing  
    **w_dir_main : str, default=os.getcwd()**  
        Working directory  
    **destination : str, default=None,**  
        Directory to create the input file(s)  
    **varfile : str, default=None**  
        Option to parse the variables using a yaml file (specify the filename)  
    **qm_input : str, default=''**  
        Keywords line for new input files (i.e. 'B3LYP/6-31G opt freq')  
    **qm_end : str, default=''**  
        Final line(s) in the new input files  
    **charge : int, default=None**  
        Charge of the calculations used in the following input files. If charge isn't defined, it defaults to 0  
    **mult : int, default=None**  
        Multiplicity of the calculations used in the following input files. If mult isn't defined, it defaults to 1  
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
    **bs_nogen : str, default=''**  
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
    **isom_type : str, default=None**  
        Check for isomerization from the initial input file to the resulting output files. It requires the extension of the initial input files (i.e. isom_type='com' or 'gjf') and the folder of the input files must be added in the isom_inputs option  
    **isom_inputs : str, default=os.getcwd()**  
        Folder containing the initial input files to check for isomerization  
    **vdwfrac : float, default=0.50**  
        Fraction of the summed VDW radii that constitutes a bond between two atoms in the isomerization filter  
    **covfrac : float, default=1.10**  
        Fraction of the summed covalent radii that constitutes a bond between two atoms in the isomerization filter  

    *-- Options related to file generation to fix issues found by QCORR --*  
    New input files are generated through the QPREP module and, therefore, all QPREP arguments can be used when calling QCORR and will overwrite default options. For example, if the user specifies qm_input='wb97xd/def2svp', all the new input files generated to fix issues will contain this keywords line. See examples in the 'Example_workflows' folder for more information.  

- [ ] QDESCP arguments:  
    **program : str, default=None**  
        Program required to create the new descriptors. Current options: 'xtb', 'nmr'  
    **w_dir_main : str, default=os.getcwd()**  
        Working directory  
    **destination : str, default=None,**  
        Directory to create the JSON file(s)  

    *-- Options for xTB descriptor generation (program='xtb') --*  
    **files : list of str, default=''**  
        Filenames of SDF/PDB/XYZ files to calculate xTB descriptors. If \*.sdf (or other strings that are not lists such as \*.pdb) are specified, the program will look for all the SDF files in the working directory through glob.glob(\*.sdf)  
    **charge : int, default=None**  
        Charge of the calculations used in the following input files. If charge isn't defined, it defaults to 0  
    **mult : int, default=None**  
        Multiplicity of the calculations used in the following input files. If mult isn't defined, it defaults to 1
    **qdescp_solvent : str, default=None**  
        Solvent used in the xTB property calculations (ALPB model)  
    **qdescp_temp : float, default=300**  
        Temperature required for the xTB property calculations  
    **qdescp_acc : float, default=0.2**  
        Accuracy required for the xTB property calculations  
    **boltz : bool, default=True**  
        Calculate Boltzmann averaged xTB properties and include RDKit molecular features  
  
    *-- Options for NMR spectra simulation (program='nmr') --*  
    **files : list of str, default=''**  
        Filenames of LOG files to retrieve NMR shifts from Gaussian calculations  
    **boltz : bool, default=True**  
        Calculate Boltzmann averaged NMR shifts  
    **nmr_atoms : list of str, default=[6, 1]**  
        List containing the atom types to consider. For example, if the user wants to retrieve NMR shifts from C and H atoms nmr_atoms=[6, 1]  
    **nmr_slope : list of float, default=[-1.0537, -1.0784]**  
        List containing the slope to apply for the raw NMR shifts calculated with Gaussian. A slope needs to be provided for each atom type in the analysis (i.e., for C and H atoms, the nmr_slope=[-1.0537, -1.0784]). These values can be adjusted using the CHESHIRE repository.  
    **nmr_intercept : list of float, default=[181.7815, 31.8723]**  
        List containing the intercept to apply for the raw NMR shifts calculated with Gaussian. An intercept needs to be provided for each atom type in the analysis (i.e., for C and H atoms, the nmr_slope=[-1.0537, -1.0784]). These values can be adjusted using the CHESHIRE repository.  
    **nmr_experim : str, default=None**  
        Filename of a CSV containing the experimental shifts. Two columnds are needed: A) 'atom_idx' should contain the indexes of the atoms to study as seen in GaussView or other molecular visualizers (i.e., the first atom of the coordinates has index 1); B) 'experimental_ppm' should contain the experimental NMR shifts in ppm observed for the atoms.  

- [ ] VISMOL arguments:  
    **files : list of str, default=''**  
        Filenames of SDF/PDB/XYZ to visualize conformers. If \*.sdf (or other strings that are not lists such as \*.pdb) are specified, the program will look for all the SDF files in the working directory through glob.glob(\*.sdf). Internal options of "line", "stick", "sphere" incorporated. Code reference from: [https://iwatobipen.wordpress.com]


## Developers and help desk
List of main developers and contact emails:  
  - [ ] [Juan V. Alegre-Requena](https://orcid.org/0000-0002-0769-7168), main developer of the CSEARCH, QCORR, QPREP and QDESCP modules. Contact: [jv.alegre@csic.es](mailto:jv.alegre@csic.es)  
  - [ ] [Shree Sowndarya S. V.](https://orcid.org/0000-0002-4568-5854), main developer of the CSEARCH, CMIN, QDESCP and VIZMOL modules. Contact: [svss@colostate.edu](mailto:svss@colostate.edu)  
  - [ ] [Turki Alturaifi](https://www.chem.pitt.edu/person/turki-alturaifi), worked in benchmarking the parameters for RDKit-based conformer generation. Contact: [tma53@pitt.edu](mailto:tma53@pitt.edu)  
  - [ ] [Raúl Pérez-Soto](https://orcid.org/0000-0002-6237-2155), worked in refactoring the code and creating the documentation. Contact: [Raul.Perez_Soto@colostate.edu](mailto:Raul.Perez_Soto@colostate.edu)  
  - [ ] [Robert S. Paton](https://orcid.org/0000-0002-0104-4166), research group supervisor and code advisor. Contact: [robert.paton@colostate.edu](mailto:robert.paton@colostate.edu)  

For suggestions and improvements of the code (greatly appreciated!), please reach out through the issues and pull requests options of Github.  

## License
AQME is freely available under an [MIT](https://opensource.org/licenses/MIT) License  

## Reference
If you use any of the AQME modules, please include this citation:  
AQME v1.4, Alegre-Requena, J. V.; Sowndarya, S.; Alturaifi, T.; Pérez-Soto, R.; Paton, R. ChemRxiv 2022, DOI: 10.26434/chemrxiv-2022-dnc48.  
  
Additionally, please include the corresponding references for the following programs:  
  * If you used CSEARCH with RDKit methods or from SMILES: [RDKit](https://www.rdkit.org)  
  * If you used CSEARCH with CREST methods: [CREST](https://crest-lab.github.io/crest-docs)  
  * If you used CMIN with xTB: [xTB](https://xtb-docs.readthedocs.io/en/latest/contents.html)  
  * If you used CMIN with ANI: [ANI](https://github.com/isayev/ASE_ANI)  
  * If you used QCORR: [cclib](https://cclib.github.io/)  
  * If you used QDESCP with xTB: [xTB](https://xtb-docs.readthedocs.io/en/latest/contents.html)