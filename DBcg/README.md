# DBcg
Conformer generator followed by generation of .com files for Gaussian starting from smiles, sdf, csv or cdx files.
The program allows for two round of optimizations if error terminated in the First round.
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

To Do list:
1. Add ENSO conformer generation
2. Automatic the work flow including the job running on the cluster.
3. Check how runtime scales with number of atoms and rotatable bonds. Provide some examples.

Limitations:
1. Salts or complexes with more than one molecule wont work. (RDkit doesn't know how to handle multiple molecules, need to figure out this!)
2. Transition states don't work (need to figure how to generalise templates)

Possible methods of invoking the script:
1.  python -m DBGEN --varfile db_gen_variables.py (from a file of .py format)
2. python -m DBGEN args (command line arguments)

N.B. this requires the location of the DBGEN directory to be added to the $PYTHONPATH environment variable
