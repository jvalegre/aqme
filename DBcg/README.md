# DBcg
Conformer generator followed by generation of .com files for Gaussian starting from smiles, sdf, csv or cdx files.
The program allows for two round of optimizations if error terminated in the First round.
The thermodynamics of these conformers for each molecule are compiled and returned as .csv files

As of now, all files are created and analyzed by the program (user only required to run the files)

Allows for creation of different input files by varying the following:
1. Level of Theory
2. Basis Set (accounts for transition metals by invoking genecp)
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

Steps involved in running the db_gen.py script

#pass the respective file as the first argument
1. Running the compute job for either csv, smi, sdf, cdx
python /DBcg/db_gen.py --compute --input example.smi

#PWD : location in which the program is run
#now after the above step all jobs would have been submitted.
2. Check for any non normal termination using the output analysis
python /DBcg/db_gen.py --analysis --path $(PWD)/gaussian/

##resubmit the gaussian jobs if the didnt normally terminate
3. Resubmit gaussian jobs if they didnt terminate normally
python /DBcg/db_gen.py --resubmit --path $(PWD)/gaussian/

##once resubmission is completed analyse the NEWLY CREATED Files and moves the normally terminated to the finished folder
4. sheck for any non normal termination using the output analysis after the second run
python /DBcg/db_gen.py --analysis --secondrun --path $(PWD)/gaussian/

##creating the NMR input files from finished log files i.e., for doing NMR after DFT optimization
5. Carry out NMR of DFT optimised files
python /DBcg/db_gen.py --nmr --path $(PWD)/gaussian/

##carrying out the boltzmann from goodvibes
##note change the name of the log files accoring to your description in the reading of log files part
4. Get the boltzmann contribution for each molecule in the database
python /DBcg/db_gen.py --boltz --path $(PWD)/gaussian/

##combing all the energies for all molecules and printing the output in one csv_file
5. Combine all the results for molecules in the data base
python //DBcg/db_gen.py --combine --path $(PWD)/gaussian/
