"""

* In this file, the variables to helper programs are collected.
* You must make sure that all the variables are correct before launching db_gen.py.

* OTHER variables USED THROUGHOUT THE PROGRAM ARE ALSO SET HERE.
"""

" INPUT FILE"
input ='smi.smi'
path = ''

"GENERAL OPTIONS FOR COMMANDLINE"
verbose = False
compute = True
analysis = False
resubmit = False
secondrun = False
nmr = False
boltz = True
combine = False

"TYPE OF OPTIMIZATION"
# Options: xTB, AN1  Default : RDKIT optimizaiton
ANI1ccx = False
xtb = True
#enso = False need to add in

" SINGLE POINTS vs  FULL OPTIMIZATION"
single_point = False


" DEFAULT PARAMETERS FOR RDKIT GENERATION AND FILTERS"
max_torsions = 5 #Skip any molecules with more than this many torsions (default 5)
max_MolWt = 1000
heavyonly = True
sample = 100 #number of conformers to sample to get non-torsional differences (default 100)
nodihedrals = True #turn to TRUE if no dihydral scan is needed.

" DEFAULT PARAMETERS FOR RDKIT OPTIMIZATION "
ff = "MMFF" #can use MMFF ro UFF
etkdg = False #use new ETKDG knowledge-based method instead of distance geometry also needs to be present in RDKIT ENV
seed = int("062609") #random seed (default 062609) for ETKDG
degree = 30 #Amount, in degrees, to enumerate torsions by (default 30.0)

" DEFAULT PARAMETERS FOR ANI1ccx OPTIMIZATION "
constraints = None

"DEFAULT PARAMTERS FOR UNIQUE CONFORMER SELECTION"
rms_threshold = 0.25 #cutoff for considering sampled conformers the same (default 0.25)
energy_threshold = 0.05 #energy difference between unique conformers
ewin = 40 #energy window to print conformers
convergence = 1.0 #Adjust convergence criteria of ANI and xtb optimizations (set at 0.005)
time = False #request run time

" ONLY LOWEST ENERGY CONFORMER REQUIRED"
lowest_only = False
lowest_n  = False # for a given threshold of energy_threshold_for_gaussian
energy_threshold_for_gaussian = 2.0 #in kJ/ mol

" DEFINITION OF ATOMS"
genecp_atoms = ['I','Ir']

"DEFINTION OF BASIS SET AND LEVEL OF THEORY AND SOLVENT"
basis_set = ['6-31g*']
basis_set_genecp_atoms = ['LANL2DZ']
level_of_theory = ['M062X','wB97XD']

#dispersion correction to be added or not
dispersion_correction = False
empirical_dispersion = 'D3BJ'

# Specify the solvation model. Options: gas_phase or any solvation model (i.e. SMD, IEFPCM, CPCM)
solvent_model = 'gas_phase'
solvent_name = 'Acetonitrile'

"DEFAULT PARAMTERS FOR GAUSSIAN OPTIMIZATION"
chk = False
nprocs=24
mem='96GB'

"TURN ON SUBMISSION OF JOBS"
qsub = False
submission_command = 'qsub_summit'

" MOLECULES now, for eg., molecule list, for later can use as total no. of molecules it is need in the boltz part to read in specific molecules"
maxnumber = 100 #max number in your database
