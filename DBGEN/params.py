"""

* In this file, the variables to helper programs are collected.
* You must make sure that all the variables are correct before launching db_gen.py.

* OTHER variables USED THROUGHOUT THE PROGRAM ARE ALSO SET HERE.
"""

" INPUT FILE"
input ='smi.smi' #input files
path ='' #path to guassian folder when we do analysis

" EXP RULES"
exp_rules = False # aply some experimental rules to the outputs
angle_off = 30


" IF METAL COMPLEX"
metal_complex= True # specify  true if metal complex
metal = 'Ir' # specify the metal
complex_coord = 6 # specify the complex coordination number
complex_type = 'octahedral' # specify the following square planar, square pyrimidal (otheriwse defaults to octahedral, Td)
m_oxi = 3 # will have to be changed for respective metals
complex_spin = 1 # will have to be changed for respective metals
charge = 0 #will automaticallu change for metal comples. default is for orgaincs.

" IF NCI COMPLEX"
nci_complex = False # specify  true if NCI complex

"GENERAL OPTIONS FOR COMMANDLINE"
verbose = True
compute = True
write_gauss = False
analysis = False
resubmit = False
nmr = False
dup = False
boltz = False
combine = False
prefix = None

"TYPE OF OPTIMIZATION"
# Options: xTB, AN1  Default : RDKIT optimizaiton
ANI1ccx = False
xtb = True

" SINGLE POINTS vs FULL OPTIMIZATION WITH or WITHOUT FREQUENCIES"
single_point = False
frequencies = True

" DEFAULT PARAMETERS FOR RDKIT GENERATION AND FILTERS"
max_torsions = 20 #Skip any molecules with more than this many torsions (default 5)
max_MolWt = 1000
heavyonly = True
sample = 10 #number of conformers to sample to get non-torsional differences (default 100)
nodihedrals = False #turn to TRUE if no dihydral scan is needed.

" DEFAULT PARAMETERS FOR RDKIT OPTIMIZATION "
ff = "MMFF" #can use MMFF ro UFF
etkdg = False #use new ETKDG knowledge-based method instead of distance geometry also needs to be present in RDKIT ENV
seed = int("062609") #random seed (default 062609) for ETKDG
degree = 30 #Amount, in degrees, to enumerate torsions by (default 30.0)

" DEFAULT PARAMETERS FOR ANI1ccx OPTIMIZATION "
constraints = None

" DEFAULT PARAMETERS FOR xTB OPTIMIZATION "
large_sys = False
STACKSIZE = '1G' #set for large system

"DEFAULT PARAMTERS FOR UNIQUE CONFORMER SELECTION"
rms_threshold = 0.25 #cutoff for considering sampled conformers the same (default 0.25)
energy_threshold = 0.05 #energy difference between unique conformers
ewin = 10 #energy window to print conformers
convergence = 1.0 #Adjust convergence criteria of ANI and xtb optimizations (set at 0.005)
time = False #request run time

" ONLY LOWEST ENERGY CONFORMER REQUIRED"
lowest_only = False
lowest_n  = False # for a given threshold of energy_threshold_for_gaussian
energy_threshold_for_gaussian = 10  #in kJ/ mol

" DEFINITION OF ATOMS"
genecp_atoms = ['I','Ir']

"DEFINTION OF BASIS SET AND LEVEL OF THEORY AND SOLVENT"
basis_set = ['def2tzvp']
basis_set_genecp_atoms = ['LANL2DZ']
level_of_theory = ['wB97XD']
max_cycle_opt = 100 #eefault is 300

#dispersion correction to be added or not
dispersion_correction = False
empirical_dispersion = 'GD3BJ'

# Specify the solvation model. Options: gas_phase or any solvation model (i.e. SMD, IEFPCM, CPCM)
solvent_model = 'gas_phase'
solvent_name = 'Acetonitrile'

"DEFAULT PARAMTERS FOR GAUSSIAN OPTIMIZATION"
chk = False
nprocs=36
mem='60GB'

"TURN ON SUBMISSION OF JOBS"
qsub = False
submission_command = 'qsub_summit'

" MOLECULES now, for eg., molecule list, for later can use as total no. of molecules it is need in the boltz part to read in specific molecules"
maxnumber = 103 #max number in your database
