"""

* In this file, the variables used by the different python scripts are collected.
* You must make sure that all the variables are correct before launching db_gen.py.

"""

" INPUT FILE "
input ='smi.smi' # input files
path ='' # path to guassian folder when we do analysis

" GENERAL OPTIONS FOR COMMANDLINE "
verbose = True
compute = True
write_gauss = False
analysis = False
resubmit = False
sp = False # write with nmr input line.
dup = False
boltz = False
combine = False
prefix = None

" EXP RULES "
exp_rules = False # apply some experimental rules to discard some outputs
angle_off = 30 # margin of error to determine angles (i.e. if angle_off is 30, and the angle is 180, angles from
		# 150 to 210 degrees will be discarded)

" CHARGE FOR XTB OPTIMIZATION AND COM FILES "
charge = 0 # final charge of the molecule (used in xTB optimization and input in the final COM input files)
		# If metal_complex = True, the script will recalculate the charge

" OPTIMIZATION PARAMETERS "

" DEFAULT PARAMTERS FOR UNIQUE CONFORMER SELECTION FOR RDKIT, XTB AND ANI1 "

rms_threshold = 0.5 #cutoff for considering sampled conformers the same (default 0.25) for RDKit and xTB duplicate filters
energy_threshold = 1 #energy difference between unique conformers for RDKit and xTB duplicate filters
ewin = 1000 #energy window to print conformers for RDKit and xTB duplicate filters
time = False #request run time

" TYPE OF OPTIMIZATION "
# Options: xTB, ANI1ccx (if True is selected).  Default : RDKit optimizaiton
ANI1ccx = False
xtb = True

" SINGLE POINTS vs FULL OPTIMIZATION WITH or WITHOUT FREQUENCIES "
single_point = False
frequencies = True

" FILTERS FOR RDKIT OPTIMIZATION "
max_torsions = 20 # Skip any molecules with more than this many torsions (default 5)
max_MolWt = 10000 # Skip any molecules with molecular weights higher than this number

" DIHEDRAL PROTOCOL FOR RDKIT OPTIMIZATION (SLOW SINCE IT SCANS MANY DIHEDRALS) "
nodihedrals = True # turn to True if no dihedral scan is needed
degree = 30 # Amount, in degrees, to enumerate torsions if nodihedrals is False

" PARAMETERS FOR RDKIT OPTIMIZATION "
ff = "MMFF" # force field used in the RDKit optimization. Options: MMFF or UFF
etkdg = False # use new ETKDG knowledge-based method instead of distance geometry also needs to be present in RDKIT ENV
seed = int("062609") #random seed (default 062609) for ETKDG
opt_steps_RDKit = 1000
heavyonly = False # If True, H from OH, NH, etc. will not be used to generate conformers (recommended: False with molecules that contain OH groups)
sample = 100 # number of conformers to sample to get non-torsional differences (default 100)

" DEFAULT PARAMETERS FOR ANI1 and xTB OPTIMIZATION "
opt_steps = 1000 # max number of cycles during optimization
opt_fmax = 0.05 # fmax value to achieve optimization

" DEFAULT PARAMETERS ONLY FOR ANI1 OPTIMIZATION "
constraints = None

" DEFAULT PARAMETERS ONLY FOR xTB OPTIMIZATION "
large_sys = False
STACKSIZE = '1G' #set for large system


" OPTIONS FOR METALS, ATOMS WITH UNCOMMON HYBRIDIZATIONS AND NCI COMPLEXES "

" IF A METAL OR AN ATOM WITH UNCOMMON HYBRIDIZATION (i.e. pentacoordinated phosphorus) IS USED"
metal_complex= True # specify True to activate this option
metal = 'Sn' # specify the metal(s) or atom(s) with uncommon hybridization, in the format 'A','B','C'...
complex_coord = 5 # specify the coordination number of the atom
complex_type = '' # specify the following: square planar, square pyrimidal (otheriwse defaults to octahedral, Td)
m_oxi = 4 # oxidation number of the atom (it is used to calculate the charge of the molecule)
complex_spin = 1 # final spin of the molecule (the code does not calculate spin, it must be defined by the user)

" IF A NCI COMPLEX IS USED "
nci_complex = False # specify  true if NCI complex


" OPTIONS FOR COM FILE GENERATION "

" ONLY LOWEST ENERGY CONFORMER REQUIRED"
lowest_only = False
lowest_n  = False # for a given threshold of energy_threshold_for_gaussian
energy_threshold_for_gaussian = 100  #in kJ/ mol, from all the conformers generated after xTB optimization
                                    # lowest_n must be True to apply this energy threshold

" DEFINITION OF A SECOND CATEGORY OF ATOMS SEPARATED IN GENECP "
genecp_atoms = []

" DEFINTION OF BASIS SET AND LEVEL OF THEORY AND SOLVENT "
basis_set = ['def2svp']
basis_set_genecp_atoms = ['LANL2DZ']
level_of_theory = ['wb97xd']
max_cycle_opt = 100 #eefault is 300

" DISPERSION CORRECTION FOR COM FILES " 
dispersion_correction = False
empirical_dispersion = 'GD3BJ'

" SOLVATION MODEL. Options: gas_phase or any solvation model (i.e. SMD, IEFPCM, CPCM)"
solvent_model = 'IEFPCM'
solvent_name = 'Chloroform'

"DEFAULT PARAMTERS FOR GAUSSIAN OPTIMIZATION "
chk = False
nprocs = 36
mem='60GB'


" OPTIONS FOR THE AUTOMATED WORKFLOW "

" TURN ON AUTOMATED SUBMISSION OF JOBS (AUTOMATED WORKFLOW) "
qsub = False
submission_command = 'qsub_summit'

" MOLECULES now, for eg., molecule list, for later can use as total no. of molecules it is need in the boltz part to read in specific molecules"
maxnumber = 103 #max number in your database
