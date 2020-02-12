"""

* In this file, the variables to helper programs are collected.
* You must make sure that all the variables are correct before launching db_gen.py.

* OTHER variables USED THROUGHOUT THE PROGRAM ARE ALSO SET HERE.
"""
##add imports for xTB

"TYPE OF OPTIMIZATION"
# Options: xTB, AN1  Default : RDKIT optimizaiton
ANI1ccx = False
xtb = False
#enso = False need to add in

" OPTIMIZATION REQUIRED OR NOT"
opt_ax = True # switch to off for single point only
opt_precision_ax = 1E-3 # toggle for optimization convergence

" SINGLE POINTS vs  FULL OPTIMIZATION"
single_point = False

" ONLY LOWEST ENERGY CONFORMER REQUIRED"
lowest_only = False
lowest_n  = False # for a given threshold of energy_threshold_for_gaussian
energy_threshold_for_gaussian = 2.0 #in kJ/ mol

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

" DEFINITION OF ATOMS"
possible_atoms = ['N', 'P', 'As', 'C', 'Si', 'Ge', 'B', 'H', 'S', 'O', 'Se', 'F', 'Br', 'Cl', 'I', 'Ir']
genecp_atoms = ['I','Ir']

"DEFINTION OF BASIS SET AND LEVEL OF THEORY AND SOLVENT"
basis_set = ['LANL2DZ', 'LANL2TZ', '6-31g*']
basis_set_genecp_atoms = ['LANL2DZ','LANL2TZ','LANL2DZ']
level_of_theory = ['wb97xd']

#dispersion correction to be added or not
dispersion_correction = False
empirical_dispersion = 'D3BJ'

# Specify the solvation model. Options: gas_phase or any solvation model (i.e. SMD, IEFPCM, CPCM)
solvent_model = 'SMD'
solvent_name = 'Acetonitrile'

#definition of input lines
if dispersion_correction == True:
    if solvent_model == 'gas_phase':
        input = 'opt freq=noraman EmpiricalDispersion=G{0}'.format(empirical_dispersion)
        input_sp = 'nmr=giao EmpiricalDispersion=G{0}'.format(empirical_dispersion)  #input for single point nmr
    else :
        input = 'opt freq=noraman SCRF=({0},Solvent={1}) EmpiricalDispersion=G{2}'.format(solvent_model, solvent_name,empirical_dispersion ) #add solvent if needed
        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao EmpiricalDispersion=G{2}'.format(solvent_model, solvent_name, empirical_dispersion)  ##add solvent if needed
else:
    if solvent_model == 'gas_phase':
        input = 'opt freq=noraman '
        input_sp = 'nmr=giao ' #input for single point nmr
    else :
        input = 'opt freq=noraman SCRF=({0},Solvent={1})'.format(solvent_model, solvent_name) #add solvent if needed
        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao'.format(solvent_model, solvent_name)  ##add solvent if needed


"DEFAULT PARAMTERS FOR GAUSSIAN OPTIMIZATION"
chk = False
nprocs=24
mem='96GB'

"TURN ON SUBMISSION OF JOBS"
QSUB = False
submission_command = 'qsub_summit'

" MOLECULES now, for eg., molecule list, for later can use as total no. of molecules it is need in the boltz part to read in specific molecules"
maxnumber = 100 #max number in your database
prefix = 'RE' #name syntax = prefix_maxnumber_confs_confsnumber.sdf

"THERMODYNAMIC DATA CALCULATED FROM GOODVIBES"
columns = ['Structure', 'E', 'ZPE', 'H', 'T.S', 'T.qh-S', 'G(T)', 'qh-G(T)']
