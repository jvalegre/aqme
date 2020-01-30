"""

* In this file, the paths to helper programs are collected.
* You must make sure that all the paths are correct before launching db_gen.py.

* OTHER variables USED THROUGHOUT THE PROGRAM ARE ALSO SET HERE.
"""
##add imports for xTB

"TYPE OF OPTIMIZATION"
# Options: xTB, AN1  Default : RDKIT optimizaiton
ANI1ccx = False
xtb = False

" OPTIMIZATION REQUIRED OR NOT"
opt_ax = True # switch to off for single point only
opt_precision_ax = 1E-3 # toggle for optimization convergence

" DEFAULT PARAMETERS FOR RDKIT GENERATION AND FILTERS"
max_torsions = 5 #Skip any molecules with more than this many torsions (default 5)
max_MolWt = 500
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
possible_atoms = ['N', 'P', 'As', 'C', 'Si', 'Ge', 'B', 'H', 'S', 'O', 'Se', 'F', 'Br', 'Cl', 'I']
genecp_atoms = ['I']

"DEFINTION OF BASIS SET AND LEVEL OF THEORY"
#dict_bs = {'M062X': ['6-31g**','6-31+g**','def2tzvp'], 'wb97xd': ['6-31g**','6-31+g**','def2tzvp'], 'B3LYP': ['6-31g** EmpiricalDispersion=GD3BJ', '6-31+g** EmpiricalDispersion=GD3BJ','def2tzvp EmpiricalDispersion=GD3BJ'  ]}
### from dictionary to pandas
#dict_bs_genecp = {'M062X-D3': ['LANL2DZ','LANL2DZ','LANL2DZ'],'wb97xd': ['LANL2TZ'], 'B3LYP': ['LANL2DZ']}
# basis_set_genecp_atoms = 'LANL2DZ'
basis_set = ['6-31g**', '6-31+g**', 'def2tzvp']
basis_set_I = 'LANL2DZ'
level_of_theory = [ 'M062X', 'wb97xd', 'b3lyp-d3']
d3bj = True #now only set for b3lyp
input = 'opt freq=noraman SCRF=(Solvent=Chloroform)' #add solvent if needed

"DEFAULT PARAMTERS FOR GAUSSIAN OPTIMIZATION"
chk = False
nprocs=24
mem='96GB'

" MOLECULES now, for eg., molecule list, for later can use as total no. of molecules"
molecules = [9,13,14,25,55]

"THERMODYNAMIC DATA CALCULATED FROM GOODVIBES"
columns = ['Structure', 'E', 'ZPE', 'H', 'T.S', 'T.qh-S', 'G(T)', 'qh-G(T)']
