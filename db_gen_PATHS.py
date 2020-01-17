"""

* In this file, the paths to helper programs are collected.
* You must make sure that all the paths are correct before launching db_gen.py.

* OTHER CONSTANTS USED THROUGHOUT THE PROGRAM ARE ALSO SET HERE.

* ALL THE TEMPLATES ARE SET HERE AS WELL. (EX: template for gaussian input)

"""

# Since there are many imports, we could save some CPU time if we moved
# some of these import inside the functions that need them
import glob,os
import numpy as np
import pandas as pd

from rdkit import Chem,DataStructs
from rdkit.Chem import rdChemReactions,AllChem,Lipinski,Descriptors
import openbabel as ob
import subprocess

from confgen_noargs import *

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
basis_set = ['6-31g**','6-31+g**','def2tzvp']
basis_set_genecp_atoms = 'LANL2DZ'
level_of_theory = ['M062X']

"DEFAULT PARAMTERS FOR GAUSSIAN OPTIMIZATION"
chk = False
nprocs=24
mem='96GB'
input = 'opt freq=noraman scrf=(smd, solvent=chloroform)'

" MAIN FUCNTION WORKING WITH MOL OBJECT TO CREATE CONFORMERS"
def conformer_generation(mol,name,args):
    valid_structure = filters(mol)
    if valid_structure:
        if args.verbose: print("o  Input Molecule:", name)

        try:
            # the conformational search
            summ_search(mol, name,args)

            # the multiple minimization returns a list of separate mol objects
            conformers, energies = mult_min(mol, name, args)
            print(energies)

            #only print if within predefined energy window
            if len(conformers) > 0:
                # list in energy order
                cids = list(range(len(conformers)))
                sortedcids = sorted(cids, key = lambda cid: energies[cid])

                #print(name, *sorted(energies), sep = ", ")
                print("o ", name, ":", np.array(sorted(energies)))
                #print(["{0:0.2f}".format(i) for i in sorted(energies)])

                sdwriter = Chem.SDWriter(name+final_output)
                glob_min = min(energies)

                write_confs = 0
                for cid in sortedcids:
                    if ANI1ccx == True:
                        if (energies[cid] - glob_min) < ewin / 2625.5:
                            sdwriter.write(conformers[cid])
                            write_confs += 1

                    elif xtb == True:
                        if (energies[cid] - glob_min) < ewin / 2625.5:
                            sdwriter.write(conformers[cid])
                            write_confs += 1

                    else:
                        if (energies[cid] - glob_min) < ewin:
                            sdwriter.write(conformers[cid])
                            write_confs += 1


                if args.verbose == True: print("o  Wrote", write_confs, "conformers to file", name+final_output, "\n")
                sdwriter.close()
            else: print("x  No conformers found!\n")

        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e: print(traceback.print_exc())
    else: print("ERROR: The structure is not valid")

" FILTER TO BE APPLIED FOR SMILES"
def filters(mol):
    valid_structure = True
    # First filter: number of rotatable bonds
    if Lipinski.NumRotatableBonds(mol) < max_torsions:
        # Second filter: molecular weight
        if Descriptors.MolWt(mol) < max_MolWt:
            # Third filter: this filters salts off (2 separated components)
            if len(Chem.MolToSmiles(mol).split('.')) == 1:
                for atom in mol.GetAtoms():
                    #Fourth filter: atoms outside the scope chosen in 'possible_atoms'
                    if atom.GetSymbol() not in possible_atoms:
                        valid_structure = False
                    elif atom.GetSymbol() in genecp_atoms:
                        genecp = True
            else: valid_structure = False
        else: valid_structure = False
    else: valid_structure = False
    return valid_structure


" MAIN FUNCTION TO CREATE GAUSSIAN JOBS"
def write_gaussian_input_file(file, name,lot, bs):

    #pathto change to
    path_write_gjf_files = 'gaussian/' + str(lot) + '-' + str(bs)
    os.chdir(path_write_gjf_files)

    path_for_file = '../../'

    gjf = '{0}_.gjf'.format(name)

    #chk option
    if chk == True:
        header = [
            '%chk={}.chk'.format(name),
            '%MEM={}'.format(mem),
            '%nprocshared={}'.format(nprocs),
            '# {0}/{1} '.format(lot, bs) + input]
    else:
        header = [
            '%MEM={}'.format(mem),
            '%nprocshared={}'.format(nprocs),
            '# {0}/{1} '.format(lot, bs) + input]


    subprocess.run(
          ['obabel', '-isdf', path_for_file+file, '-ogjf', '-O'+gjf,'-m', '-xk', '\n'.join(header)])

    os.chdir(path_for_file)

#################################needes to be edited###################
" GENECP added to the INPUT FILES "
def genecp_for_files(file):

    ecp_list,ecp_I = [],False
    read_lines = open(file,"r").readlines()
    # Detect if there are I atoms to use genecp or not (to use gen)
    for i in range(4,len(read_lines)):
        if read_lines[i].split(' ')[0] not in ecp_list and read_lines[i].split(' ')[0] in possible_atoms:
            ecp_list.append(read_lines[i].split(' ')[0])
        if read_lines[i].split(' ')[0] == 'I':
           ecp_I = True
    if ecp_I == False:
        genecp = 'gen'
    if ecp_I == True:
        genecp = 'genecp'

    with open(file, 'a+') as file:
        file.append('\nI     0\n')
        fileout.append(basis_set_genecp_atoms+'\n')
        fileout.append('****\n\n')
        fileout.append('I 0\n')
        fileout.append(basis_set_genecp_atoms+'\n\n')

##############################################################################

" WRITING TO EXCEL FILE "
