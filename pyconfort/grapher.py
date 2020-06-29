#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	    	  used for genrating a graph      	    #
#####################################################.

from pyconfort.confgen_functions import rdkit_sdf_read
from rdkit.Chem import AllChem as Chem
import cclib

def get_energy(inmols_min):
    energy_min = []
    for i,mol in enumerate(inmols_min):
        energy_min.append(["_".join(mol.GetProp('_Name').split(' ')),mol.GetProp('Energy')])
    return energy_min

def rename_name(energy,type):
    for i,_ in enumerate(energy):
        energy[i][0] = energy[i][0].split('_'+type)[0]
    return energy

def graph(sdf_rdkit,sdf_xtb,sdf_ani,log_files,args,log):

    inmols_rdkit = Chem.SDMolSupplier(sdf_rdkit, removeHs=False)
    #get the energy from sdf
    energy_rdkit = get_energy(inmols_rdkit)

    if args.xtb:
        #get energy list for all conformers from sdfs of rdkit and minimize
        inmols_xtb =  Chem.SDMolSupplier(sdf_xtb, removeHs=False)
        energy_xtb = get_energy(inmols_xtb)
        energy_xtb = rename_name(energy_xtb,'xtb')
    if args.ANI1ccx:
        #get energy list for all conformers from sdfs of rdkit and minimize
        inmols_ani = Chem.SDMolSupplier(sdf_ani, removeHs=False)
        energy_ani = get_energy(inmols_ani)
        energy_ani = rename_name(energy_ani,'ani')

    energy_xtb_dft,energy_ani_dft = [],[]
    #get energy from log FILES
    for file in log_files:
        data = cclib.io.ccread(file)
        if len(file.split('_ani.log')) == 2:
            name = file.split('_ani.log')[0]
            energy_ani_dft.append([name,data.scfenergies[0]])
        if len(file.split('_xtb.log')) == 2:
            name = file.split('_xtb.log')[0]
            energy_xtb_dft.append([name,data.scfenergies[0]])


    print(energy_rdkit)
    print(energy_xtb)
    print(energy_ani)
    print(energy_xtb_dft)
    print(energy_ani_dft)
