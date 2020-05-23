#!/usr/bin/env python

"""#####################################################.
# 		   This file stores all the functions 		    #
# 	    	  used in conformer generation			    #
######################################################"""

import math
import os
import sys
import subprocess
import glob
import shutil
import time
import yaml
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, rdMolAlign, Lipinski, Descriptors
from rdkit.Geometry import Point3D
from progress.bar import IncrementalBar
from pyconfort.writers import write_confs

# imports for xTB and ANI1
try:
	import ase
	import ase.optimize
	from ase.units import Hartree
	import torch
	os.environ['KMP_DUPLICATE_LIB_OK']='True'
	device = torch.device('cpu')
except:
	print('0')
try:
	from lib.xtb import GFN2
except:
	print('1')
try:
	import torchani
	model = torchani.models.ANI1ccx()
except:
	print('5')

hartree_to_kcal = 627.509

possible_atoms = ["", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
				 "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
				 "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
				 "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
				 "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
				 "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
				 "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
				 "Rg", "Uub", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"]

columns = ['Structure', 'E', 'ZPE', 'H', 'T.S', 'T.qh-S', 'G(T)', 'qh-G(T)']

# SUBSTITUTION WITH I
def substituted_mol(smi,args,log):
	mol = Chem.MolFromSmiles(smi)
	for atom in mol.GetAtoms():
		if atom.GetSymbol() in args.metal:
			args.metal_sym.append(atom.GetSymbol() )
			atom.SetAtomicNum(53)
			if len(atom.GetNeighbors()) == 2:
				atom.SetFormalCharge(-3)
			if len(atom.GetNeighbors()) == 3:
				atom.SetFormalCharge(-2)
			if len(atom.GetNeighbors()) == 4:
				atom.SetFormalCharge(-1)
			if len(atom.GetNeighbors()) == 5:
				atom.SetFormalCharge(0)
			if len(atom.GetNeighbors()) == 6:
				atom.SetFormalCharge(1)
			if len(atom.GetNeighbors()) == 7:
				atom.SetFormalCharge(2)
			if len(atom.GetNeighbors()) == 8:
				atom.SetFormalCharge(3)
			args.metal_idx.append(atom.GetIdx())
			args.complex_coord.append(len(atom.GetNeighbors()))

	return mol,args.metal_idx,args.complex_coord,args.metal_sym

def clean_args(args,ori_ff,smi):
	mol = Chem.MolFromSmiles(smi)
	for atom in mol.GetAtoms():
		if atom.GetSymbol() in args.metal:
			args.metal_complex= True
			break
	else:
		args.metal_complex = False
	args.ff = ori_ff
	args.metal_idx = []
	args.complex_coord = []
	args.metal_sym = []

def compute_confs(smi, name,args,log,dup_data,counter_for_template,i,start_time):
	#taking largest component for salts
	pieces = smi.split('.')
	if len(pieces) > 1:
		# take largest component by length
		smi = max(pieces, key=len)

	# Converts each line to a rdkit mol object
	if args.verbose:
		log.write("   -> Input Molecule {} is {}".format(i, smi))

	if args.metal_complex:
		mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(smi,args,log)
	else:
		mol = Chem.MolFromSmiles(smi)

	if args.metal_complex:
		# get manually for square planar and squarepyramidal
		if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyramidal':
			mol_objects = []
			if len(args.metal_idx) == 1:
				file_template = os.path.dirname(os.path.abspath(__file__)) +'/template/template-4-and-5.sdf'
				temp = Chem.SDMolSupplier(file_template)
				mol_objects_from_template,name, coord_Map, alg_Map, mol_template = template_embed_sp(mol,temp,name,args,log)
				for i,_ in enumerate(mol_objects_from_template):
					mol_objects.append([mol_objects_from_template[i],name[i],coord_Map[i],alg_Map[i],mol_template[i]])
				for [mol, name, coord_Map,alg_Map,mol_template] in mol_objects:
					conformer_generation(mol,name,start_time,args,log,dup_data,counter_for_template,coord_Map,alg_Map,mol_template)
					counter_for_template += 1
			else:
				log.write("x  Cannot use templates for complexes involving more than 1 metal or for organic molecueles.")
		else:
			conformer_generation(mol,name,start_time,args,log,dup_data,i)
	else:
		conformer_generation(mol,name,start_time,args,log,dup_data,i)

# TEMPLATE GENERATION FOR SQUAREPLANAR AND squarepyramidal
def template_embed_sp(molecule,temp,name_input,args,log):
	mol_objects = [] # a list of mol objects that will be populated
	name_return = []
	coord_Map = []

	alg_Map = []
	mol_template = []

	for atom in molecule.GetAtoms():
		if atom.GetIdx() in args.metal_idx:
			if len(atom.GetBonds()) == 5:
				atom.SetAtomicNum(14)
				atom.SetFormalCharge(1)
			if len(atom.GetBonds()) == 4:
				atom.SetAtomicNum(14)
			center_idx = atom.GetIdx()
			neighbours = atom.GetNeighbors()

	number_of_neighbours = len(neighbours)

	if number_of_neighbours == 4:
		#three cases for square planar
		for name in range(3):
			#assigning neighbours
			for atom in molecule.GetAtoms():
				if atom.GetIdx() == center_idx:
					neighbours = atom.GetNeighbors()

			#assigning order of replacement
			if name == 0:
				j = [1,2,3]
			elif name == 1:
				j = [2,3,1]
			elif name == 2:
				j = [3,1,2]

			#checking for same atom neighbours and assigning in the templates for all mols in suppl!
			for mol_1 in temp:
				for atom in mol_1.GetAtoms():
					if atom.GetSymbol() == 'F':
						mol_1 = Chem.RWMol(mol_1)
						idx = atom.GetIdx()
						mol_1.RemoveAtom(idx)
						mol_1 = mol_1.GetMol()

				site_1,site_2,site_3,site_4,metal_site  = 0,0,0,0,0
				for atom in mol_1.GetAtoms():
					if atom.GetIdx() == 4 and metal_site == 0:
						atom.SetAtomicNum(14)
						metal_site = 1
					if atom.GetIdx() == 0 and site_1 == 0:
						atom.SetAtomicNum(neighbours[0].GetAtomicNum())
						site_1 = 1
					if atom.GetIdx() == 3 and site_2 == 0:
						atom.SetAtomicNum(neighbours[j[0]].GetAtomicNum())
						site_2 = 1
					if atom.GetIdx() == 2 and site_3 == 0:
						atom.SetAtomicNum(neighbours[j[1]].GetAtomicNum())
						site_3 = 1
					if atom.GetIdx() == 1 and site_4 == 0 :
						atom.SetAtomicNum(neighbours[j[2]].GetAtomicNum())
						site_4 = 1

				#embedding of the molecule onto the core
				molecule_new, coordMap, algMap = template_embed_optimize(molecule,mol_1,args,log)

				#writing to mol_object file
				name_final = name_input + str(name)
				mol_objects.append(molecule_new)
				name_return.append(name_final)
				coord_Map.append(coordMap)
				alg_Map.append(algMap)
				mol_template.append(mol_1)

	if number_of_neighbours == 5:
		#fifteen cases for square pyrimidal
		for name_1 in range(5):
			for name_2 in range(3):
				#assigning neighbours
				for atom in molecule.GetAtoms():
					if atom.GetIdx() == center_idx:
						neighbours = atom.GetNeighbors()

				# assigning order of replacement for the top
				if name_1 == 0:
					k = 4
				elif name_1== 1:
					k = 3
				elif name_1 == 2:
					k = 2
				elif name_1== 3:
					k = 1
				elif name_1 == 4:
					k = 0

				# assigning order of replacement for the plane
				if name_2 == 0 and k == 4:
					j = [1,2,3]
				elif name_2 == 1 and k == 4:
					j = [2,3,1]
				elif name_2 == 2 and k == 4:
					j = [3,1,2]

				# assigning order of replacement for the plane
				if name_2 == 0 and k == 3:
					j = [1,2,4]
				elif name_2 == 1 and k == 3:
					j = [2,4,1]
				elif name_2 == 2 and k == 3:
					j = [4,1,2]

				# assigning order of replacement for the plane
				if name_2 == 0 and k == 2:
					j = [1,4,3]
				elif name_2 == 1 and k == 2:
					j = [4,3,1]
				elif name_2 == 2 and k == 2:
					j = [4,1,3]

				# assigning order of replacement for the plane
				if name_2 == 0 and k == 1:
					j = [4,2,3]
				elif name_2 == 1 and k == 1:
					j = [2,3,4]
				elif name_2 == 2 and k == 1:
					j = [3,4,2]

				# assigning order of replacement for the plane
				if name_2 == 0 and k == 0:
					j = [1,2,3]
				elif name_2 == 1 and k == 0:
					j = [2,3,1]
				elif name_2 == 2 and k == 0:
					j = [3,1,2]

				#checking for same atom neighbours and assigning in the templates for all mols in suppl!
				for mol_1 in temp:
					site_1,site_2,site_3,site_4,site_5,metal_site  = 0,0,0,0,0,0
					for atom in mol_1.GetAtoms():
						if atom.GetIdx()  == 5 and metal_site == 0:
							atom.SetAtomicNum(14)
							atom.SetFormalCharge(1)
							metal_site = 1
						if k!= 0:
							if atom.GetIdx()  == 1 and site_1 == 0:
								atom.SetAtomicNum(neighbours[0].GetAtomicNum())
								site_1 = 1
							elif atom.GetIdx()  == 2 and site_2 == 0:
								atom.SetAtomicNum(neighbours[j[0]].GetAtomicNum())
								site_2 = 1
							elif atom.GetIdx()  == 3 and site_3 == 0:
								atom.SetAtomicNum(neighbours[j[1]].GetAtomicNum())
								site_3 = 1
							elif atom.GetIdx()  == 4 and site_4 == 0:
								atom.SetAtomicNum(neighbours[j[2]].GetAtomicNum())
								site_4 = 1
							elif atom.GetIdx()  == 0 and site_5 == 0:
								atom.SetAtomicNum(neighbours[k].GetAtomicNum())
								site_5 = 1
						elif k == 0:
							if atom.GetIdx()  == 1 and site_1 == 0:
								atom.SetAtomicNum(neighbours[4].GetAtomicNum())
								site_1 = 1
							elif atom.GetIdx()  == 2 and site_2 == 0:
								atom.SetAtomicNum(neighbours[j[0]].GetAtomicNum())
								site_2 = 1
							elif atom.GetIdx()  == 3 and site_3 == 0:
								atom.SetAtomicNum(neighbours[j[1]].GetAtomicNum())
								site_3 = 1
							elif atom.GetIdx()  == 4 and site_4 == 0:
								atom.SetAtomicNum(neighbours[j[2]].GetAtomicNum())
								site_4 = 1
							elif atom.GetIdx() == 0 and site_5 == 0:
								atom.SetAtomicNum(neighbours[0].GetAtomicNum())
								site_5 = 1

					#assigning and embedding onto the core
					molecule_new, coordMap, algMap = template_embed_optimize(molecule,mol_1,args,log)

					#writing to mol_object file
					name_final = name_input + str(name_1)+ str(name_2)
					mol_objects.append(molecule_new)
					name_return.append(name_final)
					coord_Map.append(coordMap)
					alg_Map.append(algMap)
					mol_template.append(mol_1)

	return mol_objects, name_return, coord_Map, alg_Map, mol_template

# TEMPLATE EMBED OPTIMIZE
def template_embed_optimize(molecule_embed,mol_1,args,log):

	#assigning and embedding onto the core
	num_atom_match = molecule_embed.GetSubstructMatch(mol_1)

	#add H's to molecule
	molecule_embed = Chem.AddHs(molecule_embed)

	#definition of coordmap, the coreconfID(the firstone =-1)
	coordMap = {}
	coreConfId=-1
	randomseed=-1
	force_constant=10000

	# This part selects which atoms from molecule are the atoms of the core
	try:
		coreConf = mol_1.GetConformer(coreConfId)
	except:
		pass
	for k, idxI in enumerate(num_atom_match):
		core_mol_1 = coreConf.GetAtomPosition(k)
		coordMap[idxI] = core_mol_1

	ci = rdDistGeom.EmbedMolecule(molecule_embed, coordMap=coordMap, randomSeed=randomseed)
	if ci < 0:
		log.write('Could not embed molecule.')


	GetFF = Chem.UFFGetMoleculeForceField(molecule_embed,confId=-1)

	#algin molecule to the core
	algMap = [(k, l) for l, k in enumerate(num_atom_match)]

	for k, idxI in enumerate(num_atom_match):
		for l in range(k + 1, len(num_atom_match)):
			idxJ = num_atom_match[l]
			d = coordMap[idxI].Distance(coordMap[idxJ])
			GetFF.AddDistanceConstraint(idxI, idxJ, d, d, force_constant)
	GetFF.Initialize()
	GetFF.Minimize(maxIts=args.opt_steps_RDKit)
	# rotate the embedded conformation onto the core_mol:
	rdMolAlign.AlignMol(molecule_embed, mol_1, atomMap=algMap,reflect=True,maxIters=100)

	return molecule_embed, coordMap, algMap

# FUCNTION WORKING WITH MOL OBJECT TO CREATE CONFORMERS
def conformer_generation(mol,name,start_time,args,log,dup_data,dup_data_idx,coord_Map=None,alg_Map=None,mol_template=None):
	valid_structure = filters(mol, args,log)
	if valid_structure:
		if args.verbose:
			log.write("\n   ----- {} -----".format(name))

		try:
			# the conformational search
			gen = summ_search(mol, name,args,log,dup_data,dup_data_idx,coord_Map,alg_Map,mol_template)
			if gen != -1:
				if args.nodihedrals:
					if args.ANI1ccx:
						mult_min(name+'_'+'rdkit', args, 'ani',log,dup_data,dup_data_idx)
					if args.xtb:
						mult_min(name+'_'+'rdkit', args, 'xtb',log,dup_data,dup_data_idx)
				else:
					if args.ANI1ccx:
						if gen != 0:
							mult_min(name+'_'+'rdkit'+'_'+'rotated', args, 'ani',log,dup_data,dup_data_idx)
						else:
							log.write('\nx   No rotable dihydrals found. Using the non-rotated SDF for ANI')
							mult_min(name+'_'+'rdkit', args, 'ani',log,dup_data,dup_data_idx)
					if args.xtb:
						if gen !=0:
							mult_min(name+'_'+'rdkit'+'_'+'rotated', args, 'xtb',log,dup_data,dup_data_idx)
						else:
							log.write('\nx   No rotable dihydrals found. Using the non-rotated SDF for xTB')
							mult_min(name+'_'+'rdkit', args, 'xtb',log,dup_data,dup_data_idx)
			else:
				pass
		except (KeyboardInterrupt, SystemExit):
			raise
	else:
		log.write("ERROR: The structure is not valid")

	# removing temporary files
	temp_files = ['gfn2.out', 'xTB_opt.traj', 'ANI1_opt.traj', 'wbo', 'xtbrestart']
	for file in temp_files:
		if os.path.exists(file):
			os.remove(file)

	if args.time:
		log.write("\n Execution time: %s seconds" % (round(time.time() - start_time,2)))
		dup_data.at[dup_data_idx, 'time (seconds)'] = round(time.time() - start_time,2)

# RULES TO GET EXPERIMENTAL CONFORMERS
def exp_rules_output(mol, args,log):
	passing = True
	ligand_links = []
	atom_indexes = []
	for atom in mol.GetAtoms():
		# Finds the Ir atom and gets the atom types and indexes of all its neighbours
		if atom.GetSymbol() in args.metal:
			atomic_number = possible_atoms.index(atom.GetSymbol())
			atom.SetAtomicNum(atomic_number)
	for atom in mol.GetAtoms():
		if atom.GetAtomicNum() == atomic_number:
			metal_idx = atom.GetIdx()
			for x in atom.GetNeighbors():
				ligand_links.append(x.GetSymbol())
				atom_indexes.append(x.GetIdx())
	# I need to get the only 3D conformer generated in that mol object for rdMolTransforms
	mol_conf = mol.GetConformer(0)
	# This part will identify the pairs of C and N atoms that are part of the same Ph_Py ligand.
	# The shape of the atom pairs is '[[C1_ATOM_NUMBER, N1_ATOM_NUMBER],[C2, N2],...]'.
	# This information is required for the subsequent filtering process based on angles
	if len(atom_indexes) == args.complex_coord:
		ligand_atoms = []

		for i,_ in enumerate(atom_indexes):
			# This is a filter that excludes molecules that fell apart during DFT geometry
			# optimization (i.e. a N atom from one of the ligands separated from Ir). The
			# max distance allowed can be tuned in length_filter
			bond_length = rdMolTransforms.GetBondLength(mol_conf,metal_idx,atom_indexes[i])
			if ligand_links[i] == 'P':
				length_filter = 2.60
			else:
				length_filter = 2.25
			if bond_length > length_filter:
				passing = False
				break
			for j,_ in enumerate(atom_indexes):
				# Avoid combinations of the same atom with itself
				if atom_indexes[i] != atom_indexes[j]:
					# We know that the ligands never have 2 carbon atoms bonding the Ir atom. We
					# only use atom_indexes[i] for C atoms, and atom_indexes[j] for the potential
					# N atoms that are part of the same Ph_Py ligand
					if ligand_links[i] == 'C':
						# This part detects the Ir-C bond and breaks it, breaking the Ph_Py ring
						bond = mol.GetBondBetweenAtoms(atom_indexes[i], metal_idx)
						new_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()],addDummies=True, dummyLabels=[(atom_indexes[i], metal_idx)])
						if new_mol.GetAtomWithIdx(atom_indexes[i]).IsInRingSize(5):
							five_mem = True
						else:
							five_mem = False
						# Now, identify whether or not the initial 5-membered ring formed between
						# [-Ir-C-C-C-N-] is broken when we break the Ir-C bond. This works
						# because Ph_Py units bind Ir in the same way always, through 1 C and 1 N
						# that are in the same position, forming a 5-membered ring.
						# If this ring is broken, atom_indexes[j] will not be part of a
						# 5-membered ring (atom.IsInRingSize(5) == False) which means that
						# this atom was initially inside the same ligand as the
						# parent C of atom_indexes[i])
						if not five_mem:
							if not new_mol.GetAtomWithIdx(atom_indexes[j]).IsInRingSize(5):
								bond_2 = mol.GetBondBetweenAtoms(atom_indexes[j], metal_idx)
								new_mol_2 = Chem.FragmentOnBonds(mol, [bond_2.GetIdx()],addDummies=True, dummyLabels=[(atom_indexes[j], metal_idx)])
								#doing backwards as well eg. Ir N bond
								if not new_mol_2.GetAtomWithIdx(atom_indexes[i]).IsInRingSize(5):
									ligand_atoms.append([atom_indexes[i],atom_indexes[j]])
									break
						else:
							if not new_mol.GetAtomWithIdx(atom_indexes[j]).IsInRingSize(5):
								ligand_atoms.append([atom_indexes[i],atom_indexes[j]])
								break
		if passing:
			# This stop variable and the breaks inside the inner loops will make that if there
			# is one angle that does not meet the criteria for valid conformers, the outter (i)
			# and inner (j) loops will stop simultaneously (saves time since the molecule is
			# already an invalid geometry, it does not make sense to keep iterating)
			stop = False
			# For complexes with 3 Ph_Py ligands:
			if len(ligand_atoms) == 3:
				for i,_ in enumerate(ligand_atoms):
					if not stop:
						for j,_ in enumerate(ligand_atoms):
							# the i<=j part avoids repeating atoms, the i != j part avoid angles
							# containing the same number twice (i.e. 4-16-4, this angle will fail)
							if i <= j and i != j:
								# Calculate the angle between 2 N atoms from different Ph_Py ligands.
								# When there are 3 Ph_Py ligands, no 2 N atoms must be in 180 degrees
								angle = rdMolTransforms.GetAngleDeg(mol_conf,ligand_atoms[i][1],metal_idx,ligand_atoms[j][1])
								if (180 - args.angle_off) <= angle <= (180 + args.angle_off):
									passing = False
									break
			# For complexes with 2 Ph_Py ligands + 1 ligand that is not Ph_Py
			if len(ligand_atoms) == 2:
				# Since there are only 2 N atoms, we do not need to include a nested loop
					angle = rdMolTransforms.GetAngleDeg(mol_conf,ligand_atoms[0][1],metal_idx,ligand_atoms[1][1])
					# Calculate the angle between 2 N atoms from different Ph_Py ligands.
					# When there are 2 Ph_Py ligands, the 2 N atoms from the 2 Ph_Py ligands
					# must be in 180 degrees
					if (180 - args.angle_off) <= angle <= (180 + args.angle_off):
						pass
					else:
						passing = False
	# This is a second filter that excludes molecules that fell apart during DFT geometry
	# optimization (i.e. a N atom from one of the ligands separated from Ir). In this case,
	# it filters off molecules that the SDF only detects 5 Ir neighbours
	else:
		passing = False
	return passing

# FILTER TO BE APPLIED FOR SMILES
def filters(mol,args,log):
	valid_structure = True
	# Second filter: molecular weight
	if Descriptors.MolWt(mol) < args.max_MolWt:
		# Third filter: this filters salts off (2 separated components)
		#if len(Chem.MolToSmiles(mol).split('.')) == 1:
		for atom in mol.GetAtoms():
			#Fourth filter: atoms outside the scope chosen in 'possible_atoms'
			if atom.GetSymbol() not in possible_atoms:
				valid_structure = False
				if args.verbose:
					log.write(" Exiting as atom isn't in atoms in the periodic table")
	else:
		valid_structure = False
		if args.verbose:
			log.write(" Exiting as total molar mass > {0}".format(args.max_MolWt))
	return valid_structure

# PARSES THE ENERGIES FROM SDF FILES
def read_energies(file,log): # parses the energies from sdf files - then used to filter conformers
	energies = []
	f = open(file,"r")
	readlines = f.readlines()
	for i,_ in enumerate(readlines):
		if readlines[i].find('>  <Energy>') > -1:
			energies.append(float(readlines[i+1].split()[0]))
	f.close()
	return energies


# CALCULATES RMSD between two molecules
def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_RMSD,log):
	if heavy:
		 mol1 = Chem.RemoveHs(mol1)
		 mol2 = Chem.RemoveHs(mol2)
	rms = Chem.GetBestRMS(mol1,mol2,c1,c2,maxMatches=max_matches_RMSD)
	return rms

# DETECTS INITIAL NUMBER OF SAMPLES AUTOMATICALLY
def auto_sampling(mult_factor,mol,args,log):
	if args.metal_complex:
		if len(args.metal_idx) > 0:
			mult_factor = mult_factor*3*len(args.metal_idx) # this accounts for possible trans/cis isomers in metal complexes
	auto_samples = 0
	auto_samples += 3*(Lipinski.NumRotatableBonds(mol)) # x3, for C3 rotations
	auto_samples += 3*(Lipinski.NHOHCount(mol)) # x3, for OH/NH rotations
	auto_samples += 3*(Lipinski.NumSaturatedRings(mol)) # x3, for boat/chair/envelope confs
	if auto_samples == 0:
		auto_samples = mult_factor
	else:
		auto_samples = mult_factor*auto_samples
	return auto_samples

# DETECTS DIHEDRALS IN THE MOLECULE
def getDihedralMatches(mol, heavy,log):
	#this is rdkit's "strict" pattern
	pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
	qmol = Chem.MolFromSmarts(pattern)
	matches = mol.GetSubstructMatches(qmol)

	#these are all sets of 4 atoms, uniquify by middle two
	uniqmatches = []
	seen = set()
	for (a,b,c,d) in matches:
		if (b,c) not in seen and (c,b) not in seen:
			if heavy:
				if mol.GetAtomWithIdx(a).GetSymbol() != 'H' and mol.GetAtomWithIdx(d).GetSymbol() != 'H':
					seen.add((b,c))
					uniqmatches.append((a,b,c,d))
			if not heavy:
				if mol.GetAtomWithIdx(c).GetSymbol() == 'C' and mol.GetAtomWithIdx(d).GetSymbol() == 'H':
					pass
				else:
					seen.add((b,c))
					uniqmatches.append((a,b,c,d))
	return uniqmatches

# IF NOT USING DIHEDRALS, THIS REPLACES I BACK TO THE METAL WHEN METAL = TRUE
# AND WRITES THE RDKIT SDF FILES. WITH DIHEDRALS, IT OPTIMIZES THE ROTAMERS
def genConformer_r(mol, conf, i, matches, degree, sdwriter,args,name,log):
	if i >= len(matches): # base case, torsions should be set in conf
		#setting the metal back instead of I
		if args.metal_complex and args.nodihedrals:
			for atom in mol.GetAtoms():
				if atom.GetIdx() in args.metal_idx:
					re_symbol = args.metal_sym[args.metal_idx.index(atom.GetIdx())]
					atomic_number = possible_atoms.index(re_symbol)
					atom.SetAtomicNum(atomic_number)
		sdwriter.write(mol,conf)
		return 1
	else:
		total = 0
		deg = 0
		while deg < 360.0:
			rad = math.pi*deg / 180.0
			rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
			#recalculating energies after rotation
			if args.ff == "MMFF":
				GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
			elif args.ff == "UFF":
				GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
			else:
				log.write('   Force field {} not supported!'.format(args.ff))
				sys.exit()
			GetFF.Initialize()
			GetFF.Minimize(maxIts=args.opt_steps_RDKit)
			energy = GetFF.CalcEnergy()
			mol.SetProp("Energy",energy)
			mol.SetProp('_Name',name)
			total += genConformer_r(mol, conf, i+1, matches, degree, sdwriter,args,name,log)
			deg += degree
		return total

# AUTOMATICALLY SETS THE CHARGE FOR METAL COMPLEXES
def rules_get_charge(mol,args,log):
	C_group = ['C', 'Se', 'Ge']
	N_group = ['N', 'P', 'As']
	O_group = ['O', 'S', 'Se']
	Cl_group = ['Cl', 'Br', 'I']

	charge = np.empty(len(args.metal_idx), dtype=int)
	neighbours = []
	#get the neighbours of metal atom
	for atom in mol.GetAtoms():
		if atom.GetIdx() in args.metal_idx:
			charge_idx = args.metal_idx.index(atom.GetIdx())
			neighbours = atom.GetNeighbors()
			charge[charge_idx] = args.m_oxi[charge_idx]

			for atom in neighbours:
				#Carbon list
				if atom.GetSymbol() in C_group:
					if atom.GetTotalValence()== 4:
						charge[charge_idx] = charge[charge_idx] - 1
					if atom.GetTotalValence()== 3:
						charge[charge_idx] = charge[charge_idx] - 0
				#Nitrogen list
				if atom.GetSymbol() in N_group:
					if atom.GetTotalValence() == 3:
						charge[charge_idx] = charge[charge_idx] - 1
					if atom.GetTotalValence() == 4:
						charge[charge_idx] = charge[charge_idx] - 0
				#Oxygen list
				if atom.GetSymbol() in O_group:
					if atom.GetTotalValence() == 2:
						charge[charge_idx] = charge[charge_idx] - 1
					if atom.GetTotalValence() == 3:
						charge[charge_idx] = charge[charge_idx] - 0
				#Halogen list
				if atom.GetSymbol() in Cl_group:
					if atom.GetTotalValence() == 1:
						charge[charge_idx] = charge[charge_idx] - 1
					if atom.GetTotalValence() == 2:
						charge[charge_idx] = charge[charge_idx] - 0

	if len(neighbours) == 0:
		#no update in charge as it is an organic molecule
		return args.charge_default
	else:
		return charge

def embed_conf(mol,initial_confs,args,log,coord_Map,alg_Map, mol_template):
	if coord_Map is None and alg_Map is None and mol_template is None:
		cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs,ignoreSmoothingFailures=True, randomSeed=args.seed,numThreads = 0)
		if len(cids) == 0 or len(cids) == 1 and initial_confs != 1:
			log.write("o  Normal RDKit embeding process failed, trying to generate conformers with random coordinates (with "+str(initial_confs)+" possibilities)")
			cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0,ignoreSmoothingFailures=True, numZeroFail=1000, numThreads = 0)
		if args.verbose:
			log.write("o  "+ str(len(cids))+" conformers initially generated")
	# case of embed for templates
	else:
		cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
		if len(cids) == 0 or len(cids) == 1 and initial_confs != 1:
			log.write("o  Normal RDKit embeding process failed, trying to generate conformers with random coordinates (with "+str(initial_confs)+" possibilities)")
			cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0, numZeroFail=1000,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
		if args.verbose:
			log.write("o  "+ str(len(cids))+" conformers initially generated")

	return cids

def min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,coord_Map,alg_Map, mol_template):

	cenergy,outmols = [],[]
	bar = IncrementalBar('o  Minimizing', max = len(cids))
	for i, conf in enumerate(cids):
		if coord_Map is None and alg_Map is None and mol_template is None:
			if args.ff == "MMFF":
				GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
			elif args.ff == "UFF":
				GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
			else:
				log.write('   Force field {} not supported!'.format(args.ff))
				sys.exit()

			GetFF.Initialize()
			GetFF.Minimize(maxIts=args.opt_steps_RDKit)
			energy = GetFF.CalcEnergy()
			cenergy.append(GetFF.CalcEnergy())

		# id template realign before doing calculations
		else:
			num_atom_match = mol.GetSubstructMatch(mol_template)
			# Force field parameters
			if args.ff == "MMFF":
				GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
			elif args.ff == "UFF":
				GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
			else:
				log.write('   Force field {} not supported!'.format(args.ff))
				sys.exit()

			# clean up the conformation
			for k, idxI in enumerate(num_atom_match):
				for l in range(k + 1, len(num_atom_match)):
					idxJ = num_atom_match[l]
					d = coord_Map[idxI].Distance(coord_Map[idxJ])
					GetFF.AddDistanceConstraint(idxI, idxJ, d, d, 10000)
			GetFF.Initialize()
			GetFF.Minimize(maxIts=args.opt_steps_RDKit)
			energy = GetFF.CalcEnergy()
			# rotate the embedded conformation onto the core_mol:
			rdMolAlign.AlignMol(mol, mol_template, prbCid=conf,refCid=-1,atomMap=alg_Map,reflect=True,maxIters=100)
			cenergy.append(energy)

		# outmols is gonna be a list containing "initial_confs" mol objects with "initial_confs"
		# conformers. We do this to SetProp (Name and Energy) to the different conformers
		# and log.write in the SDF file. At the end, since all the mol objects has the same
		# conformers, but the energies are different, we can log.write conformers to SDF files
		# with the energies of the parent mol objects. We measured the computing time and
		# it's the same as using only 1 parent mol object with 10 conformers, but we couldn'temp
		# SetProp correctly
		pmol = PropertyMol.PropertyMol(mol)
		outmols.append(pmol)
		bar.next()

	bar.finish()

	for i, cid in enumerate(cids):
		outmols[cid].SetProp('_Name', name + ' conformer ' + str(i+1))
		outmols[cid].SetProp('Energy', cenergy[cid])

	cids = list(range(len(outmols)))
	sorted_all_cids = sorted(cids,key = lambda cid: cenergy[cid])

	sortedcids,nhigh_rdkit=[],0
	for i,cid in enumerate(sorted_all_cids):
		if i == 0:
			cenergy_min = cenergy[cid]
		if abs(cenergy[cid] - cenergy_min) < args.ewin_rdkit:
			sortedcids.append(cid)
		else:
			nhigh_rdkit +=1
	if args.verbose:
		log.write("o  "+str(nhigh_rdkit)+ "  Conformers rejected based on energy (E > "+str(args.ewin_rdkit)+" kcal/mol)")

	dup_data.at[dup_data_idx, 'RDKit-energy-window'] = nhigh_rdkit

	log.write("\n\no  Filters after intial embedding of "+str(initial_confs)+" conformers")
	selectedcids,selectedcids_initial, eng_dup,eng_rms_dup =[],[],-1,-1
	bar = IncrementalBar('o  Filtering based on energy (pre-filter)', max = len(sortedcids))
	for i, conf in enumerate(sortedcids):
		# This keeps track of whether or not your conformer is unique
		excluded_conf = False
		# include the first conformer in the list to start the filtering process
		if i == 0:
			selectedcids_initial.append(conf)
		# check rmsd
		for seenconf in selectedcids_initial:
			E_diff = abs(cenergy[conf] - cenergy[seenconf]) # in kcal/mol
			if E_diff < args.initial_energy_threshold:
				eng_dup += 1
				excluded_conf = True
				break
		if not excluded_conf:
			if conf not in selectedcids_initial:
				selectedcids_initial.append(conf)
		bar.next()
	bar.finish()


	if args.verbose:
		log.write("o  "+str(eng_dup)+ " Duplicates removed  pre-energy filter (E < "+str(args.initial_energy_threshold)+" kcal/mol)")
	#reduce to unique set
	if args.verbose:
		log.write("o  Removing duplicate conformers (RMSD < "+ str(args.rms_threshold)+ " and E difference < "+str(args.energy_threshold)+" kcal/mol)")

	bar = IncrementalBar('o  Filtering based on energy and RMSD', max = len(selectedcids_initial))
	#check rmsd
	for i, conf in enumerate(selectedcids_initial):

		# #set torsions to same value
		for m in rotmatches:
			rdMolTransforms.SetDihedralDeg(outmols[conf].GetConformer(conf),*m,180.0)

		# This keeps track of whether or not your conformer is unique
		excluded_conf = False
		# include the first conformer in the list to start the filtering process
		if i == 0:
			selectedcids.append(conf)
		# check rmsd
		for seenconf in selectedcids:
			E_diff = abs(cenergy[conf] - cenergy[seenconf]) # in kcal/mol
			if  E_diff < args.energy_threshold:
				rms = get_conf_RMS(outmols[conf],outmols[conf],seenconf,conf, args.heavyonly, args.max_matches_RMSD,log)
				if rms < args.rms_threshold:
					excluded_conf = True
					eng_rms_dup += 1
					break
		if not excluded_conf:
			if conf not in selectedcids:
				selectedcids.append(conf)
		bar.next()
	bar.finish()

	if args.verbose:
		log.write("o  "+str(eng_rms_dup)+ " Duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol) after rotation")
	if args.verbose:
		log.write("o  "+ str(len(selectedcids))+" unique conformers remain")

	dup_data.at[dup_data_idx, 'RDKit-energy-duplicates'] = eng_dup
	dup_data.at[dup_data_idx, 'RDKit-RMSD-and-energy-duplicates'] = eng_rms_dup
	dup_data.at[dup_data_idx, 'RDKIT-Unique-conformers'] = len(selectedcids)

	# writing charges after RDKIT
	args.charge = rules_get_charge(mol,args,log)
	dup_data.at[dup_data_idx, 'Overall charge'] = np.sum(args.charge)

	# now exhaustively drive torsions of selected conformers
	n_confs = int(len(selectedcids) * (360 / args.degree) ** len(rotmatches))
	if args.verbose and len(rotmatches) != 0:
		log.write("\n\no  Systematic generation of "+ str(n_confs)+ " confomers")
		bar = IncrementalBar('o  Generating conformations based on dihedral rotation', max = len(selectedcids))
	else:
		bar = IncrementalBar('o  Generating conformations', max = len(selectedcids))

	total = 0
	for conf in selectedcids:
		total += genConformer_r(outmols[conf], conf, 0, rotmatches, args.degree, sdwriter ,args,outmols[conf].GetProp('_Name'),log)
		bar.next()
	bar.finish()
	if args.verbose and len(rotmatches) != 0:
		log.write("o  %d total conformations generated"%total)
	status = 1

	if not args.nodihedrals:
		dup_data.at[dup_data_idx, 'RDKIT-Rotated-conformers'] = total

	return status

def filter_after_rotation(args,name,log,dup_data,dup_data_idx):
	rdmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)
	if rdmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)

	writer_mol_objects = []
	bar = IncrementalBar('o  Filtering based on energy and rms after rotation of dihedrals', max = len(rdmols))

	rd_count = 0
	rd_selectedcids,rd_dup_energy,rd_dup_rms_eng =[],-1,0
	for i, rd_mol_i in enumerate(rdmols):
		mol_rd = Chem.RWMol(rd_mol_i)
		mol_rd.SetProp('_Name',rd_mol_i.GetProp('_Name')+' '+str(i))
		# This keeps track of whether or not your conformer is unique
		excluded_conf = False
		# include the first conformer in the list to start the filtering process
		if rd_count == 0:
			rd_selectedcids.append(rd_mol_i)
			if args.metal_complex:
				for atom in mol_rd.GetAtoms():
					if atom.GetIdx() in args.metal_idx:
						re_symbol = args.metal_sym[args.metal_idx.index(atom.GetIdx())]
						atomic_number = possible_atoms.index(re_symbol)
						atom.SetAtomicNum(atomic_number)
			writer_mol_objects.append([mol_rd,float(mol_rd.GetProp('Energy'))])
		# Only the first ID gets included
		rd_count = 1
		# check rmsd
		for rd_mol_j in rd_selectedcids:
			if abs(float(rd_mol_i.GetProp('Energy')) - float(rd_mol_j.GetProp('Energy'))) < args.initial_energy_threshold: # comparison in kcal/mol
				excluded_conf = True
				rd_dup_energy += 1
				break
			if abs(float(rd_mol_i.GetProp('Energy')) - float(rd_mol_j.GetProp('Energy'))) < args.energy_threshold: # in kcal/mol
				rms = get_conf_RMS(mol_rd,rd_mol_j,-1,-1, args.heavyonly, args.max_matches_RMSD,log)
				if rms < args.rms_threshold:
					excluded_conf = True
					rd_dup_rms_eng += 1
					break
		if not excluded_conf:
			if rd_mol_i not in rd_selectedcids:
				rd_selectedcids.append(rd_mol_i)
				if args.metal_complex:
					for atom in mol_rd.GetAtoms():
						if atom.GetIdx() in args.metal_idx:
							re_symbol = args.metal_sym[args.metal_idx.index(atom.GetIdx())]
							atomic_number = possible_atoms.index(re_symbol)
							atom.SetAtomicNum(atomic_number)
				writer_mol_objects.append([mol_rd,float(mol_rd.GetProp('Energy'))] )
		bar.next()
	bar.finish()

	# writing sorted mol objects
	sdwriter_rd = Chem.SDWriter(name+'_'+'rdkit'+'_'+'rotated'+args.output)
	sortedmols = sorted(writer_mol_objects,key=lambda x: x[1])
	for i, write_mol in enumerate(sortedmols):
		sdwriter_rd.write(write_mol[0])

	sdwriter_rd.close()

	if args.verbose:
		log.write("o  "+str(rd_dup_energy)+ " Duplicates removed initial energy (E < "+str(args.initial_energy_threshold)+" kcal/mol)")
	if args.verbose:
		log.write("o  "+str(rd_dup_rms_eng)+ " Duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol) after rotation")
	if args.verbose:
		log.write("o  "+str(len(rd_selectedcids) )+ " unique conformers remain")

	# filtering process after rotations
	dup_data.at[dup_data_idx, 'RDKIT-Rotated-Unique-conformers'] = len(rd_selectedcids)

	status = 1

	return status

# EMBEDS, OPTIMIZES AND FILTERS RDKIT CONFORMERS
def summ_search(mol, name,args,log,dup_data,dup_data_idx, coord_Map = None,alg_Map=None,mol_template=None):
	sdwriter = Chem.SDWriter(name+'_'+'rdkit'+args.output)

	Chem.SanitizeMol(mol)
	if coord_Map is None and alg_Map is None and mol_template is None:
		mol = Chem.AddHs(mol)
	mol.SetProp("_Name",name)

	# detects and applies auto-detection of initial number of conformers
	if args.sample == 'auto':
		initial_confs = int(auto_sampling(args.auto_sample,mol,args,log))

	else:
		initial_confs = int(args.sample)

	dup_data.at[dup_data_idx, 'Molecule'] = name

	rotmatches = getDihedralMatches(mol, args.heavyonly,log)
	if len(rotmatches) > args.max_torsions:
		log.write("x  Too many torsions (%d). Skipping %s" %(len(rotmatches),(name+args.output)))
		status = -1

	else:
		dup_data.at[dup_data_idx, 'RDKIT-Initial-samples'] = initial_confs
		if args.nodihedrals:
			rotmatches =[]
		cids = embed_conf(mol,initial_confs,args,log,coord_Map,alg_Map, mol_template)
		#energy minimize all to get more realistic results
		#identify the atoms and decide Force Field

		for atom in mol.GetAtoms():
			if atom.GetAtomicNum() > 36: #up to Kr for MMFF, if not the code will use UFF
				args.ff = "UFF"
		if args.verbose:
			log.write("o  Optimizing "+ str(len(cids))+ " initial conformers with "+ args.ff)
		if args.verbose:
			if not args.nodihedrals:
				log.write("o  Found "+ str(len(rotmatches))+ " rotatable torsions")
			else:
				log.write("o  Systematic torsion rotation is set to OFF")

		status = min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,coord_Map,alg_Map, mol_template)
	sdwriter.close()

	if status != -1:
		#getting the energy from and mols after rotations
		if not args.nodihedrals and len(rotmatches) != 0:
			status = filter_after_rotation(args,name,log,dup_data,dup_data_idx)
		elif not args.nodihedrals and len(rotmatches) ==0:
			status = 0

	return status

# xTB AND ANI1 OPTIMIZATIONS
def optimize(mol, args, program,log,dup_data,dup_data_idx):
	# if large system increase stck size
	if args.large_sys:
		os.environ['OMP_STACKSIZE'] = args.STACKSIZE

	# removing the Ba atom if NCI complexes
	if args.nci_complex:
		for atom in mol.GetAtoms():
			if atom.GetSymbol() =='I':
				atom.SetAtomicNum(1)

	if args.metal_complex and not args.nodihedrals:
		for atom in mol.GetAtoms():
			if atom.GetIdx() in args.metal_idx:
				re_symbol = args.metal_sym[args.metal_idx.index(atom.GetIdx())]
				atomic_number = possible_atoms.index(re_symbol)
				atom.SetAtomicNum(atomic_number)

	elements = ''
	ase_metal = []
	ase_metal_idx = []
	for i,atom in enumerate(mol.GetAtoms()):
		if atom.GetIdx() in args.metal_idx:
			ase_metal.append(i)
			ase_metal_idx.append(atom.GetIdx())
		elements += atom.GetSymbol()

	args.charge = rules_get_charge(mol,args,log)
	dup_data.at[dup_data_idx, 'Overall charge'] = np.sum(args.charge)

	cartesians = mol.GetConformers()[0].GetPositions()
	coordinates = torch.tensor([cartesians.tolist()], requires_grad=True, device=device)

	if program == 'ani':
		species = model.species_to_tensor(elements).to(device).unsqueeze(0)
		_, ani_energy = model((species, coordinates))

		ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0], calculator=model.ase())
		### make a function for constraints and optimization
		if args.constraints is not None:
			fb = ase.constraints.FixBondLength(0, 1)
			ase_molecule.set_distance(0,1,2.0)
			ase_molecule.set_constraint(fb)

		optimizer = ase.optimize.BFGS(ase_molecule, trajectory='ANI1_opt.traj')
		optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
		if len(ase.io.Trajectory('xTB_opt.traj', mode='r')) != (args.opt_steps+1):
			species_coords = ase_molecule.get_positions().tolist()
			coordinates = torch.tensor([species_coords], requires_grad=True, device=device)
			converged = 0
		# Now let's compute energy:
		_, ani_energy = model((species, coordinates))
		sqm_energy = ani_energy.item() * hartree_to_kcal # Hartree to kcal/mol

	elif program == 'xtb':
		if args.metal_complex:
			# passing charges metal present
			ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=GFN2()) #define ase molecule using GFN2 Calculator
			if os.path.splitext(args.input)[1] == '.csv' or os.path.splitext(args.input)[1] == '.cdx' or os.path.splitext(args.input)[1] == '.smi':
				for i,atom in enumerate(ase_molecule):
					if i in ase_metal:
						ase_charge = args.charge[args.metal_idx.index(ase_metal_idx[ase_metal.index(i)])]

						# will update only for cdx, smi, and csv formats.
						atom.charge = ase_charge

			else:
				atom.charge = args.charge_default
				if args.verbose:
					log.write('o  The Overall charge is read from the .com file ')
		else:
			ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=GFN2()) #define ase molecule using GFN2 Calculator
		optimizer = ase.optimize.BFGS(ase_molecule, trajectory='xTB_opt.traj',logfile='xtb.opt')
		optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
		if len(ase.io.Trajectory('xTB_opt.traj', mode='r')) != (args.opt_steps+1):
			species_coords = ase_molecule.get_positions().tolist()
			coordinates = torch.tensor([species_coords], requires_grad=True, device=device)
			converged = 0
		# Now let's compute energy:
		xtb_energy = ase_molecule.get_potential_energy()
		sqm_energy = (xtb_energy / Hartree)* hartree_to_kcal

	else:
		log.write('program not defined!')

	energy, converged, cartesians = sqm_energy, converged, np.array(coordinates.tolist()[0])
	# update coordinates of mol object
	for j in range(mol.GetNumAtoms()):
		[x,y,z] = cartesians[j]
		mol.GetConformer().SetAtomPosition(j,Point3D(x,y,z))

	return mol, converged, energy

# xTB AND ANI1 OPTIMIZATION, FILTER AND WRITING SDF FILES
def mult_min(name, args, program,log,dup_data,dup_data_idx):
	inmols = Chem.SDMolSupplier(name+args.output, removeHs=False)
	if inmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)

	globmin, n_high,n_dup_energy, n_dup_rms_eng  = None, 0, 0, 0
	c_converged, c_energy, outmols = [], [], []

	if args.verbose:
		log.write("\n\no  Multiple minimization of "+ name+args.output+ " with "+ program)
	bar = IncrementalBar('o  Minimizing', max = len(inmols))

	for i,mol in enumerate(inmols):
		bar.next()
		conf = 1
		if mol is not None:
			# optimize this structure and record the energy
			mol, converged, energy = optimize(mol, args, program,log,dup_data,dup_data_idx)

			if globmin is None:
				globmin = energy
			if energy < globmin:
				globmin = energy

			if converged == 0 and abs(energy - globmin) < args.ewin_min: # comparison in kcal/mol
				unique = 0

				# compare against all previous conformers located
				for j,seenmol in enumerate(outmols):
					if abs(energy - c_energy[j]) < args.initial_energy_threshold: # comparison in kcal/mol
						unique += 1
						n_dup_energy += 1
						break

					if abs(energy - c_energy[j]) < args.energy_threshold: # comparison in kcal/mol
						rms = get_conf_RMS(mol, seenmol, 0, 0, args.heavyonly, args.max_matches_RMSD,log)
						if rms < args.rms_threshold:
							unique += 1
							n_dup_rms_eng += 1
							break

				if unique == 0:
					pmol = PropertyMol.PropertyMol(mol)
					outmols.append(pmol)
					c_converged.append(converged)
					c_energy.append(energy)
					conf += 1
			else:
				n_high += 1
		else:
			pass #log.write("No molecules to optimize")

	bar.finish()

	if args.verbose:
		log.write("o  "+str( n_dup_energy)+ " Duplicates removed initial energy (E < "+str(args.initial_energy_threshold)+" kcal/mol)")
	if args.verbose:
		log.write("o  "+str( n_dup_rms_eng)+ " Duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol)")
	if args.verbose:
		log.write("o  "+str( n_high)+ " Conformers rejected based on energy (E > "+str(args.ewin_min)+" kcal/mol)")

	# if SQM energy exists, overwrite RDKIT energies and geometries
	cids = list(range(len(outmols)))
	sortedcids = sorted(cids, key = lambda cid: c_energy[cid])

	name_mol = name.split('_rdkit')[0]

	for i, cid in enumerate(sortedcids):
		outmols[cid].SetProp('_Name', name_mol +' conformer ' + str(i+1))
		outmols[cid].SetProp('Energy', c_energy[cid])

	if program == 'xtb':
		dup_data.at[dup_data_idx, 'xTB-Initial-samples'] = len(inmols)
		dup_data.at[dup_data_idx, 'xTB-energy-window'] = n_high
		dup_data.at[dup_data_idx, 'xTB-initial_energy_threshold'] = n_dup_energy
		dup_data.at[dup_data_idx, 'xTB-RMSD-and-energy-duplicates'] = n_dup_rms_eng
		dup_data.at[dup_data_idx, 'xTB-Unique-conformers'] = len(sortedcids)

	if program == 'ani':
		dup_data.at[dup_data_idx, 'ANI1ccx-Initial-samples'] = len(inmols)
		dup_data.at[dup_data_idx, 'ANI1ccx-energy-window'] = n_high
		dup_data.at[dup_data_idx, 'ANI1ccx-initial_energy_threshold'] = n_dup_energy
		dup_data.at[dup_data_idx, 'ANI1ccx-RMSD-and-energy-duplicates'] = n_dup_rms_eng
		dup_data.at[dup_data_idx, 'ANI1ccx-Unique-conformers'] = len(sortedcids)

	# write the filtered, ordered conformers to external file
	write_confs(outmols, c_energy, name, args, program,log)
