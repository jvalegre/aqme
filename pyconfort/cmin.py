#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	    	  used in conformer minimization	    #
#####################################################.

import os
import sys
import subprocess
import time
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, rdMolAlign, Lipinski
from rdkit.Geometry import Point3D
from progress.bar import IncrementalBar
from pyconfort.argument_parser import possible_atoms
from pyconfort.qprep_gaussian import write_confs
from pyconfort.filter import filters,set_metal_atomic_number,ewin_filter,pre_E_filter,RMSD_and_E_filter

# imports for xTB and ANI1
ase_installed = True
torch_installed = True
try:
	import ase
	import ase.optimize
	from ase.units import Hartree
except (ModuleNotFoundError,AttributeError):
	ase_installed = False
	print('ASE is not installed correctly - xTB and ANI1ccx are not available')
if ase_installed:
	try:
		import torch
		os.environ['KMP_DUPLICATE_LIB_OK']='True'
		device = torch.device('cpu')
	except (ModuleNotFoundError,AttributeError):
		torch_installed = False
		print('TORCH is not installed correctly - xTB and ANI1ccx are not available')
	if torch_installed:
		try:
			from xtb.ase.calculator import XTB
		except (ModuleNotFoundError,AttributeError):
			print('xTB is not installed correctly - xTB is not available')
		try:
			import torchani
			model = torchani.models.ANI1ccx()
		except (ModuleNotFoundError,AttributeError):
			print('Torchani is not installed correctly - ANI1ccx is not available')

hartree_to_kcal = 627.509
possible_atoms = possible_atoms()

#definition of atom groups for rules to get charge
def atom_groups():
	C_group = ['C', 'Se', 'Ge']
	N_group = ['N', 'P', 'As']
	O_group = ['O', 'S', 'Se']
	F_group = ['Cl', 'Br', 'I']

	return C_group,N_group,O_group,F_group

# AUTOMATICALLY SETS THE CHARGE FOR METAL COMPLEXES
def rules_get_charge(mol,args,log):
	C_group,N_group,O_group,F_group = atom_groups()

	neighbours = []
	charge = np.empty(len(args.metal_idx), dtype=int)
	#get the neighbours of metal atom
	for atom in mol.GetAtoms():
		if atom.GetIdx() in args.metal_idx:
			charge_idx = args.metal_idx.index(atom.GetIdx())
			neighbours = atom.GetNeighbors()
			charge[charge_idx] = args.m_oxi[charge_idx]
			for neighbour in neighbours:
				if neighbour.GetTotalValence()== 4:
					if neighbour.GetSymbol() in C_group:
						charge[charge_idx] = charge[charge_idx] - 1
				elif neighbour.GetTotalValence()== 3:
					if neighbour.GetSymbol() in N_group:
						charge[charge_idx] = charge[charge_idx] - 1
				elif neighbour.GetTotalValence() == 2:
					if neighbour.GetSymbol() in O_group:
						charge[charge_idx] = charge[charge_idx] - 1
				elif neighbour.GetTotalValence() == 1:
					if neighbour.GetSymbol() in F_group:
						charge[charge_idx] = charge[charge_idx] - 1
	if len(neighbours) == 0:
		#no update in charge as it is an organic molecule
		return args.charge_default
	else:
		return charge

# ANI1 OPTIMIZER AND ENERGY CALCULATOR
def ani_calc(elements,cartesians,coordinates,args,log):
	species = model.species_to_tensor(elements).to(device).unsqueeze(0)
	_, ani_energy = model((species, coordinates))


	ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0], calculator=model.ase())

	optimizer = ase.optimize.BFGS(ase_molecule, trajectory='ANI1_opt.traj', logfile='ase.opt' )
	optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
	if len(ase.io.Trajectory('ANI1_opt.traj', mode='r')) != (args.opt_steps+1):
		species_coords = ase_molecule.get_positions().tolist()
		coordinates = torch.tensor([species_coords], requires_grad=True, device=device)

	# Now let's compute energy:
	_, ani_energy = model((species, coordinates))
	sqm_energy = ani_energy.item() * hartree_to_kcal # Hartree to kcal/mol

	return sqm_energy, coordinates

# xTB OPTIMIZER AND ENERGY CALCULATOR
def xtb_calc(elements,cartesians,coordinates,args,log,ase_metal,ase_metal_idx):
	if args.metal_complex:
		# passing charges metal present
		ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=XTB(method=args.xtb_method,accuracy=args.xtb_accuracy,electronic_temperature=args.xtb_electronic_temperature,max_iterations=args.xtb_max_iterations,solvent=args.xtb_solvent)) #define ase molecule using GFN2 Calculator
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
		ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=XTB(method=args.xtb_method,accuracy=args.xtb_accuracy,electronic_temperature=args.xtb_electronic_temperature,max_iterations=args.xtb_max_iterations,solvent=args.xtb_solvent)) #define ase molecule using GFN2 Calculator
	optimizer = ase.optimize.BFGS(ase_molecule, trajectory='xTB_opt.traj',logfile='xtb.opt')
	optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
	if len(ase.io.Trajectory('xTB_opt.traj', mode='r')) != (args.opt_steps+1):
		species_coords = ase_molecule.get_positions().tolist()
		coordinates = torch.tensor([species_coords], requires_grad=True, device=device)
	# Now let's compute energy:
	xtb_energy = ase_molecule.get_potential_energy()
	sqm_energy = (xtb_energy / Hartree)* hartree_to_kcal

	return sqm_energy, coordinates

# xTB AND ANI1 MAIN OPTIMIZATION PROCESS
def optimize(mol, args, program,log,dup_data,dup_data_idx):
	# if large system increase stack size
	if args.STACKSIZE != '1G':
		os.environ['OMP_STACKSIZE'] = args.STACKSIZE

	# removing the Ba atom if NCI complexes
	if args.nci_complex:
		for atom in mol.GetAtoms():
			if atom.GetSymbol() =='I':
				atom.SetAtomicNum(1)

	if args.metal_complex and not args.CSEARCH=='rdkit-dihedral':
		set_metal_atomic_number(mol,args)

	elements = ''
	ase_metal = []
	ase_metal_idx = []
	for i,atom in enumerate(mol.GetAtoms()):
		if atom.GetIdx() in args.metal_idx:
			ase_metal.append(i)
			ase_metal_idx.append(atom.GetIdx())
		elements += atom.GetSymbol()

	if os.path.splitext(args.input)[1] == '.cdx' or os.path.splitext(args.input)[1] == '.smi' or os.path.splitext(args.input)[1] == '.csv':
		args.charge = rules_get_charge(mol,args,log)
		dup_data.at[dup_data_idx, 'Overall charge'] = np.sum(args.charge)
	else:
		dup_data.at[dup_data_idx, 'Overall charge'] = args.charge_default

	cartesians = mol.GetConformers()[0].GetPositions()
	coordinates = torch.tensor([cartesians.tolist()], requires_grad=True, device=device)

	if program == 'ani':
		sqm_energy, coordinates = ani_calc(elements,cartesians,coordinates,args,log)

	elif program == 'xtb':
		sqm_energy, coordinates = xtb_calc(elements,cartesians,coordinates,args,log,ase_metal,ase_metal_idx)

	else:
		log.write('program not defined!')

	energy, cartesians = sqm_energy, np.array(coordinates.tolist()[0])
	# update coordinates of mol object
	for j in range(mol.GetNumAtoms()):
		[x,y,z] = cartesians[j]
		mol.GetConformer().SetAtomPosition(j,Point3D(x,y,z))

	return mol, energy

# read SDF files from RDKit optimization
def rdkit_sdf_read(name, args, log):
	inmols = Chem.SDMolSupplier(name+args.output, removeHs=False)
	if inmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)
	return inmols

# READ FILES FOR xTB AND ANI1 OPTIMIZATION, FILTER AND WRITING SDF FILES
def mult_min(name, args, program,log,dup_data,dup_data_idx):
	# read SDF files from RDKit optimization
	inmols = rdkit_sdf_read(name, args, log)

	cenergy, outmols = [],[]
	if args.verbose:
		log.write("\n\no  Multiple minimization of "+ name+args.output+ " with "+ program)
	bar = IncrementalBar('o  Minimizing', max = len(inmols))

	for i,mol in enumerate(inmols):
		bar.next()
		if mol is not None:
			# optimize this structure and record the energy
			mol,energy = optimize(mol, args, program,log,dup_data,dup_data_idx)
			pmol = PropertyMol.PropertyMol(mol)
			outmols.append(pmol)
			cenergy.append(energy)

	# if SQM energy exists, overwrite RDKIT energies and geometries
	cids = list(range(len(outmols)))
	sorted_all_cids = sorted(cids, key = lambda cid: cenergy[cid])

	name_mol = name.split('_rdkit')[0]
	for i, cid in enumerate(sorted_all_cids):
		outmols[cid].SetProp('_Name', outmols[cid].GetProp('_Name') +' '+program)
		outmols[cid].SetProp('Energy', cenergy[cid])

	#writing all conformers to files after minimization
	sdwriter = Chem.SDWriter(name.split('_rdkit')[0]+'_'+program+'_all_confs'+args.output)

	write_all_confs = 0
	for cid in sorted_all_cids:
		sdwriter.write(outmols[cid])
		write_all_confs += 1

	if args.verbose:
		log.write("\no  Writing "+str(write_all_confs )+ " conformers to file " + name+'_'+program+args.output)
	sdwriter.close()

	log.write("\n\no  Applying filters to intial conformers after "+program+" minimization")
	# filter based on energy window ewin_rdkit
	sortedcids = ewin_filter(sorted_all_cids,cenergy,args,dup_data,dup_data_idx,log,program)
	# pre-filter based on energy only
	selectedcids_initial = pre_E_filter(sortedcids,cenergy,args,dup_data,dup_data_idx,log,program)
	# filter based on energy and RMSD
	selectedcids = RMSD_and_E_filter(outmols,selectedcids_initial,cenergy,args,dup_data,dup_data_idx,log,program)

	if program == 'xtb':
		dup_data.at[dup_data_idx, 'xTB-Initial-samples'] = len(inmols)
	if program == 'ani':
		dup_data.at[dup_data_idx, 'ANI1ccx-Initial-samples'] = len(inmols)

	# write the filtered, ordered conformers to external file
	write_confs(outmols, cenergy,selectedcids, name, args, program,log)
