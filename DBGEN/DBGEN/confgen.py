#!/usr/bin/python
from __future__ import print_function
import argparse, math, os, sys, traceback,subprocess
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdchem, rdDistGeom, rdMolAlign
from rdkit.Geometry import Point3D
from periodictable import elements as elementspt

import DBGEN.db_gen_functions



### TORCHANI IMPORTS
import ase
import ase.optimize
import torch, torchani
os.environ['KMP_DUPLICATE_LIB_OK']='True'
device = torch.device('cpu')
model = torchani.models.ANI1ccx()
from ase.units import kJ,mol,Hartree,kcal

import xtb
from xtb import GFN2
from xtb.solvation import GBSA

def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_RMSD):
	'''generate RMS distance between two molecules (ignoring hydrogens)'''
	if heavy == True: mol = Chem.RemoveHs(mol)
	rms = Chem.GetBestRMS(mol1,mol2,c1,c2,maxMatches=max_matches_RMSD)
	return rms


def getDihedralMatches(mol, heavy):
	'''return list of atom indices of dihedrals'''
	#this is rdkit's "strict" pattern
	pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
	qmol = Chem.MolFromSmarts(pattern)
	matches = mol.GetSubstructMatches(qmol);

	#these are all sets of 4 atoms, uniquify by middle two
	uniqmatches = []
	seen = set()
	for (a,b,c,d) in matches:
		if (b,c) not in seen:
			if heavy == True:
				if mol.GetAtomWithIdx(a).GetSymbol() != 'H' and mol.GetAtomWithIdx(d).GetSymbol() != 'H':
					seen.add((b,c))
					uniqmatches.append((a,b,c,d))
			if heavy == False:
				seen.add((b,c))
				uniqmatches.append((a,b,c,d))
	return uniqmatches

def genConformer_r(mol, conf, i, matches, degree, sdwriter,args,name):
	'''recursively enumerate all angles for rotatable dihedrals.  i is
	which dihedral we are enumerating by degree to output conformers to out'''
	rotation_count = 1
	if i >= len(matches): #base case, torsions should be set in conf
		sdwriter.write(mol,conf)
		return 1
	else:
		#print(str(i)+'starting new else writing')
		#incr = math.pi*degree / 180.0
		total = 0
		deg = 0
		while deg < 360.0:
			#print(matches[i])
			rad = math.pi*deg / 180.0
			#print(matches[i],rad)
			rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
			#recalculating energies after rotation
			if args.ff == "MMFF":
				GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
			elif args.ff == "UFF":
				GetFF = Chem.UFFGetMoleculeForceField(mol)
			else: print('   Force field {} not supported!'.format(args.ff)); sys.exit()
			GetFF.Initialize()
			converged = GetFF.Minimize(maxIts=args.opt_steps_RDKit)
			energy = GetFF.CalcEnergy()
			mol.SetProp("Energy",energy)
			mol.SetProp('_Name',name+' - conformer from rotation - ' + str(rotation_count))
			# print(energy)
			# print('_Name',name+' - conformer from rotation - ' + str(rotation_count))
			rotation_count +=1
			#print(str(i)+'i')
			total += genConformer_r(mol, conf, i+1, matches, degree, sdwriter,args,name)
			#print(str(total)+'total')
			deg += degree
		return total

def rules_get_charge(mol, args):
	C_group = ['C', 'Se', 'Ge']
	N_group = ['N', 'P', 'As']
	O_group = ['O', 'S', 'Se']
	Cl_group = ['Cl', 'Br', 'I']

	neighbours = []
	#get the neighbours of metal atom
	for atom in mol.GetAtoms():
		if atom.GetSymbol() == args.metal:
			neighbours = atom.GetNeighbors()

	if len(neighbours) == 0 :
		if args.verbose: print("x Metal Not found! It is an organic molecule.")
		#no update in charge as it is an organic molecule
	else:
		args.charge = args.m_oxi
		for atom in neighbours:
			#Carbon list
			if atom.GetSymbol() in C_group:
				if atom.GetTotalValence()== 4:
					args.charge = args.charge - 1
				if atom.GetTotalValence()== 3:
					args.charge = args.charge - 0
			#Nitrogen list
			if atom.GetSymbol() in N_group:
				if atom.GetTotalValence() == 3:
					args.charge = args.charge - 1
				if atom.GetTotalValence() == 4:
					args.charge = args.charge - 0
			#Oxygen list
			if atom.GetSymbol() in O_group:
				if atom.GetTotalValence() == 2:
					args.charge = args.charge - 1
				if atom.GetTotalValence() == 3:
					args.charge = args.charge - 0
			#Halogen list
			if atom.GetSymbol() in Cl_group:
				if atom.GetTotalValence() == 1:
					args.charge = args.charge - 1
				if atom.GetTotalValence() == 2:
					args.charge = args.charge - 0

			print('The neighbour atoms are {0}, with valence {1}, and total charge is {2}'.format(atom.GetSymbol(),atom.GetTotalValence(),args.charge))

	return args.charge, neighbours

def summ_search(mol, name,args, coord_Map = None,alg_Map=None,mol_template=None):
	'''embeds core conformers, then optimizes and filters based on RMSD. Finally the rotatable torsions are systematically rotated'''

	Chem.SanitizeMol(mol)
	mol = Chem.AddHs(mol)
	mol.SetProp("_Name",name)

	if args.nodihedrals == False: rotmatches = getDihedralMatches(mol, args.heavyonly)
	else: rotmatches = []

	if len(rotmatches) > args.max_torsions:
		print("x  Too many torsions (%d). Skipping %s" %(len(rotmatches),(name+args.rdkit_output)))
	else:
		if coord_Map == None and alg_Map == None and mol_template == None:
			if args.etkdg:
				ps = Chem.ETKDG()
				ps.randomSeed = args.seed
				ps.ignoreSmoothingFailures=True
				ps.numThreads = 0
				cids = rdDistGeom.EmbedMultipleConfs(mol, args.sample, params=ps)
			else:
				cids = rdDistGeom.EmbedMultipleConfs(mol, args.sample,ignoreSmoothingFailures=True, rando
				mSeed=args.seed,numThreads = 0)
			if len(cids) == 0 or len(cids) == 1 and args.sample != 1:
				print("o  conformers initially sampled with random coordinates")
				cids = rdDistGeom.EmbedMultipleConfs(mol, args.sample, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0,ignoreSmoothingFailures=True, numZeroFail=1000, numThreads = 0)
			if args.verbose:
				print("o ", len(cids),"conformers initially sampled")
		# case of embed for templates
		else:
			if args.etkdg:
				ps = Chem.ETKDG()
				ps.randomSeed = args.seed
				ps.coordMap = coord_Map
				ps.ignoreSmoothingFailures=True
				ps.numThreads = 0
				cids = rdDistGeom.EmbedMultipleConfs(mol, args.sample, params=ps)
			else:
				cids = rdDistGeom.EmbedMultipleConfs(mol, args.sample, randomSeed=args.seed,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
			if len(cids) == 0 or len(cids) == 1 and args.sample != 1 :
				print("o  conformers initially sampled with random coordinates")
				cids = rdDistGeom.EmbedMultipleConfs(mol, args.sample, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0, numZeroFail=1000,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
			if args.verbose:
				print("o ", len(cids),"conformers initially sampled")

		#energy minimize all to get more realistic results
		#identify the atoms and decide Force Field
		for atom in mol.GetAtoms():
			if atom.GetAtomicNum() > 36: #upto Kr for MMFF, if not use UFF
				args.ff = "UFF"
				#print("UFF is used because there are atoms that MMFF doesn't recognise")
		if args.verbose: print("o  Optimizing", len(cids), "initial conformers with", args.ff)
		if args.verbose:
			if args.nodihedrals == False: print("o  Found", len(rotmatches), "rotatable torsions")
			else: print("o  Systematic torsion rotation is set to OFF")

		cenergy,outmols = [],[]
		for i, conf in enumerate(cids):
			if coord_Map == None and alg_Map == None and mol_template == None:
				if args.ff == "MMFF":
					GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
				elif args.ff == "UFF":
					GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
				else: print('   Force field {} not supported!'.format(args.ff)); sys.exit()

				GetFF.Initialize()
				converged = GetFF.Minimize(maxIts=args.opt_steps_RDKit)
				energy = GetFF.CalcEnergy()
				cenergy.append(GetFF.CalcEnergy())

				#if args.verbose:
				#    print("-   conformer", (i+1), "optimized: ", args.ff, "energy", GetFF.CalcEnergy())
			#id template realign before doing calculations
			else:
				num_atom_match = mol.GetSubstructMatch(mol_template)
				# Force field parameters
				if args.ff == "MMFF":
					GetFF = lambda mol,confId=conf:Chem.MMFFGetMoleculeForceField(mol,Chem.MMFFGetMoleculeProperties(mol),confId=conf)
				elif args.ff == "UFF":
					GetFF = lambda mol,confId=conf:Chem.UFFGetMoleculeForceField(mol,confId=conf)
				else: print('   Force field {} not supported!'.format(options.ff)); sys.exit()
				getForceField=GetFF

				# clean up the conformation
				ff_temp = getForceField(mol, confId=conf)
				for k, idxI in enumerate(num_atom_match):
					for l in range(k + 1, len(num_atom_match)):
						idxJ = num_atom_match[l]
						d = coord_Map[idxI].Distance(coord_Map[idxJ])
						ff_temp.AddDistanceConstraint(idxI, idxJ, d, d, 10000)
				ff_temp.Initialize()
				#reassignned n from 4 to 10 for better embed and minimzation
				n = 10
				more = ff_temp.Minimize()
				while more and n:
					more = ff_temp.Minimize()
					n -= 1
				energy = ff_temp.CalcEnergy()
				# rotate the embedded conformation onto the core_mol:
				rms = rdMolAlign.AlignMol(mol, mol_template,prbCid=conf, atomMap=alg_Map,reflect=True,maxIters=100)
				# elif len(num_atom_match) == 5:
				# 	ff_temp = GetFF(mol, confId=conf)
				# 	conf_temp = mol_template.GetConformer()
				# 	for k in range(mol_template.GetNumAtoms()):
				# 		p = conf_temp.GetAtomPosition(k)
				# 		q = mol.GetConformer(conf).GetAtomPosition(k)
				# 		pIdx = ff_temp.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
				# 		ff_temp.AddDistanceConstraint(pIdx, num_atom_match[k], 0, 0, 10000)
				# 	ff_temp.Initialize()
				# 	n = 10
				# 	more = ff_temp.Minimize(energyTol=1e-6, forceTol=1e-5)
				# 	while more and n:
				# 		more = ff_temp.Minimize(energyTol=1e-6, forceTol=1e-5)
				# 		n -= 1
				# 	# realign
				# 	energy = ff_temp.CalcEnergy()
				# 	rms = rdMolAlign.AlignMol(mol, mol_template,prbCid=conf, atomMap=alg_Map,reflect=True,maxIters=50)
				cenergy.append(energy)

			# outmols is gonna be a list containing "args.sample" mol objects with "args.sample"
			# conformers. We do this to SetProp (Name and Energy) to the different conformers
			# and print in the SDF file. At the end, since all the mol objects has the same
			# conformers, but the energies are different, we can print conformers to SDF files
			# with the energies of the parent mol objects. We measured the computing time and
			# it's the same as using only 1 parent mol object with 10 conformers, but we couldn'temp
			# SetProp correctly
			pmol = PropertyMol.PropertyMol(mol)
			outmols.append(pmol)

		for i, cid in enumerate(cids):
			outmols[cid].SetProp('_Name', name + ' conformer ' + str(i+1))
			outmols[cid].SetProp('Energy', cenergy[cid])

		#reduce to unique set
		if args.verbose: print("o  Removing duplicate conformers ( RMSD <", args.rms_threshold, "and E difference <",args.energy_threshold,"kcal/mol)")
		cids = list(range(len(outmols)))
		sortedcids = sorted(cids,key = lambda cid: cenergy[cid])

		selectedcids_initial, selectedcids = [],[]
		# First RDKit filter, only based on E
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
					excluded_conf = True
					break
			if excluded_conf == False:
				if conf not in selectedcids_initial:
					selectedcids_initial.append(conf)
		# As the RDKit filter, the get_conf_RMS is
		for i, conf in enumerate(selectedcids_initial):
			# This keeps track of whether or not your conformer is unique
			excluded_conf = False
			# include the first conformer in the list to start the filtering process
			if i == 0:
				selectedcids.append(conf)
			# check rmsd
			for seenconf in selectedcids:
				rms = get_conf_RMS(outmols[conf],outmols[conf],seenconf,conf, args.heavyonly, args.max_matches_RMSD)
				E_diff = abs(cenergy[conf] - cenergy[seenconf]) # in kcal/mol
				if rms < args.rms_threshold and E_diff < (args.energy_threshold ):
					excluded_conf = True
					break
			if excluded_conf == False:
				if conf not in selectedcids:
					selectedcids.append(conf)

		# unique_mols, unique_energies = [],[]
		# for id in selectedcids:
		# 	unique_mols.append(outmols[id])
		# 	unique_energies.append(cenergy[id])

		# print(unique_mols[0:2].GetConformers()[0].GetPositions())
		if args.verbose: print("o ", len(selectedcids),"unique conformers remain")

		#if use the rotations constraint, then this block would write the corresponding SDF file.
		if len(rotmatches) != 0:
			#print(rotmatches)
			sdwriter_rot = Chem.SDWriter(name+'_rotation_dihydral.sdf')
			# now exhaustively drive torsions of selected conformers
			n_confs = int(len(selectedcids) * (360 / args.degree) ** len(rotmatches))
			if args.verbose and len(rotmatches) != 0: print("o  Systematic generation of", n_confs, "confomers")

			total = 0
			for conf in selectedcids:
				#print(outmols[conf])
				total += genConformer_r(outmols[conf], conf, 0, rotmatches, args.degree, sdwriter_rot,args,outmols[conf].GetProp('_Name'))
			if args.verbose and len(rotmatches) != 0: print("o  %d total conformations generated"%total)
			sdwriter_rot.close()

			#doing a filter for the written molecule after rotations
			rdmols = Chem.SDMolSupplier(name+'_rotation_dihydral.sdf', removeHs=False)

			if rdmols is None:
				print("Could not open after rotations ", name+'_temp.sdf')
				sys.exit(-1)

			for mol in rdmols:
				#setting the metal back instead of I
				if args.metal_complex == True:
					for atom in mol.GetAtoms():
						if atom.GetSymbol() == 'I' and (len(atom.GetBonds()) == 6 or len(atom.GetBonds()) == 5 or len(atom.GetBonds()) == 4 or len(atom.GetBonds()) == 3 or len(atom.GetBonds()) == 2):
							for el in elementspt:
								if el.symbol == args.metal:
									atomic_number = el.number
							atom.SetAtomicNum(atomic_number)


			sdwriter = Chem.SDWriter(name+args.rdkit_output)
			rd_count = 0
			rd_selectedcids =[]
			for i in range(len(rdmols)):
				# This keeps track of whether or not your conformer is unique
				excluded_conf = False
				# include the first conformer in the list to start the filtering process
				if rd_count == 0:
					rd_selectedcids.append(i)
					sdwriter.write(rdmols[i])
				# Only the first ID gets included
				rd_count = 1
				# check rmsd
				for j in rd_selectedcids:
					rms = get_conf_RMS(rdmols[i],rdmols[j],-1,-1, args.heavyonly, args.max_matches_RMSD)
					E_diff = abs(float(rdmols[i].GetProp('Energy')) - float(rdmols[j].GetProp('Energy'))) # in kcal/mol
					if rms < args.rms_threshold and E_diff < (args.energy_threshold ):
						excluded_conf = True
						break
				if excluded_conf == False:
					sdwriter.write(rdmols[i])
					if i not in rd_selectedcids:
						rd_selectedcids.append(i)
			sdwriter.close()

		else:
			sdwriter = Chem.SDWriter(name+args.rdkit_output)
			for id in selectedcids:
				#setting the metal back instead of I
				if args.metal_complex == True:
					for atom in outmols[id].GetAtoms():
						if atom.GetSymbol() == 'I' and (len(atom.GetBonds()) == 6 or len(atom.GetBonds()) == 5 or len(atom.GetBonds()) == 4 or len(atom.GetBonds()) == 3 or len(atom.GetBonds()) == 2):
							for el in elementspt:
								if el.symbol == args.metal:
									atomic_number = el.number
							atom.SetAtomicNum(atomic_number)
				sdwriter.write(outmols[id],id)
			sdwriter.close()

		# Calculates geometries and energies with xTB or ANI1
		if args.ANI1ccx == True or args.xtb == True:

			ani_energy,xtb_energy = 0,0

			# lists containing energies, coordinates and mols generated by xTB and ANI1
			SQM_energy, SQM_cartesians, SQM_mols = [], [], []

			globmin = None

			# if large system increase stck size
			if args.large_sys == True:
				os.environ['OMP_STACKSIZE'] = args.STACKSIZE
				if args.verbose == True: print('---- The Stack size has been updated for larger molecules to {0} ----'.format(os.environ['OMP_STACKSIZE']))

			inmols = Chem.SDMolSupplier(name+args.rdkit_output, removeHs=False)

			if inmols is None:
				print("Could not open ", name+output)
				sys.exit(-1)

			for mol in inmols:
				#setting the metal back instead of I
				if args.metal_complex == True:
					for atom in mol.GetAtoms():
						if atom.GetSymbol() == 'I' and (len(atom.GetBonds()) == 6 or len(atom.GetBonds()) == 5 or len(atom.GetBonds()) == 4 or len(atom.GetBonds()) == 3 or len(atom.GetBonds()) == 2):
							for el in elementspt:
								if el.symbol == args.metal:
									atomic_number = el.number
							atom.SetAtomicNum(atomic_number)

				# removing the Ba atom if NCI complexes
				if args.nci_complex == True:
					for atom in mol.GetAtoms():
						if atom.GetSymbol() =='I':
							atom.SetAtomicNum(1)

				elements = ''
				for atom in mol.GetAtoms():
					elements += atom.GetSymbol()

				cartesians = mol.GetConformers()[0].GetPositions()
				# print(cartesians[0:2])
				if args.verbose == True: print('---- The elements are the following {0} ----'.format(elements))

				coordinates = torch.tensor([cartesians.tolist()], requires_grad=True, device=device)

				if args.ANI1ccx == True:
					final_output = args.ani1_output
					species = model.species_to_tensor(elements).to(device).unsqueeze(0)
					_, ani_energy = model((species, coordinates))
					if args.verbose: print("ANI Initial E:",ani_energy.item()*627.509,'kcal/mol')

					ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0], calculator=model.ase())
					### make a function for constraints and optimization
					if constraints != None:
						fb = ase.constraints.FixBondLength(0, 1)
						ase_molecule.set_distance(0,1,2.0)
						ase_molecule.set_constraint(fb)

					optimizer = ase.optimize.BFGS(ase_molecule, trajectory='ANI1_opt.traj')
					optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
					if len(ase.io.Trajectory('xTB_opt.traj', mode='r')) != (args.opt_steps+1):
						species_coords = ase_molecule.get_positions().tolist()
						coordinates = torch.tensor([species_coords], requires_grad=True, device=device)
						SQM_energy.append(ani_energy.item())
						cartesians = np.array(coordinates.tolist()[0])
						SQM_cartesians.append(cartesians)

						###############################################################################
						# Now let's compute energy:
						_, ani_energy = model((species, coordinates))
						aniE = ani_energy.item() #Hartree
						if args.verbose: print("ANI Final E:", aniE*627.509,'kcal/mol')
						###############################################################################
					else:
						print('ANI1 optimization could not converge (discarding conformer).')
### INCLUDE THE OPTIONS TO STORE MOLECULAR Descriptors
### CHECK THIS WEBPAGE: https://github.com/grimme-lab/xtb/tree/master/python
				elif args.xtb == True:
					final_output = args.xtb_output
					if args.metal_complex == True:
						#passing charges metal present
						if args.solvent_xtb != '':
							ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0], calculator=GFN2()) #define ase molecule using GFN2 Calculator
							#ase_molecule.set_calculator(GBSA(solvent=args.solvent_xtb, calculator=GFN2()))
						else:
							ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=GFN2()) #define ase molecule using GFN2 Calculator
						for atom in ase_molecule:
							if atom.symbol == args.metal:
								#will update only for cdx, smi, and csv formats.
								if os.path.splitext(args.input)[1] == '.csv' or os.path.splitext(args.input)[1] == '.cdx' or os.path.splitext(args.input)[1] == '.smi':
									atom.charge, neighbours = rules_get_charge(mol,args)
									if len(neighbours) != 0:
										if args.verbose == True: print('---- The Overall charge is reworked with rules for .smi, .csv, .cdx ----')
								else:
									atom.charge = args.charge
									if args.verbose == True: print('---- The Overall charge is read from the .com file ----')
								if args.verbose == True: print('---- The Overall charge considered is  {0} ----'.format(atom.charge))
					else:
						if args.solvent_xtb != '':
							print('using solvent')
							ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=GFN2()) #define ase molecule using GFN2 Calculator
							#ase_molecule.set_calculator(GBSA(solvent=args.solvent_xtb, calculator=GFN2()))
						else:
							ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=GFN2()) #define ase molecule using GFN2 Calculator

					if args.verbose: print("Initial XTB energy", (ase_molecule.get_potential_energy()/Hartree)*627.509,'kcal/mol')
					optimizer = ase.optimize.BFGS(ase_molecule, trajectory='xTB_opt.traj')
					optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
					if len(ase.io.Trajectory('xTB_opt.traj', mode='r')) != (args.opt_steps+1):
						species_coords = ase_molecule.get_positions().tolist()
						coordinates = torch.tensor([species_coords], requires_grad=True, device=device)

						###############################################################################
						# Now let's compute energy:
						xtb_energy = ase_molecule.get_potential_energy()
						SQM_energy.append((xtb_energy/Hartree)*627.509) # gets the energy in kcal/mol
						cartesians = np.array(coordinates.tolist()[0])
						SQM_cartesians.append(cartesians)
						if args.verbose: print("Final XTB E:",(xtb_energy/Hartree)*627.509,'kcal/mol')
						###############################################################################
					else:
						print('xTB optimization could not converge (discarding conformer).')

				pmol = PropertyMol.PropertyMol(mol)
				SQM_mols.append(pmol)

			if len(SQM_energy) == 0:
				if args.xtb == True:
					print('\n WARNING! No conformers converged during xTB optimization.')

				elif args.ANI1ccx == True:
					print('\n WARNING! No conformers converged during ANI1 optimization.')


			else:
				if args.verbose: print("o  Removing duplicate conformers ( RMSD <", args.rms_threshold, "and E difference <",args.energy_threshold,"kcal/mol)")
				cids = list(range(len(SQM_mols)))
				SQM_sortedcids = sorted(cids,key = lambda cid: SQM_energy[cid])
				SQM_selectedcids = []

				SQM_count = 0
				for i in SQM_sortedcids:
					# This keeps track of whether or not your conformer is unique
					excluded_conf = False
					# include the first conformer in the list to start the filtering process
					if SQM_count == 0:
						SQM_selectedcids.append(i)
					# Only the first ID gets included
					SQM_count = 1
					# check rmsd
					for j in SQM_selectedcids:
						rms = get_conf_RMS(SQM_mols[i],SQM_mols[j],-1,-1, args.heavyonly, args.max_matches_RMSD)
						E_diff = abs(SQM_energy[i] - SQM_energy[j]) # in kcal/mol
						if rms < args.rms_threshold and E_diff < (args.energy_threshold ):
							excluded_conf = True
							break
					if excluded_conf == False:
						if i not in SQM_selectedcids:
							SQM_selectedcids.append(i)

				post_outmols, post_energy = [],[]
				for id in SQM_selectedcids:
					# change mol objects with optimized x,y,z from xTB
					c = SQM_mols[id].GetConformer()
					for j in range(SQM_mols[id].GetNumAtoms()):
						[x,y,z] = SQM_cartesians[id][j]
						c.SetAtomPosition(j,Point3D(x,y,z))
					post_outmols.append(SQM_mols[id])
					post_energy.append(SQM_energy[id])

				if args.verbose: print("o ", len(SQM_selectedcids),"unique conformers remain")

				for i, cid in enumerate(SQM_selectedcids):
					SQM_mols[cid].SetProp('_Name', name + ' conformer ' + str(i+1))
					SQM_mols[cid].SetProp('Energy', SQM_energy[cid])

				sdwriter = Chem.SDWriter(name+final_output)
				for mol in post_outmols:
					sdwriter.write(mol)
				sdwriter.close()

				if len(post_outmols) == 0: print("\nx  WARNING! No conformers found!\n")
