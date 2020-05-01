#!/usr/bin/python
from __future__ import print_function
import argparse, math, os, sys, traceback, subprocess, glob, shutil, time
import os.path as path
import numpy as np
import pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdchem, rdDistGeom, rdMolAlign, PropertyMol, rdChemReactions, Lipinski, Descriptors
from rdkit.Geometry import Point3D
from periodictable import elements as elementspt
from openbabel import openbabel as ob

import progress
from progress.bar import IncrementalBar
import warnings
import yaml

from DBGEN.db_gen_functions import *

### TORCHANI IMPORTS
try:
	import ase
	import ase.optimize
	import torch, torchani
	os.environ['KMP_DUPLICATE_LIB_OK']='True'
	device = torch.device('cpu')
	model = torchani.models.ANI1ccx()
	from ase.units import Hartree,kcal,kJ
	from ase.units import mol as mol_unit
except:
	Hartree, kcal, kJ, mol_unit = 27.211386024367243, 2.611447418269555e+22, 6.241509125883258e+21, 6.022140857e+23

try:
	import xtb
	from xtb import GFN2
except: pass


def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_RMSD,log):
	'''generate RMS distance between two molecules (ignoring hydrogens)'''
	if heavy == True:
		 mol1 = Chem.RemoveHs(mol1)
		 mol2 = Chem.RemoveHs(mol2)
	rms = Chem.GetBestRMS(mol1,mol2,c1,c2,maxMatches=max_matches_RMSD)
	return rms

def get_PMIDIFF(mol1, mol2, c1, c2, heavy):
	if heavy == True:
	   mol1 = Chem.RemoveHs(mol1)
	   mol2 = Chem.RemoveHs(mol2)
	pmi1 = Chem.CalcPMI1(mol1,confId=c1) - Chem.CalcPMI1(mol2,confId=c2)
	pmi2 = Chem.CalcPMI2(mol1,confId=c1) - Chem.CalcPMI2(mol2,confId=c2)
	pmi3 = Chem.CalcPMI3(mol1,confId=c1) - Chem.CalcPMI3(mol2,confId=c2)
	diff = (pmi1 * pmi1 + pmi2 * pmi2 + pmi3 * pmi3) ** 0.5
	return diff

def get_TFD(mol1, mol2, c1, c2, heavy):
	if heavy == True:
	   mol1 = Chem.RemoveHs(mol1)
	   mol2 = Chem.RemoveHs(mol2)
	tors_list, ring_tors_list = TorsionFingerprints.CalculateTorsionLists(mol1)
	weights = TorsionFingerprints.CalculateTorsionWeights(mol1)
	torsions1 = TorsionFingerprints.CalculateTorsionAngles(mol1, tors_list, ring_tors_list, confId=c1)
	torsions2 = TorsionFingerprints.CalculateTorsionAngles(mol2, tors_list, ring_tors_list, confId=c2)
	tfd = TorsionFingerprints.CalculateTFD(torsions1, torsions2, weights=weights)
	return tfd

"AUTO-SAMPLING DETECTS INITIAL NUMBER OF SAMPLES"
def auto_sampling(mult_factor,mol,log):
	auto_samples = 0
	auto_samples += 3*(Lipinski.NumRotatableBonds(mol)) # x3, for C3 rotations
	auto_samples += 3*(Lipinski.NHOHCount(mol)) # x3, for OH/NH rotations
	auto_samples += 3*(Lipinski.NumSaturatedRings(mol)) # x3, for boat/chair/envelope confs
	auto_samples = mult_factor*auto_samples
	return auto_samples

def getDihedralMatches(mol, heavy,log):
	'''return list of atom indices of dihedrals'''
	#this is rdkit's "strict" pattern
	pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
	qmol = Chem.MolFromSmarts(pattern)
	matches = mol.GetSubstructMatches(qmol);

	#these are all sets of 4 atoms, uniquify by middle two
	uniqmatches = []
	seen = set()
	for (a,b,c,d) in matches:
		if (b,c) not in seen and (c,b) not in seen:
			if heavy == True:
				if mol.GetAtomWithIdx(a).GetSymbol() != 'H' and mol.GetAtomWithIdx(d).GetSymbol() != 'H':
					seen.add((b,c))
					uniqmatches.append((a,b,c,d))
			if heavy != True:
				if mol.GetAtomWithIdx(c).GetSymbol() == 'C' and mol.GetAtomWithIdx(d).GetSymbol() == 'H':
					pass
				else:
					seen.add((b,c))
					uniqmatches.append((a,b,c,d))
	return uniqmatches

def genConformer_r(mol, conf, i, matches, degree, sdwriter,args,name,log):
	'''recursively enumerate all angles for rotatable dihedrals.  i is
	which dihedral we are enumerating by degree to output conformers to out'''
	rotation_count = 1
	if i >= len(matches): #base case, torsions should be set in conf
		#setting the metal back instead of I
		if args.metal_complex == True and args.nodihedrals == True:
			for atom in mol.GetAtoms():
				if atom.GetSymbol() == 'I' and (len(atom.GetBonds()) == 6 or len(atom.GetBonds()) == 5 or len(atom.GetBonds()) == 4 or len(atom.GetBonds()) == 3 or len(atom.GetBonds()) == 2):
					for el in elementspt:
						if el.symbol == args.metal:
							atomic_number = el.number
					atom.SetAtomicNum(atomic_number)
		sdwriter.write(mol,conf)
		return 1
	else:
		#log.write(str(i)+'starting new else writing')
		#incr = math.pi*degree / 180.0
		total = 0
		deg = 0
		while deg < 360.0:
			#log.write(matches[i])
			rad = math.pi*deg / 180.0
			#log.write(matches[i],rad)
			rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
			#recalculating energies after rotation
			if args.ff == "MMFF":
				GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
			elif args.ff == "UFF":
				GetFF = Chem.UFFGetMoleculeForceField(mol)
			else: log.write('   Force field {} not supported!'.format(args.ff)); sys.exit()
			GetFF.Initialize()
			converged = GetFF.Minimize(maxIts=args.opt_steps_RDKit)
			energy = GetFF.CalcEnergy()
			mol.SetProp("Energy",energy)
			mol.SetProp('_Name',name+' - conformer from rotation - ' + str(rotation_count))
			rotation_count +=1
			total += genConformer_r(mol, conf, i+1, matches, degree, sdwriter,args,name,log)

			deg += degree
		return total

def rules_get_charge(mol, args,log):
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
		if args.verbose: log.write("x Metal Not found! It is an organic molecule.")
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

			# log.write('The neighbour atoms are {0}, with valence {1}, and total charge is {2}'.format(atom.GetSymbol(),atom.GetTotalValence(),args.charge))

	return args.charge, neighbours

def summ_search(mol, name,args,log,dup_data,dup_data_idx, coord_Map = None,alg_Map=None,mol_template=None):
	'''embeds core conformers, then optimizes and filters based on RMSD. Finally the rotatable torsions are systematically rotated'''

	sdwriter = Chem.SDWriter(name+'_'+'rdkit'+args.output)

	Chem.SanitizeMol(mol)
	mol = Chem.AddHs(mol)
	mol.SetProp("_Name",name)

	# detects and applies auto-detection of initial number of conformers
	if args.sample == 'auto':
		initial_confs = int(auto_sampling(args.auto_sample,mol,log))

	else:
		initial_confs = int(args.sample)

	#
	dup_data.at[dup_data_idx, 'Molecule'] = name
	dup_data.at[dup_data_idx, 'RDKIT-Initial-samples'] = initial_confs

	if args.nodihedrals == False: rotmatches = getDihedralMatches(mol, args.heavyonly,log)
	else: rotmatches = []

	if len(rotmatches) > args.max_torsions:
		log.write("x  Too many torsions (%d). Skipping %s" %(len(rotmatches),(name+args.output)))
		status = -1
	else:
		if coord_Map == None and alg_Map == None and mol_template == None:
			if args.etkdg:
				ps = Chem.ETKDG()
				ps.randomSeed = args.seed
				ps.ignoreSmoothingFailures=True
				ps.numThreads = 0
				cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, params=ps)
			else:
				cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs,ignoreSmoothingFailures=True, randomSeed=args.seed,numThreads = 0)
			if len(cids) == 0 or len(cids) == 1 and initial_confs != 1:
				log.write("o  conformers initially sampled with random coordinates")
				cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0,ignoreSmoothingFailures=True, numZeroFail=1000, numThreads = 0)
			if args.verbose:
				log.write("o  "+ str(len(cids))+" conformers initially sampled")
		# case of embed for templates
		else:
			if args.etkdg:
				ps = Chem.ETKDG()
				ps.randomSeed = args.seed
				ps.coordMap = coord_Map
				ps.ignoreSmoothingFailures=True
				ps.numThreads = 0
				cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, params=ps)
			else:
				cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
			if len(cids) == 0 or len(cids) == 1 and initial_confs != 1 :
				log.write("o  conformers initially sampled with random coordinates")
				cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0, numZeroFail=1000,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
			if args.verbose:
				log.write("o  "+ str(len(cids))+" conformers initially sampled")

		#energy minimize all to get more realistic results
		#identify the atoms and decide Force Field
		for atom in mol.GetAtoms():
			if atom.GetAtomicNum() > 36: #upto Kr for MMFF, if not use UFF
				args.ff = "UFF"
				#log.write("UFF is used because there are atoms that MMFF doesn't recognise")
		if args.verbose: log.write("o  Optimizing "+ str(len(cids))+ " initial conformers with"+ args.ff)
		if args.verbose:
			if args.nodihedrals == False:
				log.write("o  Found "+ str(len(rotmatches))+ " rotatable torsions")
				# for [a,b,c,d] in rotmatches:
				# 	log.write('  '+mol.GetAtomWithIdx(a).GetSymbol()+str(a+1)+ mol.GetAtomWithIdx(b).GetSymbol()+str(b+1)+ mol.GetAtomWithIdx(c).GetSymbol()+str(c+1)+mol.GetAtomWithIdx(d).GetSymbol()+str(d+1))
			else: log.write("o  Systematic torsion rotation is set to OFF")

		cenergy,outmols = [],[]
		bar = IncrementalBar('o  Minimizing', max = len(cids))
		for i, conf in enumerate(cids):
			if coord_Map == None and alg_Map == None and mol_template == None:
				if args.ff == "MMFF":
					GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
				elif args.ff == "UFF":
					GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
				else: log.write('   Force field {} not supported!'.format(args.ff)); sys.exit()

				GetFF.Initialize()
				converged = GetFF.Minimize(maxIts=args.opt_steps_RDKit)
				energy = GetFF.CalcEnergy()
				cenergy.append(GetFF.CalcEnergy()* 4.184)

				#if args.verbose:
				#    log.write("-   conformer", (i+1), "optimized: ", args.ff, "energy", GetFF.CalcEnergy())
			#id template realign before doing calculations
			else:
				num_atom_match = mol.GetSubstructMatch(mol_template)
				# Force field parameters
				if args.ff == "MMFF":
					GetFF = lambda mol,confId=conf:Chem.MMFFGetMoleculeForceField(mol,Chem.MMFFGetMoleculeProperties(mol),confId=conf)
				elif args.ff == "UFF":
					GetFF = lambda mol,confId=conf:Chem.UFFGetMoleculeForceField(mol,confId=conf)
				else: log.write('   Force field {} not supported!'.format(options.ff)); sys.exit()
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
				#     ff_temp = GetFF(mol, confId=conf)
				#     conf_temp = mol_template.GetConformer()
				#     for k in range(mol_template.GetNumAtoms()):
				#         p = conf_temp.GetAtomPosition(k)
				#         q = mol.GetConformer(conf).GetAtomPosition(k)
				#         pIdx = ff_temp.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
				#         ff_temp.AddDistanceConstraint(pIdx, num_atom_match[k], 0, 0, 10000)
				#     ff_temp.Initialize()
				#     n = 10
				#     more = ff_temp.Minimize(energyTol=1e-6, forceTol=1e-5)
				#     while more and n:
				#         more = ff_temp.Minimize(energyTol=1e-6, forceTol=1e-5)
				#         n -= 1
				#     # realign
				#     energy = ff_temp.CalcEnergy()
				#     rms = rdMolAlign.AlignMol(mol, mol_template,prbCid=conf, atomMap=alg_Map,reflect=True,maxIters=50)
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
		sortedcids = sorted(cids,key = lambda cid: cenergy[cid])

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
			if excluded_conf == False:
				if conf not in selectedcids_initial:
					selectedcids_initial.append(conf)
			bar.next()
		bar.finish()


		if args.verbose == True: log.write("o  "+str(eng_dup)+ " Duplicates removed  pre-energy filter (E < "+str(args.initial_energy_threshold)+" kcal/mol )")


		#reduce to unique set
		if args.verbose: log.write("o  Removing duplicate conformers ( RMSD < "+ str(args.rms_threshold)+ " and E difference < "+str(args.energy_threshold)+" kcal/mol)")

		bar = IncrementalBar('o  Filtering based on energy and rms', max = len(selectedcids_initial))
		#check rmsd
		for i, conf in enumerate(selectedcids_initial):
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
			if excluded_conf == False:
				if conf not in selectedcids:
					selectedcids.append(conf)
			bar.next()
		bar.finish()


		# unique_mols, unique_energies = [],[]
		# for id in selectedcids:
		#     unique_mols.append(outmols[id])
		#     unique_energies.append(cenergy[id])

		# log.write(unique_mols[0:2].GetConformers()[0].GetPositions())

		if args.verbose == True: log.write("o  "+str(eng_rms_dup)+ " Duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol) after rotation")
		if args.verbose: log.write("o  "+ str(len(selectedcids))+" unique (ignoring torsions) starting conformers remain")

		dup_data.at[dup_data_idx, 'RDKit-energy-duplicates'] = eng_dup
		dup_data.at[dup_data_idx, 'RDKit-RMS-and-energy-duplicates'] = eng_rms_dup
		dup_data.at[dup_data_idx, 'RDKIT-Unique-conformers'] = len(selectedcids)

		# now exhaustively drive torsions of selected conformers
		n_confs = int(len(selectedcids) * (360 / args.degree) ** len(rotmatches))
		if args.verbose and len(rotmatches) != 0:
			log.write("\n\no  Systematic generation of "+ str(n_confs)+ " confomers")
			bar = IncrementalBar('o  Generating conformations based on dihedral rotation', max = len(selectedcids))
		else:
			bar = IncrementalBar('o  Generating conformations', max = len(selectedcids))

		total = 0
		for conf in selectedcids:
			#log.write(outmols[conf])
			total += genConformer_r(outmols[conf], conf, 0, rotmatches, args.degree, sdwriter ,args,outmols[conf].GetProp('_Name'),log)
			bar.next()
		bar.finish()
		if args.verbose and len(rotmatches) != 0: log.write("o  %d total conformations generated"%total)
		status = 1
	sdwriter.close()

	#getting the energy from and mols after rotations
	if len(rotmatches) != 0:
		rdmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)
		if rdmols is None:
			log.write("Could not open "+ name+args.output)
			sys.exit(-1)


		bar = IncrementalBar('o  Filtering based on energy and rms after rotation of dihedrals', max = len(rdmols))
		sdwriter = Chem.SDWriter(name+'_'+'rdkit'+'_'+'rotated'+args.output)

		rd_count = 0
		rd_selectedcids,rd_dup_energy,rd_dup_rms_eng =[],-1,0
		for i in range(len(rdmols)):
			# This keeps track of whether or not your conformer is unique
			excluded_conf = False
			# include the first conformer in the list to start the filtering process
			if rd_count == 0:
				rd_selectedcids.append(i)
				if args.metal_complex == True:
					for atom in rdmols[i].GetAtoms():
						if atom.GetSymbol() == 'I' and (len(atom.GetBonds()) == 6 or len(atom.GetBonds()) == 5 or len(atom.GetBonds()) == 4 or len(atom.GetBonds()) == 3 or len(atom.GetBonds()) == 2):
							for el in elementspt:
								if el.symbol == args.metal:
									atomic_number = el.number
							atom.SetAtomicNum(atomic_number)
				sdwriter.write(rdmols[i])
			# Only the first ID gets included
			rd_count = 1
			# check rmsd
			for j in rd_selectedcids:
				if abs(float(rdmols[i].GetProp('Energy')) - float(rdmols[j].GetProp('Energy'))) < args.initial_energy_threshold: # comparison in kcal/mol
					excluded_conf = True
					rd_dup_energy += 1
					break
				if abs(float(rdmols[i].GetProp('Energy')) - float(rdmols[j].GetProp('Energy'))) < args.energy_threshold: # in kcal/mol
					rms = get_conf_RMS(rdmols[i],rdmols[j],-1,-1, args.heavyonly, args.max_matches_RMSD,log)
					if rms < args.rms_threshold:
						excluded_conf = True
						rd_dup_rms_eng += 1
						break
			if excluded_conf == False:
				if args.metal_complex == True:
					for atom in rdmols[i].GetAtoms():
						if atom.GetSymbol() == 'I' and (len(atom.GetBonds()) == 6 or len(atom.GetBonds()) == 5 or len(atom.GetBonds()) == 4 or len(atom.GetBonds()) == 3 or len(atom.GetBonds()) == 2):
							for el in elementspt:
								if el.symbol == args.metal:
									atomic_number = el.number
							atom.SetAtomicNum(atomic_number)
				sdwriter.write(rdmols[i])
				if i not in rd_selectedcids:
					rd_selectedcids.append(i)
			bar.next()
		bar.finish()
		sdwriter.close()

		if args.verbose == True: log.write("o  "+str(rd_dup_energy)+ " Duplicates removed initial energy ( E < "+str(args.initial_energy_threshold)+" kcal/mol )")
		if args.verbose == True: log.write("o  "+str(rd_dup_rms_eng)+ " Duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol) after rotation")
		if args.verbose == True: log.write("o  "+str(len(rd_selectedcids) )+ " unique (after torsions) conformers remain")


		#filtering process after rotations
		dup_data.at[dup_data_idx, 'RDKIT-Rotated-conformers'] = total
		dup_data.at[dup_data_idx, 'RDKIT-Rotated-Unique-conformers'] = len(rd_selectedcids)

	return status

def optimize(mol, args, program,log):

	#setup non rdkit energy calculations

	# if large system increase stck size
	if args.large_sys == True:
		os.environ['OMP_STACKSIZE'] = args.STACKSIZE
		# if args.verbose == True: log.write('---- The Stack size has been updated for larger molecules to {0} ----'.format(os.environ['OMP_STACKSIZE']))

	# removing the Ba atom if NCI complexes
	if args.nci_complex == True:
		for atom in mol.GetAtoms():
			if atom.GetSymbol() =='I':
				atom.SetAtomicNum(1)

	if args.metal_complex == True and args.nodihedrals == False:
		for atom in mol.GetAtoms():
			if atom.GetSymbol() == 'I' and (len(atom.GetBonds()) == 6 or len(atom.GetBonds()) == 5 or len(atom.GetBonds()) == 4 or len(atom.GetBonds()) == 3 or len(atom.GetBonds()) == 2):
				for el in elementspt:
					if el.symbol == args.metal:
						atomic_number = el.number
				atom.SetAtomicNum(atomic_number)

	elements = ''
	for atom in mol.GetAtoms():
		elements += atom.GetSymbol()

	cartesians = mol.GetConformers()[0].GetPositions()
	# log.write(cartesians[0:2])
	if args.verbose == True: log.write('\no  The elements are the following {0} '.format(elements))

	coordinates = torch.tensor([cartesians.tolist()], requires_grad=True, device=device)

	if program == 'ani':
		species = model.species_to_tensor(elements).to(device).unsqueeze(0)
		_, ani_energy = model((species, coordinates))

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
			converged = 0
		###############################################################################
		# Now let's compute energy:
		_, ani_energy = model((species, coordinates))
		sqm_energy = ani_energy.item() * 627.509 # Hartree to kcal/mol
		#if args.verbose: log.write("o  ANI Final E:", ani_energy.item(),'eH', ase_molecule.get_potential_energy(),'eV') #Hartree, eV
		###############################################################################

	elif program == 'xtb':

		if args.metal_complex == True:
			#passing charges metal present
			ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=GFN2()) #define ase molecule using GFN2 Calculator
			for atom in ase_molecule:
				if atom.symbol == args.metal:
					#will update only for cdx, smi, and csv formats.
					if os.path.splitext(args.input)[1] == '.csv' or os.path.splitext(args.input)[1] == '.cdx' or os.path.splitext(args.input)[1] == '.smi':
						atom.charge, neighbours = rules_get_charge(mol,args,log)
						if len(neighbours) != 0:
							if args.verbose == True: log.write('o  The Overall charge is reworked with rules for .smi, .csv, .cdx ')
					else:
						atom.charge = args.charge
						if args.verbose == True: log.write('o  The Overall charge is read from the .com file ')
					if args.verbose == True: log.write('o  The Overall charge considered is  {0} '.format(atom.charge))
		else:
			ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0],calculator=GFN2()) #define ase molecule using GFN2 Calculator
		optimizer = ase.optimize.BFGS(ase_molecule, trajectory='xTB_opt.traj',logfile='xtb.opt')
		optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
		if len(ase.io.Trajectory('xTB_opt.traj', mode='r')) != (args.opt_steps+1):
			species_coords = ase_molecule.get_positions().tolist()
			coordinates = torch.tensor([species_coords], requires_grad=True, device=device)
			converged = 0
		###############################################################################
		# Now let's compute energy:
		xtb_energy = ase_molecule.get_potential_energy()
		sqm_energy = (xtb_energy / Hartree)* 627.509
		#if args.verbose: log.write("o  Final XTB E:",xtb_energy/Hartree,'Eh',xtb_energy,'eV') #Hartree, eV
		###############################################################################

	else:
		log.write('program not defined!')

	energy, converged, cartesians = sqm_energy, converged, np.array(coordinates.tolist()[0])
	# update coordinates of mol object
	for j in range(mol.GetNumAtoms()):
		[x,y,z] = cartesians[j]
		mol.GetConformer().SetAtomPosition(j,Point3D(x,y,z))

	return mol, converged, energy

def write_confs(conformers, energies, name, args, program,log):
	n_filter = 0
	if len(conformers) > 0:
		# list in energy order
		cids = list(range(len(conformers)))
		sortedcids = sorted(cids, key = lambda cid: energies[cid])

		#if args.verbose: log.write("o ", name, "Energies:", np.array(sorted(energies)))
		name = name.split('_rdkit')[0]# a bit hacky
		sdwriter = Chem.SDWriter(name+'_'+program+args.output)
		glob_min = min(energies)

		write_confs = 0
		for cid in sortedcids:
			sdwriter.write(conformers[cid])
			write_confs += 1

		if args.verbose == True: log.write("o  Writing "+str(write_confs)+ " conformers to file " + name+'_'+program+args.output)
		sdwriter.close()
	else: log.write("x  No conformers found!")

def mult_min(name, args, program,log,dup_data,dup_data_idx):
	'''optimizes a bunch of molecules and then checks for unique conformers and then puts in order of energy'''

	inmols = Chem.SDMolSupplier(name+args.output, removeHs=False)
	if inmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)

	globmin, n_high,n_dup_energy, n_dup_rms_eng  = None, 0, 0, 0
	c_converged, c_energy, outmols = [], [], []


	if args.verbose: log.write("\n\no  Multiple minimization of "+ name+args.output+ " with "+ program)
	bar = IncrementalBar('o  Minimizing', max = len(inmols))

	for i,mol in enumerate(inmols):
		bar.next()
		conf = 1
		if mol is not None:
			# optimize this structure and record the energy
			mol, converged, energy = optimize(mol, args, program,log)
			#if args.verbose: log.write("   conformer", (i+1), energy)

			if globmin == None: globmin = energy
			if energy < globmin: globmin = energy

			if converged == 0 and (energy - globmin) < args.ewin: # comparison in kcal/mol

				#if args.verbose: log.write('   minimization converged!')
				unique,dup_id = 0, None

				# compare against all previous conformers located
				for j,seenmol in enumerate(outmols):
					if abs(energy - c_energy[j]) < args.initial_energy_threshold: # comparison in kcal/mol
						unique += 1
						dup_id = (j+1)
						n_dup_energy += 1
						break
					#pmi_diff = get_PMIDIFF(mol, seenmol, 0, 0, args.heavyonly)
					#tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol, seenmol, useWeights=False)
					#rms = get_RMS(mol, seenmol, 0, 0, args.heavyonly)
					#log.write(rms, tfd, pmi_diff)
					if abs(energy - c_energy[j]) < args.energy_threshold: # comparison in kcal/mol
						rms = get_conf_RMS(mol, seenmol, 0, 0, args.heavyonly, args.max_matches_RMSD,log)
						if rms < args.rms_threshold:
							#log.write("o  Conformer", (i+1), "matches conformer", (j+1))
							unique += 1
							dup_id = (j+1)
							n_dup_rms_eng += 1
							break

				if unique == 0:
					#if args.verbose == True: log.write("-  Conformer", (i+1), "is unique")
					pmol = PropertyMol.PropertyMol(mol)
					outmols.append(pmol); c_converged.append(converged); c_energy.append(energy)
					conf += 1
					# if args.verbose == True:log.write("x  Conformer", (i+1), "is a duplicate of", dup_id)
			else:
				#if args.verbose == True: log.write("x  Minimization of conformer", (i+1), " not converged / energy too high!", converged, (energy - globmin), args.ewin)
				n_high += 1
		else:
			pass #log.write("No molecules to optimize")

	bar.finish()
	if args.verbose == True: log.write("o  "+str( n_dup_energy)+ " Duplicates removed initial energy ( E < "+str(args.initial_energy_threshold)+" kcal/mol )")
	if args.verbose == True: log.write("o  "+str( n_dup_rms_eng)+ " Duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol)")
	if args.verbose == True: log.write("o  "+str( n_high)+ " Conformers rejected based on energy ( E > "+str(args.ewin)+" kcal/mol)")

	# if SQM energy exists, overwrite RDKIT energies and geometries
	cids = list(range(len(outmols)))
	sortedcids = sorted(cids, key = lambda cid: c_energy[cid])

	for i, cid in enumerate(sortedcids):
		outmols[cid].SetProp('_Name', name +' conformer ' + str(i+1))
		outmols[cid].SetProp('Energy', c_energy[cid])

	if program == 'xtb':
		dup_data.at[dup_data_idx, 'xTB-Initial-samples'] = len(inmols)
		dup_data.at[dup_data_idx, 'xTB-initial_energy_threshold'] = n_dup_energy
		dup_data.at[dup_data_idx, 'xTB-RMS-and-energy-duplicates'] = n_dup_rms_eng
		dup_data.at[dup_data_idx, 'xTB-Unique-conformers'] = len(sortedcids)

	if program == 'ani':
		dup_data.at[dup_data_idx, 'ANI1ccx-Initial-samples'] = len(inmols)
		dup_data.at[dup_data_idx, 'ANI1ccx-initial_energy_threshold'] = n_dup_energy
		dup_data.at[dup_data_idx, 'ANI1ccx-RMS-and-energy-duplicates'] = n_dup_rms_eng
		dup_data.at[dup_data_idx, 'ANI1ccx-Unique-conformers'] = len(sortedcids)


	#bar.finish()
	# write the filtered, ordered conformers to external file
	write_confs(outmols, c_energy, name, args, program,log)
	return outmols, c_energy
