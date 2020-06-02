#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	    	  used in conformer generation		    #
#####################################################.

import math
import os
import sys
import subprocess
import time
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, rdMolAlign, Lipinski
from rdkit.Geometry import Point3D
from progress.bar import IncrementalBar
import pyconfort
from pyconfort.writer_functions import write_confs
from pyconfort.filter_functions import filters,get_conf_RMS,set_metal_atomic_number,ewin_filter,pre_E_filter,RMSD_and_E_filter
from pyconfort.argument_parser import possible_atoms
from pyconfort.analyzer_functions import check_for_final_folder
from pyconfort.template_functions import template_embed

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
	from xtb-python import GFN2
except:
	print('1')
try:
	import torchani
	model = torchani.models.ANI1ccx()
except:
	print('5')

hartree_to_kcal = 627.509
possible_atoms = possible_atoms()

#com to xyz to sdf for obabel
def com_2_xyz_2_sdf(args):
	comfile = open(args.input,"r")
	comlines = comfile.readlines()

	emptylines=[]

	for i, line in enumerate(comlines):
		if len(comlines[i].strip()) == 0:
			emptylines.append(i)

	#assigning the charges
	charge_com = comlines[(emptylines[1]+1)].split(' ')[0]

	xyzfile = open(os.path.splitext(args.input)[0]+'.xyz',"w")
	xyzfile.write(str(emptylines[2]- (emptylines[1]+2)))
	xyzfile.write('\n')
	xyzfile.write(os.path.splitext(args.input)[0])
	xyzfile.write('\n')
	for i in range((emptylines[1]+2), emptylines[2]):
		xyzfile.write(comlines[i])

	xyzfile.close()
	comfile.close()

	subprocess.run(['obabel', '-ixyz', os.path.splitext(args.input)[0]+'.xyz', '-osdf', '-O', os.path.splitext(args.input)[0]+'.sdf','--gen3D'])

	return charge_com

# SUBSTITUTION WITH I
def substituted_mol(mol,args,log):
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

#mol from sdf
def mol_from_sdf(args):
	suppl = Chem.SDMolSupplier(args.input)
	IDs,charges = [],[]
	f = open(args.input,"r")
	readlines = f.readlines()
	for i, line in enumerate(readlines):
		if line.find('>  <ID>') > -1:
			ID = readlines[i+1].split()[0]
			IDs.append(ID)
		if line.find('M  CHG') > -1:
			charge_line =  line.split('  ')
			charge = 0
			for i in range(4,len(charge_line)):
				if (i % 2) == 0:
					if i == len(charge_line) - 1:
						charge_line[i] = charge_line[i].split('\n')[0]
					charge += int(charge_line[i])
			print(charge)
			charges.append(charge)
	if IDs == []:
		for i,_ in enumerate(suppl):
			IDs.append(os.path.splitext(args.input)[0]+'_'+str(i))
	if charges == []:
		for i,_ in enumerate(suppl):
			charges.append(0)
	return suppl, IDs, charges

def clean_args(args,ori_ff,mol):
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

def check_charge_smi(smi):
	count_plus,count_minus = 0,0
	for i in smi:
		if i == '+':
			count_plus = count_plus + 1
		if i == '-':
			count_minus = count_minus + 1
	charge = count_plus + count_minus
	return charge

def check_for_pieces(smi):
	#taking largest component for salts
	pieces = smi.split('.')
	if len(pieces) > 1:
		# take largest component by length
		smi = max(pieces, key=len)
	return smi

def compute_confs(w_dir_initial,mol, name,args,log,dup_data,counter_for_template,i,start_time):
	# Converts each line to a rdkit mol object
	if args.verbose:
		log.write("   -> Input Molecule {} is {}".format(i, Chem.MolToSmiles(mol)))

	if args.metal_complex:
		mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(mol,args,log)

	if args.metal_complex:
		# get manually for square planar and squarepyramidal
		if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyramidal' or args.complex_type == 'linear' or args.complex_type == 'trigonalplanar':
			mol_objects = []
			if len(args.metal_idx) == 1:
				try:
					os.chdir(os.path.join(pyconfort.__path__[0])+'/templates/')
				except FileNotFoundError:
					os.chdir(os.path.join(pyconfort.__path__[0])+'/templates/')
				if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyramidal':
					file_template = 'template-4-and-5.sdf'
				if args.complex_type =='linear':
					file_template = 'template-2.sdf'
				if args.complex_type =='trigonalplanar':
					file_template = 'template-3.sdf'
				temp = Chem.SDMolSupplier(file_template)
				os.chdir(w_dir_initial)
				mol_objects_from_template,name, coord_Map, alg_Map, mol_template = template_embed(mol,temp,name,args,log)
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
			set_metal_atomic_number(mol,args)
		sdwriter.write(mol,conf)
		return 1
	else:
		total = 0
		deg = 0
		while deg < 360.0:
			rad = math.pi*deg / 180.0

			rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
			#recalculating energies after rotation
			GetFF = minimize_rdkit_energy(mol,conf,args,log)
			mol.SetProp("Energy",GetFF.CalcEnergy())
			mol.SetProp('_Name',name)
			total += genConformer_r(mol, conf, i+1, matches, degree, sdwriter,args,name,log)
			deg += degree
		return total

def minimize_rdkit_energy(mol,conf,args,log):
	if args.ff == "MMFF":
		GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
	elif args.ff == "UFF":
		GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
	else:
		log.write(' Force field {} not supported!'.format(args.ff))
		sys.exit()
	GetFF.Initialize()
	GetFF.Minimize(maxIts=args.opt_steps_RDKit)
	return GetFF

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
			for atom in neighbours:
				if atom.GetTotalValence()== 4:
					if atom.GetSymbol() in C_group:
						charge[charge_idx] = charge[charge_idx] - 1
					elif atom.GetSymbol() in N_group:
						charge[charge_idx] = charge[charge_idx] - 0
				elif atom.GetTotalValence()== 3:
					if atom.GetSymbol() in C_group or atom.GetSymbol() in O_group:
						charge[charge_idx] = charge[charge_idx] - 0
					elif atom.GetSymbol() in N_group:
						charge[charge_idx] = charge[charge_idx] - 1
				elif atom.GetTotalValence() == 2:
					if atom.GetSymbol() in O_group:
						charge[charge_idx] = charge[charge_idx] - 1
					elif atom.GetSymbol() in F_group:
						charge[charge_idx] = charge[charge_idx] - 0
				elif atom.GetTotalValence() == 1:
					if atom.GetSymbol() in F_group:
						charge[charge_idx] = charge[charge_idx] - 1
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

# minimization and E calculation with RDKit after embeding
def min_and_E_calc(mol,cids,args,log,coord_Map,alg_Map,mol_template):
	cenergy,outmols = [],[]
	bar = IncrementalBar('o  Minimizing', max = len(cids))
	for i, conf in enumerate(cids):
		if coord_Map is None and alg_Map is None and mol_template is None:
			GetFF = minimize_rdkit_energy(mol,conf,args,log)
			cenergy.append(GetFF.CalcEnergy())

		# id template realign before doing calculations
		else:
			num_atom_match = mol.GetSubstructMatch(mol_template)
			GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
			for k, idxI in enumerate(num_atom_match):
				for l in range(k + 1, len(num_atom_match)):
					idxJ = num_atom_match[l]
					d = coord_Map[idxI].Distance(coord_Map[idxJ])
					GetFF.AddDistanceConstraint(idxI, idxJ, d, d, 10000)
			GetFF.Initialize()
			GetFF.Minimize(maxIts=args.opt_steps_RDKit)
			# rotate the embedded conformation onto the core_mol:
			rdMolAlign.AlignMol(mol, mol_template, prbCid=conf,refCid=-1,atomMap=alg_Map,reflect=True,maxIters=100)
			cenergy.append(GetFF.CalcEnergy())
		# outmols is gonna be a list containing "initial_confs" mol objects with "initial_confs" conformers. We do this to SetProp (Name and Energy) to the different conformers
		# and log.write in the SDF file. At the end, since all the mol objects has the same conformers, but the energies are different, we can log.write conformers to SDF files
		# with the energies of the parent mol objects. We measured the computing time and it's the same as using only 1 parent mol object with 10 conformers, but we couldn'temp SetProp correctly
		pmol = PropertyMol.PropertyMol(mol)
		outmols.append(pmol)
		bar.next()
	bar.finish()
	return outmols,cenergy

# minimizes, gets the energy and filters RDKit conformers after embeding
def min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,coord_Map,alg_Map, mol_template):
	# gets optimized mol objects and energies
	outmols,cenergy = min_and_E_calc(mol,cids,args,log,coord_Map,alg_Map,mol_template)

	for i, cid in enumerate(cids):
		outmols[cid].SetProp('_Name', name + ' conformer ' + str(i+1))
		outmols[cid].SetProp('Energy', cenergy[cid])

	# sorts the energies
	cids = list(range(len(outmols)))
	sorted_all_cids = sorted(cids,key = lambda cid: cenergy[cid])

	log.write("\n\no  Applying filters to intial conformers")

	# filter based on energy window ewin_rdkit
	sortedcids_rdkit = ewin_filter(sorted_all_cids,cenergy,args,dup_data,dup_data_idx,log,'rdkit')

	# pre-filter based on energy only
	selectedcids_initial_rdkit = pre_E_filter(sortedcids_rdkit,cenergy,args,dup_data,dup_data_idx,log,'rdkit')

	# filter based on energy and RMSD
	selectedcids_rdkit = RMSD_and_E_filter(outmols,selectedcids_initial_rdkit,cenergy,args,dup_data,dup_data_idx,log,'rdkit')

	# writing charges after RDKIT
	if os.path.splitext(args.input)[1] == '.cdx' or os.path.splitext(args.input)[1] == '.smi' or os.path.splitext(args.input)[1] == '.csv':
		args.charge = rules_get_charge(mol,args,log)
		dup_data.at[dup_data_idx, 'Overall charge'] = np.sum(args.charge)
	else:
		dup_data.at[dup_data_idx, 'Overall charge'] = args.charge_default

	# now exhaustively drive torsions of selected conformers
	n_confs = int(len(selectedcids_rdkit) * (360 / args.degree) ** len(rotmatches))
	if args.verbose and len(rotmatches) != 0:
		log.write("\n\no  Systematic generation of "+ str(n_confs)+ " confomers")
		bar = IncrementalBar('o  Generating conformations based on dihedral rotation', max = len(selectedcids_rdkit))
	else:
		bar = IncrementalBar('o  Generating conformations', max = len(selectedcids_rdkit))

	total = 0
	for conf in selectedcids_rdkit:
		total += genConformer_r(outmols[conf], conf, 0, rotmatches, args.degree, sdwriter ,args,outmols[conf].GetProp('_Name'),log)
		bar.next()
	bar.finish()
	if args.verbose and len(rotmatches) != 0:
		log.write("o  %d total conformations generated"%total)
	status = 1

	if not args.nodihedrals:
		dup_data.at[dup_data_idx, 'RDKIT-Rotated-conformers'] = total

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
			if not args.nodihedrals:
				log.write("o  Found "+ str(len(rotmatches))+ " rotatable torsions")
			else:
				log.write("o  Systematic torsion rotation is set to OFF")

		status = min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,coord_Map,alg_Map, mol_template)
	sdwriter.close()

	if status != -1:
		#getting the energy from and mols after rotations
		if not args.nodihedrals and len(rotmatches) != 0:
			rotated_energy = []
			rdmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)
			if rdmols is None:
				log.write("Could not open "+ name+args.output)
				sys.exit(-1)

			for i, rd_mol_i in enumerate(rdmols):
				rotated_energy.append(float(rd_mol_i.GetProp('Energy')))

			rotated_cids = list(range(len(rdmols)))
			sorted_rotated_cids = sorted(rotated_cids, key = lambda cid: rotated_energy[cid])

			# filter based on energy window ewin_rdkit
			sortedcids_rotated = ewin_filter(sorted_rotated_cids,rotated_energy,args,dup_data,dup_data_idx,log,'rotated_rdkit')
			# pre-filter based on energy only
			selectedcids_initial_rotated = pre_E_filter(sortedcids_rotated,rotated_energy,args,dup_data,dup_data_idx,log,'rotated_rdkit')
			# filter based on energy and RMSD
			selectedcids_rotated = RMSD_and_E_filter(rdmols,selectedcids_initial_rotated,rotated_energy,args,dup_data,dup_data_idx,log,'rotated_rdkit')

			sdwriter_rd = Chem.SDWriter(name+'_'+'rdkit'+'_'+'rotated'+args.output)
			for i, cid in enumerate(selectedcids_rotated):
				mol_rd = Chem.RWMol(rdmols[cid])
				mol_rd.SetProp('_Name',rdmols[cid].GetProp('_Name')+' '+str(i))
				if args.metal_complex:
					set_metal_atomic_number(mol_rd,args)
				sdwriter_rd.write(mol_rd)
			sdwriter_rd.close()
			status = 1
			#status=filter_after_rotation(args,name,log,dup_data,dup_data_idx)
		elif not args.nodihedrals and len(rotmatches) ==0:
			status = 0
	return status

# ANI1 OPTIMIZER AND ENERGY CALCULATOR
def ani_calc(elements,cartesians,coordinates,args,log):
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
	# Now let's compute energy:
	xtb_energy = ase_molecule.get_potential_energy()
	sqm_energy = (xtb_energy / Hartree)* hartree_to_kcal

	return sqm_energy, coordinates

# xTB AND ANI1 MAIN OPTIMIZATION PROCESS
def optimize(mol, args, program,log,dup_data,dup_data_idx):
	# if large system increase stack size
	if args.large_sys:
		os.environ['OMP_STACKSIZE'] = args.STACKSIZE

	# removing the Ba atom if NCI complexes
	if args.nci_complex:
		for atom in mol.GetAtoms():
			if atom.GetSymbol() =='I':
				atom.SetAtomicNum(1)

	if args.metal_complex and not args.nodihedrals:
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
		outmols[cid].SetProp('_Name', name_mol +' conformer ' + str(i+1))
		outmols[cid].SetProp('Energy', cenergy[cid])

	log.write("\n\no  Applying filters to intial conformers")
	# filter based on energy window ewin_rdkit
	sortedcids = ewin_filter(sorted_all_cids,cenergy,args,dup_data,dup_data_idx,log,'xtb_ani')
	# pre-filter based on energy only
	selectedcids_initial = pre_E_filter(sortedcids,cenergy,args,dup_data,dup_data_idx,log,'xtb_ani')
	# filter based on energy and RMSD
	selectedcids = RMSD_and_E_filter(outmols,selectedcids_initial,cenergy,args,dup_data,dup_data_idx,log,'xtb_ani')

	if program == 'xtb':
		dup_data.at[dup_data_idx, 'xTB-Initial-samples'] = len(inmols)
	if program == 'ani':
		dup_data.at[dup_data_idx, 'ANI1ccx-Initial-samples'] = len(inmols)

	# write the filtered, ordered conformers to external file
	write_confs(outmols, cenergy,selectedcids, name, args, program,log)
