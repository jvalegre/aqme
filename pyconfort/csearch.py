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
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, rdMolAlign, Lipinski
from rdkit.Geometry import Point3D
from progress.bar import IncrementalBar
import pyconfort
from pyconfort.qprep_gaussian import write_confs
from pyconfort.filter import filters,set_metal_atomic_number,ewin_filter,pre_E_filter,RMSD_and_E_filter
from pyconfort.argument_parser import possible_atoms
from pyconfort.tmbuild import template_embed
from pyconfort.cmin import mult_min, rules_get_charge, atom_groups
from pyconfort.fullmonte import generating_conformations_fullmonte, minimize_rdkit_energy

hartree_to_kcal = 627.509
possible_atoms = possible_atoms()

#com to xyz to sdf for obabel
def com_2_xyz_2_sdf(args):

	if os.path.splitext(args.input)[1] =='.com' or os.path.splitext(args.input)[1] =='.gjf' :

		comfile = open(args.input,"r")
		comlines = comfile.readlines()

		emptylines=[]

		for i, line in enumerate(comlines):
			if len(line.strip()) == 0:
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

	cmd_obabel = ['obabel', '-ixyz', os.path.splitext(args.input)[0]+'.xyz', '-osdf', '-O', os.path.splitext(args.input)[0]+'.sdf','--gen3D']
	subprocess.run(cmd_obabel)
	if os.path.splitext(args.input)[1] =='.com' or os.path.splitext(args.input)[1] =='.gjf':
		return charge_com
	else:
		return args.charge_default

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
def mol_from_sdf_or_mol_or_mol2(args):
	if os.path.splitext(args.input)[1] =='.sdf':
		suppl = Chem.SDMolSupplier(args.input)
	elif os.path.splitext(args.input)[1] =='.mol':
		suppl = Chem.MolFromMolFile(args.input)
	elif os.path.splitext(args.input)[1] =='.mol2':
		suppl = Chem.MolFromMol2File(args.input)

	IDs,charges = [],[]
	readlines = open(args.input,"r").readlines()

	for i, line in enumerate(readlines):
		if line.find('>  <ID>') > -1:
			ID = readlines[i+1].split()[0]
			IDs.append(ID)
		if line.find('M  CHG') > -1:
			charge_line =  line.split('  ')
			charge = 0
			for j in range(4,len(charge_line)):
				if (j % 2) == 0:
					if j == len(charge_line) - 1:
						charge_line[j] = charge_line[j].split('\n')[0]
					charge += int(charge_line[j])
			charges.append(charge)
	if IDs == []:
		if os.path.splitext(args.input)[1] =='.sdf':
			for i,_ in enumerate(suppl):
				IDs.append(os.path.splitext(args.input)[0]+'_'+str(i))
		else:
			IDs.append(os.path.splitext(args.input)[0])
	if charges == []:
		if os.path.splitext(args.input)[1] =='.sdf':
			for i,_ in enumerate(suppl):
				charges.append(0)
		else:
			charges.append(0)
	return suppl, IDs, charges

# returns the arguments to their original value after each calculation
def clean_args(args,ori_ff,mol,ori_charge):
	for atom in mol.GetAtoms():
		if atom.GetSymbol() in args.metal:
			args.metal_complex= True
			break
	else:
		args.metal_complex = False
	args.ff = ori_ff
	args.charge_default = ori_charge
	args.metal_idx = []
	args.complex_coord = []
	args.metal_sym = []

#checks the charge on the smi string
def check_charge_smi(smi):
	charge = 0
	for i,smi_letter in enumerate(smi):
		if smi_letter == '+':
			if smi[i+1] == ']':
				charge += 1
			else:
				charge += int(smi[i+1])+1
		elif smi_letter == '-':
			if smi[i+1] == ']':
				charge -= 1
			else:
				charge -= int(smi[i+1])+1
	return charge

#checks for salts
def check_for_pieces(smi):
	#taking largest component for salts
	pieces = smi.split('.')
	if len(pieces) > 1:
		# take largest component by length
		smi = max(pieces, key=len)
	return smi

#if template activated, loads it
def load_template(args):
	try:
		os.chdir(os.path.join(pyconfort.__path__[0])+'/templates/')
	except FileNotFoundError:
		print('x The templates folder was not found, probably due to a problem while installing pyCONFORT')
		sys.exit()
	if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyramidal':
		file_template = 'template-4-and-5.sdf'
	if args.complex_type =='linear':
		file_template = 'template-2.sdf'
	if args.complex_type =='trigonalplanar':
		file_template = 'template-3.sdf'

	return file_template

#function to start conf generation
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
				file_template = load_template(args)
				temp = Chem.SDMolSupplier(file_template)
				os.chdir(w_dir_initial)
				mol_objects_from_template, name_mol, coord_Map, alg_Map, mol_template = template_embed(mol,temp,name,args,log)
				for j,_ in enumerate(mol_objects_from_template):
					mol_objects.append([mol_objects_from_template[j],name_mol[j],coord_Map[j],alg_Map[j],mol_template[j]])
				for [mol_object, name_mol, coord_Map, alg_Map, mol_template] in mol_objects:
					status = conformer_generation(mol_object,name_mol,start_time,args,log,dup_data,counter_for_template,coord_Map,alg_Map,mol_template)
					counter_for_template += 1
			else:
				log.write("x  Cannot use templates for complexes involving more than 1 metal or for organic molecueles.")
		else:
			status = conformer_generation(mol,name,start_time,args,log,dup_data,i)
	else:
		 status = conformer_generation(mol,name,start_time,args,log,dup_data,i)

	return status

# FUCNTION WORKING WITH MOL OBJECT TO CREATE CONFORMERS
def conformer_generation(mol,name,start_time,args,log,dup_data,dup_data_idx,coord_Map=None,alg_Map=None,mol_template=None):
	valid_structure = filters(mol, args,log)
	if valid_structure:
		if args.verbose:
			log.write("\n   ----- {} -----".format(name))
		try:
			# the conformational search for RDKit
			status = summ_search(mol, name,args,log,dup_data,dup_data_idx,coord_Map,alg_Map,mol_template)
			if args.CMIN=='ANI1ccx' or args.CMIN=='xtb':
				if status != -1:
					if args.CMIN=='ANI1ccx' and status != 0:
						min_suffix = 'ani'
						if args.CSEARCH=='rdkit':
							mult_min(name+'_'+'rdkit', args, min_suffix, log, dup_data, dup_data_idx)
						elif args.CSEARCH=='summ':
							mult_min(name+'_'+'summ', args, min_suffix, log, dup_data, dup_data_idx)
						elif args.CSEARCH=='fullmonte':
							mult_min(name+'_'+'fullmonte', args, min_suffix, log, dup_data, dup_data_idx)
					elif args.CMIN=='xtb' and status != 0:
						min_suffix = 'xtb'
						if args.CSEARCH=='rdkit':
							mult_min(name+'_'+'rdkit', args, min_suffix, log, dup_data, dup_data_idx)
						elif args.CSEARCH=='summ':
							mult_min(name+'_'+'summ', args, min_suffix, log, dup_data, dup_data_idx)
						elif args.CSEARCH=='fullmonte':
							mult_min(name+'_'+'fullmonte', args, min_suffix, log, dup_data, dup_data_idx)
					elif status == 0:
						os.remove(name+'_'+'rdkit'+args.output)
						log.write('\nx  No rotatable dihedral found. Run again with nodihedral set to TRUE')
			else:
				if status == 0:
					os.remove(name+'_'+'rdkit'+args.output)
					log.write('\nx  No rotatable dihedral found. Run again with nodihedral set to TRUE')

		except (KeyboardInterrupt, SystemExit):
			raise
	else:
		log.write("\nx  ERROR: The structure is not valid")

	# removing temporary files
	temp_files = ['gfn2.out', 'xTB_opt.traj', 'ANI1_opt.traj', 'wbo', 'xtbrestart','ase.opt','xtb.opt','gfnff_topo']
	for file in temp_files:
		if os.path.exists(file):
			os.remove(file)

	if args.time:
		log.write("\n Execution time: %s seconds" % (round(time.time() - start_time,2)))
		dup_data.at[dup_data_idx, 'time (seconds)'] = round(time.time() - start_time,2)

	return status

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
		if args.metal_complex and args.CSEARCH=='rdkit':
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

#function to embed conformers
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
	for _, conf in enumerate(cids):
		if coord_Map is None and alg_Map is None and mol_template is None:
			GetFF = minimize_rdkit_energy(mol,conf,args,log)
			cenergy.append(GetFF.CalcEnergy())

		# id template realign before doing calculations
		else:
			mol,GetFF = realign_mol(mol,conf,coord_Map, alg_Map, mol_template,args,log)
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
		outmols[cid].SetProp('_Name', name +' '+ str(i+1))
		outmols[cid].SetProp('Energy', cenergy[cid])

	# sorts the energies
	cids = list(range(len(outmols)))
	sorted_all_cids = sorted(cids,key = lambda cid: cenergy[cid])

	log.write("\n\no  Applying filters to intial conformers")

	# filter based on energy window ewin_csearch
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

	if args.CSEARCH=='summ' or args.CSEARCH=='rdkit':
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

	if args.CSEARCH=='summ':
		dup_data.at[dup_data_idx, 'summ-conformers'] = total

	if args.CSEARCH=='fullmonte':
		status = generating_conformations_fullmonte(name,args,rotmatches,log,selectedcids_rdkit,outmols,sdwriter,dup_data,dup_data_idx,coord_Map,alg_Map, mol_template)

	return status

#conversion from rdkit to sdf
def rdkit_to_sdf(mol, name,args,log,dup_data,dup_data_idx, coord_Map, alg_Map, mol_template):
	if args.CSEARCH =='rdkit' or args.CSEARCH =='summ':
		sdwriter = Chem.SDWriter(name+'_'+'rdkit'+args.output)
	elif args.CSEARCH =='fullmonte':
		sdwriter = Chem.SDWriter(name+'_'+'fullmonte'+args.output)
	Chem.SanitizeMol(mol)
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
	elif args.CSEARCH=='summ' and len(rotmatches) == 0:
		status = 0
	elif args.CSEARCH=='fullmonte' and len(rotmatches) == 0:
		status = 0
	else:
		dup_data.at[dup_data_idx, 'RDKIT-Initial-samples'] = initial_confs
		if args.CSEARCH=='rdkit':
			rotmatches =[]
		cids = embed_conf(mol,initial_confs,args,log,coord_Map,alg_Map, mol_template)
		#energy minimize all to get more realistic results
		#identify the atoms and decide Force Field
		for atom in mol.GetAtoms():
			if atom.GetAtomicNum() > 36: #up to Kr for MMFF, if not the code will use UFF
				args.ff = "UFF"
		if args.verbose:
			log.write("o  Optimizing "+ str(len(cids))+ " initial conformers with "+ args.ff)
			if args.CSEARCH=='summ':
				log.write("o  Found "+ str(len(rotmatches))+ " rotatable torsions")
			elif args.CSEARCH=='fullmonte':
				log.write("o  Found "+ str(len(rotmatches))+ " rotatable torsions")
			else:
				log.write("o  Systematic torsion rotation is set to OFF")

		status = min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,coord_Map,alg_Map, mol_template)
	sdwriter.close()

	return status,rotmatches

#steps to realign mol
def realign_mol(mol,conf,coord_Map, alg_Map, mol_template,args,log):
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
	return mol,GetFF

#filtering after dihydral scan to sdf
def dihedral_filter_and_sdf(name,args,log,dup_data,dup_data_idx,coord_Map, alg_Map, mol_template):
	rotated_energy = []
	rdmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)
	if rdmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)

	for i, rd_mol_i in enumerate(rdmols):
		rotated_energy.append(float(rd_mol_i.GetProp('Energy')))

	rotated_cids = list(range(len(rdmols)))
	sorted_rotated_cids = sorted(rotated_cids, key = lambda cid: rotated_energy[cid])

	# filter based on energy window ewin_csearch
	sortedcids_rotated = ewin_filter(sorted_rotated_cids,rotated_energy,args,dup_data,dup_data_idx,log,'summ')
	# pre-filter based on energy only
	selectedcids_initial_rotated = pre_E_filter(sortedcids_rotated,rotated_energy,args,dup_data,dup_data_idx,log,'summ')
	# filter based on energy and RMSD
	selectedcids_rotated = RMSD_and_E_filter(rdmols,selectedcids_initial_rotated,rotated_energy,args,dup_data,dup_data_idx,log,'summ')

	sdwriter_rd = Chem.SDWriter(name+'_'+'summ'+args.output)
	for i, cid in enumerate(selectedcids_rotated):
		mol_rd = Chem.RWMol(rdmols[cid])
		mol_rd.SetProp('_Name',rdmols[cid].GetProp('_Name')+' '+str(i))
		if coord_Map is None and alg_Map is None and mol_template is None:
			if args.metal_complex:
				set_metal_atomic_number(mol_rd,args)
			sdwriter_rd.write(mol_rd)
		else:
			mol_rd_realigned,_ = realign_mol(mol_rd,-1,coord_Map, alg_Map, mol_template,args,log)
			if args.metal_complex:
				set_metal_atomic_number(mol_rd_realigned,args)
			sdwriter_rd.write(mol_rd_realigned)

	sdwriter_rd.close()
	status = 1

	#removes the rdkit file
	os.remove(name+'_'+'rdkit'+args.output)

	return status

# EMBEDS, OPTIMIZES AND FILTERS RDKIT CONFORMERS
def summ_search(mol, name,args,log,dup_data,dup_data_idx, coord_Map = None, alg_Map=None, mol_template=None):
	# writes sdf for the first RDKit conformer generation
	status,rotmatches = rdkit_to_sdf(mol, name,args,log,dup_data,dup_data_idx, coord_Map, alg_Map, mol_template)

	# reads the initial SDF files from RDKit and uses dihedral scan if selected
	if status != -1 or status != 0:
		# getting the energy and mols after rotations
		if args.CSEARCH=='summ' and len(rotmatches) != 0:
			status = dihedral_filter_and_sdf(name,args,log,dup_data,dup_data_idx,coord_Map, alg_Map, mol_template)


	return status
