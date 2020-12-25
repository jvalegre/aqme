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
from pyconfort.qprep_gaussian import write_confs
from pyconfort.filter import filters,set_metal_atomic_number,ewin_filter,pre_E_filter,RMSD_and_E_filter
from pyconfort.argument_parser import possible_atoms
from pyconfort.tmbuild import template_embed
from pyconfort.cmin import mult_min, rules_get_charge, atom_groups,substituted_mol
from pyconfort.fullmonte import generating_conformations_fullmonte, minimize_rdkit_energy,realign_mol


hartree_to_kcal = 627.509
possible_atoms = possible_atoms()


#class for logging
class Logger:
	# Class Logger to writargs.input.split('.')[0] output to a file
	def __init__(self, filein, append):
		# Logger to write the output to a file
		suffix = 'dat'
		self.log = open('{0}_{1}.{2}'.format(filein, append, suffix), 'w')

	def write(self, message):
		#print(message, end='\n')
		self.log.write(message+ "\n")

	def fatal(self, message):
		#print(message, end='\n')
		self.log.write(message + "\n")
		self.finalize()
		sys.exit(1)

	def finalize(self):
		self.log.close()

#creation of csv for csearch
def creation_of_dup_csv(args):
	# writing the list of DUPLICATES
	if args.CSEARCH=='rdkit':
		if not args.CMIN=='xtb' and not args.CMIN=='ani':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','CSEARCH time (seconds)','Overall charge'])
		elif args.CMIN=='xtb' and not args.CMIN=='ani':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','CSEARCH time (seconds)','CMIN time (seconds)','Overall charge'])
		elif args.CMIN=='ani' and not args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','ANI-Initial-samples','ANI-energy-window','ANI-initial_energy_threshold','ANI-RMSD-and-energy-duplicates','ANI-Unique-conformers','CSERACH time (seconds)','CMIN time (seconds)','Overall charge'])
		elif args.CMIN=='ani' and args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','ANI-Initial-samples','ANI-energy-window','ANI-initial_energy_threshold','ANI-RMSD-and-energy-duplicates','ANI-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','CSEARCH time (seconds)','CMIN time (seconds)','Overall charge'])
	elif args.CSEARCH=='fullmonte':
		if not args.CMIN=='xtb' and not args.CMIN=='ani':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-Unique-conformers','CSEARCH time (seconds)','Overall charge'])# ,'FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates',
		elif args.CMIN=='xtb' and not args.CMIN=='ani':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','CSEARCH time (seconds)','CMIN time (seconds)','Overall charge']) #'FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates',
		elif args.CMIN=='ani' and not args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-Unique-conformers','ANI-Initial-samples','ANI-energy-window','ANI-initial_energy_threshold','ANI-RMSD-and-energy-duplicates','ANI-Unique-conformers','CSEARCH time (seconds)','CMIN time (seconds)','Overall charge'])#'FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates',
		elif args.CMIN=='ani' and args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-Unique-conformers','ANI-Initial-samples','ANI-energy-window','ANI-initial_energy_threshold','ANI-RMSD-and-energy-duplicates','ANI-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','CSEARCH time (seconds)','CMIN time (seconds)','Overall charge']) #'FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates',
	elif args.CSEARCH=='summ':
		if not args.CMIN=='xtb' and not args.CMIN=='ani':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','CSEARCH time (seconds)','Overall charge'])
		elif args.CMIN=='xtb' and not args.CMIN=='ani':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window','RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','CSEARCH time (seconds)','CMIN time (seconds)','Overall charge'])
		elif args.CMIN=='ani' and not args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','ANI-Initial-samples','ANI-energy-window','ANI-initial_energy_threshold','ANI-RMSD-and-energy-duplicates','ANI-Unique-conformers','CSEARCH time (seconds)','CMIN time (seconds)','Overall charge'])
		elif args.CMIN=='ani' and args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','ANI-Initial-samples','ANI-energy-window','ANI-initial_energy_threshold','ANI-RMSD-and-energy-duplicates','ANI-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','CSEARCH time (seconds)',' CMIN time (seconds)','Overall charge'])
	return dup_data

#com to xyz to sdf for obabel
def com_2_xyz_2_sdf(args,start_point=None):

	if start_point is None:
		if os.path.splitext(args.input)[1] =='.com' or os.path.splitext(args.input)[1] =='.gjf' or os.path.splitext(args.input)[1] =='.xyz':
			file = args.input

	elif start_point is not None:
		file = start_point

	if os.path.splitext(args.input)[1] !='.xyz':
		comfile = open(file,"r")
		comlines = comfile.readlines()

		emptylines=[]

		for i, line in enumerate(comlines):
			if len(line.strip()) == 0:
				emptylines.append(i)

		#assigning the charges
		charge_com = comlines[(emptylines[1]+1)].split(' ')[0]

		xyzfile = open(os.path.splitext(file)[0]+'.xyz',"w")
		xyzfile.write(str(emptylines[2]- (emptylines[1]+2)))
		xyzfile.write('\n')
		xyzfile.write(os.path.splitext(file)[0])
		xyzfile.write('\n')
		for i in range((emptylines[1]+2), emptylines[2]):
			xyzfile.write(comlines[i])

		xyzfile.close()
		comfile.close()

	cmd_obabel = ['obabel', '-ixyz', os.path.splitext(file)[0]+'.xyz', '-osdf', '-O', os.path.splitext(file)[0]+'.sdf']
	subprocess.run(cmd_obabel)

	if start_point is None:
		if os.path.splitext(args.input)[1] =='.com' or os.path.splitext(args.input)[1] =='.gjf':
			return charge_com
		else:
			return args.charge_default


#mol from sdf
def mol_from_sdf_or_mol_or_mol2(input):
	if os.path.splitext(input)[1] =='.sdf':
		suppl = Chem.SDMolSupplier(input, removeHs=False)
	elif os.path.splitext(input)[1] =='.mol':
		suppl = Chem.MolFromMolFile(input, removeHs=False)
	elif os.path.splitext(input)[1] =='.mol2':
		suppl = Chem.MolFromMol2File(input, removeHs=False)

	IDs,charges = [],[]

	readlines = open(input,"r").readlines()

	molecule_count = 0
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
		if line.find('$$$$') > -1:
			molecule_count += 1
			if molecule_count != len(charges):
				charges.append(0)

	if IDs == []:
		if os.path.splitext(input)[1] =='.sdf':
			for i,_ in enumerate(suppl):
				IDs.append(os.path.splitext(input)[0]+'_'+str(i))
		else:
			IDs.append(os.path.splitext(input)[0])
	if charges == []:
		if os.path.splitext(input)[1] =='.sdf':
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
				charge += int(smi[i+1])
		elif smi_letter == '-':
			if smi[i+1] == ']':
				charge -= 1
			else:
				charge -= int(smi[i+1])
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
		log.write('x The templates folder was not found, probably due to a problem while installing pyCONFORT')
		sys.exit()
	if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyramidal':
		file_template = 'template-4-and-5.sdf'
	if args.complex_type =='linear':
		file_template = 'template-2.sdf'
	if args.complex_type =='trigonalplanar':
		file_template = 'template-3.sdf'

	return file_template

#function to start conf generation
def compute_confs(w_dir_initial, mol, name, args,i):
	try:
		os.makedirs(w_dir_initial+'/CSEARCH/dat_files')
	except OSError:
		if os.path.isdir(w_dir_initial+'/CSEARCH/dat_files'):
			pass

	log = Logger(w_dir_initial+'/CSEARCH/dat_files/'+name, args.output_name)
	# Converts each line to a rdkit mol object
	if args.verbose:
		log.write("   -> Input Molecule {} is {}".format(i, Chem.MolToSmiles(mol)))

	if args.metal_complex:
		for i,_ in enumerate(args.metal):
			args.metal_idx.append(None)
			args.complex_coord.append(None)
			args.metal_sym.append(None)

		mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(mol,args,log)

		# get pre-determined geometries for metal complexes
		if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyramidal' or args.complex_type == 'linear' or args.complex_type == 'trigonalplanar':
			mol_objects, count_metals = [],0
			for i,metal_idx_ind in enumerate(args.metal_idx):
				if metal_idx_ind is not None:
					count_metals += 1
			if count_metals == 1:
				file_template = load_template(args)
				temp = Chem.SDMolSupplier(file_template)
				os.chdir(w_dir_initial)
				mol_objects_from_template, name_mol, coord_Map, alg_Map, mol_template = template_embed(mol,temp,name,args,log)
				for j,_ in enumerate(mol_objects_from_template):
					mol_objects.append([mol_objects_from_template[j],name_mol[j],coord_Map[j],alg_Map[j],mol_template[j]])
				total_data = creation_of_dup_csv(args)
				for [mol_object, name_mol, coord_Map, alg_Map, mol_template] in mol_objects:
					data = conformer_generation(mol_object,name_mol,args,log,coord_Map,alg_Map,mol_template)
					frames = [total_data, data]
					total_data = pd.concat(frames,sort=True)
			else:
				log.write("x  Cannot use templates for complexes involving more than 1 metal or for organic molecueles.")
				total_data = None
		else:
			total_data = conformer_generation(mol,name,args,log)
	else:
		total_data = conformer_generation(mol,name,args,log)
	return total_data

# FUCNTION WORKING WITH MOL OBJECT TO CREATE CONFORMERS
def conformer_generation(mol,name,args,log,coord_Map=None,alg_Map=None,mol_template=None):
	dup_data = creation_of_dup_csv(args)
	dup_data_idx = 0
	start_time = time.time()
	valid_structure = filters(mol, args,log)
	if valid_structure:
		if args.verbose:
			log.write("\n   ----- {} -----".format(name))
		try:
			# the conformational search for RDKit
			status,update_to_rdkit = summ_search(mol, name,args,log,dup_data,dup_data_idx,coord_Map,alg_Map,mol_template)
			dup_data.at[dup_data_idx, 'status'] = status
			dup_data.at[dup_data_idx, 'update_to_rdkit'] = update_to_rdkit
		except (KeyboardInterrupt, SystemExit):
			raise
	else:
		log.write("\nx  ERROR: The structure is not valid")

	if args.time:
		log.write("\n Execution time CSEARCH: %s seconds" % (round(time.time() - start_time,2)))
		dup_data.at[dup_data_idx, 'CSEARCH time (seconds)'] = round(time.time() - start_time,2)
	return dup_data

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
def genConformer_r(mol, conf, i, matches, degree, sdwriter,args,name,log,update_to_rdkit,coord_Map,alg_Map, mol_template):
	if i >= len(matches): # base case, torsions should be set in conf
		#setting the metal back instead of I
		if args.metal_complex and (args.CSEARCH=='rdkit' or update_to_rdkit):
			if coord_Map is None and alg_Map is None and mol_template is None:
				GetFF = minimize_rdkit_energy(mol,conf,args,log)
				mol.SetProp('Energy',str(GetFF.CalcEnergy()))
			else:
				mol,GetFF = realign_mol(mol,conf,coord_Map, alg_Map, mol_template,args,log)
				mol.SetProp('Energy',str(GetFF.CalcEnergy()))
			set_metal_atomic_number(mol,args)
		sdwriter.write(mol,conf)
		return 1
	else:
		total = 0
		deg = 0
		while deg < 360.0:
			rad = math.pi*deg / 180.0
			rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
			mol.SetProp('_Name',name)
			total += genConformer_r(mol, conf, i+1, matches, degree, sdwriter,args,name,log,update_to_rdkit,coord_Map,alg_Map, mol_template)
			deg += degree
		return total

#function to embed conformers
def embed_conf(mol,initial_confs,args,log,coord_Map,alg_Map, mol_template):
	if os.path.splitext(args.input)[1] == '.sdf' or os.path.splitext(args.input)[1] == '.mol' or os.path.splitext(args.input)[1] == '.mol2':
		Chem.AssignStereochemistryFrom3D(mol)
	if coord_Map is None and alg_Map is None and mol_template is None:
		cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs,ignoreSmoothingFailures=True, randomSeed=args.seed,numThreads = 0)
		if len(cids) == 0 or len(cids) == 1 and initial_confs != 1:
			log.write("o  Normal RDKit embeding process failed, trying to generate conformers with random coordinates (with "+str(initial_confs)+" possibilities)")
			cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0,ignoreSmoothingFailures=True, numZeroFail=1000, numThreads = 0)
		if args.verbose:
			log.write("o  "+ str(len(cids))+" conformers initially generated with RDKit")
	# case of embed for templates
	else:
		cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
		if len(cids) == 0 or len(cids) == 1 and initial_confs != 1:
			log.write("o  Normal RDKit embeding process failed, trying to generate conformers with random coordinates (with "+str(initial_confs)+" possibilities)")
			cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, randomSeed=args.seed, useRandomCoords=True, boxSizeMult=10.0, numZeroFail=1000,ignoreSmoothingFailures=True, coordMap = coord_Map,numThreads = 0)
		if args.verbose:
			log.write("o  "+ str(len(cids))+" conformers initially generated with RDKit")

	if os.path.splitext(args.input)[1] == '.sdf' or os.path.splitext(args.input)[1] == '.mol' or os.path.splitext(args.input)[1] == '.mol2':
		#preserving AssignStereochemistryFrom3D
		for id in cids:
			Chem.AssignAtomChiralTagsFromStructure(mol,confId=id)

	return cids

# minimization and E calculation with RDKit after embeding
def min_and_E_calc(mol,cids,args,log,coord_Map,alg_Map,mol_template):
	cenergy,outmols = [],[]
	#bar = IncrementalBar('o  Minimizing', max = len(cids))
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
		#bar.next()
	#bar.finish()
	return outmols,cenergy

# minimizes, gets the energy and filters RDKit conformers after embeding
def min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,update_to_rdkit,coord_Map,alg_Map, mol_template):
	# gets optimized mol objects and energies
	outmols,cenergy = min_and_E_calc(mol,cids,args,log,coord_Map,alg_Map,mol_template)

	# writing charges after RDKIT
	if os.path.splitext(args.input)[1] == '.cdx' or os.path.splitext(args.input)[1] == '.smi' or os.path.splitext(args.input)[1] == '.csv':
		args.charge = rules_get_charge(mol,args,log)
		dup_data.at[dup_data_idx, 'Overall charge'] = np.sum(args.charge)
	else:
		dup_data.at[dup_data_idx, 'Overall charge'] = args.charge_default

	for i, cid in enumerate(cids):
		outmols[cid].SetProp('_Name', name +' '+ str(i+1))
		outmols[cid].SetProp('Energy', str(cenergy[cid]))
		outmols[cid].SetProp('Real charge', str(dup_data.at[dup_data_idx, 'Overall charge']))

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

	if args.CSEARCH=='summ' or args.CSEARCH=='rdkit':
		# now exhaustively drive torsions of selected conformers
		n_confs = int(len(selectedcids_rdkit) * (360 / args.degree) ** len(rotmatches))
		if args.verbose and len(rotmatches) != 0:
			log.write("\n\no  Systematic generation of "+ str(n_confs)+ " confomers")
			#bar = IncrementalBar('o  Generating conformations based on dihedral rotation', max = len(selectedcids_rdkit))
		# else:
		# 	bar = IncrementalBar('o  Writing unique conformers into an sdf file', max = len(selectedcids_rdkit))

		total = 0
		for conf in selectedcids_rdkit:
			if args.CSEARCH=='summ' and not update_to_rdkit:
				sdwriter.write(outmols[conf],conf)
				for m in rotmatches:
					rdMolTransforms.SetDihedralDeg(outmols[conf].GetConformer(conf),*m,180.0)
			total += genConformer_r(outmols[conf], conf, 0, rotmatches, args.degree, sdwriter ,args,outmols[conf].GetProp('_Name'),log,update_to_rdkit,coord_Map,alg_Map, mol_template)
			# bar.next()
		# bar.finish()
		if args.verbose and len(rotmatches) != 0:
			log.write("o  %d total conformations generated"%total)
		status = 1

	if args.CSEARCH=='summ':
		dup_data.at[dup_data_idx, 'summ-conformers'] = total

	if args.CSEARCH=='fullmonte':
		status = generating_conformations_fullmonte(name,args,rotmatches,log,selectedcids_rdkit,outmols,sdwriter,dup_data,dup_data_idx,coord_Map,alg_Map, mol_template)
		#removes the rdkit file
		os.remove(name+'_'+'rdkit'+args.output)

	return status

#conversion from rdkit to sdf
def rdkit_to_sdf(mol, name,args,log,dup_data,dup_data_idx, coord_Map, alg_Map, mol_template):

	Chem.SanitizeMol(mol)

	mol = Chem.AddHs(mol)

	mol.SetProp("_Name",name)

	# detects and applies auto-detection of initial number of conformers
	if args.sample == 'auto':
		initial_confs = int(auto_sampling(args.auto_sample,mol,args,log))
	else:
		initial_confs = int(args.sample)

	dup_data.at[dup_data_idx, 'Molecule'] = name

	update_to_rdkit=False

	rotmatches = getDihedralMatches(mol, args.heavyonly,log)
	if len(rotmatches) > args.max_torsions:
		log.write("x  Too many torsions (%d). Skipping %s" %(len(rotmatches),(name+args.output)))
		status = -1
	elif args.CSEARCH=='summ' and len(rotmatches) == 0:
		update_to_rdkit = True
		log.write('\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to SUMM SDF')
	elif args.CSEARCH=='fullmonte' and len(rotmatches) == 0:
		update_to_rdkit = True
		log.write('\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to FULLMONTE SDF')

	if update_to_rdkit and args.CSEARCH =='summ' :
		sdwriter = Chem.SDWriter(name+'_'+'summ'+args.output)
	elif update_to_rdkit and args.CSEARCH =='fullmonte':
		sdwriter = Chem.SDWriter(name+'_'+'fullmonte'+args.output)
	elif args.CSEARCH =='fullmonte':
		sdwriter = Chem.SDWriter(name+'_'+'fullmonte'+args.output)
	else:
		sdwriter = Chem.SDWriter(name+'_'+'rdkit'+args.output)

	dup_data.at[dup_data_idx, 'RDKIT-Initial-samples'] = initial_confs
	if args.CSEARCH=='rdkit':
		rotmatches =[]
	cids = embed_conf(mol,initial_confs,args,log,coord_Map,alg_Map, mol_template)

	#energy minimize all to get more realistic results
	#identify the atoms and decide Force Field
	for atom in mol.GetAtoms():
		if atom.GetAtomicNum() > 36: #up to Kr for MMFF, if not the code will use UFF
			log.write("x  "+args.ff+" is not compatible with the molecule, changing to UFF")
			args.ff = "UFF"
	if args.verbose:
		log.write("o  Optimizing "+ str(len(cids))+ " initial conformers with "+ args.ff)
		if args.CSEARCH=='summ':
			log.write("o  Found "+ str(len(rotmatches))+ " rotatable torsions")
		elif args.CSEARCH=='fullmonte':
			log.write("o  Found "+ str(len(rotmatches))+ " rotatable torsions")
		else:
			log.write("o  Systematic torsion rotation is set to OFF")

	status = min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,update_to_rdkit,coord_Map,alg_Map, mol_template)
	sdwriter.close()

	return status,rotmatches,update_to_rdkit

#filtering after dihydral scan to sdf
def dihedral_filter_and_sdf(name,args,log,dup_data,dup_data_idx,coord_Map, alg_Map, mol_template):
	rotated_energy = []

	rdmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)

	if rdmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)

	for i, rd_mol_i in enumerate(rdmols):
		if coord_Map is None and alg_Map is None and mol_template is None:
			GetFF = minimize_rdkit_energy(rd_mol_i,-1,args,log)
			rotated_energy.append(float(GetFF.CalcEnergy()))
		else:
			rd_mol_i,GetFF = realign_mol(rd_mol_i,-1,coord_Map, alg_Map, mol_template,args,log)
			rotated_energy.append(float(GetFF.CalcEnergy()))

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
		mol_rd.SetProp('Energy',str(rotated_energy[cid]))
		if args.metal_complex:
			set_metal_atomic_number(mol_rd,args)
		sdwriter_rd.write(mol_rd)
	sdwriter_rd.close()
	status = 1
	return status

# EMBEDS, OPTIMIZES AND FILTERS RDKIT CONFORMERS
def summ_search(mol, name,args,log,dup_data,dup_data_idx, coord_Map = None, alg_Map=None, mol_template=None):
	# writes sdf for the first RDKit conformer generation
	status,rotmatches,update_to_rdkit = rdkit_to_sdf(mol, name,args,log,dup_data,dup_data_idx, coord_Map, alg_Map, mol_template)

	# reads the initial SDF files from RDKit and uses dihedral scan if selected
	if status != -1 or status != 0:
		# getting the energy and mols after rotations
		if args.CSEARCH=='summ' and len(rotmatches) != 0:
			status = dihedral_filter_and_sdf(name,args,log,dup_data,dup_data_idx,coord_Map, alg_Map, mol_template)

			# removes the rdkit file
			os.remove(name+'_'+'rdkit'+args.output)

	return status,update_to_rdkit
