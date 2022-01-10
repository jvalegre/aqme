#!/usr/bin/env python

#########################################################.
# 		   This file stores all the functions 		    #
# 	       used in writing SDF and COM files,           #
#              as well as the logger and                #
#                 yaml file importer		            #
#########################################################.

import subprocess
import glob

from pathlib import Path
import pandas as pd
from aqme.utils import BasisSet

from .templates import TurbomoleTemplate,OrcaTemplate,qprep

try:
	import pybel
except ImportError:
	from openbabel import pybel # for openbabel>=3.0.0



# Aux Functions for QM input generation	
def get_molecule_list(filepath,lowest_only=False,
						lowest_n=False, energy_threshold=0.0):
	out_molecules = []
	
	molecules = [mol for mol in pybel.readfile('sdf',filepath)]
	energies = [mol.energy for mol in molecules]
	min_energy = energies[0]
	for mol,energy in zip(molecules,energies):
		is_in_threshold = energy - min_energy < energy_threshold
		title,i = mol.title.strip().rsplit(maxsplit=1) 
		mol.title = f"{title} {int(i):03d}"
		if lowest_n and is_in_threshold:
			out_molecules.append(mol)
		elif lowest_n:
			break
		elif lowest_only:
			out_molecules.append(mol)
			break
		else:
			out_molecules.append(mol)

	return out_molecules

def get_basisset_list(args,option='opt'): 
	if option=='opt':
		if args.basis_sets: 
			return [BasisSet({'all':bs}) for bs in args.basis_sets]
		elif args.basisfiles: 
			return [BasisSet.from_file(file) for file in args.basisfiles]
	elif option=='sp': 
		if args.basis_sets_sp: 
			return [BasisSet({'all':bs}) for bs in args.basis_sets_sp]
		elif args.basisfiles_sp: 
			return [BasisSet.from_file(file) for file in args.basisfiles_sp]
	return []

# MAIN FUNCTION TO CREATE QM jobs
def write_qm_input_files(destination, template, molecules, charge_data, mult='auto'):
	qm_files = []
	for mol in molecules:

		molname = mol.title.strip().rsplit(maxsplit=1)

		# Get its charge and set it
		found_charges = charge_data.query(f'Molecule=="{molname}"')['Overall charge'].tolist()
		if not found_charges: # not in the dataframe
			charge = mol.data['Real charge']
		else:
			charge = found_charges[0]
		
		if charge == 'Invalid':
			continue
		
		mol.OBMol.SetTotalCharge(int(charge))
		# Set multiplicity
		if mult != 'auto':
			mol.OBMol.SetTotalSpinMultiplicity(int(mult))

		newfile = template.write(destination,mol)
		
		qm_files.append(newfile)

def load_charge_data(filepath,backup_files):
	#read in dup_data to get the overall charge of MOLECULES
	invalid_files = []
	try:
		charge_data = pd.read_csv(filepath, usecols=['Molecule','Overall charge'])
	except:
		charge_data = pd.DataFrame()
		for i,sdf_file in enumerate(backup_files):
			if not(Path(sdf_file).exists()):
				invalid_files.append(sdf_file)
				maxsplit = 1
				if 'filter' in sdf_file: 
					maxsplit += 1
				name = sdf_file.rsplit('_',maxsplit[0])
				charge = 'Invalid'
			else:
				mol = next(pybel.readfile(sdf_file))
				name = mol.title.split(maxsplit=1)[0]
				charge = mol.data['Real charge']
			charge_data.at[i,'Molecule'] = name
			charge_data.at[i,'Overall charge'] = charge
	return charge_data,invalid_files

# MAIN QPREP FUNCTION
def main(w_dir_initial,args,log):

	if len(args.geom_rules) >= 1:
		conf_files =  glob.glob('*_rules.sdf')
	# define the SDF files to convert to COM Gaussian files
	elif args.CMIN == 'xtb': 
		conf_files =  glob.glob('*_xtb.sdf')
	elif args.CMIN=='ani':
		conf_files =  glob.glob('*_ani.sdf')
	elif args.CSEARCH=='rdkit':
		conf_files =  glob.glob('*_rdkit.sdf')
	elif args.CSEARCH=='summ':
		conf_files =  glob.glob('*_summ.sdf')
	elif args.CSEARCH=='fullmonte':
		conf_files =  glob.glob('*_fullmonte.sdf')
	else:
		conf_files =  glob.glob('*.sdf')

	if args.com_from_xyz:
		xyz_files =  glob.glob('*.xyz')
		for file in xyz_files:
			mol = next(pybel.readfile('xyz',file))
			stem = Path(file).stem
			mol.write('sdf',f'{stem}.sdf')
		conf_files =  glob.glob('*.sdf')

	if not conf_files: 
		log.write('\nx  No SDF files detected to convert to gaussian COM files')
		return 

	# names for directories created
	if args.QPREP == 'gaussian':
		qm_folder = Path(f'{w_dir_initial}/QMCALC/G16')
		template = qprep.from_args(args)
	elif args.QPREP == 'orca':
		qm_folder = Path(f'{w_dir_initial}/QMCALC/ORCA')
		template = OrcaTemplate.from_args(args)
	elif args.QPREP == 'turbomole':
		qm_folder = Path(f'{w_dir_initial}/QMCALC/TURBOMOLE')
		template = TurbomoleTemplate.from_args(args)

	csv_name = args.input.split('.')[0]
	csv_file = f"{w_dir_initial}/CSEARCH/csv_files/{csv_name}-CSEARCH-Data.csv"
	charge_data, invalid_files = load_charge_data(csv_file,conf_files)

	# remove the invalid files and non-existing files
	accept_file = lambda x: x not in invalid_files and Path(x).exists()
	conf_files = [file for file in conf_files if accept_file(file) ]
	
	# Prepare the list of molecules that are to be written 
	molecules = []
	for file in conf_files:
		filepath = f'{w_dir_initial}/{file}'
		new_mols = get_molecule_list(filepath,
							lowest_only=args.lowest_only,lowest_n=args.lowest_n,
							energy_threshold=args.energy_threshold_for_gaussian)
		molecules.extend(new_mols)
	
	basisset_list = get_basisset_list(args,'opt')

	# Update each basis set to include the elements of all the molecules
	for bs in basisset_list: 
		bs.update_atomtypes(molecules)

	# Ensure a matching length of theory and basis set
	if len(args.level_of_theory) == 1:
		items = [(args.level_of_theory[0],bs) for bs in basisset_list]
	elif len(basisset_list) == 1:
		items = [(func,basisset_list[0]) for func in args.level_of_theory]
	elif len(basisset_list) == len(args.level_of_theory): 
		items = [(func,bs) for func,bs in zip(args.level_of_theory,basisset_list)]
	else:
		raise ValueError('Inconsistent size of basis sets and theory level')
	
	for functional,basisset in items:
		folder = qm_folder/f'{functional}-{basisset.name}'
		log.write(f"\no  Preparing QM input files in {folder}")

		template.functional = functional
		template.basisset = basisset

		# this variable keeps track of folder creation
		folder.mkdir(parents=True,exist_ok=True)

		# writing the com files
		# check conf_file exists, parse energies and then write DFT input
		qm_files = write_qm_input_files(folder, template, molecules,
										charge_data, mult=args.mult)
		
		# submitting the input file on a HPC
		if args.qsub:
			for qm_file in qm_files: 
				cmd_qsub = [args.submission_command, qm_file]
				subprocess.call(cmd_qsub)
