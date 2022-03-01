######################################################.
#        This file stores functions related          #
#               to the QCORR module                  #
######################################################.

import os
import glob
import pandas as pd
import json
from pathlib import Path
from aqme.utils import (
	move_file,
)
import numpy as np

# Bondi VDW radii in Angstrom
bondi = {"H": 1.09,"He": 1.40,"Li": 1.81,"Be": 1.53,"B": 1.92,"C": 1.70,"N": 1.55,"O": 1.52,"F": 1.47,
"Ne": 1.54,"Na": 2.27,"Mg": 1.73,"Al": 1.84,"Si": 2.10,"P": 1.80,"S": 1.80,"Cl": 1.75,"Ar": 1.88,"K": 2.75,
"Ca": 2.31,"Ni": 1.63,"Cu": 1.40,"Zn": 1.39,"Ga": 1.87,"Ge": 2.11,"As": 1.85,"Se": 1.90,"Br": 1.83,"Kr": 2.02,
"Rb": 3.03,"Sr": 2.49,"Pd": 1.63,"Ag": 1.72,"Cd": 1.58,"In": 1.93,"Sn": 2.17,"Sb": 2.06,"Te": 2.06,"I": 1.98,
"Xe": 2.16,"Cs": 3.43,"Ba": 2.68,"Pt": 1.72,"Au": 1.66,"Hg": 1.55,"Tl": 1.96,"Pb": 2.02,"Bi": 2.07,"Po": 1.97,
"At": 2.02,"Rn": 2.20,"Fr": 3.48,"Ra": 2.83,"U": 1.86}

# covalent radii in Angstrom (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
rcov = {"H": 0.32,"He": 0.46,"Li": 1.33,"Be": 1.02,"B": 0.85,"C": 0.75,"N": 0.71,"O": 0.63,"F": 0.64,
"Ne": 0.67,"Na": 1.55,"Mg": 1.39,"Al": 1.26,"Si": 1.16,"P": 1.11,"S": 1.03,"Cl": 0.99,"Ar": 0.96,"K": 1.96,
"Ca": 1.71,"Sc": 1.48,"Ti": 1.36,"V": 1.34,"Cr": 1.22,"Mn": 1.19,"Fe": 1.16,"Co": 1.11,"Ni": 1.10,"Zn": 1.18,
"Ga": 1.24,"Ge": 1.21,"As": 1.21,"Se": 1.16,"Br": 1.14,"Kr": 1.17,"Rb": 2.10,"Sr": 1.85,"Y": 1.63,"Zr": 1.54,
"Nb": 1.47,"Mo": 1.38,"Tc": 1.28,"Ru": 1.25,"Rh": 1.25,"Pd": 1.20,"Ag": 1.28,"Cd": 1.36,"In": 1.42,"Sn": 1.40,
"Sb": 1.40,"Te": 1.36,"I": 1.33,"Xe": 1.31}


def detect_linear(errortype,atom_types,cclib_data):
	'''
	Check whether a linear molecule contain the right number of frequencies
	'''

	linear_options_3 = [['I', 'I', 'I'],['N', 'N', 'N'],['N', 'C', 'H'],['O', 'C', 'O'],['O', 'C', 'S'],['S', 'C', 'S'],['F','Be','F'],['F','Xe','F'],['O','C','N'],['S','C','N']]
	linear_options_4 = [['C', 'C', 'H', 'H']]

	if len(atom_types) == 3:
		for linear_3 in linear_options_3:
			if sorted(atom_types) == sorted(linear_3) and len(cclib_data['vibrations']['frequencies']) != 4:
				errortype = 'linear_mol_wrong'
				break
	elif len(atom_types) == 4:
		for linear_4 in linear_options_4:
			if sorted(atom_types) == sorted(linear_4) and len(cclib_data['vibrations']['frequencies']) != 7:
				errortype = 'linear_mol_wrong'
				break

	return errortype


def full_check(w_dir_main=os.getcwd(),destination_fullcheck='',json_files='*.json'):
	"""
	Checks that multiple calculations were done following the same protocols, including
	program and version, grid size, level of theory, dispersion and solvation model.

	Parameters
	----------
	w_dir_main : str
		Working directory
	destination_fullcheck : str
		Destination to create the file with the full check
	json_files : list of str
		json files to compare (glob.glob('*.json') and '*.json are both valid inputs to
		include all the json files from a folder)
	"""

	initial_dir = os.getcwd()
	w_dir_main = Path(w_dir_main)
	os.chdir(w_dir_main)

	if json_files == '*.json' or json_files == '*json':
		json_files=glob.glob('*.json')

	df_fullcheck = pd.DataFrame(columns=['file', 'program', 'grid_type', 'level_of_theory', 'dispersion', 'solvation'])

	for file in json_files:
		file_name = file.split('.')[0]
		with open(file) as json_file:
			cclib_data = json.load(json_file)

		program = cclib_data['metadata']['QM program']
		solvation = cclib_data['metadata']['solvation']
		dispersion = cclib_data['metadata']['dispersion model']
		grid_type = cclib_data['metadata']['grid type']
		functional = cclib_data['metadata']['functional']
		bs = cclib_data['metadata']['basis set']
		if functional != '' or bs != '':
			level_of_theory = '/'.join([functional, bs])
		else:
			level_of_theory = ''
		# designed to detect G4 calcs
		if level_of_theory == 'HF/GFHFB2':
			level_of_theory = 'G4'
		df_fullcheck.loc[len(df_fullcheck.index)] = [file_name, program, grid_type, level_of_theory, dispersion, solvation]
	
	fullcheck_file = '--QCORR_Fullcheck_Analysis--.dat'
	fullcheck_txt = '-- Full check analysis --'

	for prop in df_fullcheck.columns:
		if prop != 'file':
			unique_props = df_fullcheck[prop].unique()
			if len(unique_props) > 1:
				fullcheck_txt += f'\nx  Different {prop} used in the calculations:'
				for unique_prop in unique_props:
					file_names = df_fullcheck["file"].loc[df_fullcheck[prop] == unique_prop]
					fullcheck_txt += f'\n     * {unique_prop} in:'
					for file_name in file_names:
						adapted_name = file_name.replace('/','\\').split("\\")[-1]
						fullcheck_txt += f'\n       - {adapted_name}'
			else:
				fullcheck_txt += f'\no  Same {prop} ({unique_props[0]}) used in all the calculations'
			
	fullcheck_analysis = open(fullcheck_file, 'w')
	print(fullcheck_txt)
	fullcheck_analysis.write(fullcheck_txt)
	fullcheck_analysis.close()

	if destination_fullcheck == '':
		destination_fullcheck = w_dir_main.joinpath('successful_QM_outputs/json_files/')
	else:
		destination_fullcheck = Path(destination_fullcheck)
	move_file(destination_fullcheck, w_dir_main, fullcheck_file)
	
	os.chdir(initial_dir)


def check_isomerization(isom_data,file):
	"""
	Inputs two molecules with the atoms in the same order and checks if any bond
	is too different between them.

	Bonds are considered when the distance between two atoms is smaller than
	either the sum of their adjusted VDW radii or covalent radii (dist = n*R1 + n*R2, where
	n is a user defined parameter).

	Bonds forming part of TSs are removed.

	Parameters
	----------
	isom_data : dict
		Contains data related to coordinates, atoms and connectivity for input and output files
	file : string
		Filename (with extension)

	Returns
	-------
	isomerized : bool
		True if there is a clearly distorted bond within the geometries
	"""

	# load connectivity matrix from the starting points and convert string into matrix
	if not isom_data['Initial csv'].empty:
		filename = file.replace("_" + file.split("_")[-1], "")
		init_connectivity_string = isom_data['Initial csv'][isom_data['Initial csv']["code_name"] == filename]["initial_connectiv"][0]
		init_connectivity = json.loads(init_connectivity_string.replace(".", ",").replace(",]", "],").replace("],]", "]]"))
		isom_data['Atoms input'] = init_connectivity[0]

	else:
		init_connectivity = gen_connectivity(isom_data,isom_data['Atoms input'],isom_data['Coords input'])

	# in case the systems are not the same
	if len(isom_data['Atoms output']) != len(isom_data['Atoms input']):
		isomerized = True
	else:
		final_connectivity = gen_connectivity(isom_data,isom_data['Atoms output'],isom_data['Coords output'])

		# check connectivity differences from initial structure
		diff = final_connectivity - init_connectivity

		# remove bonds involved in TSs from connectivity matrixes
		if not isom_data['Initial csv'].empty:
			if "TS_atom_idx" in isom_data['Initial csv'].columns:
				ts_atoms = isom_data['Initial csv'][isom_data['Initial csv']["code_name"] == filename]["TS_atom_idx"][0].split(",")
				for i, ts_idx in enumerate(ts_atoms):
					for j, ts_idx_2 in enumerate(ts_atoms):
						if j > i:
							diff[int(ts_idx)][int(ts_idx_2)] = 0
							diff[int(ts_idx_2)][int(ts_idx)] = 0

		isomerized = np.any(diff)

	return isomerized


def gen_connectivity(isom_data,atom_types_conn,COORDINATES_conn):
	"""
	Use VDW radii to infer a connectivity matrix
	"""

	conn_mat = np.zeros((len(atom_types_conn), len(atom_types_conn)))
	for i, elem_i in enumerate(atom_types_conn):
		for j, elem_j in enumerate(atom_types_conn):
			if j > i:
				vdw_ij = bondi[elem_i] + bondi[elem_j]
				rcov_ij = rcov[elem_i] + rcov[elem_j]
				dist_ij = np.linalg.norm(np.array(COORDINATES_conn[i]) - np.array(COORDINATES_conn[j]))
				if dist_ij / vdw_ij < isom_data['VdW radii fraction'] or dist_ij / rcov_ij < isom_data['Covalent radii fraction']:
					conn_mat[i][j] = 1
				else:
					pass

	return conn_mat