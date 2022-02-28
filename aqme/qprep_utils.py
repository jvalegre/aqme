######################################################.
#        This file stores functions related          #
#               to the QPREP module                  #
######################################################.

import json
from rdkit.Chem import Mol
from aqme.utils import (
	cclib_atoms_coords,
	QM_coords,
	read_file)


def qprep_coords(w_dir_main,file,found_coords,qprep_self):
	'''
	Retrieve atom types and coordinates from multiple formats (LOG, OUT, JSON, MOL)
	'''
	if isinstance(file, Mol):
		atom_types = [atom.GetSymbol() for _, atom in enumerate(file.GetAtoms())]
		cartesians = file.GetConformers()[0].GetPositions()

	elif file.split('.')[1] in ['log','out']:
		# detect QM program and number of atoms
		outlines = read_file(w_dir_main,file)
		n_atoms = 0
		resume_line = 0

		for i in range(0,len(outlines)):
			if outlines[i].find('Gaussian, Inc.'):
				program = 'gaussian'
				resume_line = i
				break
			elif outlines[i].find('O   R   C   A'):
				program = 'orca'
				resume_line = i
				break

		for i in range(resume_line,len(outlines)):	
			if program == 'gaussian':
				# get charge and mult
				if outlines[i].find("Charge = ") > -1:
					charge = int(outlines[i].split()[2])
					mult = int(outlines[i].split()[5].rstrip('\n'))
				# get number of atoms
				elif outlines[i].find('Symbolic Z-matrix:') > -1:
					symbol_z_line = i
				elif outlines[i].find('GradGrad') > -1:
					gradgrad_line = i
					n_atoms = gradgrad_line-symbol_z_line-4
					break

		atom_types,cartesians = QM_coords(outlines,0,n_atoms,program)

	elif file.split('.')[1] == 'json':
		with open(file) as json_file:
			cclib_data = json.load(json_file)
		try:
			atom_types,cartesians = cclib_atoms_coords(cclib_data)
			charge = cclib_data['properties']['charge']
			mult = cclib_data['properties']['multiplicity']
		except (AttributeError,KeyError):
			print(f'x  {file} does not contain coordinates and/or atom type information')
			atom_types,cartesians = [],[]
			
	# overwrite user-defined charge and multiplicity (if any)
	if qprep_self.args.charge == None:
		qprep_self.charge = charge
	else:
		qprep_self.charge = qprep_self.args.charge
	if qprep_self.args.mult == None:
		qprep_self.mult = mult
	else:
		qprep_self.mult = qprep_self.args.mult
	if atom_types != [] and cartesians != []:
		found_coords = True

	return atom_types,cartesians,found_coords,qprep_self