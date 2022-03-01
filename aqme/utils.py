######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import sys
import getopt
import subprocess
import numpy as np
import glob
from pathlib import Path
from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem.rdmolops import RemoveHs
from openbabel import pybel
import os
from rdkit import Geometry
from rdkit.Chem import rdmolfiles
import yaml
import pandas as pd
from rdkit.Chem import AllChem as Chem
from aqme.argument_parser import set_options,var_dict


def periodic_table():
	items = """X
			H                                                                                                  He
			Li Be  B                                                                             C   N   O   F  Ne
			Na Mg Al                                                                            Si   P   S  Cl  Ar
			K Ca Sc                                           Ti  V Cr Mn Fe Co Ni Cu  Zn  Ga  Ge  As  Se  Br  Kr
			Rb Sr  Y                                           Zr Nb Mo Tc Ru Rh Pd Ag  Cd  In  Sn  Sb  Te   I  Xe
			Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta  W Re Os Ir Pt Au  Hg  Tl  Pb  Bi  Po  At  Rn
			Fr Ra Ac Th Pa  U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo
			"""
	periodic_table = items.replace("\n", " ").strip().split()
	periodic_table[0] = ""

	return periodic_table


# load paramters from yaml file
def load_from_yaml(args_, log):
	"""
	Loads the parameters for the calculation from a yaml if specified. Otherwise
	does nothing.

	Parameters
	----------
	args : argparse.args
		Dataclass
	log : Logger
		Where to log the program progress
	"""
	# Variables will be updated from YAML file
	try:
		if args_.varfile is not None:
			if os.path.exists(args_.varfile):
				if os.path.splitext(args_.varfile)[1] == ".yaml":
					log.write("\no  Importing AQME parameters from " + args_.varfile)
					with open(args_.varfile, "r") as file:
						try:
							param_list = yaml.load(file, Loader=yaml.SafeLoader)
						except yaml.scanner.ScannerError:
							log.write(
								"\nx  Error while reading "
								+ args_.varfile
								+ ". Edit the yaml file and try again (i.e. use ':' instead of '=' to specify variables)"
							)
							sys.exit(
								"\nx  Error while reading "
								+ args_.varfile
								+ ". Edit the yaml file and try again (i.e. use ':' instead of '=' to specify variables)"
							)
			for param in param_list:
				if hasattr(args_, param):
					if getattr(args_, param) != param_list[param]:
						log.write(
							"o  RESET "
							+ param
							+ " from "
							+ str(getattr(args_, param))
							+ " to "
							+ str(param_list[param])
						)
						setattr(args_, param, param_list[param])
					else:
						log.write(
							"o  DEFAULT " + param + " : " + str(getattr(args_, param))
						)
	except UnboundLocalError:  # RAUL: Is it just me or this only happens when the file exists, and ens in .yaml and is empty or does not end in .yaml?
		log.write(
			"\no  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code.\n"
		)

	return args_, log


# class for logging
class Logger:
	"""
	Class that wraps a file object to abstract the logging.
	"""

	# Class Logger to writargs.input.split('.')[0] output to a file
	def __init__(self, filein, append, suffix="dat"):
		self.log = open(f"{filein}_{append}.{suffix}", "w")

	def write(self, message):
		"""
		Appends a newline character to the message and writes it into the file.

		Parameters
		----------
		message : str
		   text to be written in the log file.
		"""
		self.log.write(f"{message}\n")

	def fatal(self, message):
		"""
		Writes the message to the file. Closes the file and raises an error exit

		Parameters
		----------
		message : str
		   text to be written in the log file.
		"""
		self.write(message)
		self.finalize()
		raise SystemExit(1)

	def finalize(self):
		"""
		Closes the file
		"""
		self.log.close()

def get_mult(mol):
	radical_index = 0
	for _,atom in enumerate(mol.GetAtoms()):
		if atom.GetNumRadicalElectrons() != 0:
			radical_index +=1
	mult = 2*(radical_index/2) + 1
	return mult

# OS utils
def creation_of_dup_csv_csearch(csearch):

	"""
	Generates a pandas.DataFrame object with the appropiate columns for the
	conformational search and the minimization.

	Parameters
	----------
	csearch : str
		Conformational search method. Current valid methods are:
		['rdkit','fullmonte','summ']

	Returns
	-------
	pandas.DataFrame
	"""
	# Boolean aliases from args
	is_rdkit = csearch == "rdkit"
	is_fullmonte = csearch == "fullmonte"
	is_crest = csearch == "crest"
	is_summ = csearch == "summ"

	# column blocks definitions
	base_columns = [
		"Molecule",
		"RDKit-Initial-samples",
		"RDKit-energy-window",
		"RDKit-initial_energy_threshold",
		"RDKit-RMSD-and-energy-duplicates",
		"RDKit-Unique-conformers",
	]
	end_columns_no_min = ["CSEARCH time (seconds)", "Overall charge"]
	fullmonte_columns = [
		"FullMonte-Unique-conformers",
	]
	#'FullMonte-conformers',
	#'FullMonte-energy-window',
	#'FullMonte-initial_energy_threshold',
	#'FullMonte-RMSD-and-energy-duplicates']
	summ_columns = [
		"summ-conformers",
		"summ-energy-window",
		"summ-initial_energy_threshold",
		"summ-RMSD-and-energy-duplicates",
		"summ-Unique-conformers",
	]
	crest_columns = ["Molecule", "crest-conformers"]

	# Check Conformer Search method
	if is_rdkit:
		columns = base_columns
	elif is_fullmonte:
		columns = base_columns + fullmonte_columns
	elif is_summ:
		columns = base_columns + summ_columns
	elif is_crest:
		columns = crest_columns
	else:
		return None
	columns += end_columns_no_min
	return pd.DataFrame(columns=columns)


def creation_of_dup_csv_cmin(cmin):

	"""
	Generates a pandas.DataFrame object with the appropiate columns for the
	conformational search and the minimization.

	Parameters
	----------
	cmin : str
		Minimization method. Current valid methods are:
		['xtb','ani']

	Returns
	-------
	pandas.DataFrame
	"""
	# Boolean aliases from args
	is_xtb = cmin == "xtb"
	is_ani = cmin == "ani"

	# column blocks definitions

	xtb_columns = [
		"xTB-Initial-samples",
		"xTB-energy-window",
		"xTB-initial_energy_threshold",
		"xTB-RMSD-and-energy-duplicates",
		"xTB-Unique-conformers",
	]
	ANI_columns = [
		"ANI-Initial-samples",
		"ANI-energy-window",
		"ANI-initial_energy_threshold",
		"ANI-RMSD-and-energy-duplicates",
		"ANI-Unique-conformers",
	]
	end_columns = ["CMIN time (seconds)", "Overall charge"]

	# Check Minimization Method
	if is_ani:
		columns = ANI_columns
	if is_xtb:  # is_ani and is_xtb will not happen, but this is what was written
		columns = xtb_columns

	columns += end_columns
	return pd.DataFrame(columns=columns)


def move_file(destination, source, file):
	"""
	Moves files from the source folder to the destination folder and creates
	the destination folders when needed.

	Parameters
	----------
	destination : str
		path to the destination folder
	src : str
		path to the source folder
	file : str
		full name of the file (file + extension)
	"""

	destination.mkdir(exist_ok=True, parents=True)
	filepath = source / file
	try:
		filepath.rename(destination / file)
	except FileExistsError:
		filepath.replace(destination / file)


# openbabel utils
# com to xyz to sdf for obabel
def com_2_xyz_2_sdf(input_file, default_charge, start_point=None):
	"""
	com to xyz to sdf for obabel

	Parameters
	----------
	input_file : str
		path to the file to convert
	start_point : str, optional
		file(path/name?) to the starting point, by default None

	Returns
	-------
	int?
		charge or None?
	"""
	extension = Path(input_file).suffix

	if start_point is None:
		if extension in ["com", "gjf", "xyz"]:
			file = Path(input_file)

	else:
		file = Path(start_point)

	filename = Path.stem

	# Create the 'xyz' file and/or get the total charge
	if (
		extension != "xyz"
	):  #  RAUL: Originally this pointed towards args.input_file, shouldn't it be to args.file?
		xyz, charge = get_info_input(file)
		xyz_txt = "\n".join(xyz)
		with open(f"{filename}.xyz", "w") as F:
			F.write(f"{len(xyz)}\n{filename}\n{xyz_txt}\n")
	else:
		charge = default_charge

	xyz_2_sdf(f"{filename}.xyz")

	return charge


def xyz_2_sdf(file, parent_dir=None):
	"""
	Creates a .sdf file from a .xyz in the specified directory. If no directory
	is specified then the files are created in the current directory.

	Parameters
	----------
	file : str
		filename and extension of an existing .xyz file
	dir : str or pathlib.Path, optional
		a path to the directory where the .xyz file is located
	"""
	if parent_dir is None:
		parent_dir = Path("")
	else:
		parent_dir = Path(parent_dir)
	mol = next(pybel.readfile("xyz", parent_dir / file))
	ofile = Path(file).stem + ".sdf"
	mol.write("sdf", parent_dir / ofile)


def get_info_input(file):
	"""
	Takes an input file and retrieves the coordinates of the atoms and the
	total charge.

	Parameters
	----------
	file : str or pathlib.Path
		A path pointing to a valid .com or .gjf file

	Returns
	-------
	coordinates : list
		A list of strings (without \\n) that contain the xyz coordinates of the
		.gjf or .com file
	charge : str
		A str with the number corresponding to the total charge of the .com or
		.gjf file
	"""

	with open(file, "r") as input_file:
		input_lines = input_file.readlines()

	_iter = input_lines.__iter__()

	line = ''
	# input for Gaussian calculations
	if file.split('.')[1] in ['com','gjf']:

		# Find the command line
		while "#" not in line:
			line = next(_iter)

		# in case the keywords are distributed in multiple lines
		while len(line.split()) > 0:
			line = next(_iter)

		# pass the title lines
		_ = next(_iter)
		while line:
			line = next(_iter).strip()

		# Read charge and multiplicity
		charge, mult = next(_iter).strip().split()

		# Store the atom types and coordinates until next empty line.
		atoms_and_coords = []
		line = next(_iter).strip()
		while line:
			atoms_and_coords.append(line.strip())
			line = next(_iter).strip()

	# input for ORCA calculations
	if file.split('.')[1] == '.inp':

		# Find the line with charge and multiplicity
		while "* xyz" not in line or "* int" not in line:
			line = next(_iter)

		# Read charge and multiplicity
		charge = line.strip().split()[-2]

		# Store the coordinates until next *
		atoms_and_coords = []
		line = next(_iter).strip()
		while len(line.split()) > 1:
			atoms_and_coords.append(line.strip())
			line = next(_iter).strip()

	return atoms_and_coords, charge


# RDKit Utils
def nci_ts_mol(smi, args, constraints_dist, constraints_angle, constraints_dihedral, name):

	if constraints_dist is not None:
		constraints_dist = [[float(y) for y in x] for x in constraints_dist]
		constraints_dist = np.array(constraints_dist)
	if constraints_angle is not None:
		constraints_angle = [[float(y) for y in x] for x in constraints_angle]
		constraints_angle = np.array(constraints_angle)
	if constraints_dihedral is not None:
		constraints_dihedral = [[float(y) for y in x] for x in constraints_dihedral]
		constraints_dihedral = np.array(constraints_dihedral)

	molsH = []
	mols = []
	for m in smi:
		mols.append(Chem.MolFromSmiles(m))
		molsH.append(Chem.AddHs(Chem.MolFromSmiles(m)))

	for m in molsH:
		Chem.EmbedMultipleConfs(m,numConfs=1)
	for m in mols:
		Chem.EmbedMultipleConfs(m,numConfs=1)

	coord = [0.0,0.0,5.0]
	molH = molsH[0]
	for fragment in molsH[1:]:
		offset_3d = Geometry.Point3D(coord[0],coord[1],coord[2])
		molH = Chem.CombineMols(molH, fragment, offset_3d)
		coord[1] += 5
		Chem.SanitizeMol(molH)

	coord = [0.0,0.0,5.0]
	mol = mols[0]
	for fragment in mols[1:]:
		offset_3d = Geometry.Point3D(coord[0],coord[1],coord[2])
		mol = Chem.CombineMols(mol, fragment, offset_3d)
		coord[1] += 5
	mol = Chem.AddHs(mol)
	Chem.SanitizeMol(mol)

	atom_map = []
	for atom in mol.GetAtoms():
		atom_map.append(atom.GetAtomMapNum())

	max_map = max(atom_map)
	for a in mol.GetAtoms():
		if a.GetSymbol() =='H':
			max_map +=1
			a.SetAtomMapNum(int(max_map))

	Chem.ConstrainedEmbed(mol,molH)
	rdmolfiles.MolToXYZFile(mol, name+'.xyz')

	# for atom in mol.GetAtoms():
		# print(atom.GetAtomMapNum(),atom.GetIdx(),atom.GetSymbol())
	nconstraints_dist = []
	if constraints_dist is not None:

		for _,r in enumerate(constraints_dist):
			# print('r:',r[:2])
			nr = []
			for _,ele in enumerate(r[:2]):
				# print('ele:',  ele)
				for atom in mol.GetAtoms():
					if ele == atom.GetAtomMapNum():
						nr.append(float(atom.GetIdx())+1)
						# print('nr:',  nr)
			nr.append(r[-1])
			# print('nr-final:',  nr)
			nconstraints_dist.append(nr)
		nconstraints_dist  = np.array(nconstraints_dist)
		# print(nconstraints_dist)
		# print('----')
	nconstraints_angle = []
	if constraints_angle is not None:

		for _,r in enumerate(constraints_angle):
			# print('r:',r[:2])
			nr = []
			for _,ele in enumerate(r[:3]):
				# print('ele:',  ele)
				for atom in mol.GetAtoms():
					if ele == atom.GetAtomMapNum():
						nr.append(float(atom.GetIdx())+1)
						# print('nr:',  nr)
			nr.append(r[-1])
			# print('nr-final:',  nr)
			nconstraints_angle.append(nr)
		nconstraints_angle  = np.array(nconstraints_angle)
		# print(nconstraints_angle)
		# print('----')
	nconstraints_dihedral = []
	if constraints_dihedral is not None:
		for _,r in enumerate(constraints_dihedral):
			# print('r:',r[:2])
			nr = []
			for _,ele in enumerate(r[:4]):
				# print('ele:',  ele)
				for atom in mol.GetAtoms():
					if ele == atom.GetAtomMapNum():
						nr.append(float(atom.GetIdx())+1)
						# print('nr:',  nr)
			nr.append(r[-1])
			# print('nr-final:',  nr)
			nconstraints_dihedral.append(nr)
		nconstraints_dihedral  = np.array(nconstraints_dihedral)
		# print(nconstraints_dihedral)
		# print('----')

	# print(constraints_dist[:2])
	# for atom in mol.GetAtoms():
	#     print(atom.GetAtomMapNum(),atom.GetIdx(),atom.GetSymbol())
	#     if constraints_dist is not None:
	#         constraints_dist = np.where(constraints_dist==float(atom.GetAtomMapNum()),float(atom.GetIdx())+1,constraints_dist)
	#     if constraints_angle is not None:
	#         constraints_angle = np.where(constraints_angle==float(atom.GetAtomMapNum()),float(atom.GetIdx())+1,constraints_angle)
	#     if constraints_dihedral is not None:
	#         constraints_dihedral = np.where(constraints_dihedral==float(atom.GetAtomMapNum()),float(atom.GetIdx())+1,constraints_dihedral)
	# rdmolfiles.MolToXYZFile(mol, name+'.xyz')
	# print(nconstraints_dist, nconstraints_angle, nconstraints_dihedral)
	return mol, nconstraints_dist, nconstraints_angle, nconstraints_dihedral


def rules_get_charge(mol, args):
	"""
	Automatically sets the charge for metal complexes

	Parameters
	----------
	mol : [type]
		[description]
	args : [type]
		[description]

	Returns
	-------
	charge : int
		Charge of the system
	"""

	C_group = ["C", "Se", "Ge"]
	N_group = ["N", "P", "As"]
	O_group = ["O", "S", "Se"]
	F_group = ["Cl", "Br", "I"]

	M_ligands, N_carbenes, bridge_atoms, neighbours = [], [], [], []
	charge_rules = np.zeros(len(args.metal_idx), dtype=int)
	neighbours = []
	for atom in mol.GetAtoms():
		# get the neighbours of metal atom and calculate the charge of metal center + ligands
		if atom.GetIdx() in args.metal_idx:
			charge_idx = args.metal_idx.index(atom.GetIdx())
			neighbours = atom.GetNeighbors()
			charge_rules[charge_idx] = args.m_oxi[charge_idx]
			for neighbour in neighbours:
				M_ligands.append(neighbour.GetIdx())
				if neighbour.GetTotalValence() == 4:
					if neighbour.GetSymbol() in C_group:
						carbene_like = False
						bridge_ligand = False
						for inside_neighbour in neighbour.GetNeighbors():
							if inside_neighbour.GetSymbol() in N_group:
								if inside_neighbour.GetTotalValence() == 4:
									for N_neighbour in inside_neighbour.GetNeighbors():
										# this option detects bridge ligands that connect two metals such as M--CN--M
										# we use I since the M is still represented as I at this point
										if N_neighbour.GetSymbol() == "I":
											bridge_ligand = True
											bridge_atoms.append(
												inside_neighbour.GetIdx()
											)
									if not bridge_ligand:
										carbene_like = True
										N_carbenes.append(inside_neighbour.GetIdx())
						if not carbene_like:
							charge_rules[charge_idx] = charge_rules[charge_idx] - 1
				elif neighbour.GetTotalValence() == 3:
					if neighbour.GetSymbol() in N_group:
						charge_rules[charge_idx] = charge_rules[charge_idx] - 1
				elif neighbour.GetTotalValence() == 2:
					if neighbour.GetSymbol() in O_group:
						nitrone_like = False
						for inside_neighbour in neighbour.GetNeighbors():
							if inside_neighbour.GetSymbol() in N_group:
								nitrone_like = True
						if not nitrone_like:
							charge_rules[charge_idx] = charge_rules[charge_idx] - 1

				elif neighbour.GetTotalValence() == 1:
					if neighbour.GetSymbol() in F_group:
						charge_rules[charge_idx] = charge_rules[charge_idx] - 1

	# recognizes charged N and O atoms in metal ligands (added to the first metal of the list as default)
	# this group contains atoms that do not count as separate charge groups (i.e. N from Py ligands)
	if len(neighbours) > 0:
		invalid_charged_atoms = M_ligands + N_carbenes + bridge_atoms
		for atom in mol.GetAtoms():
			if atom.GetIdx() not in invalid_charged_atoms:
				if atom.GetSymbol() in N_group:
					if atom.GetTotalValence() == 4:
						charge_rules[0] = charge_rules[0] + 1
				if atom.GetSymbol() in O_group:
					if atom.GetTotalValence() == 1:
						charge_rules[0] = charge_rules[0] - 1
		return charge_rules
		
	# no update in charge as it is an organic molecule
	return [args.charge_default]


def substituted_mol(mol, args):
	"""
	Returns a molecule object in which all metal atoms specified in args.metal
	are replaced by Iodine and the charge is set depending on the number of
	neighbors.

	Parameters
	----------
	mol : rdkit.Chem.rdchem.Mol
		A molecule object whose metal is to be substituted
	args : argparse.args
		[description]

	Returns
	-------
	tuple
		mol,args.metal_idx, args.complex_coord, args.metal_sym
	"""
	Neighbors2FormalCharge = dict()
	for i, j in zip(range(2, 9), range(-3, 4)):
		Neighbors2FormalCharge[i] = j

	for atom in mol.GetAtoms():
		symbol = atom.GetSymbol()
		if symbol in args.metal:
			args.metal_sym[args.metal.index(symbol)] = symbol
			args.metal_idx[args.metal.index(symbol)] = atom.GetIdx()
			args.complex_coord[args.metal.index(symbol)] = len(atom.GetNeighbors())
			atom.SetAtomicNum(53)
			n_neighbors = len(atom.GetNeighbors())
			if n_neighbors > 1:
				formal_charge = Neighbors2FormalCharge[n_neighbors]
				atom.SetFormalCharge(formal_charge)

	return mol, args.metal_idx, args.complex_coord, args.metal_sym


# DETECTS DIHEDRALS IN THE MOLECULE
def getDihedralMatches(mol, heavy, log):
	# this is rdkit's "strict" pattern
	pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
	qmol = Chem.MolFromSmarts(pattern)
	matches = mol.GetSubstructMatches(qmol)

	# these are all sets of 4 atoms, uniquify by middle two
	uniqmatches = []
	seen = set()
	for (a, b, c, d) in matches:
		if (b, c) not in seen and (c, b) not in seen:
			if heavy:
				if (
					mol.GetAtomWithIdx(a).GetSymbol() != "H"
					and mol.GetAtomWithIdx(d).GetSymbol() != "H"
				):
					seen.add((b, c))
					uniqmatches.append((a, b, c, d))
			if not heavy:
				if (
					mol.GetAtomWithIdx(c).GetSymbol() == "C"
					and mol.GetAtomWithIdx(d).GetSymbol() == "H"
				):
					pass
				else:
					seen.add((b, c))
					uniqmatches.append((a, b, c, d))
	return uniqmatches


def getDihedralMatches_v2(mol, heavy, log):  # New version using openbabel

	AtomInTripleBond = "$(*#*)"
	TerminalAtom = "D1"
	CF3 = "$(C(F)(F)F)"
	CCl3 = "$(C(Cl)(Cl)Cl)"
	CBr3 = "$(C(Br)(Br)Br)"
	tBut = "$(C([CH3])([CH3])[CH3])"
	# A 3-bonded C with a double bond to (N, O or S)
	# singlgy bonded to not ring bonded to a non-terminal N,O or S.
	CD3_1d = "$([CD3](=[N,O,S])-!@[#7,O,S!D1])"
	CD3_1r = "$([#7,O,S!D1]-!@[CD3]=[N,O,S])"  # Backwards version
	# A 3-bonded C with a double bond to (N+)
	# singlgy bonded to not ring bonded to Any non-terminal N
	CD3_2d = "$([CD3](=[N+])-!@[#7!D1])"
	CD3_2r = "$([#7!D1]-!@[CD3]=[N+])"  # Backwards version
	Atom1 = "*"
	Atom2 = f"!{AtomInTripleBond}&!{TerminalAtom}"
	Atom2 += f"&!{CF3}&!{CCl3}&!{CBr3}&!{tBut}"
	Atom2 += f"&!{CD3_1d}&!{CD3_1r}&!{CD3_2d}&!{CD3_2r}"
	Atom3 = f"!{AtomInTripleBond}&!{TerminalAtom}"
	Atom3 += f"&!{CF3}&!{CCl3}&!{CBr3}&!{tBut}"
	Atom4 = "*"
	pattern = f"{Atom1}~[{Atom2}]-!@[{Atom3}]~{Atom4}"
	smarts = pybel.Smarts(pattern)
	matches = smarts.findall(mol)

	# these are all sets of 4 atoms, uniquify by middle two
	H_atoms = set(pybel.Smarts("[#1]").findall(mol))
	C_atoms = set(pybel.Smarts("[#6]").findall(mol))
	uniqmatches = []
	seen = set()
	for (a, b, c, d) in matches:
		if (b, c) not in seen and (c, b) not in seen:
			if heavy:
				if a not in H_atoms and d not in H_atoms:
					seen.add((b, c))
					uniqmatches.append((a, b, c, d))
			if not heavy:
				# So what if a == 'H' and b == 'C'? is that valid ¿?
				if c not in C_atoms or d not in H_atoms:
					seen.add((b, c))
					uniqmatches.append((a, b, c, d))


# checks for salts
def check_for_pieces(smi):
	# taking largest component for salts
	pieces = smi.split(".")
	# if len(pieces) > 1:
	#     # take largest component by length
	#     smi = max(pieces, key=len)
	return pieces


def mol_from_sdf_or_mol_or_mol2(input_file):
	"""
	mol from sdf

	Parameters
	----------
	input_file : str
		path to a .sdf .mol or .mol2 file

	Returns
	-------
	tuple of lists?
		suppl, IDs, charges
	"""
	filename = os.path.splitext(input_file)[0]
	extension = os.path.splitext(input_file)[1]

	if extension == ".sdf":
		suppl = Chem.SDMolSupplier(input_file, removeHs=False)
	elif extension == ".mol":
		suppl = Chem.MolFromMolFile(input_file, removeHs=False)
	elif extension == ".mol2":
		suppl = Chem.MolFromMol2File(input_file, removeHs=False)

	IDs, charges = [], []

	with open(input_file, "r") as F:
		lines = F.readlines()

	molecule_count = 0
	for i, line in enumerate(lines):
		if line.find(">  <ID>") > -1:
			ID = lines[i + 1].split()[0]
			IDs.append(ID)
		if line.find("M  CHG") > -1:
			charge_line = line.split("  ")
			charge = 0
			for j in range(4, len(charge_line)):
				if (j % 2) == 0:
					if j == len(charge_line) - 1:
						charge_line[j] = charge_line[j].split("\n")[0]
					charge += int(charge_line[j])
			charges.append(charge)
		if line.find("$$$$") > -1:
			molecule_count += 1
			if molecule_count != len(charges):
				charges.append(0)

	if len(IDs) == 0:
		if extension == ".sdf":
			for i in range(len(suppl)):
				IDs.append(f"{filename}_{i}")
		else:
			IDs.append(filename)
	if len(charges) == 0:
		if extension == ".sdf":
			for _ in suppl:
				charges.append(0)
		else:
			charges.append(0)
	return suppl, IDs, charges


# checks the charge on the smi string
def check_charge_smi(smiles,ts):
	for smi in smiles:
		charge = 0
		for i, smi_letter in enumerate(smi):
			if smi_letter == "+":
				if ts:
					if smi[i + 3] == "]":
						charge += 1
					else:
						charge += int(smi[i + 1])
				else:
					if smi[i + 1] == "]":
						charge += 1
					else:
						charge += int(smi[i + 1])

			elif smi_letter == "-":
				if ts:
					if smi[i + 3] == "]":
						charge -= 1
					else:
						charge -= int(smi[i + 1])
				else:
					if smi[i + 1] == "]":
						charge -= 1
					else:
						charge -= int(smi[i + 1])
	return charge


def set_metal_atomic_number(mol, metal_idx, metal_sym):
	"""
	Changes the atomic number of the metal atoms using their indices.

	Parameters
	----------
	mol : rdkit.Chem.Mol
		rdkit molecule object
	metal_idx : list
		sorted list that contains the indices of the metal atoms in the molecule
	metal_sym : list
		sorted list (same order as metal_idx) that contains the symbols of the
		metals in the molecule.
	"""

	for atom in mol.GetAtoms():
		if atom.GetIdx() in metal_idx:
			re_symbol = metal_sym[metal_idx.index(atom.GetIdx())]
			atomic_number = periodic_table().index(re_symbol)
			atom.SetAtomicNum(atomic_number)


def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_RMSD):
	"""
	Takes in two rdkit.Chem.Mol objects and calculates the RMSD between them.
	(As side efect mol1 is left in the aligned state, if heavy is specified
	the side efect will not happen)

	Parameters
	----------
	mol1 : rdkit.Chem.Mol
		Probe molecule
	mol2 : rdkit.Chem.Mol
		Target molecule. The probe is aligned to the target to compute the RMSD.
	c1 : int
		conformation of mol1 to use for the RMSD
	c2 : int
		conformation of mol2 to use for the RMSD
	heavy : bool
		If True it will ignore the H atoms when computing the RMSD.
	max_matches_RMSD : int
		the max number of matches found in a SubstructMatch()

	Returns
	-------
	float
		Returns the best RMSD found.
	"""
	if heavy:
		mol1 = RemoveHs(mol1)
		mol2 = RemoveHs(mol2)
	return GetBestRMS(mol1, mol2, c1, c2, maxMatches=max_matches_RMSD)


def get_filenames(file_type, name):
	"""
	Finding the file type to move for analysis
	"""
	
	files = []
	if file_type == "output":
		formats = ["*.log", "*.LOG", "*.out", "*.OUT", "*json"]
	elif file_type == "input":
		formats = ["*.com", "*.gjf"]
	for _, file_format in enumerate(formats):
		if name is None:
			all_files = enumerate(glob.glob(file_format))
		else:
			all_files = enumerate(glob.glob(name + file_format))
		for _, file in all_files:
			if file not in files:
				files.append(file)
	return files


def output_to_mol(file, file_format):
	"""
	Input an XYZ, LOG or OUT file from QM calculations and converts it into
	a mol object.

	Parameters
	----------
	file : string
		Filename
	file_format : string
		File file_format
	log : Logger object
		Writes data regarding the results from QCORR

	Returns
	-------
	mol
		Mol object
	"""

	# transforms output file into mol object
	# for input (from com to xyz to mol)
	if file_format == "xyz":
		cmd_obabel = [
			"obabel",
			"-ixyz",
			os.path.splitext(file)[0] + ".xyz",
			"-omol",
			"-O",
			os.path.splitext(file)[0] + ".mol",
		]
	# for output (from log to mol)
	if file_format in ["log", "out"]:
		cmd_obabel = [
			"obabel",
			"-ilog",
			os.path.splitext(file)[0] + "." + file_format,
			"-omol",
			"-O",
			os.path.splitext(file)[0] + ".mol",
		]
	subprocess.run(cmd_obabel)
	mol = Chem.MolFromMolFile(file.split(".")[0] + ".mol")

	return mol


def read_energies(file):
	"""
	Parses the energies from sdf files .
	"""
	energies = []
	f = open(file, "r")
	readlines = f.readlines()
	for i, _ in enumerate(readlines):
		if readlines[i].find(">  <Energy>") > -1:
			energies.append(float(readlines[i + 1].split()[0]))
	f.close()
	return energies


def get_name_and_charge(name, charge_data):
	"""
	Function to get name of charge from an SDF file.
	"""

	name_list = name.split("_")

	if "xtb" in name_list or "ani" in name_list:
		if "filter" in name_list:
			name_molecule = name[:-21]
		else:
			name_molecule = name[:-4]
	elif "summ" in name_list:
		if "filter" in name_list:
			name_molecule = name[:-22]
		else:
			name_molecule = name[:-5]
	elif "rdkit" in name_list:
		if "filter" in name_list:
			name_molecule = name[:-23]
		else:
			name_molecule = name[:-6]
	elif "fullmonte" in name_list:
		if "filter" in name_list:
			name_molecule = name[:-27]
		else:
			name_molecule = name[:-10]

	if charge_data is not None:
		for i in range(len(charge_data)):
			if charge_data.loc[i, "Molecule"] == name_molecule:
				charge_com = charge_data.loc[i, "Overall charge"]
			else:
				try:
					suppl = Chem.SDMolSupplier(name + ".sdf", removeHs=False)
				except OSError:
					suppl = False
				if suppl:
					mol = suppl[0]
					charge_com = mol.GetProp("Real charge")
				else:
					charge_com = "Invalid"

		return charge_com

	return name_molecule


def command_line_args():
	'''
	Load default and user-defined arguments specified through command lines. Arrguments are loaded as a dictionary
	'''
	
	# First, create dictionary with user-defined arguments
	kwargs = {}
	available_args = ['help']
	bool_args = [
		"verbose",
		"csearch",
		"cmin",
		"qprep",
		"qcorr",
		"qstat",
		"qpred",
		"metal_complex",
		"time",
		"heavyonly",
		"cregen",
		"lowest_only",
		"lowest_n",
		"chk",
		"com_from_xyz",
		"dup",
		"fullcheck",
		"rot_dihedral",
		"nmr_online",
		"qsub",
		"qsub_ana",
		"nci_complex",
		"ts_complex"]

	for arg in var_dict:
		if arg in bool_args:
			available_args.append(f'{arg}')
		else:
			available_args.append(f'{arg} =')

	try:
		opts,_ = getopt.getopt(sys.argv[1:], 'h', available_args)
	except getopt.GetoptError as err:
		print(err)
		sys.exit()
	
	for arg,value in opts:
		arg_name = arg.split('--')[1].strip()
		if arg_name in bool_args:
			value = True
		if value == 'None':
			value = None
		if arg_name in ("h", "help"):
			print('AQME is installed correctly! For more information about the available options, see the documentation in https://github.com/jvalegre/aqme')
			sys.exit()
		else:
			kwargs[arg_name] = value
	
	# Second, load all the default variables as an "add_option" object
	args = load_variables(kwargs,_)

	return args


def load_variables(kwargs,aqme_module):
	'''
	Load default and user-defined variables
	'''

	# first, load default values and options manually added to the function
	self = set_options(kwargs)
	
	self.initial_dir = Path(os.getcwd())
	self.w_dir_main = Path(self.w_dir_main)
	if self.isom is not None:
		self.isom_inputs = Path(self.args.isom_inputs)

	# go to working folder and detect files
	os.chdir(self.w_dir_main)
	if not isinstance(self.files, list):
		self.files = glob.glob(self.files)

	# start a log file to track the QCORR module
	error_setup = False
	if aqme_module == 'qcorr':
		# detects cycle of analysis (0 represents the starting point)
		self.round_num = check_run(self.w_dir_main)
		logger_1 = 'QCORR-run'
		logger_2 = f'{str(self.round_num)}'
		try:
			self.log = Logger(self.w_dir_main / logger_1,logger_2)
		except FileNotFoundError:
			print('x  The PATH specified as input in the w_dir_main option might be invalid!')
			error_setup = True

		if len(self.files) == 0 and not error_setup:
			print(f'x  There are no output files in {self.w_dir_main}.')
			self.log.write(f'x  There are no output files in {self.w_dir_main}.')
			self.log.finalize()
			error_setup = True

	if error_setup:
		# this is added to avoid path problems in jupyter notebooks
		os.chdir(self.initial_dir)
		sys.exit()

	# this part loads variables from yaml files (if varfile is used)
	elif self.varfile is not None:
		self.yaml, self.log = load_from_yaml(self.args, self.log)
		for key,value in self.yaml.iteritems():
			setattr(self, key, value)
	
	return self


def read_file(w_dir, file):
	"""
	Reads through a file and retrieves a list with all the lines.
	"""

	os.chdir(w_dir)
	outfile = open(file, "r")
	outlines = outfile.readlines()
	outfile.close()

	return outlines


def QM_coords(outlines,min_RMS,n_atoms,program):
	'''
	Retrieves atom types and coordinates from QM output files
	'''

	atom_types,cartesians = [],[]
	per_tab = periodic_table()
	count_RMS = -1

	if program == 'gaussian':
		for i,line in enumerate(outlines):
			if line.find('Standard orientation:') > -1:
				count_RMS += 1
			if count_RMS == min_RMS:
				range_lines = [i+5,i+5+n_atoms]
				break
		for i in range(range_lines[0],range_lines[1]):
			massno = int(outlines[i].split()[1])
			if massno < len(per_tab):
				atom_symbol = per_tab[massno]
			else:
				atom_symbol = "XX"
			atom_types.append(atom_symbol)
			cartesians.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

	return atom_types,cartesians


def cclib_atoms_coords(cclib_data):
	'''
	Function to convert atomic numbers and coordinate arrays from cclib into
	a format compatible with QPREP.
	'''
	atom_numbers = cclib_data['atoms']['elements']['number']
	atom_types = []
	per_tab = periodic_table()
	for atom_n in atom_numbers:
		if atom_n < len(per_tab):
			atom_symbol = per_tab[atom_n]
		else:
			atom_symbol = "XX"
		atom_types.append(atom_symbol)

	cartesians_array = cclib_data['atoms']['coords']['3d']
	cartesians = [cartesians_array[i:i + 3] for i in range(0, len(cartesians_array), 3)]

	return atom_types,cartesians


def check_run(w_dir):
	'''
	Determines the folder where input files are gonna be generated in QCORR.
	'''

	input_folder = w_dir.joinpath('unsuccessful_QM_outputs/')
	folder_count = 1

	if os.path.exists(input_folder):
		dir_list = os.listdir(input_folder)
		for folder in dir_list:
			if folder.find('run_') > -1:
				folder_count += 1

	return folder_count
