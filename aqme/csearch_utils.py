#####################################################.
#        This file stores all the functions         #
#                 used in CSEARCH                   #
#####################################################.
import os
import sys
import subprocess
from pathlib import Path
from pkg_resources import resource_filename
import pandas as pd

try:
	import pybel
except ImportError:
	from openbabel import pybel  # for openbabel>=3.0.0
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdDistGeom, rdMolAlign
from aqme.utils import get_info_input, get_conf_RMS

TEMPLATES_PATH = Path(resource_filename("aqme", "templates"))

# Auxiliar functions of this module
def load_template(complex_type, log):
	"""
	Checks if the templates are reachable and if so returns the molecule object.

	Returns
	-------
	rdkit.Chem.Mol
			The molecule file of the corresponding template.
	"""
	type2template = dict()
	type2template["squareplanar"] = "template-4-and-5.sdf"
	type2template["squarepyramidal"] = "template-4-and-5.sdf"
	type2template["linear"] = "template-2.sdf"
	type2template["trigonalplanar"] = "template-3.sdf"

	folder = TEMPLATES_PATH

	if not folder.exists():
		log.write(
			"x The templates folder was not found, probably due to a problem while installing AQME"
		)
		sys.exit()

	file_template = folder / Path(type2template[complex_type])
	templates = Chem.SDMolSupplier(str(file_template))
	template = templates[
		-1
	]  # RAUL: I'm assuming that there is only one molecule per template and in case there are several, it's the last one

	return template


def calc_neighbours(molecule, metals_idx):
	"""
	Changes the atomic number (and charge) of the first metal found
	and returns its neighbours, the number of neighbours and the idx of the
	metal.

	Parameters
	----------
	molecule : rdkit.Chem.Mol
			[description]
	metals_idx : list
			List containing the Idx of the metals. The first match is the only one
			considered.

	Returns
	-------
	list
			list of neighbour atoms
	"""
	bonds2AtNum = dict()
	bonds2AtNum[5] = 14
	bonds2AtNum[4] = 14
	bonds2AtNum[3] = 53
	bonds2AtNum[2] = 53

	for atom in molecule.GetAtoms():
		idx = atom.GetIdx()
		if idx in metals_idx:
			n_bonds = len(atom.GetBonds())
			AtNum = bonds2AtNum[n_bonds]
			atom.SetAtomicNum(AtNum)
			if n_bonds == 5:
				atom.SetFormalCharge(1)
			neighbours = atom.GetNeighbors()
			return neighbours
	return []


def get_mappings(molecule, template, conformer_id=-1):
	match = molecule.GetSubstructMatch(template)
	conformer = template.GetConformer(conformer_id)
	coordMap = {}
	algMap = []
	for i, atomidx in enumerate(match):
		algMap.append((atomidx, i))
		coordMap[atomidx] = conformer.GetAtomPosition(i)
	return coordMap, algMap


def get_distance_constrains(coordMap):
	atom_idxs = list(coordMap.keys())
	constrains = []
	for k, i in enumerate(atom_idxs):
		for j in atom_idxs[k + 1 :]:
			d = coordMap[i].Distance(coordMap[j])
			constrains.append((i, j, d))
	return constrains


def template_embed_optimize(target, template, maxsteps, log):
	"""
	Embeds a new conformation into a molecule, optimizes it using UFF and
	realigns it.

	Parameters
	----------
	target : rdkit.Chem.Mol
			Molecule where you want to embed the new conformation
	mol_template : rdkit.Chem.Mol?
			Template molecule to identify the core of the molecule that will have
			its distances frozen in the optimization.
	maxsteps : int
			Number of maximum optimization steps in RDKit.
	log : Logger
			[description]

	Returns
	-------
	molecule, coordMap, algMap, conf_id
			molecule embedded, mapping to the atom instances,
			list of tuples with position and id and int with the conformer id.
	"""

	seed = -1
	force_constant = 10000

	coord_map, alg_map = get_mappings(target, template, conformer_id=-1)

	# add H's to molecule
	molecule = Chem.AddHs(target)

	conf_id = rdDistGeom.EmbedMolecule(molecule, coordMap=coord_map, randomSeed=seed)

	if conf_id < 0:
		log.write("Could not embed molecule.")
		return molecule, None, None, conf_id

	forcefield = Chem.UFFGetMoleculeForceField(molecule, confId=conf_id)

	constraints = get_distance_constrains(coord_map)
	for i, j, d in constraints:
		forcefield.AddDistanceConstraint(i, j, d, d, force_constant)
	forcefield.Initialize()
	forcefield.Minimize(maxIts=maxsteps)
	# rotate the embedded conformation onto the core_mol:
	rdMolAlign.AlignMol(molecule, template, atomMap=alg_map, reflect=True, maxIters=100)

	return molecule, coord_map, alg_map, conf_id


def filter_template_mol(molecule_new, mol_objects, heavyonly, max_matches):
	"""
	Returns if a mol object should be kept or not.

	Parameters
	----------
	molecule_new : [type]
			[description]
	mol_objects : [type]
			[description]
	heavyonly : bool
			If True only non-H atoms are considered for the RMS calculation
	max_matches : int
			Maximum number of matches in the RMSD?

	Returns
	-------
	bool
			Returns True when the molecule should be kept and False when it should
			be discarded.
	"""

	if not mol_objects:
		return True

	# check if molecule also exixts in the mol_objects
	for mol in mol_objects:
		rms = get_conf_RMS(mol, molecule_new, -1, -1, heavyonly, max_matches)
		if rms < 0.5:
			return False
	return True


# Decorators for documentation
def doc_parameters(f):
	"""
	Decorator that adds the "Parameters" section at the end of the
	docstrings of the decorated function. Care to use this decorator 'below' the
	doc_returns decorator.

	Parameters
	----------
	f : function
			function to decorate.

	Returns
	-------
	function
			returns the same function with the docstring modified.
	"""
	description = f.__doc__
	parameters = [
		("molecule", "rdkit.Chem.Mol", "Molecule to be embedded "),
		("template", "rdkit.Chem.Mol", "Template molecule to do the embedding"),
		("neighbours", "list", "Idx of the atoms neighbouring the metal"),
		("name_input", "str", "Base name for the embedded molecules"),
		("maxsteps", "int", "Maximum number of optimization steps run with rdkit."),
		("log", "Logger", "[description]"),
	]
	item_fmt = "{} : {}\n    {}".format
	params_txt = "\n".join([item_fmt(*items) for items in parameters])
	f.__doc__ = f"{description}\nParameters\n----------\n{params_txt}\n"
	return f


def doc_returns(f):
	"""
	Decorator that adds the "Returns" section at the end of the
	docstrings of the decorated function.

	Parameters
	----------
	f : function
			function to decorate.

	Returns
	-------
	function
			returns the same function with the docstring modified.
	"""
	description = f.__doc__
	item_fmt = "{} : {}\n    {}".format
	outputs = [
		("mol_objects", "list", "Embedded molecules."),
		("name_returns", "list", "Names for the embedded molecules"),
		(
			"coord_maps",
			"list",
			"Mappings to the Idx of the atoms afected by the embedding",
		),
		("alg_maps", "list", "Mappings to the Idx of the core to do alignments"),
		("mol_templates", "list", "List of template molecules used"),
	]
	outs_txt = "\n".join([item_fmt(*items) for items in outputs])
	f.__doc__ = f"{description}\nReturns\n-------\n{outs_txt}\n"
	return f


# Embedding functions
@doc_returns
@doc_parameters
def two_embed(molecule, template, neighbours, name, maxsteps, log):
	"""
	Embedding function for linear geometries. Requires 'linear.sdf' template.
	"""
	template.GetAtomWithIdx(0).setAtomicNum(neighbours[0].GetAtomicNum())
	template.GetAtomWithIdx(1).setAtomicNum(neighbours[1].GetAtomicNum())
	template.GetAtomWithIdx(2).setAtomicNum(53)

	# assigning and embedding onto the core
	mol_obj, coord_map, alg_map, ci = template_embed_optimize(
		molecule, template, maxsteps, log
	)
	if ci >= 0:  # writing to mol_object file
		return [mol_obj], [name], [coord_map], [alg_map], [template]

	return [], [], [], [], []


@doc_returns
@doc_parameters
def three_embed(molecule, template, neighbours, name, maxsteps, log):
	"""
	Embedding function for trigonal planar geometry. Requires
	'trigonalplanar.sdf' template.
	"""
	template.GetAtomWithIdx(0).setAtomicNum(53)
	template.GetAtomWithIdx(1).setAtomicNum(neighbours[0].GetAtomicNum())
	template.GetAtomWithIdx(2).setAtomicNum(neighbours[1].GetAtomicNum())
	template.GetAtomWithIdx(3).setAtomicNum(neighbours[2].GetAtomicNum())

	# assigning and embedding onto the core
	mol_obj, coord_map, alg_map, conf_id = template_embed_optimize(
		molecule, template, maxsteps, log
	)
	if conf_id >= 0:  # writing to mol_object file
		return [mol_obj], [name], [coord_map], [alg_map], [template]

	return [], [], [], [], []


@doc_returns
@doc_parameters
def four_embed(molecule, template, neighbours, name, maxsteps, log):
	"""
	Embedding function for squareplanar geometry. Requires 'template-4-and-5.sdf'
	template. Attempts 3 embeddings.
	"""
	mol_objects = []
	name_return = []
	coord_maps = []
	alg_maps = []
	mol_templates = []

	# Remove F atoms from the template
	for atom in template.GetAtoms():
		if atom.GetSymbol() == "F":
			template = Chem.RWMol(template)
			template.RemoveAtom(atom.GetIdx())
	template = template.GetMol()

	# three cases for square planar
	atn0 = neighbours[0].GetAtomicNum()
	atn1 = neighbours[1].GetAtomicNum()
	atn2 = neighbours[2].GetAtomicNum()
	atn3 = neighbours[3].GetAtomicNum()
	replacements_list = [
		(atn0, atn1, atn2, atn3, 14),
		(atn0, atn2, atn3, atn1, 14),
		(atn0, atn3, atn1, atn2, 14),
	]

	for i, replacements in enumerate(replacements_list):

		# Create a copy of the mol object
		mol = Chem.Mol(template)

		# Assign atomic numbers to neighbour atoms
		for idx, atn in enumerate(replacements):
			mol.GetAtomWithIdx(idx).SetAtomicNum(atn)

		# embedding of the molecule onto the core
		mol_obj, coord_map, alg_map, conf_id = template_embed_optimize(
			molecule, mol, maxsteps, log
		)

		if conf_id >= 0:
			name_out = f"{name.split()[0]}_{i}"
			mol_objects.append(mol_obj)
			name_return.append(name_out)
			coord_maps.append(coord_map)
			alg_maps.append(alg_map)
			mol_templates.append(mol)

	return mol_objects, name_return, coord_maps, alg_maps, mol_templates


@doc_returns
@doc_parameters
def five_embed(molecule, mol_template, neighbours, name, maxsteps, log):
	"""
	Embedding function for squarepyramidal geometry. Requires
	'template-4-and-5.sdf' template. Attempts 15 embeddings.
	"""
	mol_objects = []
	name_return = []
	coord_maps = []
	alg_maps = []
	mol_templates = []
	counter = 0
	atomic_numbers = [mol_template.GetAtomWithIdx(i).GetAtomicNum() for i in neighbours]
	replacements = [
		[4, 0, 1, 2, 3],
		[4, 0, 2, 3, 1],
		[4, 0, 3, 1, 2],
		[3, 0, 1, 2, 4],
		[3, 0, 2, 4, 1],
		[3, 0, 4, 1, 2],
		[2, 0, 1, 4, 3],
		[2, 0, 4, 3, 1],
		[2, 0, 4, 1, 3],
		[1, 0, 4, 2, 3],
		[1, 0, 2, 3, 4],
		[1, 0, 3, 4, 2],
		[0, 4, 1, 2, 3],
		[0, 4, 2, 3, 1],
		[0, 4, 3, 1, 2],
	]
	for replacement in replacements:
		at0, at1, at2, at3, at4 = [atomic_numbers[r] for r in replacement]
		mol_template.GetAtomWithIdx(0).SetAtomicNum(at0)
		mol_template.GetAtomWithIdx(1).SetAtomicNum(at1)
		mol_template.GetAtomWithIdx(2).SetAtomicNum(at2)
		mol_template.GetAtomWithIdx(3).SetAtomicNum(at3)
		mol_template.GetAtomWithIdx(4).SetAtomicNum(at4)
		mol_template.GetAtomWithIdx(5).SetAtomicNum(14)
		mol_template.GetAtomWithIdx(5).SetFormalCharge(1)

		# assigning and embedding onto the core
		mol_obj, coord_map, alg_map, conf_id = template_embed_optimize(
			molecule, mol_template, maxsteps, log
		)
		if conf_id >= 0:
			name_out = f"{name}_{counter}"
			mol_objects.append(mol_obj)
			name_return.append(name_out)
			coord_maps.append(coord_map)
			alg_maps.append(alg_map)
			mol_templates.append(mol_template)
			counter += 1
	return mol_objects, name_return, coord_maps, alg_maps, mol_templates


def template_embed(self, mol, complex_type, metal_idx, maxsteps, heavyonly, maxmatches):
	"""
	Wrapper function to select automatically the appropiate embedding function
	depending on the number of neighbours of the metal center.

	"""
	embed_functions = dict()
	embed_functions[2] = two_embed
	embed_functions[3] = three_embed
	embed_functions[4] = four_embed
	embed_functions[5] = five_embed

	template = load_template(complex_type, self.args.log)

	# Generate the embeddings
	neighbours = calc_neighbours(mol, metal_idx)
	embed = embed_functions[len(neighbours)]
	items = embed(mol, template, neighbours, self.args.name, maxsteps, self.args.log)

	# Filter the results
	molecules = items[0]
	if len(molecules) > 1:
		ignored = []
		for i, mol in enumerate(molecules):
			has_big_rmsd = filter_template_mol(mol, molecules, heavyonly, maxmatches)
			if has_big_rmsd:
				ignored.append(i)
		items = [item for i, item in enumerate(items) if i not in ignored]
	return items


def creation_of_dup_csv_csearch(program):

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
	is_rdkit = program == "rdkit"
	is_fullmonte = program == "fullmonte"
	is_crest = program == "crest"
	is_summ = program == "summ"

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


def prepare_direct_smi(args):
	job_inputs = []
	constraints_dist = args.constraints_dist
	constraints_angle = args.constraints_angle
	constraints_dihedral = args.constraints_dihedral
	if args.prefix == "":
		name = "".join(args.name)
	else:
		name = f"{args.prefix}_{''.join(args.name)}"
	obj = (
		args.smi,
		name,
		constraints_dist,
		constraints_angle,
		constraints_dihedral,
	)
	job_inputs.append(obj)

	return job_inputs


def prepare_smiles_files(args):
	with open(args.input) as smifile:
		lines = [line for line in smifile if line.strip()]
	job_inputs = []
	for i, line in enumerate(lines):
		(
			smi,
			name,
			constraints_dist,
			constraints_angle,
			constraints_dihedral,
		) = prepare_smiles_from_line(line, i, args)
		obj = (
			smi,
			name,
			constraints_dist,
			constraints_angle,
			constraints_dihedral,
		)
		job_inputs.append(obj)

	return job_inputs


def prepare_smiles_from_line(line, i, args):

  toks = line.split()
    # editing part
    smiles = toks[0]
    if args.prefix == "":
        name = "".join(toks[1])
    else:
        name = f"{args.prefix}_{i}_{''.join(toks[1])}"
    constraints_dist = args.constraints_dist
    constraints_angle = args.constraints_angle
    constraints_dihedral = args.constraints_dihedral
    if len(toks) > 2:
        constraints_dist = toks[2]
        constraints_dist = constraints_dist.split("/")
        for j, c in enumerate(constraints_dist):
            constraints_dist[j] = c.split("-")
        if len(toks) > 3:
            constraints_angle = toks[3]
            constraints_angle = constraints_angle.split("/")
            for k, c in enumerate(constraints_angle):
                constraints_angle[k] = c.split("-")
            if len(toks) > 4:
                constraints_dihedral = toks[4]
                constraints_dihedral = constraints_dihedral.split("/")
                for l, c in enumerate(constraints_dihedral):
                    constraints_dihedral[l] = c.split("-")

    return smiles, name, constraints_dist, constraints_angle, constraints_dihedral


def prepare_csv_files(args):
	csv_smiles = pd.read_csv(args.input)
	job_inputs = []
	for i in range(len(csv_smiles)):
		obj = generate_mol_from_csv(args, csv_smiles, i)
		job_inputs.append(obj)
	return job_inputs


def generate_mol_from_csv(args, csv_smiles, index):
	# assigning names and smi i  each loop
	try:
		smiles = csv_smiles.loc[index, "SMILES"]
	except KeyError:
		try:
			smiles = csv_smiles.loc[index, "smiles"]
		except KeyError:
			print(
				"\nx  Make sure the CSV file contains a column called 'SMILES' or 'smiles' with the SMILES of the molecules!"
			)
			args.log.write(
				"\nx  Make sure the CSV file contains a column called 'SMILES' or 'smiles' with the SMILES of the molecules!"
			)
			sys.exit()

	# pruned_smi = smi.split(".")
	# mol = Chem.MolFromSmiles(pruned_smi)

	try:
		if args.prefix == "":
			name = csv_smiles.loc[index, "code_name"]
		else:
			name = f'{args.prefix}_{csv_smiles.loc[index, "code_name"]}'
	except KeyError:
		print(
			"\nx  Make sure the CSV file contains a column called 'code_name' with the names of the molecules!"
		)
		args.log.write(
			"\nx  Make sure the CSV file contains a column called 'code_name' with the names of the molecules!"
		)
		sys.exit()

	constraints_dist = args.constraints_dist
	constraints_angle = args.constraints_angle
	constraints_dihedral = args.constraints_dihedral
	if "constraints_dist" in csv_smiles.columns:
		constraints_dist = csv_smiles.loc[index, "constraints_dist"]
		constraints_dist = constraints_dist.split("/")
		for i, c in enumerate(constraints_dist):
			constraints_dist[i] = c.split("-")
	if "constraints_angle" in csv_smiles.columns:
		constraints_angle = csv_smiles.loc[index, "constraints_angle"]
		constraints_angle = constraints_angle.split("/")
		for i, c in enumerate(constraints_angle):
			constraints_angle[i] = c.split("-")
	if "constraints_dihedral" in csv_smiles.columns:
		constraints_dihedral = csv_smiles.loc[index, "constraints_dihedral"]
		constraints_dihedral = constraints_dihedral.split("/")
		for i, c in enumerate(constraints_dihedral):
			constraints_dihedral[i] = c.split("-")

	obj = (smiles, name, constraints_dist, constraints_angle, constraints_dihedral)

	return obj


def prepare_cdx_files(args):
	# converting to smiles from chemdraw
	molecules = generate_mol_from_cdx(args)
	job_inputs = []
	for i, (smiles, _) in enumerate(molecules):
		name = f"{args.input.split('.')[0]}_{str(i)}"
		constraints_dist = args.constraints_dist
		constraints_angle = args.constraints_angle
		constraints_dihedral = args.constraints_dihedral
		obj = (smiles, name, constraints_dist, constraints_angle, constraints_dihedral)
		job_inputs.append(obj)
	return job_inputs


def generate_mol_from_cdx(args):
	cmd_cdx = ["obabel", "-icdx", args.input, "-osmi", "-Ocdx.smi"]
	subprocess.call(cmd_cdx)
	with open("cdx.smi", "r") as smifile:
		smi_lines = [str(line.strip()) for line in smifile]
	os.remove("cdx.smi")
	molecules = []
	for smi in smi_lines:
		molecule = Chem.MolFromSmiles(smi)
		molecules.append((smi, molecule))
	return molecules


def prepare_xyz_files(args):
	job_inputs = []
	charge_com, mult_com = com_2_xyz_2_sdf(args.input, args.charge, args.mult)
	name = os.path.splitext(args.input)[0]
	sdffile = f"{name}.sdf"
	suppl, _, _ = mol_from_sdf_or_mol_or_mol2(sdffile)

	for _, mol in enumerate(suppl):
		# if args.charge is None:
		# 	args.charge = charge_com
		constraints_dist = args.constraints_dist
		constraints_angle = args.constraints_angle
		constraints_dihedral = args.constraints_dihedral
		obj = (mol, name, constraints_dist, constraints_angle, constraints_dihedral)
		job_inputs.append(obj)
	return job_inputs


def prepare_sdf_files(args):
	suppl, IDs, charges = mol_from_sdf_or_mol_or_mol2(args.input)
	job_inputs = []
	for _, (mol, name, _) in enumerate(zip(suppl, IDs, charges)):
		# if args.charge == None:
		# 	args.charge = charge_sdf
		constraints_dist = args.constraints_dist
		constraints_angle = args.constraints_angle
		constraints_dihedral = args.constraints_dihedral
		obj = (mol, name, constraints_dist, constraints_angle, constraints_dihedral)
		job_inputs.append(obj)
	return job_inputs


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
	# if parent_dir is None:
	# 	parent_dir = Path("")
	# else:
	# 	parent_dir = Path(parent_dir)
	#
	# mol = next(pybel.readfile("xyz", parent_dir / file))
	# ofile = Path(file).stem + ".sdf"
	# mol.write("sdf", parent_dir / ofile)

	name1 = str(file).split(".xyz")[0]
	command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name1 + ".sdf"]
	subprocess.call(command_xyz)


def com_2_xyz_2_sdf(input_file, default_charge, default_mult, start_point=None):
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
		if extension in [".com", ".gjf", ".xyz"]:
			file = Path(input_file)
	else:
		file = Path(start_point)

	filename = str(input_file).split('.')[0]
	# Create the 'xyz' file and/or get the total charge
	if extension != ".xyz":
		xyz, charge, mult = get_info_input(file)
		xyz_txt = "\n".join(xyz)
		with open(f"{filename}.xyz", "w") as F:
			F.write(f"{len(xyz)}\n{filename}\n{xyz_txt}\n")
		xyz_2_sdf(f"{filename}.xyz")
	else:
		charge = default_charge
		mult = default_mult
		xyz_2_sdf(file)

	return charge, mult


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
		suppl = [Chem.MolFromMolFile(input_file, removeHs=False)]
	elif extension == ".mol2":
		suppl = [Chem.MolFromMol2File(input_file, removeHs=False)]

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


def minimize_rdkit_energy(mol, conf, log, FF, maxsteps):
	"""
	Minimizes a conformer of a molecule and returns the final energy.
	"""

	if FF == "MMFF":
		properties = Chem.MMFFGetMoleculeProperties(mol)
		forcefield = Chem.MMFFGetMoleculeForceField(mol, properties, confId=conf)

	if FF == "UFF" or forcefield is None:
		# if forcefield is None means that MMFF will not work. Attempt UFF.
		forcefield = Chem.UFFGetMoleculeForceField(mol, confId=conf)

	if FF not in ["MMFF", "UFF"] or forcefield is None:
		log.write(f" Force field {FF} not supported!")
		sys.exit()

	forcefield.Initialize()
	forcefield.Minimize(maxIts=maxsteps)
	energy = float(forcefield.CalcEnergy())

	return energy
