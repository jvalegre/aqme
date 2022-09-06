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
import ast
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdDistGeom, rdMolAlign
from aqme.utils import (
    get_info_input,
    get_conf_RMS,
    mol_from_sdf_or_mol_or_mol2,
    read_xyz_charge_mult,
)

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
        log.finalize()
        sys.exit()

    file_template = folder / Path(type2template[complex_type])
    templates = Chem.SDMolSupplier(str(file_template))
    template = templates[-1]

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
        for i, mol_filter in enumerate(molecules):
            has_big_rmsd = filter_template_mol(
                mol_filter, molecules, heavyonly, maxmatches
            )
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


def constraint_fix(
    constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral
):

    # this avoids problems when running AQME through command lines
    if pd.isnull(constraints_atoms):
        constraints_atoms = []
    else:
        if not isinstance(constraints_atoms, list):
            constraints_atoms = ast.literal_eval(constraints_atoms)

    if pd.isnull(constraints_dist):
        constraints_dist = []
    else:
        if not isinstance(constraints_dist, list):
            constraints_dist = ast.literal_eval(constraints_dist)

    if pd.isnull(constraints_angle):
        constraints_angle = []
    else:
        if not isinstance(constraints_angle, list):
            constraints_angle = ast.literal_eval(constraints_angle)

    if pd.isnull(constraints_dihedral):
        constraints_dihedral = []
    else:
        if not isinstance(constraints_dihedral, list):
            constraints_dihedral = ast.literal_eval(constraints_dihedral)

    return constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral


def prepare_direct_smi(args):
    job_inputs = []
    if args.prefix == "":
        name = "".join(args.name)
    else:
        name = f"{args.prefix}_{''.join(args.name)}"
    (
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
    ) = constraint_fix(
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
    )
    obj = (
        args.smi,
        name,
        args.charge,
        args.mult,
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
    )
    job_inputs.append(obj)

    return job_inputs


def prepare_smiles_files(args, csearch_file):
    with open(csearch_file) as smifile:
        lines = [line for line in smifile if line.strip()]
    job_inputs = []
    for i, line in enumerate(lines):
        (
            smi,
            name,
        ) = prepare_smiles_from_line(line, i, args)
        (
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        ) = constraint_fix(
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        )
        obj = (
            smi,
            name,
            args.charge,
            args.mult,
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
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

    return smiles, name


def prepare_csv_files(args, csearch_file):
    csv_smiles = pd.read_csv(csearch_file)
    job_inputs = []
    for i in range(len(csv_smiles)):
        obj = generate_mol_from_csv(args, csv_smiles, i)
        job_inputs.append(obj)
    return job_inputs


def generate_mol_from_csv(args, csv_smiles, index):
    # assigning names and smi in each loop
    try:
        smiles = csv_smiles.loc[index, "SMILES"]
    except KeyError:
        try:
            smiles = csv_smiles.loc[index, "smiles"]
        except KeyError:
            args.log.write(
                "\nx  Make sure the CSV file contains a column called 'SMILES' or 'smiles' with the SMILES of the molecules!"
            )
            args.log.finalize()
            sys.exit()

    try:
        if args.prefix == "":
            name = csv_smiles.loc[index, "code_name"]
        else:
            name = f'{args.prefix}_{csv_smiles.loc[index, "code_name"]}'
    except KeyError:
        args.log.write(
            "\nx  Make sure the CSV file contains a column called 'code_name' with the names of the molecules!"
        )
        args.log.finalize()
        sys.exit()

    constraints_atoms = args.constraints_atoms
    constraints_dist = args.constraints_dist
    constraints_angle = args.constraints_angle
    constraints_dihedral = args.constraints_dihedral

    if "constraints_atoms" in csv_smiles.columns:
        constraints_atoms = csv_smiles.loc[index, "constraints_atoms"]

    if "constraints_dist" in csv_smiles.columns:
        constraints_dist = csv_smiles.loc[index, "constraints_dist"]

    if "constraints_angle" in csv_smiles.columns:
        constraints_angle = csv_smiles.loc[index, "constraints_angle"]

    if "constraints_dihedral" in csv_smiles.columns:
        constraints_dihedral = csv_smiles.loc[index, "constraints_dihedral"]

    (
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
    ) = constraint_fix(
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
    )
    obj = (
        smiles,
        name,
        args.charge,
        args.mult,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
    )

    return obj


def prepare_cdx_files(args, csearch_file):
    # converting to smiles from chemdraw
    molecules = generate_mol_from_cdx(csearch_file)

    job_inputs = []
    for i, (smiles, _) in enumerate(molecules):
        name = f"{csearch_file.split('.')[0]}_{str(i)}"
        (
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        ) = constraint_fix(
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        )
        obj = (
            smiles,
            name,
            args.charge,
            args.mult,
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        )
        job_inputs.append(obj)
    return job_inputs


def generate_mol_from_cdx(csearch_file):
    cmd_cdx = ["obabel", "-icdx", csearch_file, "-osmi", "-Ocdx.smi"]
    subprocess.run(cmd_cdx, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open("cdx.smi", "r") as smifile:
        smi_lines = [str(line.strip()) for line in smifile]
    os.remove("cdx.smi")
    molecules = []
    for smi in smi_lines:
        molecule = Chem.MolFromSmiles(smi)
        molecules.append((smi, molecule))
    return molecules


def prepare_com_files(args, csearch_file):
    job_inputs = []

    if csearch_file.split(".")[1] in ["gjf", "com"]:
        xyz_file, _, _ = com_2_xyz(csearch_file)
        _, charge, mult = get_info_input(csearch_file)
    else:
        xyz_file = csearch_file
        charge, mult = read_xyz_charge_mult(xyz_file)
    xyz_2_sdf(xyz_file)
    name = os.path.splitext(csearch_file)[0]
    sdffile = f"{name}.sdf"
    suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(sdffile, "csearch")

    (
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
    ) = constraint_fix(
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
    )

    obj = (
        suppl[0],
        name,
        charge,
        mult,
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
    )
    job_inputs.append(obj)

    return job_inputs


def prepare_pdb_files(args, csearch_file):
    command_pdb = [
        "obabel",
        "-ipdb",
        csearch_file,
        "-osdf",
        f'-O{csearch_file.split(".")[0]}.sdf',
    ]
    subprocess.run(command_pdb, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    job_inputs = prepare_sdf_files(args, csearch_file)
    os.remove(f'{csearch_file.split(".")[0]}.sdf')
    return job_inputs


def prepare_sdf_files(args, csearch_file):
    suppl, charges, mults, IDs = mol_from_sdf_or_mol_or_mol2(csearch_file, "csearch")
    job_inputs = []

    for mol, charge, mult, name in zip(suppl, charges, mults, IDs):
        (
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        ) = constraint_fix(
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        )
        obj = (
            mol,
            name,
            charge,
            mult,
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
        )
        job_inputs.append(obj)
    return job_inputs


def xyz_2_sdf(file):
    """
    Creates a .sdf file from a .xyz in the specified directory. If no directory
    is specified then the files are created in the current directory.

    Parameters
    ----------
    file : str
                                                                    Filename and extension of an existing .xyz file
    """

    name = str(file).split(".xyz")[0]
    command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name + ".sdf"]
    subprocess.run(command_xyz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def com_2_xyz(input_file):
    """
    COM to XYZ to SDF for obabel
    """

    filename = input_file.split(".")[0]

    # Create the 'xyz' file and/or get the total charge
    xyz, charge, mult = get_info_input(input_file)
    xyz_txt = "\n".join(xyz)
    with open(f"{filename}.xyz", "w") as F:
        F.write(f"{len(xyz)}\n{filename}\n{xyz_txt}\n")

    return f"{filename}.xyz", charge, mult


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
        log.finalize()
        sys.exit()

    forcefield.Initialize()
    forcefield.Minimize(maxIts=maxsteps)
    energy = float(forcefield.CalcEnergy())

    return energy
