#####################################################.
#       This file stores functions related to       #
#          metal templates used in CSEARCH          #
#####################################################.

import sys
from pathlib import Path
from pkg_resources import resource_filename
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdDistGeom, rdMolAlign
from aqme.utils import get_conf_RMS, load_sdf

TEMPLATES_PATH = Path(resource_filename("aqme", "templates"))


def template_embed(self, mol, complex_type, metal_idx, maxsteps, heavyonly, maxmatches, name, geom):
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
    items = embed(mol, template, metal_idx, neighbours, name, maxsteps, self.args.log, geom)
    
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


def template_embed_optimize(target, template, metal_idx, maxsteps, log, tempalte_n=None, cumulative_algMap=[]):
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
        Molecule embedded, mapping to the atom instances,
        list of tuples with position and id and int with the conformer id.
    """

    seed = -1
    force_constant = 10000

    coord_map, alg_map, cumulative_algMap = get_mappings(target, template, metal_idx, cumulative_algMap, tempalte_n, conformer_id=-1)

    # add H's to molecule
    molecule = Chem.AddHs(target)

    conf_id = rdDistGeom.EmbedMolecule(molecule, coordMap=coord_map, randomSeed=seed)

    if conf_id < 0:
        log.write("Could not embed molecule.")
        return molecule, None, None, conf_id, cumulative_algMap

    forcefield = Chem.UFFGetMoleculeForceField(molecule, confId=conf_id)

    constraints = get_distance_constrains(coord_map)
    for i, j, d in constraints:
        forcefield.AddDistanceConstraint(i, j, d, d, force_constant)
    forcefield.Initialize()
    forcefield.Minimize(maxIts=maxsteps)
    # rotate the embedded conformation onto the core_mol:
    rdMolAlign.AlignMol(molecule, template, atomMap=alg_map, reflect=True, maxIters=100)

    return molecule, coord_map, alg_map, conf_id, cumulative_algMap


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
        Returns True when the molecule should be kept and False when it should be discarded.
    """

    if not mol_objects:
        return True

    # check if molecule also exixts in the mol_objects
    for mol in mol_objects:
        rms = get_conf_RMS(mol, molecule_new, -1, -1, heavyonly, max_matches)
        if rms < 0.5:
            return False
    return True


def get_mappings(molecule, template, metal_idx, cumulative_algMap, tempalte_n, conformer_id=-1):
    match = molecule.GetSubstructMatch(template)
    conformer = template.GetConformer(conformer_id)
    coordMap = {}
    algMap = []
    for i, atomidx in enumerate(match):
        algMap.append((atomidx, i))
        coordMap[atomidx] = conformer.GetAtomPosition(i)
    # in cases where the neighbours of the metal are the same types of atoms but within different ligands
    if algMap in cumulative_algMap and tempalte_n is not None:
        new_coordMap,new_algMap = {},[]
        # create a list with the ligands and the metal separated. This is a fix that accounts for
        # matchings that put the metal atom in different order
        neigh_coordMap,neigh_algMap = {},[]
        metal_coordMap,metal_algMap = {},[]
        for idx,match in algMap:
            if idx == metal_idx[0]:
                metal_algMap.append((idx,match))
                metal_coordMap[idx] = coordMap[idx]
            else:
                neigh_algMap.append((idx,match))
                neigh_coordMap[idx] = coordMap[idx]

        if tempalte_n == 1:
            rel_idx = [2,3,1]
        if tempalte_n == 2:
            rel_idx = [3,1,2]
        coord_keys = list(neigh_coordMap.keys()) # add ligands
        new_coordMap[coord_keys[0]] = neigh_coordMap[coord_keys[0]]
        new_coordMap[coord_keys[rel_idx[0]]] = neigh_coordMap[coord_keys[1]]
        new_coordMap[coord_keys[rel_idx[1]]] = neigh_coordMap[coord_keys[2]]
        new_coordMap[coord_keys[rel_idx[2]]] = neigh_coordMap[coord_keys[3]]
        new_coordMap[metal_idx[0]] = metal_coordMap[metal_idx[0]] # add metal
        coordMap = new_coordMap
        new_algMap.append(neigh_algMap[0]) # add ligands
        new_algMap.append((neigh_algMap[rel_idx[0]][0],neigh_algMap[1][1]))
        new_algMap.append((neigh_algMap[rel_idx[1]][0],neigh_algMap[2][1]))
        new_algMap.append((neigh_algMap[rel_idx[2]][0],neigh_algMap[3][1]))
        new_algMap.append(metal_algMap[0]) # add metal
        algMap = new_algMap

    cumulative_algMap.append(algMap)

    return coordMap, algMap, cumulative_algMap


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
        log.write("x  The templates folder was not found, probably due to a problem while installing AQME")
        log.finalize()
        sys.exit()

    file_template = folder / Path(type2template[complex_type])
    templates = load_sdf(str(file_template))
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
        Mol object
    metals_idx : list
        List containing the Idx of the metals. The first match is the only one considered.

    Returns
    -------
    neighbours : list
        List of neighbour atoms
    """

    # depending on the amount of neighbours, use Si or I atoms to fit the templates
    bonds2AtNum = dict()
    bonds2AtNum[5] = 14
    bonds2AtNum[4] = 14
    bonds2AtNum[3] = 53
    bonds2AtNum[2] = 53

    for atom in molecule.GetAtoms():
        idx = atom.GetIdx()
        if idx in metals_idx:
            neighbours = atom.GetNeighbors()
            # in case metals are used with different bonds (i.e., Cu2+ and CuL2)
            n_bonds = len(atom.GetBonds())
            AtNum = bonds2AtNum[n_bonds]
            atom.SetAtomicNum(AtNum)
            if n_bonds == 5:
                atom.SetFormalCharge(1)
            return neighbours
    return []


def check_metal_neigh(mol, complex_type, metal_idx_ind, log, valid_template):
    """
    Checks if the metal and the template contain the same number of ligands.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Mol object
    complex_type : str
        Type of template to be used (i.e., squareplanar)
    metals_idx : list
        List containing the Idx of the metals

    Returns
    -------
    valid_template : bool
        Whether the complexes are compatible with the template selected
    """

    if complex_type == "linear":
        expect_neigh = 2
    elif complex_type == "trigonalplanar":
        expect_neigh = 3
    elif complex_type == "squareplanar":
        expect_neigh = 4
    elif complex_type == "squarepyramidal":
        expect_neigh = 5
    metal_atom = mol.GetAtoms()[metal_idx_ind]
    metal_neigh = metal_atom.GetNeighbors()
    if len(metal_neigh) != expect_neigh:
        log.write(f"x  The number of neighbours of the metal ({len(metal_neigh)}) does not match the number of expected neighbours for the template selected ({complex_type}). No templates will be applied to this system.")
        valid_template = False

    return valid_template


def get_distance_constrains(coordMap):
    atom_idxs = list(coordMap.keys())
    constrains = []
    for k, i in enumerate(atom_idxs):
        for j in atom_idxs[k + 1 :]:
            d = coordMap[i].Distance(coordMap[j])
            constrains.append((i, j, d))
    return constrains


# Decorators for documentation
def doc_parameters(f):
    """
    Decorator that adds the "Parameters" section at the end of the
    docstrings of the decorated function. Care to use this decorator 'below' the
    doc_returns decorator.

    Parameters
    ----------
    f : function
        Function to decorate.

    Returns
    -------
    function
        Returns the same function with the docstring modified.
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
        Function to decorate.

    Returns
    -------
    function
        Returns the same function with the docstring modified.
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
def two_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log, geom):
    """
    Embedding function for linear geometries. Requires 'linear.sdf' template.
    """

    template.GetAtomWithIdx(0).SetAtomicNum(neighbours[0].GetAtomicNum())
    template.GetAtomWithIdx(1).SetAtomicNum(neighbours[1].GetAtomicNum())
    template.GetAtomWithIdx(2).SetAtomicNum(53)

    # assigning and embedding onto the core
    mol_obj, coord_map, alg_map, ci, _ = template_embed_optimize(
        molecule, template, metal_idx, maxsteps, log
    )

    original_atn_list = [None] # only working for Ir squareplanar

    if ci >= 0:  # writing to mol_object file
        return [mol_obj], [name], [coord_map], [alg_map], [template], original_atn_list

    return [], [], [], [], [], original_atn_list


@doc_returns
@doc_parameters
def three_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log, geom):
    """
    Embedding function for trigonal planar geometry. Requires
    'trigonalplanar.sdf' template.
    """

    template.GetAtomWithIdx(0).SetAtomicNum(53)
    template.GetAtomWithIdx(1).SetAtomicNum(neighbours[0].GetAtomicNum())
    template.GetAtomWithIdx(2).SetAtomicNum(neighbours[1].GetAtomicNum())
    template.GetAtomWithIdx(3).SetAtomicNum(neighbours[2].GetAtomicNum())

    # assigning and embedding onto the core
    mol_obj, coord_map, alg_map, conf_id, _ = template_embed_optimize(
        molecule, template, metal_idx, maxsteps, log
    )
    
    original_atn_list = [None] # only working for Ir squareplanar

    if conf_id >= 0:  # writing to mol_object file
        return [mol_obj], [name], [coord_map], [alg_map], [template], original_atn_list

    return [], [], [], [], [], original_atn_list


@doc_returns
@doc_parameters
def four_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log, geom):
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

    # there are atoms that give problems in the embedding for Ir squareplanar complexes (i.e. As).
    # This replaces them for P and them restore them back to normal after the embeding (in genConformer_r() from base.py)
    original_atn = None
    if geom == ['Ir_squareplanar']:
        invalid_atn = [33]
        atn_list = [atn0,atn1,atn2,atn3]
        for i,atn in enumerate(atn_list):
            if atn in invalid_atn:
                original_atn = [atn,neighbours[i].GetIdx()]
                molecule.GetAtomWithIdx(neighbours[i].GetIdx()).SetAtomicNum(15)
                atn_list[i] = 15

        atn0 = atn_list[0]
        atn1 = atn_list[1]
        atn2 = atn_list[2]
        atn3 = atn_list[3]

    replacements_list = [
        (atn0, atn1, atn2, atn3, 14),
        (atn0, atn2, atn3, atn1, 14),
        (atn0, atn3, atn1, atn2, 14),
    ]
    
    cumulative_algMap = []
    original_atn_list = []
    for i, replacements in enumerate(replacements_list):

        # Create a copy of the mol object
        template_mol = Chem.Mol(template)

        # Assign atomic numbers to neighbour atoms
        for idx, atn in enumerate(replacements):
            template_mol.GetAtomWithIdx(idx).SetAtomicNum(atn)

        # embedding of the molecule onto the core
        mol_obj, coord_map, alg_map, conf_id, cumulative_algMap = template_embed_optimize(
            molecule, template_mol, metal_idx, maxsteps, log, tempalte_n=i, cumulative_algMap=cumulative_algMap
        )

        if conf_id >= 0:
            name_out = f"{name.split()[0]}_{i}"
            mol_objects.append(mol_obj)
            name_return.append(name_out)
            coord_maps.append(coord_map)
            alg_maps.append(alg_map)
            mol_templates.append(template_mol)
            original_atn_list.append(original_atn)

    return mol_objects, name_return, coord_maps, alg_maps, mol_templates, original_atn_list


@doc_returns
@doc_parameters
def five_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log, geom):
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
    atomic_numbers = []
    for _,atom in enumerate(neighbours):
        atomic_numbers.append(atom.GetAtomicNum())
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

    original_atn_list = []
    for replacement in replacements:
        at0, at1, at2, at3, at4 = [atomic_numbers[r] for r in replacement]
        template.GetAtomWithIdx(0).SetAtomicNum(at0)
        template.GetAtomWithIdx(1).SetAtomicNum(at1)
        template.GetAtomWithIdx(2).SetAtomicNum(at2)
        template.GetAtomWithIdx(3).SetAtomicNum(at3)
        template.GetAtomWithIdx(4).SetAtomicNum(at4)
        template.GetAtomWithIdx(5).SetAtomicNum(14)
        template.GetAtomWithIdx(5).SetFormalCharge(1)

        # assigning and embedding onto the core
        mol_obj, coord_map, alg_map, conf_id, _ = template_embed_optimize(
            molecule, template, metal_idx, maxsteps, log
        )
        if conf_id >= 0:
            name_out = f"{name}_{counter}"
            mol_objects.append(mol_obj)
            name_return.append(name_out)
            coord_maps.append(coord_map)
            alg_maps.append(alg_map)
            mol_templates.append(template)
            original_atn_list.append(None) # only working for Ir squareplanar
            counter += 1
    return mol_objects, name_return, coord_maps, alg_maps, mol_templates, original_atn_list
