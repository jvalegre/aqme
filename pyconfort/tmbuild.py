#####################################################.
#        This file stores all the functions         #
#    used in template based conformer generation    #
#####################################################.
import os
import sys
from pathlib import Path
from pkg_resources import resource_filename 

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdDistGeom, rdMolAlign

from pyconfort.utils import get_conf_RMS

TEMPLATES_PATH = Path(resource_filename('pyconfort','templates'))

def load_template(complex_type,log):
    """
    Checks if the templates are reachable and if so returns the molecule object. 

    Returns
    -------
    rdkit.Chem.Mol
        The molecule file of the corresponding template.
    """
    type2template = dict() 
    type2template['squareplanar'] = 'template-4-and-5.sdf'
    type2template['squarepyramidal'] = 'template-4-and-5.sdf'
    type2template['linear'] = 'template-2.sdf'
    type2template['trigonalplanar'] = 'template-3.sdf'

    folder = TEMPLATES_PATH
    
    if not folder.exists(): 
        log.write('x The templates folder was not found, probably due to a problem while installing pyCONFORT')
        sys.exit()

    file_template = folder/Path(type2template[complex_type])
    templates = list(Chem.SDMolSupplier(file_template))
    template = templates[-1] # RAUL: I'm assuming that there is only one molecule per template and in case there are several, it's the last one

    return template

def calc_neighbours(molecule,args):
    """
    Changes the atomic number (and charge) of the first metal found 
    and returns its neighbours, the number of neighbours and the idx of the 
    metal.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        [description]
    args : [type]
        [description]

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
    metals_idx = args.metal_idx
    # RAUL: Are we interested in the first Metal found in the molecule or do we 
    #       know beforehand the idx of the Metal atom?. If we know beforehand:  
    #idx = metals_idx[0]
    #atom = molecule.GetAtomWithIdx(idx)
    #n_bonds = len(atom.GetBonds())
    #atom.SetAtomicNum(bonds2AtNum[n_bonds])
    #if n_bonds == 5: 
    #    atom.SetFormalCharge(1)
    #neighbours = atom.GetNeighbors()
    #return idx, neighbours
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

#GET THE LINEAR GEOMETRY
def two_embed(molecule_embed, mol_template, neighbours, name_input, args, log):
    
    mol_template.GetAtomWithIdx(0).setAtomicNum(neighbours[0].GetAtomicNum())
    mol_template.GetAtomWithIdx(1).setAtomicNum(neighbours[1].GetAtomicNum())
    mol_template.GetAtomWithIdx(2).setAtomicNum(53)

    #assigning and embedding onto the core
    mol_obj, coord_map, alg_map,ci = template_embed_optimize(molecule_embed,
                                                             mol_template,
                                                             args, log)
    if ci >= 0: #writing to mol_object file
        return [mol_obj], [name_input], [coord_map], [alg_map], [mol_template]

    return [], [], [], [], []

# GET THE TRIGONAL PLANAR GEOMETRY
def three_embed(molecule_embed, mol_template,neighbours,name_input,args,log):

    mol_template.GetAtomWithIdx(0).setAtomicNum(53)
    mol_template.GetAtomWithIdx(1).setAtomicNum(neighbours[0].GetAtomicNum())
    mol_template.GetAtomWithIdx(2).setAtomicNum(neighbours[1].GetAtomicNum())
    mol_template.GetAtomWithIdx(3).setAtomicNum(neighbours[2].GetAtomicNum())

    #assigning and embedding onto the core
    mol_obj, coord_map, alg_map, conf_id = template_embed_optimize(molecule_embed,
                                                                   mol_template,
                                                                   args, log)
    if conf_id >= 0: #writing to mol_object file
        return [mol_obj], [name_input], [coord_map], [alg_map], [mol_template]

    return [], [], [], [], []

# GET THE SQUAREPLANAR GEOMETRY
def four_embed(molecule_embed, mol_1,neighbours,name_input,args,log):

    mol_objects = []
    name_return = []
    coord_maps = []
    alg_maps = []
    mol_templates = []
    
    # Remove F atoms from the template
    mol_template = Chem.RWMol(mol_1)
    for atom in mol_template.GetAtoms():
        if atom.GetSymbol() == 'F': 
            mol_template.RemoveAtom(atom.GetIdx())
    mol_template = mol_template.GetMol()

    #three cases for square planar
    atn0 = neighbours[0].GetAtomicNum()
    atn1 = neighbours[1].GetAtomicNum()
    atn2 = neighbours[2].GetAtomicNum()
    atn3 = neighbours[3].GetAtomicNum()
    replacements_list = [(atn0,atn1,atn2,atn3,14),
                         (atn0,atn2,atn3,atn1,14),
                         (atn0,atn3,atn1,atn2,14)]

    for i,replacements in enumerate(replacements_list):
        
        # Create a copy of the mol object 
        mol_1 = Chem.Mol(mol_template)

        # Assign atomic numbers to neighbour atoms
        for idx,atn in enumerate(replacements):
            mol_1.GetAtomWithIdx(idx).setAtomicNum(atn)
        
        #embedding of the molecule onto the core
        mol_obj, coord_map, alg_map, conf_id = template_embed_optimize(molecule_embed,
                                                                       mol_1,
                                                                       args, log)

        if conf_id >= 0:
            check = filter_template_mol(mol_obj, mol_objects,args)
            if check:
                #writing to mol_object file
                name = f'{name_input.split()[0]}_{i}'
                mol_objects.append(mol_obj)
                name_return.append(name)
                coord_maps.append(coord_map)
                alg_maps.append(alg_map)
                mol_templates.append(mol_1)

    return mol_objects, name_return, coord_maps, alg_maps, mol_templates

# GET THE SQUAREPYRAMIDAL GEOMETRY
def five_embed(molecule_embed, mol_template,neighbours,name_input,args,log):
    mol_objects   = []
    name_return   = []
    coord_maps    = []
    alg_maps      = []
    mol_templates = []
    counter = 0
    atomic_numbers = [mol_template.GetAtomWithIdx(i).GetAtomicNum() for i in neighbours]
    replacements = [[4,0,1,2,3],
                    [4,0,2,3,1],
                    [4,0,3,1,2],
                    [3,0,1,2,4],
                    [3,0,2,4,1],
                    [3,0,4,1,2],
                    [2,0,1,4,3],
                    [2,0,4,3,1],
                    [2,0,4,1,3],
                    [1,0,4,2,3],
                    [1,0,2,3,4],
                    [1,0,3,4,2],
                    [0,4,1,2,3],
                    [0,4,2,3,1],
                    [0,4,3,1,2]]
    for replacement in replacements:
        at0,at1,at2,at3,at4 = [atomic_numbers[r] for r in replacement]
        mol_template.GetAtomWithIdx(0).SetAtomicNum(at0)
        mol_template.GetAtomWithIdx(1).SetAtomicNum(at1)
        mol_template.GetAtomWithIdx(2).SetAtomicNum(at2)
        mol_template.GetAtomWithIdx(3).SetAtomicNum(at3)
        mol_template.GetAtomWithIdx(4).SetAtomicNum(at4)
        mol_template.GetAtomWithIdx(5).SetAtomicNum(14)
        mol_template.GetAtomWithIdx(5).SetFormalCharge(1)

        #assigning and embedding onto the core
        mol_obj, coord_map, alg_map, conf_id = template_embed_optimize(molecule_embed,
                                                                       mol_template,
                                                                       args,log)
        if conf_id >= 0:
            check = filter_template_mol(mol_obj, mol_objects, args)
            if check:
                #writing to mol_object file
                name = f'{name_input}_{counter}'
                mol_objects.append(mol_obj)
                name_return.append(name)
                coord_maps.append(coord_map)
                alg_maps.append(alg_map)
                mol_templates.append(mol_template)
            counter += 1
    return mol_objects, name_return, coord_maps, alg_maps, mol_templates

# Auxiliar function to get the mappings of the core atoms
def get_mappings(molecule,template,conformer_id=-1):
    match = molecule.GetSubstructMatch(template)
    conformer = template.GetConformer(conformer_id)
    coordMap = {}
    algMap = []
    for i,atomidx in enumerate(match):
        algMap.append((atomidx,i))
        coordMap[atomidx] = conformer.GetAtomPosition(i)
    return coordMap, algMap
def get_distance_constrains(coordMap):
    atom_idxs = list(coordMap.keys())
    constrains = []
    for k, i in enumerate(atom_idxs):
        for j in atom_idxs[k+1:]:
            d = coordMap[i].Distance(coordMap[j])
            constrains.append((i,j,d))
    return constrains

# TEMPLATE GENERATION FOR SQUAREPLANAR AND squarepyramidal
def template_embed(molecule,name_input,args,log):
    """
    Wrapper function to select automatically the appropiate embedding function
    depending on the number of neighbours of the metal center.

    Parameters
    ----------
    molecule : [type]
        [description]
    temp : [type]
        [description]
    name_input : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]

    Returns
    -------
    tuple
        mol_objects, name_return, coord_maps, alg_maps, mol_templates
    """
    embed_functions = dict()
    embed_functions[2] = two_embed
    embed_functions[3] = three_embed
    embed_functions[4] = four_embed
    embed_functions[5] = five_embed

    template = load_template(args.complex_type,log) 

    neighbours = calc_neighbours(molecule,args)
    embed = embed_functions[len(neighbours)]
    
    return embed(molecule,template,neighbours,name_input,args,log)

def template_embed_optimize(target,template,args,log):
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
    args : argparse.args
        [description]
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

    coord_map, alg_map = get_mappings(target,template,conformer_id=-1)

    #add H's to molecule
    molecule = Chem.AddHs(target)
    
    conf_id = rdDistGeom.EmbedMolecule(molecule, coordMap=coord_map, randomSeed=seed)

    if conf_id < 0:
        log.write('Could not embed molecule.')
        return molecule, None, None, conf_id

    forcefield = Chem.UFFGetMoleculeForceField(molecule,confId=conf_id)

    constraints = get_distance_constrains(coord_map)
    for i,j,d in constraints: 
        forcefield.AddDistanceConstraint(i, j, d, d, force_constant)
    forcefield.Initialize()
    forcefield.Minimize(maxIts=args.opt_steps_RDKit)
    # rotate the embedded conformation onto the core_mol:
    rdMolAlign.AlignMol(molecule, template, 
                        atomMap=alg_map,
                        reflect=True,
                        maxIters=100)

    return molecule, coord_map, alg_map, conf_id

# FILTER FOR REMOVING MOLS IF LIGANDS ARE THE SAME
def filter_template_mol(molecule_new, mol_objects,args):

    if not mol_objects:
        return True
    
    #check if molecule also exixts in the mol_objects
    for mol in mol_objects:
        rms = get_conf_RMS(mol, molecule_new, -1, -1, args.heavyonly, args.max_matches_RMSD)
        if rms < 0.5:
            return False
    return True
