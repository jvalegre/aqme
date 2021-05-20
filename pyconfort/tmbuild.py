#####################################################.
#        This file stores all the functions         #
#    used in template based conformer generation    #
#####################################################.

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdDistGeom, rdMolAlign
from pyconfort.utils import get_conf_RMS

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
def two_embed(molecule_embed, molecule, mol_template, 
              neighbours, name_input, args, log):
    
    mol_template.GetAtomWithIdx(0).setAtomicNum(neighbours[0].GetAtomicNum())
    mol_template.GetAtomWithIdx(1).setAtomicNum(neighbours[1].GetAtomicNum())
    mol_template.GetAtomWithIdx(2).setAtomicNum(53)

    #assigning and embedding onto the core
    mol_obj, coord_map, alg_map,ci = template_embed_optimize(molecule_embed,
                                                             molecule,
                                                             mol_template,
                                                             args, log)
    if ci >= 0: #writing to mol_object file
        return [mol_obj], [name_input], [coord_map], [alg_map], [mol_template]

    return [], [], [], [], []

# GET THE TRIGONAL PLANAR GEOMETRY
def three_embed(molecule_embed,molecule,mol_template,neighbours,name_input,args,log):

    mol_template.GetAtomWithIdx(0).setAtomicNum(53)
    mol_template.GetAtomWithIdx(1).setAtomicNum(neighbours[0].GetAtomicNum())
    mol_template.GetAtomWithIdx(2).setAtomicNum(neighbours[1].GetAtomicNum())
    mol_template.GetAtomWithIdx(3).setAtomicNum(neighbours[2].GetAtomicNum())

    #assigning and embedding onto the core
    mol_obj, coord_map, alg_map,ci = template_embed_optimize(molecule_embed,
                                                             molecule,
                                                             mol_template,
                                                             args, log)
    if ci>=0: #writing to mol_object file
        return [mol_obj], [name_input], [coord_map], [alg_map], [mol_template]

    return [], [], [], [], []

# GET THE SQUAREPLANAR GEOMETRY
def four_embed(molecule_embed,molecule,mol_1,neighbours,name_input,args,log):

    mol_objects = []
    name_return = []
    coord_Map = []
    alg_Map = []
    mol_template = []
    
    # Remove F atoms from the template
    mol_template = Chem.RWMol(mol_template)
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
        mol_obj, coord_map, alg_map,ci = template_embed_optimize(molecule_embed,
                                                                 molecule,
                                                                 mol_1,
                                                                 args, log)

        if ci >= 0:
            check = filter_template_mol(mol_obj, mol_objects,args,log)
            if check:
                #writing to mol_object file
                name = f'{name_input.split()[0]}_{i}'
                mol_objects.append(mol_obj)
                name_return.append(name)
                coord_Map.append(coord_map)
                alg_Map.append(alg_map)
                mol_template.append(mol_1)

    return mol_objects, name_return, coord_Map, alg_Map, mol_template

# GET THE SQUAREPYRAMIDAL GEOMETRY
def five_embed(molecule_embed,molecule,mol_1,neighbours,name_input,args,log):
    mol_objects,name_return,coord_Map,alg_Map,mol_template = [],[],[],[],[]
    
    #fifteen cases for square pyrimidal
    counter = 0
    for name_1 in range(5):
        for name_2 in range(3):

            # assigning order of replacement for the top
            if name_1 == 0:
                k = 4
            elif name_1== 1:
                k = 3
            elif name_1 == 2:
                k = 2
            elif name_1== 3:
                k = 1
            elif name_1 == 4:
                k = 0

            # assigning order of replacement for the plane
            if name_2 == 0 and k == 4:
                j = [1,2,3]
            elif name_2 == 1 and k == 4:
                j = [2,3,1]
            elif name_2 == 2 and k == 4:
                j = [3,1,2]

            # assigning order of replacement for the plane
            if name_2 == 0 and k == 3:
                j = [1,2,4]
            elif name_2 == 1 and k == 3:
                j = [2,4,1]
            elif name_2 == 2 and k == 3:
                j = [4,1,2]

            # assigning order of replacement for the plane
            if name_2 == 0 and k == 2:
                j = [1,4,3]
            elif name_2 == 1 and k == 2:
                j = [4,3,1]
            elif name_2 == 2 and k == 2:
                j = [4,1,3]

            # assigning order of replacement for the plane
            if name_2 == 0 and k == 1:
                j = [4,2,3]
            elif name_2 == 1 and k == 1:
                j = [2,3,4]
            elif name_2 == 2 and k == 1:
                j = [3,4,2]

            # assigning order of replacement for the plane
            if name_2 == 0 and k == 0:
                j = [1,2,3]
            elif name_2 == 1 and k == 0:
                j = [2,3,1]
            elif name_2 == 2 and k == 0:
                j = [3,1,2]

            #checking for same atom neighbours and assigning in the templates for all mols in suppl!
            for atom in mol_1.GetAtoms():
                if atom.GetIdx()  == 5:
                    atom.SetAtomicNum(14)
                    atom.SetFormalCharge(1)
                if atom.GetIdx()  == 1:
                    if k!= 0:
                        atom.SetAtomicNum(neighbours[0].GetAtomicNum())
                    elif k == 0:
                        atom.SetAtomicNum(neighbours[4].GetAtomicNum())
                elif atom.GetIdx()  == 2:
                    atom.SetAtomicNum(neighbours[j[0]].GetAtomicNum())
                elif atom.GetIdx()  == 3:
                    atom.SetAtomicNum(neighbours[j[1]].GetAtomicNum())
                elif atom.GetIdx()  == 4:
                    atom.SetAtomicNum(neighbours[j[2]].GetAtomicNum())
                elif atom.GetIdx()  == 0:
                    if k!= 0:
                        atom.SetAtomicNum(neighbours[k].GetAtomicNum())
                    elif k == 0:
                        atom.SetAtomicNum(neighbours[0].GetAtomicNum())

            #assigning and embedding onto the core
            molecule_new, coordMap, algMap,ci = template_embed_optimize(molecule_embed,molecule,mol_1,args,log)
            if ci>=0:
                check=filter_template_mol(molecule_new, mol_objects,args,log)
                if check:
                    name = str(counter)
                    #writing to mol_object file
                    name_final = name_input +'_'+name
                    mol_objects.append(molecule_new)
                    name_return.append(name_final)
                    coord_Map.append(coordMap)
                    alg_Map.append(algMap)
                    mol_template.append(mol_1)
            else:
                pass
            counter += 1
    return mol_objects, name_return, coord_Map, alg_Map, mol_template

# TEMPLATE GENERATION FOR SQUAREPLANAR AND squarepyramidal
def template_embed(molecule,temp,name_input,args,log):
    neighbours = calc_neighbours(molecule,args)
    number_of_neighbours = len(neighbours)

    mol_objects,name_return,coord_Map,alg_Map,mol_template = [],[],[],[],[]

    # RAUL: Is there more than one? In that case what's the purpose of 
    #       mol_objects, name_return, coord_Map, alg_Map, mol_template?
    for mol_1 in temp:
        if number_of_neighbours == 2:
            mol_objects, name_return, coord_Map, alg_Map, mol_template = two_embed(molecule,molecule,mol_1,neighbours,name_input,args,log)

        elif number_of_neighbours == 3:
            mol_objects, name_return, coord_Map, alg_Map, mol_template = three_embed(molecule,molecule,mol_1,neighbours,name_input,args,log)

        elif number_of_neighbours == 4:
            mol_objects, name_return, coord_Map, alg_Map, mol_template = four_embed(molecule,molecule,mol_1,neighbours,name_input,args,log)

        elif number_of_neighbours == 5:
            mol_objects, name_return, coord_Map, alg_Map, mol_template = five_embed(molecule,molecule,mol_1,neighbours,name_input,args,log)

    return mol_objects, name_return, coord_Map, alg_Map, mol_template

# TEMPLATE EMBED OPTIMIZE
def template_embed_optimize(molecule_embed,molecule,mol_1,args,log):

    #assigning and embedding onto the core
    num_atom_match = molecule_embed.GetSubstructMatch(mol_1)

    #add H's to molecule
    molecule_embed = Chem.AddHs(molecule_embed)

    #definition of coordmap, the coreconfID(the firstone =-1)
    coordMap = {}
    coreConfId=-1
    randomseed=-1
    force_constant=10000

    # This part selects which atoms from molecule are the atoms of the core
    try:
        coreConf = mol_1.GetConformer(coreConfId)
    except:
        pass
    for k, idxI in enumerate(num_atom_match):
        core_mol_1 = coreConf.GetAtomPosition(k)
        coordMap[idxI] = core_mol_1

    ci = rdDistGeom.EmbedMolecule(molecule_embed, coordMap=coordMap, randomSeed=randomseed)

    if ci < 0:
        log.write('Could not embed molecule.')
        coordMap = None
        algMap = None

    if ci >= 0:
        GetFF = Chem.UFFGetMoleculeForceField(molecule_embed,confId=-1)

        #algin molecule to the core
        algMap = [(k, l) for l, k in enumerate(num_atom_match)]

        for k, idxI in enumerate(num_atom_match):
            for l in range(k + 1, len(num_atom_match)):
                idxJ = num_atom_match[l]
                d = coordMap[idxI].Distance(coordMap[idxJ])
                GetFF.AddDistanceConstraint(idxI, idxJ, d, d, force_constant)
        GetFF.Initialize()
        GetFF.Minimize(maxIts=args.opt_steps_RDKit)
        # rotate the embedded conformation onto the core_mol:
        rdMolAlign.AlignMol(molecule_embed, mol_1, atomMap=algMap,reflect=True,maxIters=100)

    return molecule_embed, coordMap, algMap, ci

# FILTER FOR REMOVING MOLS IF LIGANDS ARE THE SAME
def filter_template_mol(molecule_new, mol_objects,args,log):
    if len(mol_objects) ==0:
        check = True
    else:
        check = True
        #check if molecule also exixts in the mol_objects
        for mol in mol_objects:
            rms = get_conf_RMS(mol, molecule_new, -1, -1, args.heavyonly, args.max_matches_RMSD)
            if rms < 0.5:
                check = False
                break
    return check
