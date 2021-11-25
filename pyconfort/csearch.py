#####################################################.
#            This file stores all the functions     #
#               used in conformer generation        #
#####################################################.

import math
import os
import sys
import subprocess
import time
from pathlib import Path
import concurrent.futures as futures     # RAUL: This is for the main

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem

from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, Lipinski
from progress.bar import IncrementalBar  # RAUL: This is for the main

from pyconfort.filter import (filters, ewin_filter, 
                              pre_E_filter, RMSD_and_E_filter)
from pyconfort.tmbuild import template_embed
from pyconfort.cmin import rules_get_charge, substituted_mol
from pyconfort.fullmonte import (generating_conformations_fullmonte, 
                                 minimize_rdkit_energy, realign_mol)
from pyconfort.utils import Logger, set_metal_atomic_number, com_2_xyz_2_sdf

SUPPORTED_INPUTS = ['.smi', '.sdf', '.cdx', 
                    '.csv', '.com', '.gjf', 
                    '.mol', '.mol2','.xyz',
                    '.txt','.yaml','.yml',
                    '.rtf']


def mol_from_sdf_or_mol_or_mol2(input):
    """
    mol from sdf

    Parameters
    ----------
    input : str
        path to a .sdf .mol or .mol2 file

    Returns
    -------
    tuple of lists?
        suppl, IDs, charges
    """
    filename = os.path.splitext(input)[0]
    extension = os.path.splitext(input)[1]

    if extension =='.sdf':
        suppl = Chem.SDMolSupplier(input, removeHs=False)
    elif extension =='.mol':
        suppl = Chem.MolFromMolFile(input, removeHs=False)
    elif extension =='.mol2':
        suppl = Chem.MolFromMol2File(input, removeHs=False)

    IDs,charges = [],[]

    with open(input,"r") as F: 
        lines = F.readlines()

    molecule_count = 0
    for i, line in enumerate(lines):
        if line.find('>  <ID>') > -1:
            ID = lines[i+1].split()[0]
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

    if len(IDs) == 0:
        if extension == '.sdf':
            for i in range(len(suppl)):
                IDs.append(f'{filename}_{i}')
        else:
            IDs.append(filename)
    if len(charges) == 0:
        if extension == '.sdf':
            for _ in suppl:
                charges.append(0)
        else:
            charges.append(0)
    return suppl, IDs, charges

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

def getDihedralMatches_v2(mol, heavy,log): #New version using openbabel
    # If this version is selected, bring the import to the top
    import pybel # from openbabel import pybel for openbabel>=3.0.0
    AtomInTripleBond = '$(*#*)'
    TerminalAtom = 'D1'
    CF3 = '$(C(F)(F)F)'
    CCl3 = '$(C(Cl)(Cl)Cl)'
    CBr3 = '$(C(Br)(Br)Br)'
    tBut = '$(C([CH3])([CH3])[CH3])'
    # A 3-bonded C with a double bond to (N, O or S) 
    # singlgy bonded to not ring bonded to a non-terminal N,O or S.
    CD3_1d = '$([CD3](=[N,O,S])-!@[#7,O,S!D1])' 
    CD3_1r = '$([#7,O,S!D1]-!@[CD3]=[N,O,S])' # Backwards version
    # A 3-bonded C with a double bond to (N+) 
    # singlgy bonded to not ring bonded to Any non-terminal N
    CD3_2d = '$([CD3](=[N+])-!@[#7!D1])'
    CD3_2r = '$([#7!D1]-!@[CD3]=[N+])' # Backwards version 
    Atom1 = '*'
    Atom2 =  f'!{AtomInTripleBond}&!{TerminalAtom}'
    Atom2 += f'&!{CF3}&!{CCl3}&!{CBr3}&!{tBut}'
    Atom2 += f'&!{CD3_1d}&!{CD3_1r}&!{CD3_2d}&!{CD3_2r}'
    Atom3 =  f'!{AtomInTripleBond}&!{TerminalAtom}'
    Atom3 += f'&!{CF3}&!{CCl3}&!{CBr3}&!{tBut}'
    Atom4 = '*'
    pattern = f'{Atom1}~[{Atom2}]-!@[{Atom3}]~{Atom4}'
    smarts = pybel.Smarts(pattern)
    matches = smarts.findall(mol)

    #these are all sets of 4 atoms, uniquify by middle two
    H_atoms = set(pybel.Smarts('[#1]').findall(mol))
    C_atoms = set(pybel.Smarts('[#6]').findall(mol))
    uniqmatches = []
    seen = set()
    for (a,b,c,d) in matches:
        if (b,c) not in seen and (c,b) not in seen:
            if heavy:
                if a not in H_atoms and d not in H_atoms:
                    seen.add((b,c))
                    uniqmatches.append((a,b,c,d))
            if not heavy:
                # So what if a == 'H' and b == 'C'? is that valid ¿? 
                if c not in C_atoms or d not in H_atoms:
                    seen.add((b,c))
                    uniqmatches.append((a,b,c,d))

#creation of csv for csearch
def creation_of_dup_csv(csearch,cmin):
    """
    Generates a pandas.DataFrame object with the appropiate columns for the 
    conformational search and the minimization. 

    Parameters
    ----------
    csearch : str
        Conformational search method. Current valid methods are: 
        ['rdkit','fullmonte','summ']
    cmin : str
        Minimization method. Current valid methods are: 
        ['xtb','ani']

    Returns
    -------
    pandas.DataFrame
    """
    # Boolean aliases from args
    is_rdkit = csearch=='rdkit'
    is_fullmonte = csearch=='fullmonte'
    is_summ = csearch=='summ'
    is_xtb = cmin == 'xtb'
    is_ani = cmin == 'ani'
    
    # column blocks definitions
    base_columns = ['Molecule',
                    'RDKit-Initial-samples',
                    'RDKit-energy-window',
                    'RDKit-initial_energy_threshold',
                    'RDKit-RMSD-and-energy-duplicates',
                    'RDKit-Unique-conformers']
    xtb_columns = ['xTB-Initial-samples',
                   'xTB-energy-window',
                   'xTB-initial_energy_threshold',
                   'xTB-RMSD-and-energy-duplicates',
                   'xTB-Unique-conformers']
    ANI_columns = ['ANI-Initial-samples',
                   'ANI-energy-window',
                   'ANI-initial_energy_threshold',
                   'ANI-RMSD-and-energy-duplicates',
                   'ANI-Unique-conformers']
    end_columns_no_min = ['CSEARCH time (seconds)',
                          'Overall charge']
    end_columns_min = ['CSEARCH time (seconds)',
                       'CMIN time (seconds)',
                       'Overall charge']
    fullmonte_columns = ['FullMonte-Unique-conformers',] 
                        #'FullMonte-conformers',
                        #'FullMonte-energy-window',
                        #'FullMonte-initial_energy_threshold',
                        #'FullMonte-RMSD-and-energy-duplicates']
    summ_columns = ['summ-conformers',
                    'summ-energy-window',
                    'summ-initial_energy_threshold',
                    'summ-RMSD-and-energy-duplicates',
                    'summ-Unique-conformers']
    
    # Check Conformer Search method
    if is_rdkit:
        columns = base_columns
    elif is_fullmonte:
        columns = base_columns + fullmonte_columns
    elif is_summ:
        columns = base_columns + summ_columns
    else:
        return None
    
    # Check Minimization Method
    if is_ani: 
        columns += ANI_columns
    if is_xtb:  # is_ani and is_xtb will not happen, but this is what was written
        columns += xtb_columns
    if is_ani or is_xtb: 
        columns += end_columns_min
    else:
        columns += end_columns_no_min
    
    return pd.DataFrame(columns=columns)

def clean_args(args,ori_ff,mol,ori_charge):                                # RAUL: I hope this function does not survive the clean-ups
    """
    returns the arguments to their original value after each calculation

    Parameters
    ----------
    args : argparse.args
        [description]
    ori_ff : [type]
        original forcefield
    mol : rdkit.Chem.Mol
        [description]
    ori_charge : [type]
        original charge
    """
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

def compute_confs(w_dir_initial, mol, name, args,i):
    """
    function to start conf generation

    Parameters
    ----------
    w_dir_initial : [type]
        [description]
    mol : rdkit.Chem.Mol
        [description]
    name : [type]
        [description]
    i : [type]
        [description]

    Returns
    -------
    pandas.Dataframe
        total_data
    """
    
    csearch_dir = Path(w_dir_initial) / 'CSEARCH'
    dat_dir = csearch_dir / 'dat_files'
    dat_dir.mkdir(parents=True, exist_ok=True)
    
    log = Logger(dat_dir/name, args.output_name)
    # Converts each line to a rdkit mol object
    if args.verbose:
        log.write(f"   -> Input Molecule {i} is {Chem.MolToSmiles(mol)}")

    if args.metal_complex:
        for _ in args.metal:
            args.metal_idx.append(None)
            args.complex_coord.append(None)
            args.metal_sym.append(None)

        args.mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(mol,args)

        # get pre-determined geometries for metal complexes
        accepted_complex_types = ['squareplanar',
                                  'squarepyramidal',
                                  'linear',
                                  'trigonalplanar']
        if args.complex_type in accepted_complex_types:
            count_metals = 0
            for metal_idx_ind in args.metal_idx:
                if metal_idx_ind is not None:
                    count_metals += 1
            if count_metals == 1:
                os.chdir(w_dir_initial)
                template_kwargs = dict() 
                template_kwargs['complex_type'] = args.complex_type
                template_kwargs['metal_idx'] = args.metal_idx
                template_kwargs['maxsteps'] = args.opt_steps_RDKit
                template_kwargs['heavyonly'] = args.heavyonly
                template_kwargs['maxmatches'] = args.max_matches_RMSD
                items = template_embed(mol,name,log,**template_kwargs)

                total_data = creation_of_dup_csv(args.CSEARCH,args.CMIN)
                
                for mol_obj, name_in, coord_map, alg_map, template in zip(*items):
                    data = conformer_generation(mol_obj, name_in, args, log,
                                                coord_map, alg_map, template)
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

def conformer_generation(mol,name,args,log,coord_Map=None,alg_Map=None,mol_template=None):
    """
    FUNCTION WORKING WITH MOL OBJECT TO CREATE CONFORMERS

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        [description]
    name : [type]
        [description]
    args : argparse.args
        [description]
    log : [type]
        [description]
    coord_Map : [type], optional
        [description], by default None
    alg_Map : [type], optional
        [description], by default None
    mol_template : [type], optional
        [description], by default None

    Returns
    -------
    pd.Dataframe
        dup_data
    """
    dup_data = creation_of_dup_csv(args.CSEARCH,args.CMIN)
    
    dup_data_idx = 0
    start_time = time.time()
    valid_structure = filters(mol,log,args.max_MolWt,args.verbose)
    if valid_structure:
        if args.verbose:
            log.write(f"\n   ----- {name} -----")
        try:
            # the conformational search for RDKit
            status,update_to_rdkit = summ_search(mol, name,args,log,dup_data,dup_data_idx,coord_Map,alg_Map,mol_template)
            dup_data.at[dup_data_idx, 'status'] = status
            dup_data.at[dup_data_idx, 'update_to_rdkit'] = update_to_rdkit
        except (KeyboardInterrupt, SystemExit): # RAUL: This try-except is useless.
            raise
    else:
        log.write("\nx  ERROR: The structure is not valid")

    if args.time:
        n_seconds = (round(time.time() - start_time,2))
        log.write(f"\n Execution time CSEARCH: {n_seconds} seconds")
        dup_data.at[dup_data_idx, 'CSEARCH time (seconds)'] = n_seconds
    return dup_data

def auto_sampling(mult_factor,mol,args,log):
    """
    DETECTS INITIAL NUMBER OF SAMPLES AUTOMATICALLY

    Parameters
    ----------
    mult_factor : [type]
        [description]
    mol : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]

    Returns
    -------
    int
        auto_samples
    """
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

def genConformer_r(mol, conf, i, matches, degree, sdwriter,args,name,log,update_to_rdkit,coord_Map,alg_Map, mol_template):
    """
    IF NOT USING DIHEDRALS, THIS REPLACES I BACK TO THE METAL WHEN METAL = TRUE
    AND WRITES THE RDKIT SDF FILES. WITH DIHEDRALS, IT OPTIMIZES THE ROTAMERS

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        [description]
    conf : [type]
        [description]
    i : [type]
        [description]
    matches : [type]
        [description]
    degree : [type]
        [description]
    sdwriter : [type]
        [description]
    args : [type]
        [description]
    name : [type]
        [description]
    log : [type]
        [description]
    update_to_rdkit : [type]
        [description]
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]

    Returns
    -------
    int
        total number of conformers generated
    """
    if i >= len(matches): # base case, torsions should be set in conf
        #setting the metal back instead of I
        if args.metal_complex and (args.CSEARCH=='rdkit' or update_to_rdkit):
            if coord_Map is None and alg_Map is None and mol_template is None:
                energy = minimize_rdkit_energy(mol,conf,log,args.ff,args.opt_steps_RDKit)
            else:
                mol,energy = realign_mol(mol,conf,coord_Map, alg_Map, mol_template,args.opt_steps_RDKit)
            mol.SetProp('Energy',str(energy))
            set_metal_atomic_number(mol,args.metal_idx,args.metal_sym)
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

def embed_conf(mol,initial_confs,args,log,coord_Map,alg_Map, mol_template):
    """
    function to embed conformers

    Parameters
    ----------
    mol : [type]
        [description]
    initial_confs : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]

    Returns
    -------
    list
        cids
    """
    is_sdf_mol_or_mol2 = os.path.splitext(args.input)[1] in ['.sdf','.mol','.mol2']
    
    if is_sdf_mol_or_mol2:
            Chem.AssignStereochemistryFrom3D(mol)

    embed_kwargs = dict()
    embed_kwargs['ignoreSmoothingFailures'] = True
    embed_kwargs['randomSeed'] = args.seed
    embed_kwargs['numThreads'] = 0
    
    if (coord_Map,alg_Map,mol_template) != (None,None,None):
        embed_kwargs['coordMap'] = coord_Map
    cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

    if len(cids) <= 1 and initial_confs != 1:
        log.write("o  Normal RDKit embeding process failed, trying to " \
                  "generate conformers with random coordinates " \
                  f"(with {str(initial_confs)} possibilities)")
        embed_kwargs['useRandomCoords'] = True
        embed_kwargs['boxSizeMult'] = 10.0
        embed_kwargs['numZeroFail'] = 1000
        cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)
    
    if is_sdf_mol_or_mol2:
        #preserving AssignStereochemistryFrom3D
        for cid in cids:
            Chem.AssignAtomChiralTagsFromStructure(mol,confId=cid)
    
    return cids

def min_and_E_calc(mol,cids,args,log,coord_Map,alg_Map,mol_template):
    """
    minimization and E calculation with RDKit after embeding

    Parameters
    ----------
    mol : [type]
        [description]
    cids : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]

    Returns
    -------
    outmols,cenergy
        outmols is gonna be a list containing "initial_confs" mol objects 
        with "initial_confs" conformers. We do this to SetProp 
        (Name and Energy) to the different conformers
        and log.write in the SDF file. At the end, since all the mol 
        objects has the same conformers, but the energies are different, 
        we can log.write conformers to SDF files with the energies of the 
        parent mol objects. We measured the computing time and it's the 
        same as using only 1 parent mol object with 10 conformers, but we 
        couldn'temp SetProp correctly.
    """
    
    cenergy,outmols = [],[]
    #bar = IncrementalBar('o  Minimizing', max = len(cids))
    for _, conf in enumerate(cids):
        if coord_Map is None and alg_Map is None and mol_template is None:
            energy = minimize_rdkit_energy(mol,conf,log,args.ff,args.opt_steps_RDKit)
        else: # id template realign before doing calculations
            mol,energy = realign_mol(mol,conf,coord_Map, alg_Map, mol_template,args.opt_steps_RDKit)
        cenergy.append(energy)
        pmol = PropertyMol.PropertyMol(mol)
        outmols.append(pmol)
        #bar.next()
    #bar.finish()
    return outmols,cenergy

def min_after_embed(mol,cids,name,initial_confs,rotmatches,dup_data,dup_data_idx,sdwriter,args,log,update_to_rdkit,coord_Map,alg_Map, mol_template):
    """
    minimizes, gets the energy and filters RDKit conformers after embeding

    Parameters
    ----------
    mol : [type]
        [description]
    cids : [type]
        [description]
    name : [type]
        [description]
    initial_confs : [type]
        [description]
    rotmatches : [type]
        [description]
    dup_data : [type]
        [description]
    dup_data_idx : [type]
        [description]
    sdwriter : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]
    update_to_rdkit : [type]
        [description]
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]

    Returns
    -------
    int
        status
    """
    
    # gets optimized mol objects and energies
    outmols,cenergy = min_and_E_calc(mol,cids,args,log,coord_Map,alg_Map,mol_template)
    
    # writing charges after RDKit
    if os.path.splitext(args.input)[1] == '.cdx' or os.path.splitext(args.input)[1] == '.smi' or os.path.splitext(args.input)[1] == '.csv':
        args.charge = rules_get_charge(mol,args)
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
    sortedcids_rdkit = ewin_filter(sorted_all_cids,cenergy,args,dup_data,dup_data_idx,log,'rdkit',args.ewin_csearch)

    # pre-filter based on energy only
    selectedcids_initial_rdkit = pre_E_filter(sortedcids_rdkit,cenergy,dup_data,dup_data_idx,log,'rdkit',args.initial_energy_threshold,args.verbose)

    # filter based on energy and RMSD
    selectedcids_rdkit = RMSD_and_E_filter(outmols,selectedcids_initial_rdkit,cenergy,args,dup_data,dup_data_idx,log,'rdkit')

    if args.CSEARCH=='summ' or args.CSEARCH=='rdkit':
        # now exhaustively drive torsions of selected conformers
        n_confs = int(len(selectedcids_rdkit) * (360 / args.degree) ** len(rotmatches))
        if args.verbose and len(rotmatches) != 0:
            log.write("\n\no  Systematic generation of "+ str(n_confs)+ " confomers")
            #bar = IncrementalBar('o  Generating conformations based on dihedral rotation', max = len(selectedcids_rdkit))
        # else:
        #     bar = IncrementalBar('o  Writing unique conformers into an sdf file', max = len(selectedcids_rdkit))

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

def rdkit_to_sdf(mol, name,args,log,dup_data,dup_data_idx, coord_Map, alg_Map, mol_template):
    """
    conversion from rdkit to sdf

    Parameters
    ----------
    mol : [type]
        [description]
    name : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]
    dup_data : [type]
        [description]
    dup_data_idx : [type]
        [description]
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]

    Returns
    -------
    tuple
        status,rotmatches,update_to_rdkit
    """

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

    dup_data.at[dup_data_idx, 'RDKit-Initial-samples'] = initial_confs
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

def dihedral_filter_and_sdf(name,args,log,dup_data,dup_data_idx,coord_Map, alg_Map, mol_template):
    """
    filtering after dihedral scan to sdf

    Parameters
    ----------
    name : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]
    dup_data : [type]
        [description]
    dup_data_idx : [type]
        [description]
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]

    Returns
    -------
    int
        status (job I guess?)
    """
    rotated_energy = []

    rdmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)

    if rdmols is None:
        log.write("Could not open "+ name+args.output)
        sys.exit(-1)

    for i, rd_mol_i in enumerate(rdmols):
        if coord_Map is None and alg_Map is None and mol_template is None:
            energy = minimize_rdkit_energy(rd_mol_i,-1,log,args.ff,args.opt_steps_RDKit)
        else:
            rd_mol_i,energy = realign_mol(rd_mol_i,-1,coord_Map, alg_Map, mol_template,args.opt_steps_RDKit)
        rotated_energy.append(energy)

    rotated_cids = list(range(len(rdmols)))
    sorted_rotated_cids = sorted(rotated_cids, key = lambda cid: rotated_energy[cid])

    # filter based on energy window ewin_csearch
    sortedcids_rotated = ewin_filter(sorted_rotated_cids,rotated_energy,args,dup_data,dup_data_idx,log,'summ',args.ewin_csearch)
    # pre-filter based on energy only
    selectedcids_initial_rotated = pre_E_filter(sortedcids_rotated,rotated_energy,dup_data,dup_data_idx,log,'summ',args.initial_energy_threshold,args.verbose)
    # filter based on energy and RMSD
    selectedcids_rotated = RMSD_and_E_filter(rdmols,selectedcids_initial_rotated,rotated_energy,args,dup_data,dup_data_idx,log,'summ')

    sdwriter_rd = Chem.SDWriter(name+'_'+'summ'+args.output)
    for i, cid in enumerate(selectedcids_rotated):
        mol_rd = Chem.RWMol(rdmols[cid])
        mol_rd.SetProp('_Name',rdmols[cid].GetProp('_Name')+' '+str(i))
        mol_rd.SetProp('Energy',str(rotated_energy[cid]))
        if args.metal_complex:
            set_metal_atomic_number(mol_rd,args.metal_idx,args.metal_sym)
        sdwriter_rd.write(mol_rd)
    sdwriter_rd.close()
    status = 1
    return status

def summ_search(mol, name,args,log,dup_data,dup_data_idx, coord_Map = None, alg_Map=None, mol_template=None):
    """
    EMBEDS, OPTIMIZES AND FILTERS RDKIT CONFORMERS

    Parameters
    ----------
    mol : [type]
        [description]
    name : [type]
        [description]
    args : [type]
        [description]
    log : [type]
        [description]
    dup_data : [type]
        [description]
    dup_data_idx : [type]
        [description]
    coord_Map : [type], optional
        [description], by default None
    alg_Map : [type], optional
        [description], by default None
    mol_template : [type], optional
        [description], by default None

    Returns
    -------
    tuple
        status, update_to_rdkit
    """
    
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


# MAIN FUNCTION 

# main function to generate conformers
def csearch_main_v2(w_dir_initial,args,log_overall):

    file_format = os.path.splitext(args.input)[1]

    # Checks
    if file_format not in SUPPORTED_INPUTS:
        log_overall.write("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!")
        sys.exit()
    if not os.path.exists(args.input):
        log_overall.write("\nx  INPUT FILE NOT FOUND!")
        sys.exit()

    # sets up the chosen force field (this fixes some problems in case MMFF is replaced by UFF)
    ori_ff = args.ff
    ori_charge = args.charge_default

    # if large system increase stack size
    # if args.STACKSIZE != '1G':
    #     os.environ['OMP_STACKSIZE'] = args.STACKSIZE
    smi_derivatives = ['.smi', '.txt', '.yaml', '.yml', '.rtf']
    Extension2inputgen = dict()
    for key in smi_derivatives:
        Extension2inputgen[key] = prepare_smiles_files
    Extension2inputgen['.csv'] = prepare_csv_files
    Extension2inputgen['.cdx'] = prepare_cdx_files
    Extension2inputgen['.gjf'] = prepare_gaussian_files
    Extension2inputgen['.com'] = prepare_gaussian_files
    Extension2inputgen['.xyz'] = prepare_gaussian_files
    Extension2inputgen['.sdf'] = prepare_sdf_files
    Extension2inputgen['.mol'] = prepare_mol_files
    Extension2inputgen['.mol2'] = prepare_mol_files

    with futures.ProcessPoolExecutor(max_workers=args.cpus) as executor:
        # Submit a set of asynchronous jobs
        jobs = []
        count_mol = 0
        
        # Prepare the Jobs
        prepare_function = Extension2inputgen[file_format]
        job_inputs = prepare_function(args,ori_ff,ori_charge,w_dir_initial)
        
        # Submit the Jobs
        if file_format in smi_derivatives:
            for job_input in job_inputs:
                try:
                    compute_confs,w_dir_initial,mol,name,args,i = job_input
                    job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
                    jobs.append(job)
                    count_mol += 1
                    # compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)
                except AttributeError:                                          # TODO I NEED CONFIRMATION OF WHEN DOES THE ATTRIBUTE ERROR APPEAR
                    smi = Chem.MolToSmiles(mol)
                    log_overall.write(f"\nx  Wrong SMILES string ({smi}) found (not compatible with RDKit or ANI/xTB if selected)! This compound will be omitted\n")
        else:     
            for job_input in job_inputs: 
                compute_confs,w_dir_initial,mol,name,args,i = job_input
                job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
                jobs.append(job)
                count_mol += 1

        final_dup_data = creation_of_dup_csv(args)
        bar = IncrementalBar('o  Number of finished jobs from CSEARCH', max = count_mol)
        # Process the job results (in submission order) and save the conformers.
        for i,job in enumerate(jobs):
            total_data = job.result()
            frames = [final_dup_data, total_data]
            final_dup_data = pd.concat(frames,ignore_index=True,sort=True)
            bar.next()
        bar.finish()

        # removing temporary files
        temp_files = ['gfn2.out', 'xTB_opt.traj', 'ANI1_opt.traj', 'wbo', 'xtbrestart','ase.opt','xtb.opt','gfnff_topo']
        for file in temp_files:
            if os.path.exists(file):
                os.remove(file)

    return final_dup_data

def prepare_smiles_files(args,ori_ff,ori_charge,w_dir_initial): 
    with open(args.input) as smifile:
        lines = [line for line in smifile if line.strip()]
    job_inputs = []
    for i, line in enumerate(lines):
        mol,name,args= prepare_smiles_from_line(line,i,args,ori_ff,ori_charge)
        job_inp = (compute_confs,w_dir_initial,mol,name,args,i)
        job_inputs.append(job_inp)
    return job_inputs
def prepare_smiles_from_line(line,i,args,ori_ff,ori_charge):
    toks = line.split()
    #editing part
    smi = toks[0]
    smi = check_for_pieces(smi)
    mol = Chem.MolFromSmiles(smi)
    clean_args(args,ori_ff,mol,ori_charge)                          # I assume no AttributeError 
    if args.charge_default == 'auto':                               # I assume no AttributeError 
        if not args.metal_complex:                                  # I assume no AttributeError 
            args.charge_default = check_charge_smi(smi)             # I assume no AttributeError 
    if args.prefix == 'None':                                       # I assume no AttributeError 
        name = ''.join(toks[1:])                                    # I assume no AttributeError 
    else:                                                           # I assume no AttributeError 
        name = f"{args.prefix}_{i}_{''.join(toks[1:])}"             # I assume no AttributeError
    return mol,name,args

def prepare_csv_files(args,ori_ff,ori_charge,w_dir_initial): 
    csv_smiles = pd.read_csv(args.input)
    job_inputs = []
    for i in range(len(csv_smiles)):
        job_inp = generate_mol_from_csv(args,w_dir_initial,ori_charge,ori_ff,csv_smiles,i)
        job_inputs.append(job_inp)

    return job_inputs
def generate_mol_from_csv(args,w_dir_initial,ori_charge,ori_ff,csv_smiles,index):
    #assigning names and smi i  each loop
    smi = csv_smiles.loc[index, 'SMILES']
    pruned_smi = check_for_pieces(smi)
    mol = Chem.MolFromSmiles(pruned_smi)
    clean_args(args,ori_ff,mol,ori_charge)
    if args.charge_default == 'auto':
        if not args.metal_complex:
            args.charge_default = check_charge_smi(pruned_smi)
    if args.prefix == 'None':
        name = csv_smiles.loc[index, 'code_name']
    else:
        name = 'comp_'+str(index)+'_'+csv_smiles.loc[index, 'code_name']
    return (compute_confs,w_dir_initial,mol,name,args,index)

def prepare_cdx_files(args,ori_ff,ori_charge,w_dir_initial):
    #converting to smiles from chemdraw
    molecules = generate_mol_from_cdx(args)
    job_inputs = []
    for i,(smi,mol) in enumerate(molecules):   
        clean_args(args,ori_ff,mol,ori_charge)
        name = f"{args.input.split('.')[0]}_{str(i)}"
        if args.charge_default == 'auto':
            if not args.metal_complex:
                args.charge_default = check_charge_smi(smi)
        job_inputs.append((compute_confs,w_dir_initial,mol,name,args,i))
    return job_inputs
def generate_mol_from_cdx(args):
    cmd_cdx = ['obabel', '-icdx', args.input, '-osmi', '-Ocdx.smi']
    subprocess.call(cmd_cdx)
    with open('cdx.smi',"r") as smifile:
        smi_lines = [line for line in smifile]
    os.remove('cdx.smi')
    molecules = []
    for smi in smi_lines:
        pruned_smi = check_for_pieces(smi)
        molecule = Chem.MolFromSmiles(pruned_smi)
        molecules.append((pruned_smi,molecule))
    return molecules

def prepare_gaussian_files(args,ori_ff,ori_charge,w_dir_initial):
    job_inputs = []
    charge_com = com_2_xyz_2_sdf(args.input,args.default_charge)
    name = os.path.splitext(args.input)[0]
    sdffile = f'{name}.sdf'
    suppl, _, _ = mol_from_sdf_or_mol_or_mol2(sdffile)

    for i,mol in enumerate(suppl):
        clean_args(args,ori_ff,mol,ori_charge)
        if args.charge_default == 'auto':
            args.charge_default = charge_com
        job_inputs.append((compute_confs,w_dir_initial,mol,name,args,i))

def prepare_xyz_files(args,ori_ff,ori_charge,w_dir_initial):
    job_inputs = []
    name = os.path.splitext(args.input)[0]
    sdffile = f'{name}.sdf'
    suppl, _, charge_com_list = mol_from_sdf_or_mol_or_mol2(sdffile)
    charge_com = charge_com_list[0]

    for i,mol in enumerate(suppl):
        clean_args(args,ori_ff,mol,ori_charge)
        if args.charge_default == 'auto':
            args.charge_default = charge_com
        job_inputs.append((compute_confs,w_dir_initial,mol,name,args,i))

def prepare_sdf_file(ori_ff,ori_charge,charge_sdf,w_dir_initial,mol,name,args,i):
    clean_args(args,ori_ff,mol,ori_charge)
    if args.charge_default == 'auto':
        args.charge_default = charge_sdf
    return compute_confs,w_dir_initial,mol,name,args,i
def prepare_sdf_files(args,ori_ff,ori_charge,w_dir_initial):
    suppl, IDs, charges = mol_from_sdf_or_mol_or_mol2(args.input)
    job_inputs = []
    for i,(mol,name,charge_sdf) in enumerate(zip(suppl,IDs,charges)):
        job_inp = prepare_sdf_file(ori_ff,ori_charge,charge_sdf,w_dir_initial,mol,name,args,i)
        job_inputs.append(job_inp)
    return job_inputs

def prepare_mol_file(suppl,name,charge,args,w_dir_initial):
    if args.charge_default == 'auto':
        args.charge_default = charge
    mol = suppl
    return compute_confs,w_dir_initial,mol,name,args,0
def prepare_mol_files(args,ori_ff,ori_charge,w_dir_initial): # The extra variables are for API Consistency.
    suppl, IDs, charges = mol_from_sdf_or_mol_or_mol2(args.input)
    job_inputs = []
    job_inp = prepare_mol_file(suppl,IDs[0],charges[0],args,w_dir_initial)
    job_inputs.append(job_inp)
    return job_inputs
