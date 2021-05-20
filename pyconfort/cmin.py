#####################################################.
#            This file stores all the functions     #
#               used in conformer minimization      #
#####################################################.

import os
import sys
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.PropertyMol import PropertyMol
from rdkit.Geometry import Point3D
from pyconfort.qprep_gaussian import write_confs
from pyconfort.filter import (ewin_filter, pre_E_filter, RMSD_and_E_filter)
from pyconfort.utils import set_metal_atomic_number

hartree_to_kcal = 627.509


def rules_get_charge(mol,args):
    """
    AUTOMATICALLY SETS THE CHARGE FOR METAL COMPLEXES

    Parameters
    ----------
    mol : [type]
        [description]
    args : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    C_group = ['C', 'Se', 'Ge']
    N_group = ['N', 'P', 'As']
    O_group = ['O', 'S', 'Se']
    F_group = ['Cl', 'Br', 'I']

    M_ligands, N_carbenes, bridge_atoms, neighbours = [],[],[],[]
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
                if neighbour.GetTotalValence()== 4:
                    if neighbour.GetSymbol() in C_group:
                        carbene_like = False
                        bridge_ligand = False
                        for inside_neighbour in neighbour.GetNeighbors():
                            if inside_neighbour.GetSymbol() in N_group:
                                if inside_neighbour.GetTotalValence() == 4:
                                    for N_neighbour in inside_neighbour.GetNeighbors():
                                        # this option detects bridge ligands that connect two metals such as M--CN--M
                                        # we use I since the M is still represented as I at this point
                                        if N_neighbour.GetSymbol() == 'I':
                                            bridge_ligand = True
                                            bridge_atoms.append(inside_neighbour.GetIdx())
                                    if not bridge_ligand:
                                        carbene_like = True
                                        N_carbenes.append(inside_neighbour.GetIdx())
                        if not carbene_like:
                            charge_rules[charge_idx] = charge_rules[charge_idx] - 1
                elif neighbour.GetTotalValence()== 3:
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
    else:
        # no update in charge as it is an organic molecule
        return [args.charge_default,]

def substituted_mol(mol,args):
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
    for i,j in zip(range(2,9),range(-3,4)): 
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

    return mol,args.metal_idx,args.complex_coord,args.metal_sym

def rdkit_sdf_read(name, args, log):
    """
    Reads sdf files and stops the execution if the file was not accesible. 

    Parameters
    ----------
    name : [type]
        [description]
    args : argparse.args
        [description]
    log : Logger
        [description]

    Returns
    -------
    list
        rdkit.Chem.Mol objects
    """
    inmols = Chem.SDMolSupplier(name+args.output, removeHs=False)
    if inmols is None:
        log.write(f"Could not open {name}{args.output}")
        sys.exit(-1)
    return inmols

def mult_min(name, args, program,log,dup_data,dup_data_idx):
    """
    READ FILES FOR xTB AND ANI1 OPTIMIZATION, FILTER AND WRITING SDF FILES

    Parameters
    ----------
    name : [type]
        [description]
    args : [type]
        [description]
    program : [type]
        [description]
    log : [type]
        [description]
    dup_data : [type]
        [description]
    dup_data_idx : [type]
        [description]
    """

    if args.verbose:
        if args.CMIN=='xtb':
            if args.xtb_solvent == 'none':
                method = f"xTB ({args.xtb_method}"
            else:
                method = f"xTB ({args.xtb_method} in {args.xtb_solvent})"
        if args.CMIN=='ani':
            method = f"ANI ({args.ani_method})"
        if args.CMIN in ['xtb','ani']: 
            filename = name+args.output
            log.write(f"\n\no  Multiple minimization of {filename} with {method}")

    # bar = IncrementalBar('o  Minimizing', max = len(inmols))

    # read SDF files from RDKit optimization
    inmols = rdkit_sdf_read(name, args, log)

    cenergy, outmols = [],[]

    for mol in inmols:
        # bar.next()
        if mol is not None:

            # optimize this structure and record the energy
            if args.metal_complex:
                args.metal_idx = []
                args.complex_coord = []
                args.metal_sym = []
                # fill the lists with None for every metal in the option
                for _ in args.metal:
                    args.metal_idx.append(None)
                    args.complex_coord.append(None)
                    args.metal_sym.append(None)

                mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(mol,args,log)

            mol,energy,ani_incompatible = optimize(mol, args, program,log,dup_data,dup_data_idx)
            if not ani_incompatible:
                pmol = PropertyMol(mol)
                outmols.append(pmol)
                cenergy.append(energy)

    # if SQM energy exists, overwrite RDKit energies and geometries
    cids = list(range(len(outmols)))
    sorted_all_cids = sorted(cids, key = lambda cid: cenergy[cid])

    name_mol = name.split('_'+args.CSEARCH)[0]

    for cid in sorted_all_cids:
        outmols[cid].SetProp('_Name', outmols[cid].GetProp('_Name') +' '+program)
        outmols[cid].SetProp('Energy', cenergy[cid])


    #writing all conformers to files after minimization
    sdwriter = Chem.SDWriter(f'{name_mol}_{program}_all_confs{args.output}')

    write_all_confs = 0
    for cid in sorted_all_cids:
        sdwriter.write(outmols[cid])
        write_all_confs += 1
    sdwriter.close()

    if args.verbose:
        log.write(f"\no  Writing {str(write_all_confs)} conformers to file {name_mol}_{program}_all_confs{args.output}")

    log.write(f"\n\no  Applying filters to intial conformers after {program} minimization")
    # filter based on energy window ewin_csearch
    sortedcids = ewin_filter(sorted_all_cids,cenergy,args,dup_data,dup_data_idx,log,program,args.ewin_csearch)
    # pre-filter based on energy only
    selectedcids_initial = pre_E_filter(sortedcids,cenergy,dup_data,dup_data_idx,log,program,args.initial_energy_threshold,args.verbose)
    # filter based on energy and RMSD
    selectedcids = RMSD_and_E_filter(outmols,selectedcids_initial,cenergy,args,dup_data,dup_data_idx,log,program)

    if program == 'xtb':
        dup_data.at[dup_data_idx, 'xTB-Initial-samples'] = len(inmols)
    if program == 'ani':
        dup_data.at[dup_data_idx, 'ANI-Initial-samples'] = len(inmols)

    # write the filtered, ordered conformers to external file
    write_confs(outmols, cenergy,selectedcids, name_mol, args, program,log)

############# ASE  AND TORCH AND TORCHANI DEPENDENT FUNCTIONS ##################

log = 'something' # RAUL: When the Logger is refactored This has to be changed appropiatedly 

# imports for xTB and ANI1
try:
    import ase
except (ModuleNotFoundError,AttributeError):
    log.write('\nx  ASE is not installed correctly - xTB and ANI are not available')
    raise ModuleNotFoundError
else:
    import ase.optimize
    from ase.units import Hartree

try:
    import torch
except (ModuleNotFoundError,AttributeError):
    log.write('\nx  TORCH is not installed correctly - xTB and ANI are not available')
    raise ModuleNotFoundError
else:
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
    DEVICE = torch.device('cpu')

# Attempt an XTB import, if it fails log it and mock xtb_calc to delay the 
# system exit until it is used. 
try:
    from xtb.ase.calculator import XTB
except (ModuleNotFoundError,AttributeError):
    log.write('\nx  xTB is not installed correctly - xTB is not available')
    xtb_calc = lambda *x,**y: sys.exit()
else: 
    def xtb_calc(elements,coordinates,args,log,ase_metal,ase_metal_idx):
        """
        Run an xtb optimization and return the optimized coordinates and the energy
        of the molecule. 

        Parameters
        ----------
        elements : [type]
            [description]
        coordinates : [type]
            [description]
        args : argparse.args
            [description]
        log : Logger
            [description]
        ase_metal : [type]
            [description]
        ase_metal_idx : [type]
            [description]

        Returns
        -------
        tuple
            sqm_energy, coordinates
        """
        
        xtb_calculator = XTB(method=args.xtb_method,
                        accuracy=args.xtb_accuracy,
                        electronic_temperature=args.xtb_electronic_temperature,
                        max_iterations=args.xtb_max_iterations,
                        solvent=args.xtb_solvent)

        #define ase molecule using GFN2 Calculator
        ase_molecule = ase.Atoms(elements, 
                                positions=coordinates.tolist()[0],
                                calculator=xtb_calculator)
        
        # Adjust the charge of the metal atoms
        if args.metal_complex:
            extension = os.path.splitext(args.input)[1]
            if extension in ['.csv','.cdx','.smi']:
                for i,atom in enumerate(ase_molecule):
                    if i in ase_metal:
                        ase_charge = args.charge[args.metal_idx.index(ase_metal_idx[ase_metal.index(i)])]
                        # will update only for cdx, smi, and csv formats.
                        atom.charge = ase_charge
                else:
                    atom.charge = args.charge_default
                    if args.verbose:
                        log.write('o  The Overall charge is read from the .com file ')
        
            
        optimizer = ase.optimize.BFGS(ase_molecule, 
                                    trajectory='xTB_opt.traj',
                                    logfile='xtb.opt')

        optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)

        if len(ase.io.Trajectory('xTB_opt.traj', mode='r')) != (args.opt_steps+1):
            species_coords = ase_molecule.get_positions().tolist()
            coordinates = torch.tensor([species_coords], requires_grad=True, device=DEVICE)
        
        # Now let's compute energy:
        xtb_energy = ase_molecule.get_potential_energy()
        sqm_energy = (xtb_energy / Hartree)* hartree_to_kcal

        return sqm_energy, coordinates

# Attempt a torchani import, if it fails log it and mock ani_calc function
# to raise a sys.exit() if called 
try:
    import torchani
except (ModuleNotFoundError,AttributeError):
    log.write('\nx  Torchani is not installed correctly - ANI is not available')
    ani_calc = lambda *x,**y: sys.exit()
else:
    def ani_calc(elements,coordinates,args):
        """
        Run a ANI1 optimization and return the energy and optimized coordinates. 

        Parameters
        ----------
        ase : [type]
            [description]
        torch : [type]
            [description]
        model : [type]
            [description]
        device : [type]
            [description]
        elements : [type]
            [description]
        coordinates : [type]
            [description]
        args : [type]
            [description]

        Returns
        -------
        tuple
            sqm_energy, coordinates
        """
        # Select the model
        ANI_method = args.ani_method
        if ANI_method == 'ANI1x':
            model = torchani.models.ANI1x()
        if ANI_method == 'ANI1ccx':
            model = torchani.models.ANI1ccx()
        if ANI_method == 'ANI2x':
            model = torchani.models.ANI2x()
        if ANI_method == 'ANI2ccx':
            model = torchani.models.ANI2ccx()
        if ANI_method == 'ANI3x':
            model = torchani.models.ANI3x()
        if ANI_method == 'ANI3ccx':
            model = torchani.models.ANI3ccx()
        # A bit Fancier
        #model = getattr(torchani.models,ANI_method)()
        
        species = model.species_to_tensor(elements).to(DEVICE).unsqueeze(0)
        _, ani_energy = model((species, coordinates))

        ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0], calculator=model.ase())

        optimizer = ase.optimize.BFGS(ase_molecule, trajectory='ANI1_opt.traj', logfile='ase.opt' )
        optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
        if len(ase.io.Trajectory('ANI1_opt.traj', mode='r')) != (args.opt_steps+1):
            species_coords = ase_molecule.get_positions().tolist()
            coordinates = torch.tensor([species_coords], requires_grad=True, device=DEVICE)

        # Now let's compute energy:
        _, ani_energy = model((species, coordinates))
        sqm_energy = ani_energy.item() * hartree_to_kcal # Hartree to kcal/mol

        return sqm_energy, coordinates

# xTB AND ANI1 MAIN OPTIMIZATION PROCESS
def optimize(mol, args, program,log,dup_data,dup_data_idx):

    # if large system increase stack size
    if args.STACKSIZE != '1G':
        os.environ['OMP_STACKSIZE'] = args.STACKSIZE

    # removing the Ba atom if NCI complexes
    if args.nci_complex:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'I':
                atom.SetAtomicNum(1)

    if args.metal_complex and not args.CSEARCH=='summ':
        set_metal_atomic_number(mol,args.metal_idx,args.metal_sym)

    elements = ''
    ase_metal = []
    ase_metal_idx = []
    for i,atom in enumerate(mol.GetAtoms()):
        if atom.GetIdx() in args.metal_idx:
            ase_metal.append(i)
            ase_metal_idx.append(atom.GetIdx())
        elements += atom.GetSymbol()

    extension = os.path.splitext(args.input)[1]
    if  extension in ['.cdx', '.smi', '.csv']:
        args.charge = rules_get_charge(mol,args,log)
        # replace None values if there are metals that are not used
        for i,charge in enumerate(args.charge):
            if charge is None:
                args.charge[i] = 0
        dup_data.at[dup_data_idx, 'Overall charge'] = np.sum(args.charge)
    else:
        dup_data.at[dup_data_idx, 'Overall charge'] = np.sum(args.charge)

    cartesians = mol.GetConformers()[0].GetPositions()
    coordinates = torch.tensor([cartesians.tolist()], requires_grad=True, device=DEVICE)

    ani_incompatible = False
    if program == 'ani':
        try:
            energy, coordinates = ani_calc(elements,cartesians,coordinates,args,log)
        except KeyError:
            log.write(f'\nx  {args.ani_method} could not optimize this molecule (i.e. check of atoms that are not compatible)')
            ani_incompatible = True
            energy = 0
            coordinates = np.zeros((1,3))

    elif program == 'xtb':
        energy, coordinates = xtb_calc(elements,coordinates,args,log,ase_metal,ase_metal_idx)

    else:
        log.write('x  Option not compatible with CMIN (check the available options)!')

    cartesians = np.array(coordinates.tolist()[0])
    # update coordinates of mol object
    for j in range(mol.GetNumAtoms()):
        [x,y,z] = cartesians[j]
        mol.GetConformer().SetAtomPosition(j,Point3D(x,y,z))

    return mol, energy, ani_incompatible
