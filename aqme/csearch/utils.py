#####################################################.
#            This file stores functions             #
#                 used in CSEARCH                   #
#####################################################.

import os
import sys
import subprocess
import pandas as pd
import ast
from pathlib import Path
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolAlign

from aqme.utils import (
    get_info_input,
    mol_from_sdf_or_mol_or_mol2,
    read_xyz_charge_mult,
    add_prefix_suffix,
)
from aqme.csearch.crest import nci_ts_mol


def csv_2_list(contraints):
    try:
        if pd.isnull(contraints):
            contraints = []
    except ValueError:
        pass
    if not isinstance(contraints, list):
        contraints = ast.literal_eval(contraints)
    
    return contraints


def prepare_direct_smi(args):
    job_inputs = []

    if args.name is not None:
        name = args.name
        name = add_prefix_suffix(name, args)
    else:
        args.log.write(f"\nx  Specify a name ('name' option) when using the 'smi' option!")
        args.log.finalize()
        sys.exit()

    obj = (
        args.smi,
        name,
        args.charge,
        args.mult,
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
        args.complex_type,
        args.geom,
        args.sample
    )
    job_inputs.append(obj)

    return job_inputs


def prepare_smiles_files(args, csearch_file):
    with open(csearch_file) as smifile:
        lines = [line for line in smifile if line.strip()]
    job_inputs = []
    for _, line in enumerate(lines):
        (
            smi,
            name,
        ) = prepare_smiles_from_line(line, args)
        obj = (
            smi,
            name,
            args.charge,
            args.mult,
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
            args.complex_type,
            args.geom,
            args.sample
        )
        job_inputs.append(obj)

    return job_inputs


def prepare_smiles_from_line(line, args):

    toks = line.split()
    # editing part
    smiles = toks[0]
    name = toks[1]
    name = add_prefix_suffix(name, args)

    # Check for N@@ or N@ in SMILES
    if "N@@" in smiles or "N@" in smiles:
        args.log.write(f"\nx  WARNING! AQME does not support quiral N atoms in SMILES strings (N@@ or N@). These atoms were replaced by N in the SMILES: {smiles}.")
        smiles = smiles.replace("N@@", "N").replace("N@", "N")

    return smiles, name


def prepare_csv_files(args, csearch_file):
    # Check if file is empty before reading
    if os.path.getsize(csearch_file) == 0:
        args.log.write(f"File {args.input} is empty!")
        args.log.finalize()
        sys.exit()
    
    csv_smiles = pd.read_csv(csearch_file)
    
    # Check if DataFrame is empty
    if csv_smiles.empty:
        args.log.write(f"File {args.input} is empty!")
        args.log.finalize()
        sys.exit()
    
    # Check if 'code_name' exists and is entirely empty/NaN
    elif "code_name" in csv_smiles.columns and csv_smiles["code_name"].dropna().empty:
        args.log.write(f"File {args.input} has a 'code_name' column with no values.")
        args.log.finalize()
        sys.exit()

    # avoid running calcs with special signs (i.e. *)
    for name_csv_indiv in csv_smiles['code_name']:
        if '*' in f'{name_csv_indiv}':
            args.log.write(f"\nx  WARNING! The names provided in the CSV contain * (i.e. {name_csv_indiv}). Please, remove all the * characters.")
            args.log.finalize()
            sys.exit()

    job_inputs = []
    # run conformer searches only for unique SMILES
    unique_smiles = set()
    smi_col = False
    for column_index, column in enumerate(csv_smiles.columns):
        if "SMILES" == column.upper() or "SMILES_" in column.upper():
            smi_col = True
            for i in range(len(csv_smiles)):
                obj = generate_mol_from_csv(args, csv_smiles, i, column_index)
                if obj is not None:
                    if obj[0] not in unique_smiles:
                        job_inputs.append(obj)
                        unique_smiles.add(obj[0])
                    else:
                        args.log.write(f'\nx  SMILES "{obj[0]}" used in {obj[1]} is a duplicate, it was already used with a different code_name!')
                       
    if not smi_col:
        args.log.write("\nx  Make sure the CSV file contains a column called 'SMILES', 'smiles' or 'SMILES_' with the SMILES of the molecules!")
        args.log.finalize()
        sys.exit()
    return job_inputs


def generate_mol_from_csv(args, csv_smiles, index, column_index):
    # assigning names and smi in each loop

    smi_valid = False
    for column in csv_smiles.columns:
        if column.upper() == "SMILES" or column.upper().startswith("SMILES_"):
            column_name = csv_smiles.columns[column_index]
            smiles = csv_smiles.loc[index, column_name]
            # this part avoids empty cells at the end of CSV files
            if str(smiles).lower() != 'nan':
                smi_valid = True
            break

    if smi_valid:
        # Check for N@@ or N@ in SMILES
        if "N@@" in smiles or "N@" in smiles:
            args.log.write(f"\nx  WARNING! AQME does not support quiral N atoms in SMILES strings (N@@ or N@). These atoms were replaced by N in the SMILES: {smiles}.")
            smiles = smiles.replace("N@@", "N").replace("N@", "N")

        try:
            name = str(csv_smiles.loc[index, "code_name"])
            column_name = csv_smiles.columns[column_index]
            if column_name.upper() == "SMILES" or not "_" in column_name:
                name += ""
            else:
                suffix = column_name.split("_")[-1]
                name += "_" + suffix
        except KeyError:
            args.log.write("\nx  Make sure the CSV file contains a column called 'code_name' with the names of the molecules!")
            args.log.finalize()
            sys.exit()

        constraints_atoms = args.constraints_atoms
        constraints_dist = args.constraints_dist
        constraints_angle = args.constraints_angle
        constraints_dihedral = args.constraints_dihedral

        if "constraints_atoms" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'constraints_atoms']).lower() != 'nan':
                constraints_atoms = csv_smiles.loc[index, "constraints_atoms"]

        if "constraints_dist" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'constraints_dist']).lower() != 'nan':
                constraints_dist = csv_smiles.loc[index, "constraints_dist"]

        if "constraints_angle" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'constraints_angle']).lower() != 'nan':
                constraints_angle = csv_smiles.loc[index, "constraints_angle"]

        if "constraints_dihedral" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'constraints_dihedral']).lower() != 'nan':
                constraints_dihedral = csv_smiles.loc[index, "constraints_dihedral"]

        constraints_atoms = csv_2_list(constraints_atoms)
        constraints_dist = csv_2_list(constraints_dist)
        constraints_angle = csv_2_list(constraints_angle)
        constraints_dihedral = csv_2_list(constraints_dihedral)

        charge = args.charge
        mult = args.mult
        if "charge" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'charge']).lower() != 'nan':
                charge = csv_smiles.loc[index, "charge"]
        if "mult" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'mult']).lower() != 'nan':
                mult = csv_smiles.loc[index, "mult"]

        complex_type = args.complex_type
        if "complex_type" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'complex_type']).lower() != 'nan':
                complex_type = csv_smiles.loc[index, "complex_type"]

        geom = args.geom
        if "geom" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'geom']).lower() != 'nan':
                geom = csv_smiles.loc[index, "geom"]
        geom = csv_2_list(geom)

        sample = args.sample
        if "sample" in csv_smiles.columns:
            if str(csv_smiles.loc[index, 'sample']).lower() != 'nan':
                sample = csv_smiles.loc[index, "sample"]

        obj = (
            smiles,
            name,
            charge,
            mult,
            constraints_atoms,
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
            complex_type,
            geom,
            sample
        )

        return obj
    
    else:
        return None


def prepare_cdx_files(args, csearch_file):
    # converting to smiles from chemdraw
    molecules = generate_mol_from_cdx(csearch_file)

    job_inputs = []
    for i, (smiles, _) in enumerate(molecules):
        name = f"{'.'.join(os.path.basename(Path(csearch_file)).split('.')[:-1])}_{str(i)}"
        name = add_prefix_suffix(name, args)

        obj = (
            smiles,
            name,
            args.charge,
            args.mult,
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
            args.complex_type,
            args.geom,
            args.sample
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

    filename = os.path.basename(Path(csearch_file))
    if os.path.basename(Path(filename)).split('.')[-1] in ["gjf", "com"]:
        xyz_file, _, _ = com_2_xyz(csearch_file)
        if args.charge is None:
            _, charge, _ = get_info_input(csearch_file)
        else:
            charge = args.charge
        if args.mult is None:
            _, _, mult = get_info_input(csearch_file)
        else:
            mult = args.mult
    else:
        xyz_file = csearch_file
        if args.charge is None:
            charge, _ = read_xyz_charge_mult(csearch_file)
        else:
            charge = args.charge
        if args.mult is None:
            _, mult = read_xyz_charge_mult(csearch_file)
        else:
            mult = args.mult
    _ = xyz_2_sdf(xyz_file)

    sdffile = f'{os.path.dirname(Path(csearch_file))}/{".".join(filename.split(".")[:-1])}.sdf'
    suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(sdffile, "csearch", args)

    name = f'{".".join(filename.split(".")[:-1])}'
    name = add_prefix_suffix(name, args)

    obj = (
        suppl[0],
        name,
        charge,
        mult,
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle,
        args.constraints_dihedral,
        args.complex_type,
        args.geom,
        args.sample
    )
    job_inputs.append(obj)
    if os.path.basename(Path(csearch_file)).split('.')[-1] in ["gjf", "com"]:
        os.remove(xyz_file)
    os.remove(sdffile)

    return job_inputs


def prepare_pdb_files(args, csearch_file):
    filename = os.path.basename(csearch_file)
    sdffile = f'{os.path.dirname(csearch_file)}/{".".join(filename.split(".")[:-1])}.sdf'
    command_pdb = [
        "obabel",
        "-ipdb",
        csearch_file,
        "-osdf",
        f'-O{sdffile}',
    ]
    subprocess.run(command_pdb, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    job_inputs = prepare_sdf_files(args, csearch_file)
    os.remove(sdffile)
    return job_inputs


def prepare_sdf_files(args, csearch_file):
    filename = os.path.basename(csearch_file)
    sdffile = f'{os.path.dirname(csearch_file)}/{filename}'

    suppl, charges, mults, IDs = mol_from_sdf_or_mol_or_mol2(sdffile, "csearch", args)
    IDs = [os.path.basename(x) for x in IDs]

    job_inputs = []
    for mol, charge, mult, name in zip(suppl, charges, mults, IDs):
        name = add_prefix_suffix(name, args)
        obj = (
            mol,
            name,
            charge,
            mult,
            args.constraints_atoms,
            args.constraints_dist,
            args.constraints_angle,
            args.constraints_dihedral,
            args.complex_type,
            args.geom,
            args.sample
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

    name = f'{os.path.dirname(Path(file))}/{os.path.basename(Path(file)).split(".xyz")[0]}'
    command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name + ".sdf"]
    subprocess.run(command_xyz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def check_constraints(self):
    if (
        (len(self.args.constraints_atoms) != 0)
        or (len(self.args.constraints_dist) != 0)
        or (len(self.args.constraints_angle) != 0)
        or (len(self.args.constraints_dihedral) != 0)
    ):
        complex_ts = True
    else:
        complex_ts = False

    return complex_ts


def com_2_xyz(input_file):
    """
    COM to XYZ to SDF for obabel
    """

    filename = '.'.join(os.path.basename(Path(input_file)).split('.')[:-1])
    path_xyz = f'{os.path.dirname(input_file)}/{filename}.xyz'

    # Create the 'xyz' file and/or get the total charge
    xyz, charge, mult = get_info_input(input_file)
    xyz_txt = "\n".join(xyz)
    with open(path_xyz, "w") as F:
        F.write(f"{len(xyz)}\n{filename}\n{xyz_txt}\n")

    return path_xyz, charge, mult


def realign_mol(
    mol, conf, coord_Map, alg_Map, mol_template, maxsteps
):  # RAUL: This function requires a clear separation between minimization and alignment.
    """
    Minimizes and aligns the molecule provided freezing the atoms that match the mol_template

    Parameters
    ----------
    mol : RDKit mol object
        Molecule to be minimized and aligned
    conf : int
        Number that indicates which conformation of the molecule will be minimized and aligned
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]
    maxsteps : int
        Maximum number of iterations in FF minimization

    Returns
    -------
    mol,energy
        The updated mol object and the final forcefield energy.
    """

    num_atom_match = mol.GetSubstructMatch(mol_template)
    forcefield = Chem.UFFGetMoleculeForceField(mol, confId=conf)
    for i, idxI in enumerate(num_atom_match):
        for idxJ in num_atom_match[i + 1 :]:
            d = coord_Map[idxI].Distance(coord_Map[idxJ])
            forcefield.AddDistanceConstraint(idxI, idxJ, d, d, 10000)
    forcefield.Initialize()
    forcefield.Minimize(maxIts=maxsteps)
    # rotate the embedded conformation onto the core_mol:
    rdMolAlign.AlignMol(
        mol,
        mol_template,
        prbCid=conf,
        refCid=-1,
        atomMap=alg_Map,
        reflect=True,
        maxIters=100,
    )
    energy = float(forcefield.CalcEnergy())

    return mol, energy


def minimize_rdkit_energy(mol, conf, log, FF, maxsteps):
    """
    Minimizes a conformer of a molecule and returns the final energy.
    """

    forcefield = None
    if FF.upper() == "NO FF":
        energy = 0

    else:
        if FF.upper() == "MMFF":
            properties = Chem.MMFFGetMoleculeProperties(mol)
            forcefield = Chem.MMFFGetMoleculeForceField(mol, properties, confId=conf)
            if forcefield is None:
                log.write(f"x  Force field {FF} did not work! Changing to UFF.")

        if FF.upper() == "UFF" or forcefield is None:
            # if forcefield is None means that MMFF will not work. Attempt UFF.
            forcefield = Chem.UFFGetMoleculeForceField(mol, confId=conf)

        forcefield.Initialize()
        try:
            forcefield.Minimize(maxIts=maxsteps)
        except RuntimeError:
            log.write(f"\nx  Geometry minimization failed with {FF}, using non-optimized geometry.")
        energy = float(forcefield.CalcEnergy())

    return energy


def getDihedralMatches(mol, heavy):
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


def smi_to_mol(
    smi,
    program,
    log,
    seed,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
):

    complex_ts = False
    smi = smi.split(".")
    if (
        len(smi) > 1
        or len(constraints_atoms) != 0
        or len(constraints_dist) != 0
        or len(constraints_angle) != 0
        or len(constraints_dihedral) != 0
    ):
        if program not in ["crest"]:
            log.write(f"\nx  {program} not supported for conformer generation of complexes and TSs (your SMILES has {len(smi)} parts, separated by a period)! Specify: program='crest' for complexes")
            sys.exit()

        (
            mol,
            constraints_atoms,
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
        ) = nci_ts_mol(
            smi,
            log,
            seed,
            constraints_atoms,
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
        )
        complex_ts = True

    else:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        smi = smi[0]
        try:
            # fix mapped atoms
            if ':' in smi:
                # smi = fix_mapped_atoms(smi)
                log.write(f"\nx  WARNING! The SMILES string provided ( {smi} ) contains mapped atoms, make sure you include their corresponding H atoms explicitly in the SMILES (otherwise they'll be omitted). For example, use [C:1]([H])([H])([H])C instead of [C:1]C.\n")

            mol = Chem.MolFromSmiles(smi, params)
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
        except Chem.AtomValenceException:
            log.write(f"\nx  The SMILES string provided ( {smi} ) contains errors or the molecule needs to be drawn in a different way. For example, N atoms from ligands of metal complexes should be N+ since they're drawn with four bonds in ChemDraw, same for O atoms in carbonyl ligands, etc.\n")
            mol = None

    return (
        mol,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
        complex_ts
    )

# this function was disabled to allow ROBERT users to use atom idx for atomic descriptor generation
# def fix_mapped_atoms(smi):
#     '''
#     This protocol to handle mapped SMILES. Otherwise, Hs are not added right and charges/mult
#     are not calculated correctly either.
#     '''

#     map_list = []
#     smi_map = smi.replace(']','[').split('[')
#     for piece in smi_map:
#         if ':' in piece:
#             map_list.append(f'[{piece}]')
#     for map_atom in map_list:
#         new_atom = map_atom.replace(':','[').split('[')
#         while('' in new_atom):
#             new_atom.remove('')
#         new_atom = new_atom[0]
#         smi = smi.replace(map_atom,new_atom)
    
#     return smi

def substituted_mol(mol, checkI, metal_atoms):
    """
    Returns a molecule object in which all metal atoms are replaced by Iodine
    and the charge is set depending on the number of neighbors.

    """

    metal_idx = []
    complex_coord = []
    metal_sym = []

    for _ in metal_atoms:
        metal_idx.append(None)
        complex_coord.append(None)
        metal_sym.append(None)

    Neighbors2FormalCharge = dict()
    for i, j in zip(range(2, 9), range(-3, 4)):
        Neighbors2FormalCharge[i] = j

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in metal_atoms:
            metal_sym[metal_atoms.index(symbol)] = symbol
            metal_idx[metal_atoms.index(symbol)] = atom.GetIdx()
            complex_coord[metal_atoms.index(symbol)] = len(
                atom.GetNeighbors()
            )
            if checkI == "I":
                atom.SetAtomicNum(53)
                n_neighbors = len(atom.GetNeighbors())
                if n_neighbors > 1:
                    formal_charge = Neighbors2FormalCharge[n_neighbors]
                    atom.SetFormalCharge(formal_charge)

    return metal_idx, complex_coord, metal_sym
