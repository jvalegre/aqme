#####################################################.
#        This file stores all the functions         #
#             used for filtering                    #
#####################################################.

# from functools import partial
# from rdkit.Chem import AllChem as Chem
# from rdkit.Chem import rdMolTransforms, Descriptors
from rdkit.Chem import Descriptors
from aqme.utils import periodic_table, get_conf_RMS

# Aux functions of the geometry filter
# def is_carbene_like(neighbours):
#     """
#     Takes a list of the metal neighbour atoms and returns True if they can be considered as if it were a carbene.

#     Parameters
#     ----------
#     neighbors : list
#         list of neighbour atoms of the metal.

#     Returns
#     -------
#     bool
#     """

#     carbene_like = False
#     N_group = ["N", "P", "As"]
#     for inside_neighbour in neighbours:
#         if inside_neighbour.GetSymbol() in N_group:
#             if inside_neighbour.GetTotalValence() == 4:
#                 carbene_like = True
#     return carbene_like


# def get_fragmented_versions(mol, atom_i, atom_j, metal_idx):
#     """
#     Creates two copies of the mol object where in the first one the atom_i and
#     in the second the atom_j are not bonded to the metal.

#     Parameters
#     ----------
#     mol : rdkit.Chem.Mol
#             The molecule to fragment.
#     atom_i : int
#             Idx of the first atom
#     atom_j : int
#             Idx of the second atom
#     metal_idx : int
#             Idx of the metal

#     Returns
#     -------
#     tuple
#             mol_i, mol_j
#     """

#     bond_i = mol.GetBondBetweenAtoms(atom_i, metal_idx)
#     mol_i = Chem.FragmentOnBonds(
#         mol, [bond_i.GetIdx()], addDummies=True, dummyLabels=[(atom_i, metal_idx)]
#     )
#     bond_j = mol.GetBondBetweenAtoms(atom_j, metal_idx)
#     mol_j = Chem.FragmentOnBonds(
#         mol, [bond_j.GetIdx()], addDummies=True, dummyLabels=[(atom_j, metal_idx)]
#     )
#     return mol_i, mol_j


# def is_around_angle(test, angle, offset):
#     """
#     Checks if a test angle is close to the angle or not.

#     Parameters
#     ----------
#     test : float
#             Angle to test in Degrees.
#     angle : float
#             Angle to compare in Degrees.
#     offset : float
#             Tolerance around 'angle' in degrees.

#     Returns
#     -------
#     bool
#             True if it is in the range [angle-offset,angle+offset].
#     """

#     return (angle - offset) <= test <= (angle + offset)


# def passes_Ir_bidentate_x3_rule_angle_requirements(
#     ligand_atoms, mol, metal_idx, offset
# ):
#     """
#     Checks if for complexes with Ph_Py ligands if the angles between the N
#     atoms are acceptable. Linear when there are 2 Ph_Py and not linear when
#     there are 3 Ph_Py Ligands.

#     Parameters
#     ----------
#     ligand_atoms : list
#             A list of tuples where the second value of each tuple is the Idx in the
#             provided mol of the N atom of that ligand.
#     mol : rdkit.Chem.Mol
#             The molcule object of the complex
#     metal_idx : int
#             Idx of the metal atom
#     offset : float
#             Angle error in degrees to consider or not an angle to be linear.

#     Returns
#     -------
#     bool
#             Whether it passes the angle requirements or not
#     """

#     is_lineal = partial(is_around_angle, angle=180, offset=offset)
#     get_angle = lambda i, j: rdMolTransforms.GetAngleDeg(mol, i, metal_idx, j)
#     passing = True

#     if len(ligand_atoms) == 3:  # For complexes with 3 Ph_Py ligands:
#         i, j, k = ligand_atoms
#         angles = [get_angle(a[1], b[1]) for a, b in [(i, j), (i, k), (j, k)]]
#         passing = not any([is_lineal(angle) for angle in angles])

#     if (
#         len(ligand_atoms) == 2
#     ):  # For complexes with 2 Ph_Py ligands + 1 ligand that is not Ph_Py
#         (i, i2), (j, j2) = ligand_atoms
#         angle = get_angle(i2, j2)
#         passing = is_lineal(angle)

#     return passing


# def passes_Ir_bidentate_x3_rule(mol, angle_off):
#     """
#     Checks if a mol containing an Ir metal complex passes the bidentate x3
#     rule or not.

#     Parameters
#     ----------
#     mol : rdkit.Chem.Mol
#             The molecule to be tested.
#     angle_off : float
#             Angle in degrees that acts as tolerance around 180º

#     Returns
#     -------
#     bool
#             True if it has passed the rule
#     """

#     bond_threshold = 2.60  # based on observation from DFT optimized geometries
#     ligand_links = []
#     atom_indexes = []

#     passing = True
#     # Finds the Ir atom and gets the atom types and indexes of all its neighbours
#     # the filter is compatible with molecules that do not contain Ir (always passing)
#     metal_idx = None
#     for atom in mol.GetAtoms():
#         if atom.GetAtomicNum() == 77:
#             metal_idx = atom.GetIdx()
#             for x in atom.GetNeighbors():
#                 ligand_links.append(x.GetSymbol())
#                 atom_indexes.append(x.GetIdx())

#     # I need to get the only 3D conformer generated in that mol object for rdMolTransforms
#     mol_conf = mol.GetConformer(0)
#     # This part will identify the pairs of C and N atoms that are part of the same Ph_Py ligand.
#     # The shape of the atom pairs is '[[C1_ATOM_NUMBER, N1_ATOM_NUMBER],[C2, N2],...]'.
#     # This information is required for the subsequent filtering process based on angles

#     # it filters off molecules that the SDF only detects 5 Ir neighbours
#     if len(atom_indexes) != 6:
#         if metal_idx is not None:
#             passing = False
#         return passing

#     # Filter by distance of the neighbour atoms and the metal
#     for atom_idx in atom_indexes:
#         bond_length = rdMolTransforms.GetBondLength(mol_conf, metal_idx, atom_idx)
#         if bond_length > bond_threshold:
#             return False

#     ligand_atoms = []
#     for atom_i, sym_i in zip(atom_indexes, ligand_links):
#         # This is a filter that excludes molecules that fell apart during DFT geometry
#         # optimization (i.e. a N atom from one of the ligands separated from Ir). The
#         # max distance allowed can be tuned in length_filter
#         for atom_j in atom_indexes:
#             # Avoid combinations of the same atom with itself
#             if atom_i != atom_j and sym_i == "C":
#                 # We know that the ligands never have 2 carbon atoms bonding the Ir atom except
#                 # for carbenes. We only use atom_indexes[i] for C atoms, and atom_indexes[j] for the potential
#                 # N atoms that are part of the same Ph_Py ligand
#                 neighbours = mol.GetAtoms()[atom_i].GetNeighbors()
#                 if not is_carbene_like(neighbours):
#                     # This part detects the Ir-C bond and breaks it, breaking the Ph_Py ring
#                     mol_i, mol_j = get_fragmented_versions(
#                         mol, atom_i, atom_j, metal_idx
#                     )
#                     # identify whether or not the initial 5-membered ring formed between [-Ir-C-C-C-N-] is broken when we break the Ir-C bond. This works
#                     # because Ph_Py units bind Ir in the same way always, through 1 C and 1 N that are in the same position, forming a 5-membered ring.
#                     # If this ring is broken, atom_indexes[j] will not be part of a 5-membered ring (atom.IsInRingSize(5) == False) which means that
#                     # this atom was initially inside the same ligand as the parent C of atom_indexes[i])
#                     if not mol_i.GetAtomWithIdx(atom_i).IsInRingSize(5):
#                         # doing backwards as well eg. Ir N bond
#                         if not mol_j.GetAtomWithIdx(atom_i).IsInRingSize(5):
#                             ligand_atoms.append([atom_i, atom_j])
#                             break
#                     else:
#                         if not mol_i.GetAtomWithIdx(atom_j).IsInRingSize(5):
#                             if mol.GetAtomWithIdx(atom_j).IsInRingSize(5):
#                                 ligand_atoms.append([atom_i, atom_j])
#                                 break
#     # Check Angles
#     passing = passing and passes_Ir_bidentate_x3_rule_angle_requirements(
#         ligand_atoms, mol_conf, metal_idx, angle_off
#     )
#     return passing


# # Main API of the geometry filter
# def geom_rules_output(mol, args, log, file, print_error_geom_rules):
#     """
#     returns if a mol object passes all the discarding rules or not.

#     Parameters
#     ----------
#     mol : rdkit.Chem.Mol
#             molecule to be tested.
#     args : argparse.args
#             [description]
#     log : Logger
#             [description]
#     file : str
#             Only used to write to the log
#     print_error_geom_rules : bool
#             Controls extra writing to the log.

#     Returns
#     -------
#     bool
#             If True, it means that it is in accordance with the rules
#     """

#     passing = True
#     for rule in args.geom_rules:
#         if rule == "Ir_bidentate_x3":
#             passing = passes_Ir_bidentate_x3_rule(mol, args.angle_off)
#         else:
#             var = rule.split(",")
#             if len(var) < 2:
#                 log.write("x  The geom_rules parameter(s) was not correctly defined, this filter will be turned off")
#                 return True
#             atoms_filter = var[0].split("-")
#             angle_rules = int(var[1])
#             # the elements of this initial list will be replaced by the corresponding atom id numebrs
#             atom_idx = ["ATOM1", "ATOM2", "ATOM3"]

#             find_angle = 0
#             incompatibility_found = False
#             for atom in mol.GetAtoms():
#                 # count matches
#                 neigh_count_first = 0
#                 neigh_count_second = 0
#                 sym_0, sym_1, sym_2 = atoms_filter
#                 # Finds the metal atom and gets the atom types and indexes of all its neighbours
#                 if atom.GetSymbol() == sym_1:
#                     # idx of the central atom
#                     atom_idx[1] = atom.GetIdx()
#                     for x in atom.GetNeighbors():
#                         sym = x.GetSymbol()
#                         if sym == sym_0 or sym == sym_2:
#                             # this ensures that both neighbours are used
#                             if sym == sym_0 and sym == sym_2:
#                                 if neigh_count_first <= neigh_count_second:
#                                     neigh_count_first += 1
#                                     atom_idx[0] = x.GetIdx()
#                                 else:
#                                     neigh_count_second += 1
#                                     atom_idx[2] = x.GetIdx()
#                             elif sym == sym_0:
#                                 neigh_count_first += 1
#                                 atom_idx[0] = x.GetIdx()
#                             elif sym == sym_2:
#                                 neigh_count_second += 1
#                                 atom_idx[2] = x.GetIdx()
#                             # count matches
#                             matches = neigh_count_first + neigh_count_second
#                             if matches > 2:
#                                 if not print_error_geom_rules:
#                                     log.write(f"x  There are multiple options in geom_rules for {file}, this filter will be turned off")
#                                     incompatibility_found = True
#                                     break
#                     if neigh_count_first == 1 and neigh_count_second == 1:
#                         find_angle += 1
#             if find_angle == 0 and not incompatibility_found:
#                 if not print_error_geom_rules:
#                     log.write(f"x  No angles matching the description from geom_rules in {file}, this filter will be turned off")
#             elif find_angle > 1:
#                 log.write(f"x  {file} contain more than one atom that meets the geom_rules criteria, this filter will be turned off")
#             elif find_angle == 1:
#                 # Retrieve the only 3D conformer generated in that mol object for rdMolTransforms
#                 mol_conf = mol.GetConformer(0)
#                 # Calculate the angle between the 3 elements
#                 angle = rdMolTransforms.GetAngleDeg(
#                     mol_conf, atom_idx[0], atom_idx[1], atom_idx[2]
#                 )
#                 passing = not is_around_angle(angle, angle_rules, args.angle_off)
#             if not passing:
#                 break
#     return passing

# Mol filters
def filters(mol, log, molwt_cutoff):
    """
    Applies some basic filters (molwt, salts[currently off], weird atom symbols)
    that only require SMILES data from a compound and returns if the molecule
    passes the filters or not.
    """

    # Filter 1
    molwt_cutoff = float(molwt_cutoff)
    if Descriptors.MolWt(mol) >= molwt_cutoff and molwt_cutoff > 0:
        log.write(f"x   Skipping this molecule as total molar mass > {molwt_cutoff}")
        return False

    # Filter 2
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    unknown_atoms = list(set(symbols) - set(periodic_table()))
    if unknown_atoms:
        log.write(f" Exiting as atoms [{','.join(unknown_atoms)}] are not in the periodic table")
        return False

    # Passed
    return True


def ewin_filter(
    sorted_all_cids,
    cenergy,
    dup_data,
    dup_data_idx,
    calc_type,
    energy_window,
):
    """
    Given a sorted list of Compound Ids and a sorted list of their energies
    it discards all compound Ids that have an energy higher than the
    args.ewin_csearch with respect to the lowest one.

    Parameters
    ----------
    sorted_all_cids : list
            [description]
    cenergy : list
            [description]
    args : argparse.args
            [description]
    dup_data : pd.Dataframe
            [description]
    dup_data_idx : pd.Dataframe?
            [description]
    calc_type : str
            A string that points towards the column of the dataframe that should
            be filled with the number of duplicates. The current choices are:
            ['rdkit','summ','ani','xtb']
    energy_window : float
            Minimum energy difference with respect to the lowest compound
            discard a compound.

    Returns
    -------
    list
            list of cids accepted
    """

    sortedcids = []
    count = 0
    energy_window = float(energy_window)

    cenergy_min = cenergy[sorted_all_cids[0]]
    # Filter by Energy Window
    for cid in sorted_all_cids:
        if abs(cenergy[cid] - cenergy_min) < energy_window:
            sortedcids.append(cid)
        else:
            count += 1

    if calc_type == "rdkit":
        key = "RDKit"
    elif calc_type == "summ":
        key = "summ"
    elif calc_type == "ani":
        key = "ANI"
    elif calc_type == "xtb":
        key = "xTB"

    # Write it
    dup_data.at[dup_data_idx, f"{key}-energy-window"] = count

    return sortedcids


def pre_E_filter(
    sortedcids, cenergy, dup_data, dup_data_idx, calc_type, threshold
):
    """
    This filter selects the first compound that it finds with energy an energy
    difference higher or equal to the threshold with respect to the previously
    admitted compounds. (Thought as filter for rdkit)

    Parameters
    ----------
    sortedcids : list or pd.Dataframe?
            List of compound Ids.
    cenergy : list or pd.Dataframe?
            list of compound energies
    dup_data : pd.Dataframe
            [description]
    dup_data_idx : pd.Dataframe?
            [description]
    calc_type : str
            A string that points towards the column of the dataframe that should
            be filled with the number of duplicates. The current choices are:
            ['rdkit','summ','ani','xtb']
    threshold : float
            Minimum energy difference to consider two compounds as different.
            (kcal/mol)

    Returns
    -------
    list
            list of accepted compound Ids
    """

    selectedcids_initial = []
    eng_dup = 0
    threshold = float(threshold)

    # Add the first one
    selectedcids_initial.append(sortedcids[0])
    for conf in sortedcids[1:]:
        is_unique = True
        # check rmsd
        for seenconf in selectedcids_initial:
            E_diff = abs(cenergy[conf] - cenergy[seenconf])  # in kcal/mol
            if E_diff < threshold:
                eng_dup += 1
                is_unique = False
                break
        if is_unique:
            if conf not in selectedcids_initial:
                selectedcids_initial.append(conf)

    if calc_type == "rdkit":
        column = "RDKit-initial_energy_threshold"
    elif calc_type == "summ":
        column = "summ-initial_energy_threshold"
    # if calc_type == 'fullmonte': column = 'FullMonte-initial_energy_threshold'
    elif calc_type == "ani":
        column = "ANI-initial_energy_threshold"
    elif calc_type == "xtb":
        column = "xTB-initial_energy_threshold"
    else:
        column = ""

    if column:
        dup_data.at[dup_data_idx, column] = eng_dup
    return selectedcids_initial


def RMSD_and_E_filter(
    outmols, selectedcids_initial, cenergy, args, dup_data, dup_data_idx, calc_type
):
    """
    This filter selects the first compound that it finds with energy an energy
    difference lower than the threshold with a higher than the threshold rms
    with respect to the nearest (in energy) accepted compound.
    """

    selectedcids = []
    eng_rms_dup = 0
    selectedcids.append(selectedcids_initial[0])
    energy_threshold = float(args.energy_threshold)
    rms_threshold = float(args.rms_threshold)
    max_matches_rmsd = int(args.max_matches_rmsd)

    for _,conf in enumerate(selectedcids_initial[1:]):
        # This keeps track of whether or not your conformer is unique
        excluded_conf = False

        # check energy and rmsd
        for seenconf in selectedcids:
            E_diff = abs(cenergy[conf] - cenergy[seenconf])  # in kcal/mol
            if E_diff < energy_threshold:
                if calc_type == "rdkit":
                    rms = get_conf_RMS(
                        outmols[seenconf],
                        outmols[conf],
                        seenconf,
                        conf,
                        args.heavyonly,
                        max_matches_rmsd
                    )
                # elif calc_type == 'summ' or calc_type == 'fullmonte' or calc_type =='xtb' or calc_type =='ani':
                elif calc_type == "summ" or calc_type == "xtb" or calc_type == "ani":
                    rms = get_conf_RMS(
                        outmols[conf],
                        outmols[seenconf],
                        -1,
                        -1,
                        args.heavyonly,
                        max_matches_rmsd
                    )
                if rms < rms_threshold:
                    excluded_conf = True
                    eng_rms_dup += 1
                    break

        if not excluded_conf:
            if conf not in selectedcids:
                selectedcids.append(conf)

    # Write the found duplicates:
    if calc_type == "rdkit":
        key = "RDKit"
    elif calc_type == "summ":
        key = "summ"
    elif calc_type == "ani":
        key = "ANI"
    elif calc_type == "xtb":
        key = "xTB"
    else:
        key = ""

    if key:
        duplicates_column = f"{key}-RMSD-and-energy-duplicates"
        uniques_column = f"{key}-Unique-conformers"
        dup_data.at[dup_data_idx, duplicates_column] = eng_rms_dup
        dup_data.at[dup_data_idx, uniques_column] = len(selectedcids)

    return selectedcids
