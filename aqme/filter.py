#####################################################.
#        This file stores all the functions         #
#             used for filtering                    #
#####################################################.
from functools import partial

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, Descriptors

from aqme.utils import periodic_table, get_conf_RMS

# Aux functions for passes_Ir_bidentate_x3_rule
def is_carbene_like(neighbours):
    """
    Takes a list of the metal neighbour atoms and returns True if they can be considered as if it were a carbene.

    Parameters
    ----------
    neighbors : list
        list of neighbour atoms of the metal.

    Returns
    -------
    bool
    """

    carbene_like = False
    N_group = ["N", "P", "As"]
    for inside_neighbour in neighbours:
        if inside_neighbour.GetSymbol() in N_group:
            if inside_neighbour.GetTotalValence() == 4:
                carbene_like = True
    return carbene_like


def get_fragmented_versions(mol, atom_i, atom_j, metal_idx):
    """
    Creates two copies of the mol object where in the first one the atom_i and
    in the second the atom_j are not bonded to the metal.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
            The molecule to fragment.
    atom_i : int
            Idx of the first atom
    atom_j : int
            Idx of the second atom
    metal_idx : int
            Idx of the metal

    Returns
    -------
    tuple
            mol_i, mol_j
    """
    bond_i = mol.GetBondBetweenAtoms(atom_i, metal_idx)
    mol_i = Chem.FragmentOnBonds(
        mol, [bond_i.GetIdx()], addDummies=True, dummyLabels=[(atom_i, metal_idx)]
    )
    bond_j = mol.GetBondBetweenAtoms(atom_j, metal_idx)
    mol_j = Chem.FragmentOnBonds(
        mol, [bond_j.GetIdx()], addDummies=True, dummyLabels=[(atom_j, metal_idx)]
    )
    return mol_i, mol_j


def is_around_angle(test, angle, offset):
    """
    Checks if a test angle is close to the angle or not.

    Parameters
    ----------
    test : float
            Angle to test in Degrees.
    angle : float
            Angle to compare in Degrees.
    offset : float
            Tolerance around 'angle' in degrees.

    Returns
    -------
    bool
            True if it is in the range [angle-offset,angle+offset].
    """
    return (angle - offset) <= test <= (angle + offset)


def passes_Ir_bidentate_x3_rule_angle_requirements(
    ligand_atoms, mol, metal_idx, offset
):
    """
    Checks if for complexes with Ph_Py ligands if the angles between the N
    atoms are acceptable. Linear when there are 2 Ph_Py and not linear when
    there are 3 Ph_Py Ligands.

    Parameters
    ----------
    ligand_atoms : list
            A list of tuples where the second value of each tuple is the Idx in the
            provided mol of the N atom of that ligand.
    mol : rdkit.Chem.Mol
            The molcule object of the complex
    metal_idx : int
            Idx of the metal atom
    offset : float
            Angle error in degrees to consider or not an angle to be linear.

    Returns
    -------
    bool
            Whether it passes the angle requirements or not
    """
    is_lineal = partial(is_around_angle, angle=180, offset=offset)
    get_angle = lambda i, j: rdMolTransforms.GetAngleDeg(mol, i, metal_idx, j)
    passing = True

    if len(ligand_atoms) == 3:  # For complexes with 3 Ph_Py ligands:
        i, j, k = ligand_atoms
        angles = [get_angle(a[1], b[1]) for a, b in [(i, j), (i, k), (j, k)]]
        passing = not any([is_lineal(angle) for angle in angles])

    if (
        len(ligand_atoms) == 2
    ):  # For complexes with 2 Ph_Py ligands + 1 ligand that is not Ph_Py
        (i, i2), (j, j2) = ligand_atoms
        angle = get_angle(i2, j2)
        passing = is_lineal(angle)

    return passing


def passes_Ir_bidentate_x3_rule(mol, angle_off):
    """
    Checks if a mol containing an Ir metal complex passes the bidentate x3
    rule or not.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
            The molecule to be tested.
    angle_off : float
            Angle in degrees that acts as tolerance around 180º

    Returns
    -------
    bool
            True if it has passed the rule
    """
    bond_threshold = 2.60  # based on observation from DFT optimized geometries
    ligand_links = []
    atom_indexes = []

    passing = True
    # Finds the Ir atom and gets the atom types and indexes of all its neighbours
    # the filter is compatible with molecules that do not contain Ir (always passing)
    metal_idx = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 77:
            metal_idx = atom.GetIdx()
            for x in atom.GetNeighbors():
                ligand_links.append(x.GetSymbol())
                atom_indexes.append(x.GetIdx())

    # I need to get the only 3D conformer generated in that mol object for rdMolTransforms
    mol_conf = mol.GetConformer(0)
    # This part will identify the pairs of C and N atoms that are part of the same Ph_Py ligand.
    # The shape of the atom pairs is '[[C1_ATOM_NUMBER, N1_ATOM_NUMBER],[C2, N2],...]'.
    # This information is required for the subsequent filtering process based on angles

    # it filters off molecules that the SDF only detects 5 Ir neighbours
    if len(atom_indexes) != 6:
        if metal_idx is not None:
            passing = False
        return passing

    # Filter by distance of the neighbour atoms and the metal
    for atom_idx in atom_indexes:
        bond_length = rdMolTransforms.GetBondLength(mol_conf, metal_idx, atom_idx)
        if bond_length > bond_threshold:
            return False

    ligand_atoms = []
    for atom_i, sym_i in zip(atom_indexes, ligand_links):
        # This is a filter that excludes molecules that fell apart during DFT geometry
        # optimization (i.e. a N atom from one of the ligands separated from Ir). The
        # max distance allowed can be tuned in length_filter
        for atom_j in atom_indexes:
            # Avoid combinations of the same atom with itself
            if atom_i != atom_j and sym_i == "C":
                # We know that the ligands never have 2 carbon atoms bonding the Ir atom except
                # for carbenes. We only use atom_indexes[i] for C atoms, and atom_indexes[j] for the potential
                # N atoms that are part of the same Ph_Py ligand
                neighbours = mol.GetAtoms()[atom_i].GetNeighbors()
                if not is_carbene_like(neighbours):
                    # This part detects the Ir-C bond and breaks it, breaking the Ph_Py ring
                    mol_i, mol_j = get_fragmented_versions(
                        mol, atom_i, atom_j, metal_idx
                    )
                    # identify whether or not the initial 5-membered ring formed between [-Ir-C-C-C-N-] is broken when we break the Ir-C bond. This works
                    # because Ph_Py units bind Ir in the same way always, through 1 C and 1 N that are in the same position, forming a 5-membered ring.
                    # If this ring is broken, atom_indexes[j] will not be part of a 5-membered ring (atom.IsInRingSize(5) == False) which means that
                    # this atom was initially inside the same ligand as the parent C of atom_indexes[i])
                    if not mol_i.GetAtomWithIdx(atom_i).IsInRingSize(5):
                        # doing backwards as well eg. Ir N bond
                        if not mol_j.GetAtomWithIdx(atom_i).IsInRingSize(5):
                            ligand_atoms.append([atom_i, atom_j])
                            break
                    else:
                        if not mol_i.GetAtomWithIdx(atom_j).IsInRingSize(5):
                            if mol.GetAtomWithIdx(atom_j).IsInRingSize(5):
                                ligand_atoms.append([atom_i, atom_j])
                                break
    # Check Angles
    passing = passing and passes_Ir_bidentate_x3_rule_angle_requirements(
        ligand_atoms, mol_conf, metal_idx, angle_off
    )
    return passing


# Main API of the module
def geom_rules_output(mol, args, log, file, print_error_geom_rules):
    """
    returns if a mol object passes all the discarding rules or not.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
            molecule to be tested.
    args : argparse.args
            [description]
    log : Logger
            [description]
    file : str
            Only used to write to the log
    print_error_geom_rules : bool
            Controls extra writing to the log.

    Returns
    -------
    bool
            If True, it means that it is in accordance with the rules
    """
    passing = True
    for rule in args.geom_rules:
        if rule == "Ir_bidentate_x3":
            passing = passes_Ir_bidentate_x3_rule(mol, args.angle_off)
        else:
            var = rule.split(",")
            if len(var) < 2:
                log.write(
                    "x  The geom_rules parameter(s) was not correctly defined, this filter will be turned off"
                )
                return True
            atoms_filter = var[0].split("-")
            angle_rules = int(var[1])
            # the elements of this initial list will be replaced by the corresponding atom id numebrs
            atom_idx = ["ATOM1", "ATOM2", "ATOM3"]

            find_angle = 0
            incompatibility_found = False
            for atom in mol.GetAtoms():
                # count matches
                neigh_count_first = 0
                neigh_count_second = 0
                sym_0, sym_1, sym_2 = atoms_filter
                # Finds the metal atom and gets the atom types and indexes of all its neighbours
                if atom.GetSymbol() == sym_1:
                    # idx of the central atom
                    atom_idx[1] = atom.GetIdx()
                    for x in atom.GetNeighbors():
                        sym = x.GetSymbol()
                        if sym == sym_0 or sym == sym_2:
                            # this ensures that both neighbours are used
                            if sym == sym_0 and sym == sym_2:
                                if neigh_count_first <= neigh_count_second:
                                    neigh_count_first += 1
                                    atom_idx[0] = x.GetIdx()
                                else:
                                    neigh_count_second += 1
                                    atom_idx[2] = x.GetIdx()
                            elif sym == sym_0:
                                neigh_count_first += 1
                                atom_idx[0] = x.GetIdx()
                            elif sym == sym_2:
                                neigh_count_second += 1
                                atom_idx[2] = x.GetIdx()
                            # count matches
                            matches = neigh_count_first + neigh_count_second
                            if matches > 2:
                                if not print_error_geom_rules:
                                    log.write(
                                        f"x  There are multiple options in geom_rules for {file}, this filter will be turned off"
                                    )
                                    incompatibility_found = True
                                    break
                    if neigh_count_first == 1 and neigh_count_second == 1:
                        find_angle += 1
            if find_angle == 0 and not incompatibility_found:
                if not print_error_geom_rules:
                    log.write(
                        f"x  No angles matching the description from geom_rules in {file}, this filter will be turned off"
                    )
            elif find_angle > 1:
                log.write(
                    f"x  {file} contain more than one atom that meets the geom_rules criteria, this filter will be turned off"
                )
            elif find_angle == 1:
                # I need to get the only 3D conformer generated in that mol object for rdMolTransforms
                mol_conf = mol.GetConformer(0)
                # Calculate the angle between the 3 elements
                angle = rdMolTransforms.GetAngleDeg(
                    mol_conf, atom_idx[0], atom_idx[1], atom_idx[2]
                )
                passing = not is_around_angle(angle, angle_rules, args.angle_off)
            if not passing:
                break
    return passing


def filters(mol, log, molwt_cutoff, verbose):
    """
    Applies some basic filters (molwt, salts[currently off], weird atom symbols)
    that only require SMILES data from a compound and returns if the molecule
    passes the filters or not.
    """

    # Filter 1
    if Descriptors.MolWt(mol) >= molwt_cutoff and molwt_cutoff > 0:
        if verbose:
            log.write(f"x   Exiting as total molar mass > {molwt_cutoff}")
        return False

    # Filter 2
    # if len(Chem.MolToSmiles(mol).split('.')) > 1: # Disconnected Molecule
    #    return False

    # Filter 3
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    unknown_atoms = list(set(symbols) - set(periodic_table()))
    if unknown_atoms:
        if verbose:
            log.write(
                f" Exiting as atoms [{','.join(unknown_atoms)}] are not in the periodic table"
            )
        return False

    # Passed
    return True


def ewin_filter(
    sorted_all_cids,
    cenergy,
    args,
    dup_data,
    dup_data_idx,
    log,
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
    log : Logger
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
    verbose = args.verbose

    sortedcids = []
    count = 0

    cenergy_min = cenergy[sorted_all_cids[0]]
    # Filter by Energy Window
    for cid in sorted_all_cids:
        if abs(cenergy[cid] - cenergy_min) < energy_window:
            sortedcids.append(cid)
        else:
            count += 1

    log_msg = ""
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

    # Log it
    if verbose:
        msg = "o  {} conformers rejected {}based on energy window ewin_csearch (E > {} kcal/mol)"
        log.write(msg.format(count, log_msg, energy_window))

    return sortedcids


def pre_E_filter(
    sortedcids, cenergy, dup_data, dup_data_idx, log, calc_type, threshold, verbose
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
    log : Logger
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

    if verbose:
        log.write(
            f"o  {str(eng_dup)} duplicates removed  pre-energy filter (E < {threshold} kcal/mol)"
        )

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
    outmols, selectedcids_initial, cenergy, args, dup_data, dup_data_idx, log, calc_type
):
    """
    This filter selects the first compound that it finds with energy an energy
    difference lower than the threshold with a higher than the threshold rms
    with respect to the nearest (in energy) accepted compound.
    """
    rms_threshold = args.rms_threshold
    energy_threshold = args.energy_threshold
    verbose = args.verbose

    if calc_type not in ["rdkit", "summ", "xtb", "ani"]:
        # RAUL: If this is true "ALL conformers are going to be written"
        #       I Think this should be taken care at the beggining of the
        #       function or before using the function.
        pass

    if verbose:
        log.write(
            f"o  Removing duplicate conformers (RMSD < {rms_threshold} and E difference < {energy_threshold} kcal/mol)"
        )
    # bar = IncrementalBar('o  Filtering based on energy and RMSD', max = len(selectedcids_initial))

    selectedcids = []
    eng_rms_dup = 0
    selectedcids.append(selectedcids_initial[0])

    for _,conf in enumerate(selectedcids_initial[1:]):
        # This keeps track of whether or not your conformer is unique
        excluded_conf = False

        # check energy and rmsd
        for seenconf in selectedcids:
            E_diff = abs(cenergy[conf] - cenergy[seenconf])  # in kcal/mol
            if E_diff < args.energy_threshold:
                if calc_type == "rdkit":
                    rms = get_conf_RMS(
                        outmols[seenconf],
                        outmols[conf],
                        seenconf,
                        conf,
                        args.heavyonly,
                        args.max_matches_rmsd,
                    )
                # elif calc_type == 'summ' or calc_type == 'fullmonte' or calc_type =='xtb' or calc_type =='ani':
                elif calc_type == "summ" or calc_type == "xtb" or calc_type == "ani":
                    rms = get_conf_RMS(
                        outmols[conf],
                        outmols[seenconf],
                        -1,
                        -1,
                        args.heavyonly,
                        args.max_matches_rmsd,
                    )
                if rms < args.rms_threshold:
                    excluded_conf = True
                    eng_rms_dup += 1
                    break
        if not excluded_conf:
            if conf not in selectedcids:
                selectedcids.append(conf)
        # bar.next()
    # bar.finish()

    if verbose:
        log.write(
            f"o  {eng_rms_dup} duplicates removed (RMSD < {rms_threshold} / E < {energy_threshold} kcal/mol)"
        )
        log.write(f"o  {len(selectedcids)} unique conformers remain")

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


# Base classes for the filters (Raul's refactoring part, work in progress)

# class Filter(object):
#     """
#     Base class for the definition of different types of filters. Defines the
#     basic API that any Filter object should have.

#     Parameters
#     ----------
#     function : Function
#             A single parameter function that returns true when the item should pass
#             the filter and False otherwise.

#     Attributes
#     ----------
#     dataset : list or None
#             A list containing a reference to all the elements that where filtered.
#     outcomes : list or None
#             A list of booleans with the same order as the dataset elements.
#     discarded : list or None
#             A list of discarded items from the dataset.
#     accepted : list or None
#             A list of accepted items from the dataset.
#     """

#     def __init__(self, function=None):
#         if function is not None:
#             self.function = function
#         elif getattr(self, "function", None) is None:
#             self.function = lambda x: True
#         self.dataset = None
#         self.outcomes = None
#         self._discarded = None
#         self._accepted = None

#     def add_dummy_parameter(self, function):
#         """
#         Adds a dummy parameter as the first positional parameter to a given
#         function and returns the wrapped function.
#         """

#         def f(dummy, parameter):
#             return function(parameter)

#         return f

#     def apply(self, dataset, key=None, force=False):
#         """
#         Applies the filter to an iterable. Setting the 'dataset' and 'outcomes'
#         attributes of the Filter in the process.

#         Parameters
#         ----------
#         dataset : iterable
#                 Iterable that contains the items to run the filtering
#         key : [type], optional
#                 [description], by default None
#         force : bool, optional
#                 A True value will apply the filter to the dataset forcefully,
#                 overwriting any previous usage of the filter to other dataset.

#         Raises
#         ------
#         ValueError
#                 If the dataset has been already applied to a dataset and the force
#                 keyword is set to 'false'
#         """
#         if self.dataset is not None and not force:
#             msg = "Attempting to apply a previously applied filter with force='false'"
#             raise ValueError(msg)
#         elif self.dataset is not None:
#             self.clean()

#         if key is None:
#             key = lambda x: x
#         self.dataset = list(dataset)
#         outcomes = self.calc_outcomes(dataset, key)
#         self.outcomes = tuple(outcomes)

#     def calc_outcomes(self, dataset, key):
#         """
#         Runs the filter on a dataset and returns the outcomes without storing
#         the values.

#         Parameters
#         ----------
#         dataset : iterable
#                 Iterable that contains the items to run the filtering
#         key : function
#                 [description]
#         """
#         return tuple(self.function(key(data)) for data in dataset)

#     def extend(self, dataset, key=None):
#         """
#         Extends the current dataset with the new dataset and includes the
#         outcomes of applying the filter to the new dataset.

#         Parameters
#         ----------
#         dataset : iterable
#                 Iterable that contains the items to run the filtering
#         key : [type], optional
#                 [description], by default None
#         """
#         if key is None:
#             key = lambda x: x
#         new = list(dataset)
#         outcomes = self.calc_outcomes(new, key)

#         self.dataset = self.dataset + new
#         self.outcomes = tuple(i for i in chain(self.outcomes, outcomes))
#         # And resets the cache of discarded and accepted
#         self._accepted = None
#         self._discarded = None

#     def clean(self):
#         """
#         Resets the state of the Filter removing all the information stored
#         about the last application of the of the filter to a set of data.
#         """
#         self.dataset = None
#         self.outcomes = None
#         self._discarded = None
#         self._accepted = None

#     @property
#     def discarded(self):
#         if self._discarded is None and self.dataset is not None:
#             self._discarded = [
#                 d for out, d in zip(self.outcomes, self.dataset) if not out
#             ]
#         return self._discarded

#     @property
#     def accepted(self):
#         if self._accepted is None and self.dataset is not None:
#             self._accepted = [d for out, d in zip(self.outcomes, self.dataset) if out]
#         return self._accepted


# class CompoundFilter(Filter):
#     """
#     Class used to apply several filters to a same dataset in a certain order and
#     store the information about which ones were discarded in which filter.

#     Parameters
#     ----------
#     *filters : Filter
#             An undefined number of Filter objects

#     Attributes
#     ----------
#     filters : list
#             List of Filter objects in application order.
#     dataset : list or None
#             A list containing a reference to all the elements that where filtered.
#     outcomes : list or None
#             A list of booleans with the same order as the dataset elements.
#     discarded : list or None
#             A list of discarded items from the dataset.
#     accepted : list or None
#             A list of accepted items from the dataset.
#     """

#     def __init__(self, *filters):
#         self.filters = filters
#         super().__init__()

#     def insert(self, index, item):
#         """
#         Inserts a filter before the index position.

#         Parameters
#         ----------
#         item : Filter
#                 The filter object to include.
#         """
#         self.filters.insert(index, item)

#     def pop(self, index=-1):
#         """
#         Remove and return a filter at index (default last).
#         """
#         return self.filters.pop(index)

#     def apply(self, dataset, key=None, keys=None, force=False):
#         """
#         Applies the filter to an iterable. Setting the 'dataset' and 'outcomes'
#         attributes of each Filter in the process.

#         Parameters
#         ----------
#         dataset : iterable
#                 Iterable that contains the items to run the filtering
#         key : function, optional
#                 If no function is provided it is assumed that the filter can process
#                 each item in the dataset without any change, by default None.
#         keys : iterable, optional
#                 Iterable of functions with the same length as the number of filters
#                 to provide different, by default None. Overrides the key argument.
#         force : bool, optional
#                 A True value will apply the filter to the dataset forcefully,
#                 overwriting any previous usage of the filter to other dataset.

#         Raises
#         ------
#         ValueError
#                 If the dataset has been already applied to a dataset and the force
#                 keyword is set to 'false'
#         """
#         if self.dataset is not None and not force:
#             msg = "Attempting to apply a previously applied filter with force='false'"
#             raise ValueError(msg)
#         elif any([f.dataset is not None for f in self.filters]) and not force:
#             msg = (
#                 "At least one of the filters has already been applied and force='false'"
#             )
#             raise ValueError(msg)
#         elif force:
#             self.clean()

#         # Assign the keys iterable
#         if keys is None and key is None:
#             key = lambda x: x
#             keys = [key for _ in self.filters]
#         elif key is not None:
#             keys = [key for _ in self.filters]

#         self.dataset = dataset
#         outcomes = self.calc_outcomes(dataset, keys)

#         # Set the attributes for all the Filters
#         self.filters[0].dataset = dataset
#         self.filters[0].outcomes = outcomes[0]
#         out_old = outcomes[0]
#         for f, out in zip(self.filters[1:], outcomes[1:]):
#             dataset_n = [d for i, d in zip(out_old, self.dataset) if i]
#             f.dataset = dataset_n
#             f.outcomes = out
#             out_old = out

#         # Set the attributes for the current object
#         outcomes = self.homogenize_outcomes(outcomes)

#         self.outcomes = outcomes

#     def calc_outcomes(self, dataset, keys):
#         """
#         Runs the filter on a dataset and returns the outcomes without storing
#         the values.

#         Parameters
#         ----------
#         dataset : iterable
#                 Iterable that contains the items to run the filtering
#         keys : iterable
#                 Iterable of functions that ensure proper input per each filter.
#         """
#         outcomes = []
#         dataset_old = dataset
#         # outcomes_old = (True,)*len(dataset)
#         for f, key in zip(self.filters, keys):
#             # Use the filter and get the dataset for the next filter
#             out = f.calc_outcomes(dataset_old, key=key)
#             dataset_old = [d for i, d in zip(out, dataset_old) if i]
#             outcomes.append(tuple(out))
#         return outcomes

#     def homogenize_outcomes(self, outcomes):
#         """
#         Returns a list of tuples where all tuples are the same size.

#         Parameters
#         ----------
#         outcomes : list
#                 List of tuples with the outcomes of each of the filters.
#         n : int
#                 size of the dataset

#         Returns
#         -------
#         list
#                 list of the homogenized tuples.
#         """
#         homogenized = []
#         outcomes_old = (True,) * len(outcomes[0])
#         for _, out in zip(self.filters, outcomes):
#             _iter = out.__iter__()
#             outcomes_new = [False if not i else next(_iter) for i in outcomes_old]
#             homogenized.append(tuple(outcomes_new))
#             outcomes_old = outcomes_new
#         return homogenized

#     def extend(self, dataset, key=None, keys=None):
#         """
#         Extends the current dataset with the new dataset and includes the
#         outcomes of applying the filter to the new dataset.

#         Parameters
#         ----------
#         dataset : iterable
#                 Iterable that contains the items to run the filtering
#         key : function, optional
#                 If no function is provided it is assumed that the filter can process
#                 each item in the dataset without any change, by default None.
#         keys : iterable, optional
#                 Iterable of functions with the same length as the number of filters
#                 to provide different, by default None. Overrides the key argument.
#         """
#         # Assign the keys iterable
#         if keys is None and key is None:
#             key = lambda x: x
#             keys = [key for _ in self.filters]
#         elif key is not None:
#             keys = [key for _ in self.filters]

#         new = list(dataset)

#         _outcomes = self.outcomes
#         outcomes = self.calc_outcomes(new, keys)

#         # Set the attributes for the current object
#         self.dataset = self.dataset + new

#         # reset the cache of the properties
#         self._accepted = None
#         self._discarded = None

#         # Set the attributes for all the Filters
#         self.filters[0].dataset = self.dataset
#         self.filters[0].outcomes = self.filters[0].outcomes + outcomes[0]
#         out_old = outcomes[0]
#         for f, out in zip(self.filters[1:], outcomes[1:]):
#             dataset_n = [d for i, d in zip(out_old, self.dataset) if i]
#             f.dataset = f.dataset + dataset_n
#             f.outcomes = f.outcomes + out
#             # reset the cache of the properties
#             f._accepted = None
#             f._discarded = None
#         outcomes = self.homogenize_outcomes(outcomes)
#         self.outcomes = [
#             out_old + out_new for out_old, out_new in zip(_outcomes, outcomes)
#         ]

#     def clean(self):
#         """
#         Resets the state of the CompoundFilter removing all the information
#         stored about the last application of the of the filter to a set of data.
#         It cleans all the component filters.
#         """
#         super().clean()
#         for f in self.filters:
#             f.clean()

#     def accepted_from(self, index):
#         """
#         returns the list of accepted items at the specified filter.
#         """
#         return [d for out, d in zip(self.outcomes[index], self.dataset) if out]

#     def discarded_from(self, index):
#         """
#         returns the total list of discarded items after the specified filter.
#         """
#         return [d for out, d in zip(self.outcomes[index], self.dataset) if not out]

#     @Filter.discarded.getter
#     def discarded(self):
#         if self._discarded is None and self.dataset is not None:
#             self._discarded = self.discarded_from(-1)
#         return self._discarded

#     @Filter.accepted.getter
#     def accepted(self):
#         if self._accepted is None and self.dataset is not None:
#             self._accepted = self.accepted_from(-1)
#         return self._accepted


# class RMSDFilter(Filter):
#     """
#     This filter inputs tuples of (molecule,cid). Each time a conformer that
#     passes the filter is found it is added to the pool of conformers.

#     Note: The RMSD calculation done by rdkit has the side effect of leaving the
#     target conformer aligned to the probe conformer.

#     Parameters
#     ----------
#     threshold : float
#             Minimum RMSD to accept a conformer.
#     maxmatches : int
#             maximum number of atoms should match in the alignment previous to the
#             RMSD calculation.
#     heavyonly : [type], optional
#             [description], by default True.
#     reverse : bool, optional
#             Reverses the threshold. If True, only conformers with an RMSD < threshold
#             will be accepted. by default False.
#     is_rdkit : bool
#             If the conformers to compare have been generated with rdkit and the
#             cid to use to access the conformer is the one provided instead of -1.
#             by default False.

#     """

#     def __init__(
#         self, threshold, maxmatches, heavyonly=True, reverse=False, is_rdkit=False
#     ):
#         self.threshold = threshold
#         self.maxmatches = maxmatches
#         self.heavyonly = heavyonly
#         self.pool = []
#         self.is_rdkit = is_rdkit
#         super().__init__(function=self.function)

#     def set_initial_pool(self, pool, key=None):
#         """
#         Sets the initial pool of seen conformers.

#         Parameters
#         ----------
#         pool : list
#                 A list of conformers ( or things convertible to conformers through
#                 the key function)
#         key : function, optional
#                 A function to convert each item in the pool to a conformer that
#                 can be subjected to the filter, by default None.
#                 i.e
#                 >>> myconformers[cid<-int] = conformer
#                 >>> pool = [cid1,cid2,...,cidn]
#                 >>> key = lambda x: myconformers[x]
#         """
#         if key is not None:
#             pool = list(map(key, pool))
#         self.pool = pool

#     def clean(self, also_pool=True):
#         """
#         Resets the state of the Filter removing all the information stored
#         about the last application of the of the filter to a set of data.

#         Parameters
#         ----------

#         also_pool : bool
#                 If true it will also clear the initial pool of conformers.
#         """
#         if also_pool:
#             self.pool = []
#         super().__init__()

#     def function(self, item):
#         """
#         main function of the filter. inputs a tuple (molecule,cid) aligns it to
#         all the previously stored molecules and if it passes the filter, stores
#         it and returns True, otherwise returns False.
#         """
#         probe, cid_p = item
#         if item in self.pool:
#             return True
#         for target, cid_t in self.pool:
#             c1 = c2 = -1
#             if self.is_rdkit:
#                 probe, target = target, probe
#                 c1, c2 = cid_t, cid_p
#             rms = get_conf_RMS(probe, target, c1, c2, self.heavyonly, self.maxmatches)
#             if self.reverse:
#                 reject = rms < self.threshold
#             else:
#                 reject = rms > self.threshold
#             if reject:
#                 return False
#         else:
#             self.pool.append(item)
#             return True


# class EnergyFilter(Filter):
#     """
#     This filter inputs energy values. Each time a conformer that
#     passes the filter is found it is added to the pool of conformers.

#     Note: The RMSD calculation done by rdkit has the side effect of leaving the
#     target conformer aligned to the probe conformer.

#     Parameters
#     ----------
#     threshold : float
#             Energy threshold in kcal/mol.
#     mode : str, ['difference','window']
#             'difference' accepts all the conformers whose energy difference
#             with respect to all the previous conformers found is larger than
#             the threshold.
#             'window' accepts all the conformers whose energy difference with respect
#             to the lowest in energy is smaller than the threshold.

#     Attributes
#     ----------
#     pool : list
#             In window mode, it is the list of minima used to calculate the energy
#             window. If the lowest conformer was either set as initial pool or
#             provided as the first item of the dataset the pool will remain of len=1.
#             In difference mode, it corresponds to all the conformers accepted.
#     """

#     def __init__(self, threshold, mode):
#         self.mode = mode  # difference or window
#         self.threshold = threshold
#         self.pool = []
#         super().__init__()

#     def set_initial_pool(self, pool, key=None):
#         """
#         Sets the initial pool of seen conformers.

#         Parameters
#         ----------
#         pool : list
#                 A list of conformers ( or things convertible to conformers through
#                 the key function)
#         key : function, optional
#                 A function to convert each item in the pool to a conformer that
#                 can be subjected to the filter, by default None.
#                 i.e
#                 >>> myconformers[cid<-int] = conformer
#                 >>> pool = [cid1,cid2,...,cidn]
#                 >>> key = lambda x: myconformers[x]
#         """
#         if key is not None:
#             pool = list(map(key, pool))
#         self.pool = pool

#     def clean(self, also_pool=True):
#         """
#         Resets the state of the Filter removing all the information stored
#         about the last application of the of the filter to a set of data.

#         Parameters
#         ----------

#         also_pool : bool
#                 If true it will also clear the initial pool of conformers.
#         """
#         if also_pool:
#             self.pool = []
#         super().__init__()

#     def function(self, item):
#         """
#         main function of the filter. Inputs a tuple (cid,energy) and accepts or
#         rejects the item depending on the energy threshold and filter mode.
#         """
#         assert self.mode in ["window", "difference"]
#         if self.mode == "window":
#             out = self._window(item)
#         else:
#             out = self._difference(item)
#         return out

#     def _difference(self, item):
#         cid_p, energy_p = item
#         if item in self.pool:
#             return True
#         for cid_t, energy_t in self.pool:
#             if abs(energy_p - energy_t) < self.threshold:
#                 return False
#         else:
#             self.pool.append(item)
#             return True

#     def _window(self, item):
#         if not self.pool:  # Ensure the pool has len = 1
#             self.pool.append(item)
#             return True
#         cid_p, energy_p = item
#         cid_t, energy_t = self.pool[-1]
#         reject = abs(energy_p - energy_t) < self.threshold
#         is_lowest = energy_p < energy_t
#         if reject:
#             return False
#         if is_lowest:
#             self.pool.append(item)
#         return True
