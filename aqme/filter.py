#####################################################.
#        This file stores all the functions         #
#             used for filtering                    #
#####################################################.

import os
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, Descriptors
from rdkit.ML.Cluster import Butina

from aqme.utils import periodic_table, get_conf_RMS


# Main API of the geometry filter
def geom_filter(self,mol,geom):
    """
    Returns whether a mol object passes all the geometric rules.

    Parameters
    ----------
    self : argparse.args
            Self object with the AQME arguments used
    mol : rdkit.Chem.Mol
            Molecule to be tested.

    Returns
    -------
    bool
        If True, it means that it is in accordance with the rules
    """

    # detects if the geometry rule is for atoms, bonds, angles or dihedrals

    passing = True
    if geom != []:
        passing = False
        if geom == ['Ir_squareplanar']:
            new_geom = Ir_SP_filter(mol)
            if len(new_geom) == 0:
                self.args.log.write(f"x  This molecule is not one of the supported Ir squareplanar complexes! It was discarded by the geom filter")
                passing = False
            passing = matching_fun(self,mol,'Ir_squareplanar',new_geom,passing)

        else:
            if len(geom) != 2:
                self.args.log.write(f"x  The geom option {geom} was not correctly defined, the geometric filter will be turned off! Correct format: [SMARTS,THRESHOLD], for example [CCCO,180] for a 180 degree dihedral")
                return passing
            passing = matching_fun(self,mol,'regular_rule',geom,passing)
    
    return passing


def Ir_SP_filter(mol):
    '''
    Special geometry rule designed to filter the correct conformers of Ir squareplanar complexes.
    So far, the ligands tested are those in DOI: https://doi.org/10.1039/D0SC00445F
    '''

    # get Ir and its potential neighbours
    smarts_list = ['[Ir][C-]','[Ir][N+]','[Ir][n+]','[Ir][N]','[Ir][n]','[Ir][P+]','[Ir][p+]','[Ir][As+]']
    Ir_neighs = []
    L_atom_1, L_atom_2, Ir_idx = None, None, None
    for smarts in smarts_list:
        pairs = list(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
        if len(pairs) > 0:
            for pair in pairs:
                if pair not in Ir_neighs:
                    Ir_neighs.append(pair)

    for Ir_neigh in Ir_neighs:
        for idx in Ir_neigh:
            correct_neigh = False
            # find Ir
            if mol.GetAtoms()[idx].GetAtomicNum() == 77 and Ir_idx is None:
                Ir_idx = idx
            # find right C for lugand of type A and discard all the others from types B and C (from the Chemical Science paper)
            elif mol.GetAtoms()[idx].GetAtomicNum() == 6:
                N_neigh = 0
                for C_neigh in mol.GetAtoms()[idx].GetNeighbors():
                    if C_neigh.GetAtomicNum() == 7:
                        N_neigh += 1
                if N_neigh == 2:
                    correct_neigh = True
            
            # find right N for lugand of type A and discard all the others from types B and C (from the Chemical Science paper)
            elif mol.GetAtoms()[idx].GetAtomicNum() == 7:
                C_neigh = 0
                for N_neigh in mol.GetAtoms()[idx].GetNeighbors():
                    if N_neigh.GetAtomicNum() == 6:
                        C_neigh += 1
                if C_neigh in [2,3]:
                    correct_neigh = True

            # find P and As atoms from ligands of type A
            elif mol.GetAtoms()[idx].GetAtomicNum() in [15,33]:
                correct_neigh = True

            if correct_neigh == True:
                if L_atom_1 is None:
                    L_atom_1 = idx
                elif L_atom_2 is None:
                    L_atom_2 = idx

    # rule: the two ligands of type A are positioned in trans to each other
    if None not in [L_atom_1, L_atom_2, Ir_idx]:
        new_geom = [L_atom_1, Ir_idx, L_atom_2, 180]
    else:
        new_geom = []

    return new_geom

def matching_fun(self,mol,type_match,geom,passing):
    '''
    Checks matches and analyzed if they pass the geometry rules
    '''

    # SMARTS match to detect the atoms and calculate the geometric value. Then, check if the value
    # is within the threshold

    if type_match == 'regular_rule':
        matches = []
        smarts = geom[0]
        geom_val = geom[1]
        smarts_content = ''.join(smarts.replace('[',']').split(']')) # this way both 'ATOM' and '[ATOM]' work
        try:
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        except: # I tried to make this except more specific for Boost.Python.ArgumentError, but apparently it's not as simple as it looks
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(f'[{smarts}]'))
        if len(matches) > 0:
            matches = list(matches[0])
        
    elif type_match == 'Ir_squareplanar':
        matches = geom[:3]
        geom_val = geom[3]
        smarts_content = 'Ir_squareplanar'

    mol_conf = mol.GetConformer(0) # Retrieve the only 3D conformer generated in that mol object for rdMolTransforms
    if smarts_content in periodic_table():
        if len(matches) >= 1:
            passing = True
    elif len(matches) == 2:
        mol_val = rdMolTransforms.GetBondLength(mol_conf, matches[0], matches[1])
        passing = (geom_val - self.args.bond_thres) <= mol_val <= (geom_val + self.args.bond_thres)
    elif len(matches) == 3:
        mol_val = rdMolTransforms.GetAngleDeg(mol_conf, matches[0], matches[1], matches[2])
        passing = (geom_val - self.args.angle_thres) <= mol_val <= (geom_val + self.args.angle_thres)
    elif len(matches) == 4:
        mol_val = rdMolTransforms.GetDihedralDeg(mol_conf, matches[0], matches[1], matches[2], matches[3])
        passing = (geom_val - self.args.dihedral_thres) <= mol_val <= (geom_val + self.args.dihedral_thres)

    return passing


def filters(mol, log, molwt_cutoff):
    """
    Applies some basic filters (molwt, salts[currently off], weird atom symbols)
    that only require SMILES data from a compound and returns if the molecule
    passes the filters or not.
    """

    # Filter 1
    molwt_cutoff = float(molwt_cutoff)
    if Descriptors.MolWt(mol) >= molwt_cutoff and molwt_cutoff > 0:
        log.write(f"x  Skipping this molecule as total molar mass > {molwt_cutoff}")
        return False

    # Filter 2
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    unknown_atoms = list(set(symbols) - set(periodic_table()))
    if unknown_atoms:
        log.write(f"x  Exiting as atoms [{','.join(unknown_atoms)}] are not in the periodic table")
        return False

    # Passed
    return True


def conformer_filters(self,sorted_all_cids,cenergy,outmols):
    '''
    Sequence of three filters based on energy and RMSD
    '''
    
    # filter based on energy window ewin_cmin
    sortedcids = ewin_filter(
        sorted_all_cids,
        cenergy,
        self.args.ewin_cmin,
    )
    # pre-filter based on energy only
    selectedcids_initial = pre_E_filter(
        sortedcids,
        cenergy,
        self.args.initial_energy_threshold,
    )
    # filter based on energy and RMSD
    selectedcids = RMSD_and_E_filter(
        outmols,
        selectedcids_initial,
        cenergy,
        self.args,
        self.args.program.lower(),
    )

    return selectedcids


def ewin_filter(
    sorted_all_cids,
    cenergy,
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
    energy_window : float
            Minimum energy difference with respect to the lowest compound
            discard a compound.

    Returns
    -------
    list
            list of cids accepted
    """

    sortedcids = []
    energy_window = float(energy_window)

    cenergy_min = cenergy[sorted_all_cids[0]]
    # Filter by Energy Window
    for cid in sorted_all_cids:
        if abs(cenergy[cid] - cenergy_min) < energy_window:
            sortedcids.append(cid)

    return sortedcids


def pre_E_filter(
    sortedcids, cenergy, threshold
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

    return selectedcids_initial


def RMSD_and_E_filter(
    outmols, selectedcids_initial, cenergy, args, calc_type
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
                n_mol_1 = -1
                n_mol_2 = -1
                if calc_type == "rdkit":
                    n_mol_1 = seenconf
                    n_mol_2 = conf
                try:
                    rms = get_conf_RMS(
                        outmols[seenconf],
                        outmols[conf],
                        n_mol_1,
                        n_mol_2,
                        args.heavyonly,
                        max_matches_rmsd
                    )
                except RuntimeError:
                    rms = rms_threshold + 1
                    args.log.write('\nx  The mols loaded by RDKit from the SDF file have different substructures and the RMS filter failed. The duplicate filter will only use E on some conformers!')
                if rms < rms_threshold:
                    excluded_conf = True
                    eng_rms_dup += 1
                    break

        if not excluded_conf:
            if conf not in selectedcids:
                selectedcids.append(conf)

    return selectedcids



def cluster_conformers(self, mols, program, csearch_file, name):
    '''
    Performs a Butina clustering based on the RMS differences of the conformers
    '''

    dists = []
    for i in range(len(mols)):
        for j in range(i):
            # using 100 matches only since the molecules are aligned and share the same atom numbering
            dists.append(get_conf_RMS(mols[i], mols[j], -1, -1, self.args.heavyonly, 100))

    if program.lower() == 'rdkit' or self.args.crest_runs == 1:
        # Step 1. Automatically adjust the RMS threshold to meet self.args.sample - 20% of points 
        # that will be added later with the most stable confs from RDKit
        stable_points = int(round(self.args.sample*0.2))
        cluster_points = self.args.sample - stable_points
        if program.lower() == 'rdkit':
            self.args.log.write(f'\no  Selecting {self.args.sample} conformers using a combination of energies and Butina RMS-based clustering. Users might disable this option with --auto_cluster False ({os.path.basename(Path(name))})')
        else:
            self.args.log.write(f'\no  Selecting the most stable RDKit conformer to start a CREST search')

    if program.lower() == 'crest':
        # Step 1. Automatically adjust the RMS threshold to meet self.args.crest_runs
        cluster_points = self.args.crest_runs
        self.args.log.write(f'\no  Selecting {self.args.crest_runs} RDKit conformers using Butina RMS-based clustering to start {self.args.crest_runs} different CREST searches')

    # Threshold instructions: elements within this range of each other are considered to be neighbors
    # and, therefore, the higher the threshold the fewer number of clusters will be generated
    cluster_thr = 1.5

    clusts = Butina.ClusterData(dists, len(mols), cluster_thr, isDistData=True, reordering=True)

    if len(clusts) > cluster_points:
        while len(clusts) > cluster_points:
            cluster_thr = cluster_thr * 1.02 # each iteration the threshold is increased by 2%
            clusts = Butina.ClusterData(dists, len(mols), cluster_thr, isDistData=True, reordering=True)
    elif len(clusts) < cluster_points:
        while len(clusts) < cluster_points:
            cluster_thr = cluster_thr * 0.98 # each iteration the threshold is reduced by 2%
            clusts = Butina.ClusterData(dists, len(mols), cluster_thr, isDistData=True, reordering=True)

    # get centroids (first element of each cluster) and their corresponding mols
    centroids = [x[0] for x in clusts]
    cluster_mols = [mols[x] for x in centroids]

    # adjust to the exact number of clusters to match self.args.sample at the end
    cluster_mols = cluster_mols[:cluster_points]
    allenergy = []
    for mol in cluster_mols:
        allenergy.append(float(mol.GetProp('Energy')))

    if program.lower() == 'rdkit':
        # Step 2. Fill the other positions with add the most stable conformers found initially
        for mol in mols: # already sorted by energy
            if float(mol.GetProp('Energy')) not in allenergy:
                if len(cluster_mols) < self.args.sample:
                    cluster_mols.append(mol)
                    allenergy.append(float(mol.GetProp('Energy')))
                else:
                    break

    # sort mols by energy
    cluster_mols_sorted = [mol for _, mol in sorted(zip(allenergy, cluster_mols), key=lambda pair: pair[0])]

    if program.lower() == 'rdkit':
        # Step 3. Replace the original SDF with the updated SDF
        if not self.args.pytest_testing:
            os.remove(f'{csearch_file}')
        else:
            shutil.move(f'{csearch_file}',f'{Path(str(csearch_file).replace(".sdf","_all_confs.sdf"))}')
        sdwriter = Chem.SDWriter(f'{csearch_file}')
        for mol in cluster_mols_sorted:
            sdwriter.write(mol)
        sdwriter.close()

    return cluster_mols_sorted


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

#     from functools import partial

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


# def passes_Ir_bidentate_x3_rule(mol, angle_thres):
#     """
#     Checks if a mol containing an Ir metal complex passes the bidentate x3
#     rule or not.

#     Parameters
#     ----------
#     mol : rdkit.Chem.Mol
#             The molecule to be tested.
#     angle_thres : float
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
#         ligand_atoms, mol_conf, metal_idx, angle_thres
#     )
#     return passing
