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


def geom_filter(self, mol_ensemb, mol_geom, geom):
    """Check if a molecule passes geometric filtering rules.
    
    Applies geometry-based filters to conformers, including specialized
    rules for specific metal complexes (e.g., Ir squareplanar).

    Args:
        self: AQME instance with arguments
        mol_ensemb (rdkit.Chem.Mol): Molecule ensemble with all conformers
        mol_geom (rdkit.Chem.Mol): Specific conformer to test
        geom (list): Geometry rule specification
            - Empty list [] means no filtering
            - ['Ir_squareplanar'] for Ir complex special rule
            - [SMARTS, THRESHOLD] for custom geometric constraints

    Returns:
        bool: True if molecule passes all geometric rules, False otherwise
    """
    if not geom:
        return True
    
    # Handle Ir squareplanar special case
    if geom == ['Ir_squareplanar']:
        new_geom = get_ir_squareplanar_geometry(mol_ensemb)
        if not new_geom:
            self.args.log.write(
                "x  This molecule is not one of the supported Ir squareplanar complexes! "
                "It was discarded by the geom filter"
            )
            return False
        return check_geometric_match(self, mol_ensemb, mol_geom, 'Ir_squareplanar', new_geom)
    
    # Handle custom geometry rules
    if len(geom) != 2:
        self.args.log.write(
            f"x  The geom option {geom} was not correctly defined, the geometric filter will be turned off! "
            f"Correct format: [SMARTS,THRESHOLD], for example [CCCO,180] for a 180 degree dihedral"
        )
        return True
    
    return check_geometric_match(self, mol_ensemb, mol_geom, 'regular_rule', geom)


def get_ir_squareplanar_geometry(mol):
    """Extract geometry parameters for Ir squareplanar complexes.
    
    Identifies the two ligands of type A in trans configuration
    for Ir squareplanar complexes based on the ligands from
    DOI: https://doi.org/10.1039/D0SC00445F
    
    Args:
        mol (rdkit.Chem.Mol): Molecule to analyze
    
    Returns:
        list: [L_atom_1, Ir_idx, L_atom_2, 180] if valid geometry found,
              empty list otherwise
    """
    def _is_correct_carbon_neighbor(atom):
        """Check if carbon has exactly 2 nitrogen neighbors (type A ligand)."""
        if atom.GetAtomicNum() != 6:
            return False
        n_neighbors = sum(1 for neighbor in atom.GetNeighbors() 
                         if neighbor.GetAtomicNum() == 7)
        return n_neighbors == 2
    
    def _is_correct_nitrogen_neighbor(atom):
        """Check if nitrogen has 2-3 carbon neighbors (type A ligand)."""
        if atom.GetAtomicNum() != 7:
            return False
        c_neighbors = sum(1 for neighbor in atom.GetNeighbors() 
                         if neighbor.GetAtomicNum() == 6)
        return c_neighbors in [2, 3]
    
    def _is_correct_p_or_as_neighbor(atom):
        """Check if atom is P or As (type A ligand)."""
        return atom.GetAtomicNum() in [15, 33]
    
    # SMARTS patterns to find Ir and its neighbors
    smarts_patterns = [
        '[Ir][C]', '[Ir][N+]', '[Ir][n+]', '[Ir][N]', 
        '[Ir][n]', '[Ir][P+]', '[Ir][p+]', '[Ir][As+]'
    ]
    
    # Collect all Ir-neighbor pairs
    ir_neighbor_pairs = []
    for pattern in smarts_patterns:
        pairs = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        for pair in pairs:
            if pair not in ir_neighbor_pairs:
                ir_neighbor_pairs.append(pair)
    
    # Identify Ir and valid ligand atoms
    ir_idx = None
    ligand_atoms = []
    
    for pair in ir_neighbor_pairs:
        for idx in pair:
            atom = mol.GetAtomWithIdx(idx)
            
            # Find Ir atom
            if atom.GetAtomicNum() == 77 and ir_idx is None:
                ir_idx = idx
                continue
            
            # Check if this is a valid ligand atom
            is_valid = (_is_correct_carbon_neighbor(atom) or
                       _is_correct_nitrogen_neighbor(atom) or
                       _is_correct_p_or_as_neighbor(atom))
            
            if is_valid and idx not in ligand_atoms:
                ligand_atoms.append(idx)
                if len(ligand_atoms) >= 2:
                    break
        
        if len(ligand_atoms) >= 2:
            break
    
    # Return geometry if we found Ir and two ligand atoms (trans configuration)
    if ir_idx is not None and len(ligand_atoms) >= 2:
        return [ligand_atoms[0], ir_idx, ligand_atoms[1], 180]
    
    return []


def check_geometric_match(self, mol_ensemb, mol_geom, match_type, geom):
    """Check if molecule matches geometric constraints.
    
    Validates geometric parameters (bonds, angles, dihedrals) against
    specified thresholds using SMARTS pattern matching.
    
    Args:
        self: AQME instance with threshold arguments
        mol_ensemb (rdkit.Chem.Mol): Molecule ensemble for SMARTS matching
        mol_geom (rdkit.Chem.Mol): Specific conformer for geometry calculation
        match_type (str): Type of matching ('regular_rule' or 'Ir_squareplanar')
        geom (list): Geometry specification
            For regular_rule: [SMARTS, threshold_value]
            For Ir_squareplanar: [atom1, atom2, atom3, angle_value]
    
    Returns:
        bool: True if geometry passes the constraint, False otherwise
    """
    def _find_smarts_matches(mol, smarts_pattern):
        """Find SMARTS pattern matches in molecule."""
        try:
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_pattern))
        except:
            # Handle both 'ATOM' and '[ATOM]' formats
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(f'[{smarts_pattern}]'))
        
        return list(matches[0]) if matches else []
    
    def _check_single_atom_match(smarts_content):
        """Check if SMARTS is just a single atom."""
        return smarts_content in periodic_table()
    
    def _check_bond_length(mol_geom, atom_indices, threshold, bond_thres):
        """Check if bond length is within threshold."""
        mol_val = rdMolTransforms.GetBondLength(mol_geom, atom_indices[0], atom_indices[1])
        return (threshold - bond_thres) <= mol_val <= (threshold + bond_thres)
    
    def _check_angle(mol_geom, atom_indices, threshold, angle_thres):
        """Check if angle is within threshold."""
        mol_val = rdMolTransforms.GetAngleDeg(
            mol_geom, atom_indices[0], atom_indices[1], atom_indices[2]
        )
        return (threshold - angle_thres) <= mol_val <= (threshold + angle_thres)
    
    def _check_dihedral(mol_geom, atom_indices, threshold, dihedral_thres):
        """Check if dihedral angle is within threshold."""
        mol_val = rdMolTransforms.GetDihedralDeg(
            mol_geom, atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3]
        )
        return (threshold - dihedral_thres) <= mol_val <= (threshold + dihedral_thres)
    
    # Extract atom indices and threshold based on match type
    if match_type == 'regular_rule':
        smarts, threshold = geom
        # Remove brackets to get clean atom content
        smarts_content = ''.join(smarts.replace('[', ']').split(']'))
        atom_indices = _find_smarts_matches(mol_ensemb, smarts)
    elif match_type == 'Ir_squareplanar':
        atom_indices = geom[:3]
        threshold = geom[3]
        smarts_content = 'Ir_squareplanar'
    else:
        return False
    
    # Check based on number of matched atoms
    if _check_single_atom_match(smarts_content):
        return len(atom_indices) >= 1
    elif len(atom_indices) == 2:
        return _check_bond_length(mol_geom, atom_indices, threshold, self.args.bond_thres)
    elif len(atom_indices) == 3:
        return _check_angle(mol_geom, atom_indices, threshold, self.args.angle_thres)
    elif len(atom_indices) == 4:
        return _check_dihedral(mol_geom, atom_indices, threshold, self.args.dihedral_thres)
    
    return False


def filters(mol, log, molwt_cutoff):
    """Apply basic molecular filters based on molecular weight and atom types.
    
    Filters molecules that exceed weight cutoff or contain unsupported atoms.
    
    Args:
        mol (rdkit.Chem.Mol): Molecule to filter
        log: Logger object for writing messages
        molwt_cutoff (float): Maximum allowed molecular weight (0 = no limit)
    
    Returns:
        bool: True if molecule passes all filters, False otherwise
    """
    def _check_molecular_weight(mol, cutoff, log):
        """Check if molecular weight is below cutoff."""
        cutoff = float(cutoff)
        if cutoff <= 0:
            return True
        
        mol_weight = Descriptors.MolWt(mol)
        if mol_weight >= cutoff:
            log.write(f"x  Skipping this molecule as total molar mass > {cutoff}")
            return False
        return True
    
    def _check_atom_types(mol, log):
        """Check if all atoms are in the periodic table."""
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        unknown_atoms = list(set(symbols) - set(periodic_table()))
        
        if unknown_atoms:
            log.write(f"x  Exiting as atoms [{','.join(unknown_atoms)}] are not in the periodic table")
            return False
        return True
    
    # Apply filters
    if not _check_molecular_weight(mol, molwt_cutoff, log):
        return False
    
    if not _check_atom_types(mol, log):
        return False
    
    return True


def conformer_filters(self, sorted_all_cids, cenergy, outmols):
    """Apply sequential energy and RMSD-based conformer filters.
    
    Three-stage filtering process:
    1. Energy window filter (ewin_cmin)
    2. Pre-filter based on energy differences
    3. Combined energy and RMSD filter
    
    Args:
        self: AQME instance with filter arguments
        sorted_all_cids (list): List of conformer IDs sorted by energy
        cenergy (dict/list): Conformer energies indexed by ID
        outmols (dict/list): Conformer molecule objects indexed by ID
    
    Returns:
        list: Selected conformer IDs that pass all filters
    """
    # Stage 1: Filter by energy window
    sortedcids = apply_energy_window_filter(
        sorted_all_cids,
        cenergy,
        self.args.ewin_cmin,
    )
    
    # Stage 2: Pre-filter based on energy differences only
    selectedcids_initial = apply_pre_energy_filter(
        sortedcids,
        cenergy,
        self.args.initial_energy_threshold,
    )
    
    # Stage 3: Filter based on both energy and RMSD
    selectedcids = apply_rmsd_and_energy_filter(
        outmols,
        selectedcids_initial,
        cenergy,
        self.args,
    )
    
    return selectedcids


def apply_energy_window_filter(sorted_all_cids, cenergy, energy_window):
    """Filter conformers by energy window from the lowest energy conformer.
    
    Discards conformers with energy higher than energy_window relative
    to the lowest energy conformer.
    
    Args:
        sorted_all_cids (list): Conformer IDs sorted by energy
        cenergy (dict/list): Conformer energies indexed by ID
        energy_window (float): Maximum energy difference from minimum (kcal/mol)
    
    Returns:
        list: Conformer IDs within the energy window
    """
    energy_window = float(energy_window)
    min_energy = cenergy[sorted_all_cids[0]]
    
    filtered_cids = []
    for cid in sorted_all_cids:
        if abs(cenergy[cid] - min_energy) < energy_window:
            filtered_cids.append(cid)
    
    return filtered_cids


def apply_pre_energy_filter(sortedcids, cenergy, threshold):
    """Pre-filter conformers based solely on energy differences.
    
    Selects conformers with energy differences greater than or equal to
    the threshold compared to previously selected conformers. This reduces
    the number of conformers before the more expensive RMSD filtering.
    
    Args:
        sortedcids (list): Conformer IDs sorted by energy
        cenergy (dict/list): Conformer energies indexed by ID
        threshold (float): Minimum energy difference to consider unique (kcal/mol)
    
    Returns:
        list: Selected conformer IDs passing energy pre-filter
    """
    threshold = float(threshold)
    selected_cids = [sortedcids[0]]  # Always include lowest energy
    
    for cid in sortedcids[1:]:
        is_unique = True
        
        # Check energy difference with all selected conformers
        for selected_cid in selected_cids:
            energy_diff = abs(cenergy[cid] - cenergy[selected_cid])
            if energy_diff < threshold:
                is_unique = False
                break
        
        if is_unique and cid not in selected_cids:
            selected_cids.append(cid)
    
    return selected_cids


def apply_rmsd_and_energy_filter(outmols, selectedcids_initial, cenergy, args):
    """Filter conformers based on combined energy and RMSD criteria.
    
    For each conformer, compares with already selected conformers. If energy
    difference is below threshold AND RMSD is below threshold with any selected
    conformer, the new conformer is rejected as a duplicate.
    
    Args:
        outmols (dict/list): Conformer molecule objects indexed by ID
        selectedcids_initial (list): Pre-filtered conformer IDs
        cenergy (dict/list): Conformer energies indexed by ID
        args: Arguments object with thresholds:
            - energy_threshold: Energy similarity threshold (kcal/mol)
            - rms_threshold: RMSD similarity threshold (Angstroms)
            - heavyonly: Use only heavy atoms for RMSD
            - max_matches_rmsd: Maximum atom matches for RMSD calculation
    
    Returns:
        list: Final selected conformer IDs
    """
    selected_cids = [selectedcids_initial[0]]  # Always include lowest energy
    energy_threshold = float(args.energy_threshold)
    rms_threshold = float(args.rms_threshold)
    max_matches_rmsd = int(args.max_matches_rmsd)
    
    for cid in selectedcids_initial[1:]:
        is_duplicate = False
        
        # Check against all already selected conformers
        for selected_cid in selected_cids:
            energy_diff = abs(cenergy[cid] - cenergy[selected_cid])
            
            # Only calculate RMSD if energies are similar
            if energy_diff < energy_threshold:
                try:
                    rms = get_conf_RMS(
                        outmols[selected_cid],
                        outmols[cid],
                        -1,  # n_mol_1
                        -1,  # n_mol_2
                        args.heavyonly,
                        max_matches_rmsd
                    )
                except RuntimeError:
                    # Different substructures - use only energy filter
                    rms = rms_threshold + 1
                    args.log.write(
                        '\nx  The mols loaded by RDKit from the SDF file have different substructures '
                        'and the RMS filter failed. The duplicate filter will only use E on some conformers!'
                    )
                
                if rms < rms_threshold:
                    is_duplicate = True
                    break
        
        if not is_duplicate and cid not in selected_cids:
            selected_cids.append(cid)
    
    return selected_cids



def compute_pairwise_rms_distances(self, mols):
    """Compute pairwise RMS distances for all conformers.
    
    Creates a distance matrix of RMS values between all pairs of conformers.
    Uses only 100 atom matches since molecules are aligned with same numbering.
    
    Args:
        mols (list): List of RDKit molecule objects with conformers
    
    Returns:
        list: Flattened upper triangular distance matrix
    """
    dists = []
    for i in range(len(mols)):
        for j in range(i):
            # using 100 matches only since the molecules are aligned and share the same atom numbering
            rms = get_conf_RMS(mols[i], mols[j], -1, -1, self.args.heavyonly, 100)
            dists.append(rms)
    return dists


def determine_cluster_points(self, program, sample, name):
    """Determine target number of cluster points based on program and settings.
    
    For RDKit: Uses 80% of sample size for clustering, reserves 20% for stable conformers.
    For CREST: Uses number of CREST runs as cluster points.
    
    Args:
        program (str): Program name ('rdkit' or 'crest')
        sample (int): Total number of conformers to select
        name (str): Molecule name for logging
    
    Returns:
        int: Target number of cluster points
    """
    if program.lower() == 'rdkit' or self.args.crest_runs == 1:
        # Reserve 20% of points for most stable conformers
        stable_points = int(round(sample * 0.2))
        cluster_points = sample - stable_points
        
        if program.lower() == 'rdkit':
            self.args.log.write(
                f'\no  Selecting {sample} conformers using a combination of energies and '
                f'Butina RMS-based clustering. Users might disable this option with '
                f'--auto_cluster False ({os.path.basename(Path(name))})'
            )
        else:
            self.args.log.write(
                f'\no  Selecting the most stable RDKit conformer to start a CREST search'
            )
    elif program.lower() == 'crest':
        # Use number of CREST runs as cluster points
        cluster_points = self.args.crest_runs
        self.args.log.write(
            f'\no  Selecting {self.args.crest_runs} RDKit conformers using Butina '
            f'RMS-based clustering to start {self.args.crest_runs} different CREST searches'
        )
    
    return cluster_points


def adjust_cluster_threshold(dists, num_mols, cluster_points):
    """Automatically adjust clustering threshold to reach target cluster count.
    
    Uses Butina clustering algorithm and iteratively adjusts threshold until
    the number of clusters matches the target. Higher threshold = fewer clusters.
    
    Args:
        dists (list): Pairwise RMS distances
        num_mols (int): Number of molecules to cluster
        cluster_points (int): Target number of clusters
    
    Returns:
        tuple: (clusters, final_threshold)
    """
    # Threshold: elements within this range are considered neighbors
    # Higher threshold = fewer clusters
    cluster_thr = 1.5
    
    clusts = Butina.ClusterData(dists, num_mols, cluster_thr, isDistData=True, reordering=True)
    
    # Increase threshold if too many clusters
    while len(clusts) > cluster_points:
        cluster_thr = cluster_thr * 1.02  # increase by 2%
        clusts = Butina.ClusterData(dists, num_mols, cluster_thr, isDistData=True, reordering=True)
    
    # Decrease threshold if too few clusters
    while len(clusts) < cluster_points:
        cluster_thr = cluster_thr * 0.98  # decrease by 2%
        clusts = Butina.ClusterData(dists, num_mols, cluster_thr, isDistData=True, reordering=True)
    
    return clusts


def extract_cluster_centroids(mols, clusts, cluster_points):
    """Extract centroid conformers from each cluster.
    
    Gets the first element (centroid) from each cluster and retrieves
    the corresponding molecule objects.
    
    Args:
        mols (list): List of all molecule objects
        clusts (list): List of clusters from Butina clustering
        cluster_points (int): Maximum number of centroids to extract
    
    Returns:
        tuple: (cluster_mols, allenergy) - selected molecules and their energies
    """
    # Get centroids (first element of each cluster)
    centroids = [x[0] for x in clusts]
    cluster_mols = [mols[x] for x in centroids]
    
    # Adjust to exact number of clusters
    cluster_mols = cluster_mols[:cluster_points]
    
    # Extract energies
    allenergy = [float(mol.GetProp('Energy')) for mol in cluster_mols]
    
    return cluster_mols, allenergy


def add_stable_conformers(mols, cluster_mols, allenergy, sample):
    """Fill remaining slots with most stable conformers.
    
    After clustering, add the lowest energy conformers that weren't
    selected as centroids to reach the target sample size.
    
    Args:
        mols (list): All molecules sorted by energy
        cluster_mols (list): Currently selected cluster centroids
        allenergy (list): Energies of selected molecules
        sample (int): Target total number of conformers
    
    Returns:
        tuple: (updated_cluster_mols, updated_allenergy)
    """
    for mol in mols:  # already sorted by energy
        energy = float(mol.GetProp('Energy'))
        
        if energy not in allenergy:
            if len(cluster_mols) < sample:
                cluster_mols.append(mol)
                allenergy.append(energy)
            else:
                break
    
    return cluster_mols, allenergy


def write_clustered_sdf(self, cluster_mols_sorted, csearch_file):
    """Write clustered conformers to SDF file.
    
    Replaces the original SDF file with the filtered conformers.
    In pytest mode, moves original file instead of removing it.
    
    Args:
        cluster_mols_sorted (list): Molecules sorted by energy
        csearch_file (str): Path to SDF file to write
    """
    # Handle original file
    if not self.args.pytest_testing:
        os.remove(f'{csearch_file}')
    else:
        backup_path = Path(str(csearch_file).replace(".sdf", "_all_confs.sdf"))
        shutil.move(f'{csearch_file}', f'{backup_path}')
    
    # Write filtered conformers
    sdwriter = Chem.SDWriter(f'{csearch_file}')
    for mol in cluster_mols_sorted:
        sdwriter.write(mol)
    sdwriter.close()


def cluster_conformers(self, mols, program, csearch_file, name, sample):
    """Performs Butina clustering based on RMS differences of conformers.
    
    Uses a two-step approach for RDKit:
    1. Cluster conformers using RMS-based Butina algorithm (80% of sample)
    2. Fill remaining slots with lowest energy conformers (20% of sample)
    
    For CREST, selects conformers as starting points for CREST searches.
    
    Args:
        mols (list): List of RDKit molecule objects with conformers
        program (str): Program name ('rdkit' or 'crest')
        csearch_file (str): Path to SDF file
        name (str): Molecule name for logging
        sample (int): Total number of conformers to select
    
    Returns:
        list: Selected and sorted molecule objects
    """
    # Step 1: Compute pairwise RMS distances
    dists = compute_pairwise_rms_distances(self,mols)
    
    # Step 2: Determine target number of clusters
    cluster_points = determine_cluster_points(self,program, sample, name)
    
    # Step 3: Adjust threshold to reach target cluster count
    clusts = adjust_cluster_threshold(dists, len(mols), cluster_points)
    
    # Step 4: Extract cluster centroids
    cluster_mols, allenergy = extract_cluster_centroids(mols, clusts, cluster_points)
    
    # Step 5: For RDKit, add most stable conformers
    if program.lower() == 'rdkit':
        cluster_mols, allenergy = add_stable_conformers(mols, cluster_mols, allenergy, sample)
    
    # Step 6: Sort molecules by energy
    cluster_mols_sorted = [mol for _, mol in sorted(zip(allenergy, cluster_mols), key=lambda pair: pair[0])]
    
    # Step 7: Write to SDF file for RDKit
    if program.lower() == 'rdkit':
        write_clustered_sdf(self, cluster_mols_sorted, csearch_file)
    
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
