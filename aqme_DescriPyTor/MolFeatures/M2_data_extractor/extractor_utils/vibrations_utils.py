from typing import List, Optional
import numpy as np
import numpy.typing as npt
import pandas as pd
import math
import networkx as nx

try:
    from ...utils.help_functions import *
except ImportError:
    from utils.help_functions import * 

def find_atoms_with_similar_amplitude(
    vibration_mode_dict: dict,
    coordinates_array: np.ndarray,
    atom_pair: list,
    similarity_tol: float = 0.1,   # 10% similarity by default
    amplitude_type: str = "projection"  # or "norm"
):
    """
    For each vibration mode (by frequency), find atoms with similar amplitude of motion to the given atom pair.

    Parameters
    ----------
    vibration_mode_dict: dict
        {frequency: np.ndarray of shape (n_atoms, 3)}
    coordinates_array: np.ndarray
        shape (n_atoms, 3), the geometry of the molecule
    atom_pair: list of int
        Indices of two atoms
    similarity_tol: float
        Relative tolerance for "similar" amplitude (e.g., 0.1 for ±10%)
    amplitude_type: str
        'projection' to use the bond direction, or 'norm' for total displacement

    Returns
    -------
    List of dicts with keys:
        'frequency', 'pair_amplitude', 'similar_atoms', 'similar_amps'
    """
    results = []
    bond_vec = coordinates_array[atom_pair[0]] - coordinates_array[atom_pair[1]]
    bond_vec_norm = bond_vec / np.linalg.norm(bond_vec)

    for freq, displacements in vibration_mode_dict.items():
        if amplitude_type == "projection":
            # For the atom pair, sum of absolute projections
            pair_amp = (
                abs(np.dot(displacements[atom_pair[0]], bond_vec_norm)) +
                abs(np.dot(displacements[atom_pair[1]], bond_vec_norm))
            )
            # All atom amplitudes (projection)
            atom_amps = [abs(np.dot(vec, bond_vec_norm)) for vec in displacements]
        elif amplitude_type == "norm":
            pair_amp = (
                np.linalg.norm(displacements[atom_pair[0]]) +
                np.linalg.norm(displacements[atom_pair[1]])
            )
            atom_amps = [np.linalg.norm(vec) for vec in displacements]
        else:
            raise ValueError("amplitude_type must be 'projection' or 'norm'")

        # Find atoms not in atom_pair with similar amplitude
        similar_atoms = [
            i for i, amp in enumerate(atom_amps)
            if i not in atom_pair and np.isclose(amp, pair_amp, rtol=similarity_tol)
        ]
        similar_amps = [atom_amps[i] for i in similar_atoms]

        results.append({
            "frequency": freq,
            "pair_amplitude": pair_amp,
            "similar_atoms": similar_atoms,
            "similar_amps": similar_amps
        })
    return results

def build_vibration_mode_dict(vibration_dict, info_df):
    """
    Returns:
        vibration_mode_dict: {frequency: np.array shape (n_atoms, 3)}
    """
    n_atoms = len(vibration_dict)

    # Use the minimum number of modes across all atoms (in case of mismatch)
    n_modes_per_atom = [np.array(v).shape[0] for v in vibration_dict.values()]
    min_modes = min(n_modes_per_atom)

    frequencies = info_df['Frequency'].values

    vibration_mode_dict = {}
    for i, freq in enumerate(frequencies[:min_modes]):
        mode_displacements = np.vstack([
            np.array(vibration_dict[f'vibration_atom_{atom_idx+1}'])[i]
            for atom_idx in range(n_atoms)
        ])  # shape (n_atoms, 3)
        
        vibration_mode_dict[freq] = mode_displacements
    
    return vibration_mode_dict

def calc_max_frequency_magnitude(vibration_array: npt.ArrayLike, info_df: pd.DataFrame,
                                 threshhold: int = 1500) -> pd.DataFrame:  ##add option to return ordered_info_df-like dot.prod.info
    """
    a function that gets vibration and info dataframes and returns
    the frequecy and IR with the max magnitude for the vibration.
    splits the coordinates of vibration to 3 coordinates and calculates the magnituede. takes frequencys greter than 1500
    and returns the frequency and IR corresponding to max magnitude.

    Parameters
    ----------
    vibration_array: np.array
        organized vibration file in array.
    info_df: np.dataframe
        organized info file in dataframe.
    threshhold: int
        the threshhold for the frequency. default is 1500.
        Frequency      IR
   0      20.3253  0.0008
   1      25.3713  0.0023
   2      29.0304  0.0019

    Returns
    -------
    dataframe
        max frequency and IR for specific vibration.

    Output:
                                    54
            Frequency         1689.5945
            IR intensity        6.5260
    """

    magnitude = np.linalg.norm(vibration_array, axis=1)
    df = info_df.T
    df['magnitude'] = magnitude
    outer_finger = (df['Frequency'].astype(float) > threshhold)
    index_max = df[outer_finger]['magnitude'].idxmax()
    return info_df[index_max]


def check_pair_in_bonds(pair, bonds_df):  ##help functions for gen_vibration
    """
    a function that checks that the all given atom pair exists as a bond in the bonds_df
    """
    bonds_list = (bonds_df.astype(int)).values.tolist()
    bool_check = (pair in bonds_list) or (pair[::-1] in bonds_list)
    
    return bool_check


def vibrations_dict_to_list(vibration_dict: dict, vibration_atom_nums: list[int]):
    """
    a function that gets a vibration dictionary and a list of atom numbers and returns a list of vibration arrays
    for each atom number.
    ----------
    vibration_dict : dict
        dictionary containing vibration arrays.
    vibration_atom_nums : list
        a list of chosen atoms to get vibration for.
    Returns
    -------
    vibration_array_pairs ([0]) : list
    containing the arrays from the vibration file, each variable containing two
    arrays corresponding to the bond pairs
        .
    vibration_array_list ([1]): list
        containing all vibration array coordinates.

    """
    try:
        vibration_array_list = [vibration_dict[f'vibration_atom_{num}']for num in vibration_atom_nums] 
    except Exception as e:
        return print('Error: no vibration for those atoms-pick another one')
    vibration_array_pairs = list(
        zip(vibration_array_list[::2], vibration_array_list[1::2]))  # [::2] means every second element
    return vibration_array_pairs, vibration_array_list  ##optional- do not split to pairs
    


def calc_vibration_dot_product(extended_vib_df, coordinates_vector):
    vibration_dot_product_list=[]
    for i in range (extended_vib_df.shape[0]):
        vibration_dot_product_list.append(abs(np.dot(extended_vib_df[[0,1,2]].iloc[i], coordinates_vector)) + abs(np.dot(extended_vib_df[[3,4,5]].iloc[i], coordinates_vector)))
    # vibration_dot_product=abs(np.dot(extended_vib_df[[0,1,2]], coordinates_vector)) + abs(np.dot(extended_vib_df[[3,4,5]], coordinates_vector))
    return vibration_dot_product_list

def extended_df_for_stretch(vibration_dict: dict, info_df:pd.DataFrame,atom_pair,threshhold: int = 1400, upper_threshold: int = 3500):
    vibration_array_pairs,_= vibrations_dict_to_list(vibration_dict, atom_pair)
    array=pd.DataFrame(np.hstack([vibration_array_pairs[0][0],vibration_array_pairs[0][1]]))
    df=pd.concat([array,info_df['Frequency']],axis=1)
    filter_df=df[df['Frequency']>threshhold]
    if upper_threshold:
        filter_df=filter_df[filter_df['Frequency']<upper_threshold]
    return filter_df


# def calc_vibration_dot_product_from_pairs(coordinates_array: npt.ArrayLike,
#                                           vibration_dict: dict,
#                                           atom_pair: list, info_df: pd.DataFrame, operation:str='dot',threshold=3000 , vibration_mode_dict = None) -> List[float]:
#     """
#     Calculates the dot product between a vibration mode vector and the bond vector between two atoms.

#     Parameters
#     ----------
#     coordinates_array: np.array
#         An array of x, y, z coordinates for each atom in the molecule.
#     vibration_dict: dict
#         A dictionary where the keys are atom indices and the values are numpy arrays representing the
#         x, y, z components of the vibration mode vectors for that atom.
#     atom_pairs: list
#         A list of pairs of atom indices representing the bond vectors to calculate the dot product with.

#     Returns
#     -------
#     vibration_dot_product: list
#         A list of dot products between each vibration mode vector and the bond vector between the two atoms
#         in each pair in `atom_pairs`.
#     """
    
#     atoms = adjust_indices(atom_pair)
#     extended_df=extended_df_for_stretch(vibration_dict,info_df,atom_pair,threshold)
#     print(f'extended_df shape: {extended_df}')
#     coordinates_vector=coordinates_array[atoms[0]]-coordinates_array[atoms[1]]
#     vibration_dot_product = calc_vibration_dot_product(extended_df, coordinates_vector)
#     extended_df['Amplitude']=vibration_dot_product
#     ## add the frequency to the extended_df
#     extended_df['Frequency'] = info_df['Frequency']
#     # take the frequency of the highest amplitude, got to vibration_mode_dict and find the atoms with similar amplitude, using the new bond vector
   
#     extended_df.reset_index(drop=True,inplace=True)

def check_directional_symmetry(vec1, vec2, threshold=0.0):
    """Helper to classify alignment as symmetric/asymmetric."""
    v1 = np.array(vec1)
    v2 = np.array(vec2)
    v1n = v1 / np.linalg.norm(v1) if np.linalg.norm(v1) > 0 else v1
    v2n = v2 / np.linalg.norm(v2) if np.linalg.norm(v2) > 0 else v2
    sim = np.dot(v1n, v2n)

    if sim > threshold:
        return "asymmetric", sim
    elif sim < threshold:
        return "symmetric", sim
    else:
        return "unclear", sim



def calc_vibration_dot_product_from_pairs(
    coordinates_array: np.ndarray,
    vibration_dict: dict,
    atom_pair: list,
    info_df: pd.DataFrame,
    operation: str = 'dot',
    threshold: float = 1600,
    upper_threshold: float = 3500,
    vibration_mode_dict: dict = None,
    similarity_tol: float = 0.1
) -> pd.DataFrame:
    """
    Calculates the dot product between vibration mode vectors and the bond vector between two atoms.
    For the two highest-amplitude modes, checks vector alignment and tags as symmetric/asymmetric.
    Adds original coordinates of top moving atoms for each mode for comparison/plotting.
    """

    atoms = adjust_indices(atom_pair)
    extended_df = extended_df_for_stretch(vibration_dict, info_df, atom_pair, threshold, upper_threshold)
    coordinates_vector = coordinates_array[atoms[0]] - coordinates_array[atoms[1]]
    vibration_dot_product = calc_vibration_dot_product(extended_df, coordinates_vector)
    extended_df['Amplitude'] = vibration_dot_product

    # Add frequency column (ensure alignment)
    if 'Frequency' not in extended_df.columns:
        extended_df['Frequency'] = info_df.loc[extended_df.index, 'Frequency'].values

    extended_df.reset_index(drop=True, inplace=True)
    return extended_df
    # tag_results = []

    # # Now: for the two highest amplitude modes, do symmetry tagging
    # if vibration_mode_dict is not None and not extended_df.empty:
    #     # Get indices for top 2 amplitudes
    #     top_idxs = extended_df['Amplitude'].abs().nlargest(2).index.tolist()
    #     # tag_results = []
        
    #     tag_list = []
    #     for idx in top_idxs:
    #         freq = extended_df.loc[idx, 'Frequency']
    #         amp = extended_df.loc[idx, 'Amplitude']
    #         freq_key = freq
    #         if freq_key in vibration_mode_dict:
    #             mode_vectors = vibration_mode_dict[freq_key]  # shape: (n_atoms, 3)
    #         else:
    #             mode_vectors = list(vibration_mode_dict.values())[idx]
    #         # Get amplitudes of each atom for this mode
          
    #         amplitudes = np.linalg.norm(mode_vectors, axis=1)
    #         # Sort by closeness to 'amp'
    #         sorted_idx = np.argsort(np.abs(amplitudes - amp))

    #         first_atom = sorted_idx[0]
    #         second_atom = sorted_idx[1] if len(sorted_idx) > 1 else first_atom
            
    #         min_ratio = 0.5
    #         if amplitudes[second_atom] < min_ratio * amplitudes[first_atom]:

    #             tag = "insufficient_second_movement"
    #             sim = np.nan
    #             print(f"Second atom (index {second_atom}) movement ({amplitudes[second_atom]:.3f}) is less than half of first atom ({amplitudes[first_atom]:.3f} , {first_atom}) for mode {idx} (skipping symmetry tag).")
    #             vibration_df, idx = calc_max_frequency_gen_vibration(extended_df)
    #             vibration_df.rename(index={idx: f'Stretch_{atom_pair[0]}_{atom_pair[1]}'})
    #             return vibration_df
                
    #         else:
             

    #             if first_atom > second_atom:
    #                 first_atom, second_atom = second_atom, first_atom
    #             vec1 = mode_vectors[first_atom]
    #             vec2 = mode_vectors[second_atom]
                
    #             tag, sim = check_directional_symmetry(vec1, vec2)
    #             # print atoms and amplitudes
    #             print(f"Mode {idx}: Atoms {first_atom} and {second_atom} with amplitudes {amplitudes[first_atom]:.3f} and {amplitudes[second_atom]:.3f}, similarity: {sim:.3f}")
    #             Amplitude= amp
    #             tag+=f'_Stretch_{atom_pair[0]}_{atom_pair[1]}'
    #             tag_results.append({
    #                 'Frequency': freq,
    #                 'Amplitude': Amplitude,
    #             })
    #             tag_list.append(tag)
                
    #     tag_df = pd.DataFrame(tag_results,index=tag_list) 
    #     # check tag list, if its the same tag for both , check symmetry of vec from third_atom and forth_atom
    #     if len(tag_list) == 2 and tag_list[0] == tag_list[1]:
    #         print(f"Both modes have the same tag: {tag_list[0]}. Checking symmetry with third and fourth atoms.")
    #         tag_results_backup = []
    #         tag_list_backup = []
    #         for idx in top_idxs:
    #             freq = extended_df.loc[idx, 'Frequency']
    #             amp = extended_df.loc[idx, 'Amplitude']
    #             freq_key = freq
    #             if freq_key in vibration_mode_dict:
    #                 mode_vectors = vibration_mode_dict[freq_key]  # shape: (n_atoms, 3)
    #             amplitudes= np.linalg.norm(mode_vectors, axis=1)
    #             sorted_idx = np.argsort(np.abs(amplitudes))[::-1]
    #             third_atom = sorted_idx[2] if len(sorted_idx) > 2 else first_atom
    #             fourth_atom = sorted_idx[3] if len(sorted_idx) > 3 else second_atom
    #             vec3 = mode_vectors[third_atom]
    #             vec4 = mode_vectors[fourth_atom]
                
    #             tag, sim = check_directional_symmetry(vec3, vec4)
    #             tag += f'_Stretch_{atom_pair[0]}_{atom_pair[1]}'
    #             tag_results_backup.append({
    #                 'Frequency': freq,
    #                 'Amplitude': Amplitude,
    #             })
    #             tag_list_backup.append(tag)
            
    #         tag_df = pd.DataFrame(tag_results_backup,index=tag_list_backup)

    #     print(tag_df, "\nWe strongly recommend to visualize the vibration modes to confirm the symmetry/asymmetry.")
        # return tag_df

def calc_max_frequency_gen_vibration(extended_df):  ## fix so 0 is IR and 1 is frequency
    """
    a function that gets info dataframe and vibration dot product and returns the max frequency and IR for each vibration.
    Parameters
    ----------
    vibration_dot_product: list
        list of dot product for each vibration.
    info_df:
        Frequency      IR
   0      20.3253  0.0008
   1      25.3713  0.0023
   2      29.0304  0.0019
    threshhold: int
        the threshhold for the frequency. default is 500.
    Returns
    -------
    IR            0.556928
    Frequency  1124.742600
    IR            0.403272
    Frequency  1128.584700
    """

    # df=pd.DataFrame(np.vstack([vibration_dot_product, (info_df)['Frequency']]),
    #                         index=XYZConstants.VIBRATION_INDEX.value).T
    
    index_max=extended_df['Amplitude'].idxmax()

    max_frequency_vibration=(pd.DataFrame(extended_df.iloc[index_max]).T)[XYZConstants.VIBRATION_INDEX.value]

    return max_frequency_vibration, index_max


def vibration_ring_array_list_to_vector(vibration_array_list: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """
    A function that receives a list of vibration arrays, representing a ring, and returns the sum of
    the first, third and fifth arrays as one vector and the sum of the second, fourth and sixth arrays as
    another vector. These vectors represent the ring's vibrations.
    Parameters
    ----------
    vibration_array_list: List[np.ndarray]
        A list containing six vibration arrays.
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]:
        A tuple containing two np.ndarray vectors representing the ring's vibrations.
    """
    atom_one,atom_two,atom_three=vibration_array_list[0], vibration_array_list[4],vibration_array_list[2]
    atom_four,atom_five,atom_six=vibration_array_list[1], vibration_array_list[3],vibration_array_list[5]
    vec_sum_1_3_5 = (atom_one+atom_three+atom_five) ## try
    vec_sum_2_4_6 = (atom_two+atom_four+atom_six)

  
    return vec_sum_1_3_5, vec_sum_2_4_6



def get_data_for_ring_vibration(info_df: pd.DataFrame, vibration_array_list: List[np.ndarray],
                                coordinates_vector: np.ndarray) -> pd.DataFrame:
    
    product = [np.dot(row_1, row_2) for row_1, row_2 in zip(*vibration_ring_array_list_to_vector(vibration_array_list))]
    _, vibration_array_list = vibration_ring_array_list_to_vector(vibration_array_list)

    
    sin_angle = [abs(math.sin(calc_angle(row, coordinates_vector))) for row in vibration_array_list]

    data_df = pd.DataFrame(np.vstack([product, (info_df)['Frequency'], sin_angle]),index=['Product','Frequency','Sin_angle'] )
  
    return data_df.T



def get_filter_ring_vibration_df(data_df: pd.DataFrame, prods_threshhold: float = 0.01,
                                 frequency_min_threshhold: float = 1600,
                                 frequency_max_threshhold: float = 1750) -> pd.DataFrame:
    # Filter based on product value
    filter_prods = (abs(data_df[XYZConstants.RING_VIBRATION_INDEX.value[0]]) > prods_threshhold) & \
                   (data_df[XYZConstants.RING_VIBRATION_INDEX.value[0]] != 0)
    filter_frequency = (data_df[XYZConstants.RING_VIBRATION_INDEX.value[1]] > frequency_min_threshhold) & \
                       (data_df[XYZConstants.RING_VIBRATION_INDEX.value[1]] < frequency_max_threshhold)
    # Apply combined filter
    filtered_df = data_df[filter_prods & filter_frequency].reset_index()

    if filtered_df.empty:
        print('No data within the specified thresholds. Adjust your thresholds.')
    
    return filtered_df



def indices_to_coordinates_vector(coordinates_array,indices):
    """
    a function that recives coordinates_array and indices of two atoms
    and returns the bond vector between them
    """
    try:
        if  isinstance(indices[0], (list, np.ndarray, tuple)):
            bond_vector=[(coordinates_array[index[1]]-coordinates_array[index[0]]) for index in indices]
            
        else:
            bond_vector= coordinates_array[indices[1]]-coordinates_array[indices[0]]
    except:
        if  isinstance(indices[0], tuple):
            bond_vector=[(coordinates_array[index[0]]-coordinates_array[index[1]]) for index in indices]
        else:
            bond_vector= coordinates_array[indices[0]]-coordinates_array[indices[1]]

       
    return bond_vector

def get_filtered_ring_df(info_df: pd.DataFrame,
                         coordinates_array: np.ndarray,
                         vibration_dict: dict,
                         ring_atom_indices: list,
                         debug: bool = False) -> pd.DataFrame:
    """
    A function that returns a filtered DataFrame of ring vibrations based on their product, frequency,
    and sin(angle) values.

    Parameters
    ----------
    info_df : pd.DataFrame
        DataFrame with vibrational frequencies/intensities.

    coordinates_array : np.ndarray
        (N_atoms, 3) array of x,y,z coordinates.

    vibration_dict : dict
        Vibrational modes → (freq, intensities, vectors).

    ring_atom_indices : list
        Atom indices defining atoms in the ring.

    debug : bool
        If True, prints debug info at each step.

    Returns
    -------
    filtered_df : pd.DataFrame
        A DataFrame that contains the filtered ring vibrations.
    """

    if debug:
        print(f"\n[DEBUG] --- get_filtered_ring_df called ---")
        print(f"[DEBUG] ring_atom_indices (raw): {ring_atom_indices}")

    # adjust indices
    ring_indices = adjust_indices(ring_atom_indices)
    if debug:
        print(f"[DEBUG] ring_indices (after adjust): {ring_indices}")

    # get coordinate vector
    coordinates_vector = indices_to_coordinates_vector(coordinates_array, ring_indices)[0]
    if debug:
        print(f"[DEBUG] coordinates_vector shape: {coordinates_vector.shape}")
        print(f"[DEBUG] coordinates_vector (first few): {coordinates_vector[:5]}")

    # vibration arrays
    vibration_atom_nums = flatten_list(ring_atom_indices)
    _, vibration_array_list = vibrations_dict_to_list(vibration_dict, vibration_atom_nums)
    if debug:
        print(f"[DEBUG] vibration_atom_nums: {vibration_atom_nums}")
        print(f"[DEBUG] vibration_array_list length: {len(vibration_array_list)}")
        if vibration_array_list:
            print(f"[DEBUG] vibration_array_list[0] shape: {np.shape(vibration_array_list[0])}")

    # build DataFrame
    data_df = get_data_for_ring_vibration(info_df, vibration_array_list, coordinates_vector)
    if debug:
        print(f"[DEBUG] data_df shape: {data_df.shape}")
        print(f"[DEBUG] data_df columns: {list(data_df.columns)}")
        print(f"[DEBUG] data_df head:\n{data_df}")

    # filtering
    filtered_df = get_filter_ring_vibration_df(data_df)
    if debug:
        if filtered_df is None or filtered_df.empty:
            print(f"[DEBUG] filtered_df is empty after filtering")
        else:
            print(f"[DEBUG] filtered_df shape: {filtered_df.shape}")
            print(f"[DEBUG] filtered_df head:\n{filtered_df}")

    return filtered_df



def calc_min_max_ring_vibration(filtered_df: pd.DataFrame, ring_atom_indices: list) -> pd.DataFrame:
    """
    Calculates the minimum and maximum vibration frequency and the angle between the vibration vector
    and the plane of the ring, safely handling asin domain errors.
    """
    import math
    import numpy as np
    import pandas as pd

    freq_col = XYZConstants.RING_VIBRATION_INDEX.value[2]

    max_idx = filtered_df[freq_col].idxmin()
    min_idx = filtered_df[freq_col].idxmax()

    max_vibration_frequency = filtered_df.iloc[max_idx, 2]
    min_vibration_frequency = filtered_df.iloc[min_idx, 2]

    # ---- fix: safely clip the values to [-1, 1] before asin ----
    def safe_asin(x):
        return math.degrees(math.asin(float(np.clip(x, -1.0, 1.0))))

    asin_max = safe_asin(filtered_df.iloc[max_idx, 2])
    asin_min = safe_asin(filtered_df.iloc[min_idx, 2])

    df = pd.DataFrame(
        [(max_vibration_frequency, asin_max, min_vibration_frequency, asin_min)],
        columns=XYZConstants.RING_VIBRATION_COLUMNS.value
    )

    df.columns = [f"{col}" for col in df.columns]

    return df



### bending vibration
def find_center_atom(atom1: str, atom2: str, adjacency_dict: Dict[str, List[str]]) -> bool:
    
    neighbors1 = adjacency_dict.get(atom1, [])
    neighbors2 = adjacency_dict.get(atom2, [])
 
    
    common_atoms = set(neighbors1).intersection(neighbors2)

    
    return bool(common_atoms)


def create_adjacency_dict_for_pair(bond_df: pd.DataFrame, pair: List[str]) -> Dict[str, List[str]]:
    
    def create_adjacency_list(bond_df: pd.DataFrame) -> Dict[int, List[int]]:
        adjacency_list = {i: [] for i in np.unique(bond_df.values)}
        for atom1, atom2 in bond_df.values:
            adjacency_list[atom1].append(atom2)
            adjacency_list[atom2].append(atom1)
        
        return adjacency_list

    adjacency_dict = create_adjacency_list(bond_df)
    
    # Check if pair atoms exist in the adjacency dict before filtering
    for atom in pair:
        if atom not in adjacency_dict:
            print(f"[create_adjacency_dict_for_pair] WARNING: atom {atom} not found in bond_df!")
    
    result = {atom: adjacency_dict.get(atom, []) for atom in pair}
   
    return result



def reindex_and_preserve(df, new_index_order):
        # Select rows that are in the new index order
        reindexed_part = df.loc[df.index.intersection(new_index_order)]

        # Select rows that are not in the new index order
        non_reindexed_part = df.loc[~df.index.isin(new_index_order)]

        # Concatenate the two parts
        return pd.concat([reindexed_part, non_reindexed_part])


def get_benzene_ring_indices(bonds_df, ring_atoms):
    """
    Identifies benzene ring indices from a bond dataframe and a set of ring atoms.
    Also detects and prints fused benzene rings (two rings sharing an edge).
    """
    # Read atom indices
    atom1_idx = ring_atoms[0]
    atom2_idx = ring_atoms[1] if len(ring_atoms) > 1 else None

    # Build the molecular graph
    G = nx.Graph()
    for _, row in bonds_df.iterrows():
        a1, a2 = int(row[0]), int(row[1])
        G.add_edge(a1, a2)

    # Find all simple cycles and filter 6-membered rings that include atom1_idx
    cycles = nx.cycle_basis(G)
    # print(f"Debug: Found {len(cycles)} cycles in the graph: {cycles}")
    benzene_rings = [cycle for cycle in cycles if len(cycle) == 6 and atom1_idx in cycle]

    if not benzene_rings:


        # Check if atom2_idx was provided and is expected to be in a ring
        if atom2_idx is not None:
            rings_with_atom2 = [i for i, ring in enumerate(benzene_rings) if atom2_idx in ring]
            print(f"Debug: Rings containing atom {atom2_idx}: {rings_with_atom2}")
            
            # Find rings containing both atoms (if atom2_idx was specified)
            rings_with_both = [i for i, ring in enumerate(benzene_rings) 
                    if atom1_idx in ring and atom2_idx in ring]
            if rings_with_both:
                print(f"Debug: Rings containing both atoms {atom1_idx} and {atom2_idx}: {rings_with_both}")
            else:
                print(f"Debug: No rings found containing both atoms {atom1_idx} and {atom2_idx}")
        
        return None

  

    # Detect fused rings: those sharing exactly two atoms
    if len(benzene_rings) > 1:
        print("\nDetected fused benzene ring pairs:")
        for i in range(len(benzene_rings)):
            for j in range(i + 1, len(benzene_rings)):
                shared = set(benzene_rings[i]).intersection(benzene_rings[j])
                if len(shared) == 2:
                    print(f"  Ring1: {benzene_rings[i]}\n  Ring2: {benzene_rings[j]}\n  Shared atoms: {sorted(shared)}\n")
                # take the second ring if it exists

    # fix for future
    selected_ring = benzene_rings[1] if len(benzene_rings) > 1 else benzene_rings[0]
    return (
        selected_ring[3],
        selected_ring[0],
        selected_ring[1],
        selected_ring[-1],
        selected_ring[2],
        selected_ring[4]
    )
