import numpy as np
import pandas as pd 
import numpy.typing as npt
import igraph as ig
import matplotlib.pyplot as plt
plt.ion() 

try:
    from ...utils.help_functions import *
    from ...utils.visualize import *
except ImportError:
    from utils.help_functions import * 
    from utils.visualize import *


class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    
    PYYKKO_RADII= {
            'H': 0.31, 'He': 0.28, 'Li': 1.28,
            'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
            'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
            'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
            'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
            'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
            'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
            'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
            'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
            'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
            'Am': 1.80, 'Cm': 1.69
    }
    
    BONDI_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98, 
        'N': 1.55, 'O': 1.52, 'Co': 2.00, 
        'Br': 1.83, 'Si': 2.10,'Ni': 2.00,
        'P': 1.80, 'Cl': 1.75, 
    }
    CPK_RADII = {
    'C': 1.50,
    'C3': 1.60,
    'C6/N6': 1.70,
    'H': 1.00,
    'N': 1.50,
    'N4': 1.45,
    'O': 1.35,
    'O2': 1.35,
    'P': 1.40,
    'S': 1.70,
    'S1': 1.00,
    'F': 1.35,
    'Cl': 1.80,
    'S4': 1.40,
    'Br': 1.95,
    'I': 2.15,
    'X': 1.92
}
    

def generate_circle(center_x, center_y, radius, n_points=20):
    """
    Generate circle coordinates given a center and radius.
    Returns a DataFrame with columns 'x' and 'y'.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)
    return np.column_stack((x, y))



    
def calc_new_base_atoms(coordinates_array: npt.ArrayLike, atom_indices: npt.ArrayLike):  #help function for calc_coordinates_transformation
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices[0] is list, compute the new origin as the middle of the first atoms.
    """
   
    if isinstance(atom_indices[0], list):
        new_origin=np.mean(coordinates_array[atom_indices[0]], axis=0)
    else:
        new_origin=coordinates_array[atom_indices[0]]
    new_y=(coordinates_array[atom_indices[1]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[1]]-new_origin))
    coplane=((coordinates_array[atom_indices[2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[2]]-new_origin)+0.00000001))
    
    return (new_origin,new_y,coplane)

def np_cross_and_vstack(plane_1, plane_2):
    cross_plane=np.cross(plane_1, plane_2)
    united_results=np.vstack([plane_1, plane_2, cross_plane])
    return united_results

from numpy.typing import ArrayLike

def calc_basis_vector(origin, y: ArrayLike, coplane: ArrayLike):
    """
    Calculate the new basis vector.
    
    Parameters
    ----------
    origin : array-like
        The origin of the new basis.
    y : array-like
        The new basis's y direction.
    coplane : array-like
        A vector coplanar with the new y direction.
    
    Returns
    -------
    new_basis : np.array
        The computed new basis matrix.
    """
    coef_mat = np_cross_and_vstack(coplane, y)
    angle_new_y_coplane = calc_angle(coplane, y)
    cop_ang_x = angle_new_y_coplane - (np.pi/2)
    result_vector = [np.cos(cop_ang_x), 0, 0]
    new_x, _, _, _ = np.linalg.lstsq(coef_mat, result_vector, rcond=None)
    new_basis = np_cross_and_vstack(new_x, y)
    return new_basis

def transform_row(row_array, new_basis, new_origin, round_digits):
    """
    Transform a single row of coordinates.
    
    Parameters
    ----------
    row_array : array-like
        A row of coordinates.
    new_basis : array-like
        New basis matrix.
    new_origin : array-like
        The new origin.
    round_digits : int
        Number of decimal places to round to.
    
    Returns
    -------
    np.array
        The transformed coordinates.
    """
    row_array = np.squeeze(row_array)
    translocated_row = row_array - new_origin
    ## make sure or change to shape (3,) if needed
    if translocated_row.shape != (3,):
        try:
            translocated_row = np.reshape(translocated_row, (3,))
        except Exception as e:
            print("transform_row: Error reshaping translocated_row:", e)
    
    result = np.dot(new_basis, translocated_row).round(round_digits)
    return result

def calc_coordinates_transformation(coordinates_array: ArrayLike, base_atoms_indices: ArrayLike, round_digits: int = 4, origin: ArrayLike = None) -> ArrayLike:
    """
    Transform a coordinates array to a new basis defined by the atoms in base_atoms_indices.
    
    Parameters
    ----------
    coordinates_array : np.array
        The original xyz molecule coordinates.
    base_atoms_indices : list or array-like
        Indices of atoms to define the new basis.
    round_digits : int, optional
        Number of digits to round the result.
    origin : array-like, optional
        A new origin for the transformation. If None, the first base atom is used.
    
    Returns
    -------
    transformed_coordinates : np.array
        The transformed coordinates.
    """
    # indices = adjust_indices(base_atoms_indices)
    # Calculate new basis using helper function calc_new_base_atoms (assumed to return origin, y, and coplane)
    new_basis = calc_basis_vector(*calc_new_base_atoms(coordinates_array, base_atoms_indices))
    if origin is None:
        new_origin = coordinates_array[base_atoms_indices[0]]
    else:
        new_origin = origin
   
    transformed_coordinates = np.apply_along_axis(lambda x: transform_row(x, new_basis, new_origin, round_digits), 1, coordinates_array)
  
    return transformed_coordinates


def get_center_of_mass(xyz_df):
    coordinates=np.array(xyz_df[['x','y','z']].values,dtype=float)
    atoms_symbol=np.array(xyz_df['atom'].values)
    masses=np.array([GeneralConstants.ATOMIC_WEIGHTS.value[symbol] for symbol in atoms_symbol])
    center_of_mass=np.sum(coordinates*masses[:,None],axis=0)/np.sum(masses)
    return center_of_mass

def get_closest_atom_to_center(xyz_df,center_of_mass):
    distances = np.sqrt((xyz_df['x'] - center_of_mass[0]) ** 2 + (xyz_df['y'] - center_of_mass[1]) ** 2 + (xyz_df['z'] - center_of_mass[2]) ** 2)
    idx_closest = np.argmin(distances)
    center_atom = xyz_df.loc[idx_closest]
    return center_atom

def get_sterimol_base_atoms(center_atom, bonds_df):
    
    center_atom_id=int(center_atom.name)+1
    base_atoms = [center_atom_id]
    if (any(bonds_df[0] == center_atom_id)):
        base_atoms.append(int(bonds_df[(bonds_df[0]==center_atom_id)][1].iloc[0]))
    else:
        base_atoms.append(int(bonds_df[(bonds_df[1]==center_atom_id)][0].iloc[0]))
    
    return base_atoms

def center_substructure(coordinates_array,atom_indices):
    atom_indices=adjust_indices(atom_indices)
    substructure=coordinates_array[atom_indices]
    center_substructure=np.mean(substructure,axis=0)
    return center_substructure



def get_sterimol_indices(coordinates,bonds_df):
    center=get_center_of_mass(coordinates)
    center_atom=get_closest_atom_to_center(coordinates,center)
    base_atoms=get_sterimol_base_atoms(center_atom,bonds_df)
    return base_atoms

def filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df):
    """
   
    """

    allowed_bonds_indices= pd.concat([bonded_atoms_df['index_1'],bonded_atoms_df['index_2']],axis=1).reset_index(drop=True)
    atom_filter=adjust_indices(np.unique([atom for sublist in allowed_bonds_indices.values.tolist() for atom in sublist]))
    # save current index as a list
    index_list=coordinates_df.loc[atom_filter].index.tolist()
    edited_coordinates_df=coordinates_df.loc[atom_filter].reset_index(drop=True)

    return edited_coordinates_df, index_list


def get_extended_df_for_sterimol(coordinates_df, bonds_df, radii='CPK'):
    """
    A function that adds information to the regular coordinates_df

    Parameters
    ----------
    coordinates_df : dataframe
    bond_type : str
        The bond type of the molecule
    radii : str, optional
        The type of F to use ('bondi' or 'CPK'), by default 'bondi'

    Returns
    -------
    dataframe
        The extended dataframe with additional columns

    """

    radii_map = GeneralConstants.CPK_RADII.value if radii == 'CPK' else GeneralConstants.BONDI_RADII.value
    df = coordinates_df.copy()  # make a copy of the dataframe to avoid modifying the original
    if radii == 'bondi':
        df['atype']=df['atom']
    else:
        df['atype']=nob_atype(coordinates_df, bonds_df)
    df['magnitude']  = np.linalg.norm((df[['x', 'z']].astype(float)), axis=1)
    df['radius'] = df['atype'].map(radii_map)
    df['B5'] = df['radius'] + df['magnitude']
    df['L'] = df['y'] + df['radius']

    return df


def get_transfomed_plane_for_sterimol(plane,degree):
    """
    a function that gets a plane and rotates it by a given degree
    in the case of sterimol the plane is the x,z plane.
    Parameters:
    ----------
    plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
    degree : float
    """
  
    cos_deg=np.cos(degree*(np.pi/180))
    sin_deg=np.sin(degree*(np.pi/180))
    rot_matrix=np.array([[cos_deg,-1*sin_deg],[sin_deg,cos_deg]])
    transformed_plane=np.vstack([np.matmul(rot_matrix,row) for row in plane]).round(4)

    
    ## return the inversed rotation matrix to transform the plane back

    return transformed_plane


def calc_B1(transformed_plane,avs,edited_coordinates_df,column_index):
    """
    Parameters
    ----------
    transformed_plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
            [-0.7384 -0.5135]
            [-0.3759 -0.271 ]
            [-1.1046 -0.8966]
            [ 0.6763  0.5885]
    avs : list
        the max & min of the [x,z] columns from the transformed_plane.
        example:[0.6763, -1.1046, 0.5885, -0.8966
                 ]
    edited_coordinates_df : TYPE
        DESCRIPTION.
    column_index : int
        0 or 1 depending- being used for transformed plane.
    """
    
    ## get the index of the min value in the column compared to the avs.min
    idx=np.where(np.isclose(np.abs(transformed_plane[:,column_index]),(avs.min())))[0][0]  ## .round(4)
   
    if transformed_plane[idx,column_index]<0:
       
        new_idx=np.where(np.isclose(transformed_plane[:,column_index],transformed_plane[:,column_index].min()))[0][0]
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[new_idx,column_index],
                                 transformed_plane[:,column_index]<=transformed_plane[new_idx,column_index]+1)
        
        transformed_plane[:,column_index]=-transformed_plane[:,column_index]
    else:
     
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[idx,column_index]-1,
                                 transformed_plane[:,column_index]<=transformed_plane[idx,column_index])
        
    against,against_loc=[],[]
    B1,B1_loc=[],[]

    ### return the part of the transformed plane that b1 is calculated from
    ### convert it back with the inverse matrix to get the original coordinates of b1 location
    
    for i in range(0,transformed_plane.shape[0]): 
        if bool_list[i]:
            against.append(np.array(transformed_plane[i,column_index]+edited_coordinates_df['radius'].iloc[i]))
            against_loc.append(edited_coordinates_df['L'].iloc[i])
        

        if len(against)>0:
        
            B1.append(max(against))
            B1_loc.append(against_loc[against.index(max(against))])
          
        else:
       
            B1.append(np.abs(transformed_plane[idx,column_index]+edited_coordinates_df['radius'].iloc[idx]))
            B1_loc.append(edited_coordinates_df['radius'].iloc[idx])
            
            
      
    return [B1,B1_loc] 

def b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane, b1_planes):
    """
    For each degree in degree_list, rotate the plane and compute B1 values and the B1-B5 angle.
    Instead of returning after the first iteration, this version accumulates all results
    into a DataFrame.
    
    Parameters
    ----------
    extended_df : pd.DataFrame
        DataFrame with at least columns 'x', 'z', 'radius', 'L'.
    b1s : list
        (Unused here; kept for compatibility)
    b1s_loc : list
        (Unused here; kept for compatibility)
    degree_list : list
        List of rotation angles (in degrees) to scan.
    plane : np.array
        Array of shape (n_points_total, 2) that contains the combined circle points.
    b1_planes : list
        List to store the rotated plane from each iteration.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with one row per degree, containing:
         - 'degree': the rotation angle,
         - 'B1': the minimum extreme value from the rotated plane,
         - 'B1_B5_angle': the angle between the B1 and B5 arrows in degrees,
         - 'b1_coords': the coordinate (as a tuple) chosen for B1,
         - 'b5_value': the distance of the farthest point (B5) from the origin.
    """
    results = []  # List to accumulate results
    
    for degree in degree_list:
        transformed_plane = get_transfomed_plane_for_sterimol(plane, degree)  # Rotate the plane
     
        max_x = np.max(transformed_plane[:, 0])
        min_x = np.min(transformed_plane[:, 0])
        max_y = np.max(transformed_plane[:, 1])
        min_y = np.min(transformed_plane[:, 1])
        avs = np.abs([max_x, min_x, max_y, min_y])
       
        min_val = np.min(avs)
        min_index = np.argmin(avs)
        
        # Mimic R's switch to pick the B1 coordinate:
        if min_index == 0:
            b1_coords = (max_x, 0)
        elif min_index == 1:
            b1_coords = (min_x, 0)
        elif min_index == 2:
            b1_coords = (0, max_y)
        else:
            b1_coords = (0, min_y)
        
        # Determine B5 as the farthest point from the origin.
        norms_sq = np.sum(transformed_plane**2, axis=1)
        b5_index = np.argmax(norms_sq)
        b5_point = transformed_plane[b5_index]
        b5_value = np.linalg.norm(b5_point)
        
        # Calculate angles for the arrows.
        angle_b1 = np.arctan2(b1_coords[1], b1_coords[0]) % (2*np.pi)
        angle_b5 = np.arctan2(b5_point[1], b5_point[0]) % (2*np.pi)
        angle_diff = abs(angle_b5 - angle_b1)
        if angle_diff > np.pi:
            angle_diff = 2*np.pi - angle_diff
        # Convert the angle difference to degrees.
        B1_B5 = np.degrees(angle_diff)
        B1 = min_val
        
        # Save the transformed plane for this iteration.
        b1_planes.append(transformed_plane)
        
        # Accumulate the result for this degree.
        results.append({
            'degree': degree,
            'B1': B1,
            'B1_B5_angle': B1_B5,
            'b1_coords': b1_coords,
            'b5_value': b5_value,
            'plane': transformed_plane
        })
    
    # Create a DataFrame from the accumulated results.
    sterimol_df = pd.DataFrame(results)
  
    return sterimol_df


def preform_coordination_transformation(xyz_df, indices=None, origin=None):
    """
    Perform a coordination transformation on the xyz DataFrame.
    
    Parameters
    ----------
    xyz_df : pd.DataFrame
        DataFrame containing columns 'x', 'y', 'z'.
    indices : array-like, optional
        Atom indices to use for the new basis. If None, default indices [1,2,3] are used.
    
    Returns
    -------
    xyz_copy : pd.DataFrame
        DataFrame with transformed coordinates.
    """
    xyz_copy = xyz_df.copy()
    coordinates = np.array(xyz_copy[['x', 'y', 'z']].values)
    if origin is not None:
       
        geometric_origin=np.mean(xyz_df.iloc[indices][['x','y','z']].values,axis=0)
    else:
        geometric_origin=None
    if indices is None:
        transformed = calc_coordinates_transformation(coordinates, [1,2,3], origin=geometric_origin)
    else:
   
        transformed = calc_coordinates_transformation(coordinates, indices, origin=geometric_origin)
    
    xyz_copy[['x', 'y', 'z']] = transformed
    return xyz_copy



def direction_atoms_for_sterimol(bonds_df, base_atoms) -> list:
    """
    Return [origin, direction, third] for Sterimol coordinate setup.
    - bonds_df: DataFrame with two columns (0,1), undirected bonds between integer atom indices.
    - base_atoms: [origin, direction] or [origin, direction, third].
       If third == origin, choose a neighbor of 'direction' different from origin.
    Raises ValueError with a clear message if it can't determine 'third'.
    """
    import pandas as pd

    if not isinstance(bonds_df, pd.DataFrame) or not set(bonds_df.columns[:2]) >= {0, 1}:
        raise ValueError("bonds_df must be a DataFrame with columns 0 and 1 listing bonds.")
    
    # if base atoms list of list do base_atoms=base_atoms[0]
    if isinstance(base_atoms[0], list):
        base_atoms=base_atoms[0]

    if len(base_atoms) < 2:
        raise ValueError(f"base_atoms must contain at least two indices [origin, direction]; got: {base_atoms}")

    # Normalize to ints
    origin, direction = int(base_atoms[0]), int(base_atoms[1])

    # Build neighbor sets
    # neighbors(n): all atoms bonded to n (undirected)
    col0, col1 = bonds_df[0].astype(int), bonds_df[1].astype(int)
    nbr_dir = set(col1[col0 == direction].astype(int)).union(set(col0[col1 == direction].astype(int)))
    nbr_org = set(col1[col0 == origin].astype(int)).union(set(col0[col1 == origin].astype(int)))

    # Helper to pick a 'third' atom:
    def pick_third():
        # Prefer neighbor of direction that is not origin
        cand = [a for a in nbr_dir if a != origin]
        if cand:
            return int(cand[0])
        # If direction is terminal (only bonded to origin), try neighbor of origin that is not direction
        cand2 = [a for a in nbr_org if a != direction]
        if cand2:
            return int(cand2[0])
        # Nothing usable
        return None

    # If user supplied a 3rd atom explicitly
    if len(base_atoms) >= 3:
        third = int(base_atoms[2])
        if third == origin:
            third = pick_third()
            if third is None:
                raise ValueError(
                    f"Cannot determine third atom: 'direction'={direction} has no neighbor != origin={origin}, "
                    f"and 'origin' has no neighbor != direction."
                )
            return [origin, direction, third]
        else:
            # If the provided third is bonded to direction, accept it; otherwise try to fix gracefully
            if (third in nbr_dir) and (third != origin):
                return [origin, direction, third]
            # Try automatic pick if given third is not valid
            auto = pick_third()
            if auto is not None:
                return [origin, direction, auto]
            raise ValueError(
                f"Provided third atom {third} is not bonded to 'direction'={direction}; "
                f"also failed to auto-pick a neighbor."
            )

    # Only two atoms given → pick a third automatically
    third = pick_third()
    if third is None:
        raise ValueError(
            f"Cannot determine third atom from bonds: 'direction'={direction} neighbors={sorted(nbr_dir)}, "
            f"'origin'={origin} neighbors={sorted(nbr_org)}."
        )

    return [origin, direction, third]



    

def get_molecule_connections(bonds_df, source, direction, mode='all'):
    # Ensure correct column names for igraph
    bonds_df = bonds_df.copy()
    bonds_df.columns = [0, 1]

    # Ensure integer atom indices
    bonds_df[[0, 1]] = bonds_df[[0, 1]].apply(pd.to_numeric, errors='raise').astype(int)

    # Build graph
    graph = ig.Graph.DataFrame(edges=bonds_df, directed=False)

    # If edge (source, direction) doesn’t exist → add it
    if not graph.are_connected(source, direction):
        graph.add_edge(source, direction)

    # Get all simple paths from source
    paths = graph.get_all_simple_paths(v=source, mode='all')

    # Keep only paths starting with [source, direction]
    paths_with_start = [path for path in paths 
                        if len(path) >= 2 and path[0] == source and path[1] == direction]

    if mode == 'all':
        longest_path = np.unique(flatten_list(paths_with_start))
        
        return longest_path

    elif mode == 'shortest':
        if paths_with_start:
            shortest_path = min(paths_with_start, key=len)
            return shortest_path
        else:
            return []
        
def find_longest_simple_path(bonds_df):
    """
    bonds_df : pandas.DataFrame with exactly two columns [u, v] (integer atom indices)
    
    Returns
    -------
    List[int]
        The sequence of atom indices along the longest simple path (graph diameter)
        in the undirected molecule graph.
    """
    # 1) Build an undirected igraph from your bond list
    g = ig.Graph.DataFrame(edges=bonds_df, directed=False)
    
    # 2) Compute the diameter: the furthest‐apart pair of nodes and the path between them
    diameter_path = g.get_diameter()  # returns a list of vertex indices
    
    # 3) Convert to plain Python ints and return
    return [int(v) for v in diameter_path]

def get_specific_bonded_atoms_df(bonds_df, longest_path, coordinates_df):
    """
    Returns a DataFrame of atoms bonded along the longest path.
    Includes debug prints to trace internal state.

    Parameters
    ----------
    bonds_df : pd.DataFrame
        DataFrame of bonds, with atom indices (1-based)
    longest_path : list or None
        List of atom indices to include (1-based)
    coordinates_df : pd.DataFrame
        DataFrame of atomic coordinates, must include 'atom' column

    Returns
    -------
    pd.DataFrame
        A DataFrame with bonded atom names and their original bond indices
    """

    # 1) Filter bonds based on longest_path
    if longest_path is not None:
   
        mask = bonds_df.isin(longest_path)
  
        edited_bonds_df = bonds_df[mask.all(axis=1)].dropna().reset_index(drop=True)
       
    else:
        edited_bonds_df = bonds_df.copy()

    if edited_bonds_df.empty:
        return pd.DataFrame(columns=XYZConstants.BONDED_COLUMNS.value)


    # 2) Build bonds_array (zero-based indices)
    bonds_array = (edited_bonds_df.values - 1).astype(int)
    # 3) Lookup atom names for each bond
    atom_bonds_list = []
    for i, bond in enumerate(bonds_array):
        try:
            coords = coordinates_df.iloc[bond]['atom'].values
            atom_bonds_list.append(coords)
        except Exception as e:

            print(f"Error retrieving atom names for bond {i}: {bond} - {e}")

    # 4) Stack into shape (n_bonds, 2)
    atom_bonds = np.vstack(atom_bonds_list).reshape(-1, 2)

    # 5) Concatenate with edited_bonds_df
    bonded_atoms_df = pd.concat([pd.DataFrame(atom_bonds), edited_bonds_df], axis=1)
    bonded_atoms_df.columns = XYZConstants.BONDED_COLUMNS.value
    return bonded_atoms_df
    
def get_b1s_list(extended_df, scans=90//5,plot_result=False):
    """
    Calculate B1 values by scanning over a range of rotation angles.
    Instead of using only the center points, this version generates circle points
    for each substituent from extended_df and then stacks them into one plane.
    
    Parameters
    ----------
    extended_df : pd.DataFrame
        DataFrame with at least the columns: 'x', 'z', 'radius', 'L'.
    scans : int, optional
        Degree step for the initial scan.
        
    Returns
    -------
    tuple
        (np.array of B1 values, np.array of B1 location values, list of rotated planes)
    """
    b1s, b1s_loc, b1_planes = [], [], []
    
    degree_list = list(range(18, 108, scans))
    circles_y = []
    circles = []

    for idx, row in extended_df.iterrows():
        
        circle_points = generate_circle(row['x'], row['z'], row['radius'], n_points=100)
        circles.append(circle_points)
        circle_points_y = generate_circle(row['x'], row['z'], row['radius'], n_points=100)
        circles_y.append(circle_points_y)
    plane_xz = np.vstack(circles)  # All circle points combined.
    plane_yz = np.vstack(circles_y)  # All circle points combined.

    sterimol_df=b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane_xz, b1_planes)

    b1s=sterimol_df['B1']
    
    try:
        back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
        front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans
        new_degree_list=range(back_ang,front_ang+1)
    except:
        
        back_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]-scans
        front_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]+scans
        new_degree_list=range(back_ang,front_ang+1)
    # plane=sterimol_df['plane']
    best_idx = np.where(sterimol_df['B1']==sterimol_df['B1'].min())[0][0]
    # 4) *you should* take that rotated plane:
    plane_xz = sterimol_df.loc[best_idx, 'plane']

    b1s, b1s_loc, b1_planes = [], [], []
    sterimol_df=b1s_for_loop_function(extended_df, b1s, b1s_loc, list(new_degree_list), plane_xz, b1_planes)

    b1s=sterimol_df['B1']

    b1_b5_angle=sterimol_df['B1_B5_angle']
    plane_xz=sterimol_df['plane']
    

    return [b1s, b1_b5_angle, plane_xz]




def calc_sterimol(bonded_atoms_df,extended_df,visualize_bool=False):
    edited_coordinates_df,index_list=filter_atoms_for_sterimol(bonded_atoms_df,extended_df)
    b1s,b1_b5_angle,plane=get_b1s_list(edited_coordinates_df)
    valid_indices = np.where(b1s >= 0)[0]
    best_idx = valid_indices[np.argmin(b1s[valid_indices])]
    best_b1_plane = plane[best_idx]
    max_x = np.max(best_b1_plane[:, 0])
    min_x = np.min(best_b1_plane[:, 0])
    max_y = np.max(best_b1_plane[:, 1])
    min_y = np.min(best_b1_plane[:, 1])
    avs = np.abs([max_x, min_x, max_y, min_y])
    min_val = np.min(avs)
    B1= min_val
    b1_index=np.where(b1s==B1)[0][0]
    angle=b1_b5_angle[b1_index]
    norms_sq = np.sum(best_b1_plane**2, axis=1)
    b5_index = np.argmax(norms_sq)
    b5_point = best_b1_plane[b5_index]
    b5_value = np.linalg.norm(b5_point)
    max_row=edited_coordinates_df['B5'].idxmax()
    max_row=int(max_row)
    L=max(edited_coordinates_df['L'].values)
    loc_B5 = edited_coordinates_df['y'].iloc[np.where(edited_coordinates_df['B5']==max(edited_coordinates_df['B5']))[0][0]]

    # --- B1 direction in the X-Z perpendicular plane ---
    min_index = np.argmin(avs)
    if min_index == 0:
        b1_xz = np.array([max_x, 0.0])
    elif min_index == 1:
        b1_xz = np.array([min_x, 0.0])
    elif min_index == 2:
        b1_xz = np.array([0.0, max_y])
    else:
        b1_xz = np.array([0.0, min_y])

    # 3D vectors in the Sterimol frame (Y = primary/L axis, X-Z = perpendicular plane)
    vis_params = {
        'B1_coords': [float(b1_xz[0]), 0.0, float(b1_xz[1])],
        'B5_coords': [float(b5_point[0]), 0.0, float(b5_point[1])],
        'L_coords':  [0.0, float(L), 0.0],
        'B1_value':  float(B1),
        'B5_value':  float(b5_value),
        'L_value':   float(L),
    }

    sterimol_df = pd.DataFrame([B1, b5_value, L ,loc_B5,angle], index=XYZConstants.STERIMOL_INDEX.value)
    if visualize_bool:
        # change edit coordinates_df index to original index
        edited_coordinates_df.index = [index_list[i] + 1 for i in range(len(index_list))]
        plot_b1_visualization(best_b1_plane, edited_coordinates_df=edited_coordinates_df)
        plt.show()
        print('B1 B5 Plane')
        # plot_L_B5_plane(edited_coordinates_df, sterimol_df.T)
        # plt.show()
        # print('L B5 Plane')


    return sterimol_df.T, vis_params


def get_sterimol_df(coordinates_df, bonds_df, base_atoms, connected_from_direction, 
                    radii='CPK', sub_structure=True, drop_atoms=None, visualize_bool=False, mode='all'):
    # Drop atoms and update indices
    if drop_atoms is not None:
  
        drop_atoms_adjusted = adjust_indices(drop_atoms)
        coordinates_df = coordinates_df.drop(drop_atoms_adjusted)
        coordinates_df = coordinates_df.reset_index(drop=True)
        
        # Create a mapping from old to new indices
        old_to_new = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted(coordinates_df.index))}
        # Map bonds_df to new indices, drop any bonds involving dropped atoms
        mask = ~bonds_df[0].isin(drop_atoms) & ~bonds_df[1].isin(drop_atoms)
  
        bonds_df = bonds_df[mask].copy()
        bonds_df[0] = bonds_df[0].map(old_to_new)
        bonds_df[1] = bonds_df[1].map(old_to_new)
        bonds_df = bonds_df.dropna().reset_index(drop=True)  # Drop any rows with NaN values after mapping
        # Update base_atoms as well
        base_atoms = [old_to_new[a] for a in base_atoms if a in old_to_new]
        if connected_from_direction is not None:
            connected_from_direction = [old_to_new[a] for a in connected_from_direction if a in old_to_new]
    
    bonds_direction = adjust_indices(direction_atoms_for_sterimol(bonds_df, base_atoms))
    if isinstance(base_atoms[0], list):
        base_atoms=base_atoms[0]
        # print(f'Debug! base_atoms adjusted to: {base_atoms}')
    # verify bonds direction is 3 atoms
    if len(bonds_direction) != 3:
        raise ValueError("Bonds direction must contain exactly 3 atoms, check molecule connectivity.")
    new_coordinates_df = preform_coordination_transformation(coordinates_df, bonds_direction)
    
    if sub_structure:
        if connected_from_direction is None:
            connected_from_direction = get_molecule_connections(bonds_df, base_atoms[0], base_atoms[1], mode=mode)
    else:
        connected_from_direction = None
  
    bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction, new_coordinates_df)
    extended_df = get_extended_df_for_sterimol(new_coordinates_df, bonds_df, radii)

    sterimol_df, vis_params = calc_sterimol(bonded_atoms_df, extended_df, visualize_bool)
    sterimol_df = sterimol_df.rename(index={0: str(base_atoms[0]) + '-' + str(base_atoms[1])})
    sterimol_df = sterimol_df.round(4)

    return sterimol_df