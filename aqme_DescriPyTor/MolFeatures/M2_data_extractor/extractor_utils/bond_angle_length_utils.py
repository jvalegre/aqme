from typing import List, Optional
import numpy as np
import numpy.typing as npt
import pandas as pd


try:
    from ...utils.help_functions import *
except ImportError:
    from utils.help_functions import * 




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






def get_bonds_vector_for_calc_angle(coordinates_array, atoms_indices):
    """
    Prepares bond vectors needed for calculating angles between atoms.
    This function takes atomic coordinates and indices of atoms, adjusts the indices,
    and creates pairs of indices to generate vectors between atoms. The function is 
    specifically designed to work with the calc_angle_between_atoms function.
    Parameters
    ----------
    coordinates_array : numpy.ndarray
        Array containing the 3D coordinates of all atoms in the molecule.
    atoms_indices : list
        List of atom indices (0-based or 1-based) for which bond vectors need to be calculated.
        Typically contains 3 indices for a single angle or 4 indices for a dihedral angle.
    Returns
    -------
    list
        List of bond vectors calculated from the coordinates of the specified atoms.
        For a 3-atom angle calculation, returns 2 bond vectors.
        For a 4-atom dihedral angle calculation, returns 3 bond vectors.
    Notes
    -----
    The function works by:
    1. Adjusting indices to ensure they are 0-based
    2. Creating pairs of indices representing connections between atoms
    3. Converting these pairs to actual coordinate vectors
    For 3 atoms (angle), it creates 2 vectors: from atom1→atom2 and atom2→atom3
    For 4 atoms (dihedral), it creates 3 vectors: atom1→atom2, atom2→atom3, and atom3→atom4
    """

    indices=adjust_indices(atoms_indices)#t
    augmented_indices=[indices[0],indices[1],indices[1],indices[2]]
    if len(indices)==4:
        augmented_indices.extend([indices[2],indices[3]])
    indices_pairs=list(zip(augmented_indices[::2],augmented_indices[1::2]))
  
    bond_vector=indices_to_coordinates_vector(coordinates_array,indices_pairs)
    return bond_vector
  

def calc_angle_between_atoms(coordinates_array,atoms_indices): #gets a list of atom indices
    """
    a function that gets 3/4 atom indices, and returns the angle between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atoms_indices- list of ints
        a list of atom indices to calculate the angle between- [2,3,4]
   
    Returns
    -------
    angle: float
        the bond angle between the atoms
    """
    bonds_list=get_bonds_vector_for_calc_angle(coordinates_array,atoms_indices)
    if len(atoms_indices)==3:
 
        angle=calc_angle(bonds_list[0], bonds_list[1]*(-1), degrees=True)
    else:
        first_cross=np.cross(bonds_list[0],bonds_list[1]*(-1))
        second_cross=np.cross(bonds_list[2]*(-1),bonds_list[1]*(-1)) 
        angle=calc_angle(first_cross, second_cross, degrees=True)
    return angle

def get_angle_df(coordinates_array, atom_indices):
    """
    a function that gets a list of atom indices, and returns a dataframe of angles between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
    
    atom_indices- list of lists of ints
        a list of atom indices to calculate the angle between- [[2,3,4],[2,3,4,5]]
    """
 
    if isinstance(atom_indices, list) and all(isinstance(elem, list) for elem in atom_indices):
        indices_list=['angle_{}'.format(index) if len(index)==3 else 'dihedral_{}'.format(index) for index in atom_indices]
        angle_list=[calc_angle_between_atoms(coordinates_array,index) for index in atom_indices]
        return pd.DataFrame(angle_list,index=indices_list)
    else:
        indices_list=['angle_{}'.format(atom_indices) if len(atom_indices)==3 else 'dihedral_{}'.format(atom_indices)]
        angle=[calc_angle_between_atoms(coordinates_array,atom_indices)]
        return pd.DataFrame(angle,index=indices_list)



def calc_single_bond_length(coordinates_array,atom_indices):
    """
    a function that gets 2 atom indices, and returns the distance between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atom_indices- list of ints
        a list of atom indices to calculate the distance between- [2,3]
   
    Returns
    -------
    distance: float
        the bond distance between the atoms
    """
    indices=adjust_indices(atom_indices)
    distance=np.linalg.norm(coordinates_array[indices[0]]-coordinates_array[indices[1]])
    return distance

def calc_bonds_length(coordinates_array,atom_pairs): 
    """
    a function that calculates the distance between each pair of atoms.
    help function for molecule class
    
    Parameters
    ----------
    coordinates_array: np.array
        xyz coordinates 
        
    atom_pairs : iterable
        list containing atom pairs-(([2,3],[4,5]))
        
    Returns
    -------
    pairs_df : dataframe
        distance between each pair 
        
        Output:
                               0
    bond length[2, 3]            1.525692
    bond length[4, 5]            2.881145
    
    """
    # Check the order of bond list and pairs
    bond_list=[calc_single_bond_length(coordinates_array,pair) for pair in atom_pairs]
    pairs=adjust_indices(atom_pairs)
    index=[('bond_length')+str(pair) for pair in pairs]
    pairs_df=pd.DataFrame(bond_list,index=index)
    return pairs_df



def calculate_bond_lengths_matrix(coords, connections_df):
    """
    Calculate a matrix of bond lengths between connected atoms in a molecule.
    This function computes the Euclidean distances between pairs of atoms that are 
    connected according to the connections dataframe. The result is a symmetric matrix
    where element [i,j] represents the bond length between atoms i and j.
    Parameters
    ----------
    coords : numpy.ndarray
        Array of shape (num_atoms, 3) containing the 3D coordinates of each atom in the molecule.
    connections_df : pandas.DataFrame
        DataFrame containing pairs of connected atoms, with columns [0, 1] representing
        the indices of the connected atoms.
    Returns
    -------
    numpy.ndarray
        A symmetric square matrix of shape (num_atoms, num_atoms) where each element [i,j]
        contains the bond length between atoms i and j if they are connected, and 0 otherwise.
    """
    num_atoms = coords.shape[0]
    bond_lengths = np.zeros((num_atoms, num_atoms))

    for _, row in connections_df.iterrows():
        atom1 = int(row[0])
        atom2 = int(row[1])
        length = np.linalg.norm(coords[atom1] - coords[atom2])
        bond_lengths[atom1, atom2] = length
        bond_lengths[atom2, atom1] = length  # Symmetric matrix

    return bond_lengths

def calculate_angles_matrix(coords, connections_df):
    num_atoms = coords.shape[0]
    angles = np.zeros((num_atoms, num_atoms, num_atoms))

    for i in range(num_atoms):
        for j in range(num_atoms):
            for k in range(num_atoms):
                if i != j and j != k and i != k:
                    if (connections_df[(connections_df[0] == i) & (connections_df[1] == j)].empty == False or 
                        connections_df[(connections_df[0] == j) & (connections_df[1] == i)].empty == False) and (
                        connections_df[(connections_df[0] == i) & (connections_df[1] == k)].empty == False or 
                        connections_df[(connections_df[0] == k) & (connections_df[1] == i)].empty == False):
                        
                        v1 = coords[j] - coords[i]
                        v2 = coords[k] - coords[i]
                        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                        angles[i, j, k] = np.arccos(np.clip(cos_angle, -1.0, 1.0))

    return angles