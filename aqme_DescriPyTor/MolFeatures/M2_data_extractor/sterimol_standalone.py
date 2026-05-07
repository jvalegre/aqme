import pandas as pd
import numpy as np
import os
import sys
import math
from enum import Enum
import igraph as ig

from typing import *

import warnings
from scipy.spatial.distance import pdist, squareform

from sklearn.preprocessing import MinMaxScaler
from morfeus import Sterimol, read_xyz
warnings.filterwarnings("ignore", category=RuntimeWarning)


def flatten_list(nested_list_arg: List[list]) -> List:
    """
    Flatten a nested list.
    turn [[1,2],[3,4]] to [1,2,3,4]
    """
    flat_list=[item for sublist in nested_list_arg for item in sublist]
    return flat_list

def split_strings(strings_list):
    split_list = []
    for string in strings_list:
        split_list.extend(string.split())
    return split_list

def get_df_from_file(filename,columns=['atom','x','y','z'],index=None):
    """
    Parameters
    ----------
    filename : str
        full file name to read.
    columns : str , optional
        list of column names for DataFrame. The default is None.
    splitter : str, optional
        input for [.split().] , for csv-',' for txt leave empty. The default is None.
    dtype : type, optional
        type of variables for dataframe. The default is None.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    with open(filename, 'r') as f:
        lines=f.readlines()[2:]
    splitted_lines=split_strings(lines)
    df=pd.DataFrame(np.array(splitted_lines).reshape(-1,4),columns=columns,index=index)
    df[['x','y','z']]=df[['x','y','z']].astype(float)
    return df


class XYZConstants(Enum):
    """
    Constants related to XYZ file processing
    """
    DF_COLUMNS=['atom','x','y','z']
    STERIMOL_INDEX = ['B1', 'B5', 'L', 'loc_B1','loc_B5']
    DIPOLE_COLUMNS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole']
    RING_VIBRATION_COLUMNS = ['cross', 'cross_angle', 'para', 'para_angle']
    RING_VIBRATION_INDEX=['Product','Frequency','Sin_angle']
    VIBRATION_INDEX = ['Frequency', 'Amplitude']
    BONDED_COLUMNS = ['atom_1', 'atom_2', 'index_1', 'index_2']
    NOF_ATOMS = ['N', 'O', 'F']
    STERIC_PARAMETERS = ['B1', 'B5', 'L', 'loc_B1', 'loc_B5','RMSD']
    ELECTROSTATIC_PARAMETERS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole','energy']

class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    COVALENT_RADII= {
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
        'N': 1.55, 'O': 1.52,
        'Br': 1.83, 'Si': 2.10,
        'P': 1.80, 'Cl': 1.75,
        # alkali / alkaline-earth (Alvarez 2013, DOI: 10.1039/c3dt50599e)
        'Li': 1.82, 'Be': 1.53, 'Na': 2.27, 'Mg': 1.73,
        'K': 2.75, 'Ca': 2.31, 'Rb': 3.03, 'Sr': 2.49,
        'Cs': 3.43, 'Ba': 2.68,
        # p-block metals / heavier metalloids
        'Al': 1.84, 'Ga': 1.87, 'Ge': 2.11, 'As': 1.85,
        'In': 1.93, 'Sn': 2.17, 'Sb': 2.06, 'Te': 2.06,
        'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07,
        # first-row transition metals
        'Sc': 2.18, 'Ti': 2.11, 'V': 2.07, 'Cr': 2.06,
        'Mn': 2.05, 'Fe': 2.04, 'Co': 2.00, 'Ni': 1.97,
        'Cu': 1.96, 'Zn': 2.01,
        # second-row transition metals
        'Y': 2.32, 'Zr': 2.23, 'Nb': 2.18, 'Mo': 2.17,
        'Tc': 2.16, 'Ru': 2.13, 'Rh': 2.10, 'Pd': 2.10,
        'Ag': 2.11, 'Cd': 2.18,
        # third-row transition metals
        'Hf': 2.23, 'Ta': 2.22, 'W': 2.18, 'Re': 2.16,
        'Os': 2.16, 'Ir': 2.13, 'Pt': 2.09, 'Au': 2.14,
        'Hg': 2.23,
        # lanthanides
        'La': 2.43, 'Ce': 2.42, 'Pr': 2.40, 'Nd': 2.39,
        'Pm': 2.38, 'Sm': 2.36, 'Eu': 2.35, 'Gd': 2.34,
        'Tb': 2.33, 'Dy': 2.31, 'Ho': 2.30, 'Er': 2.29,
        'Tm': 2.27, 'Yb': 2.26, 'Lu': 2.24,
        # actinides
        'Ac': 2.47, 'Th': 2.45, 'Pa': 2.43, 'U': 2.41,
        'Np': 2.39, 'Pu': 2.37, 'Am': 2.35, 'Cm': 2.35,
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
    'X': 1.92,
    # metals — element symbols used directly (Alvarez 2013 VdW values)
    'Li': 1.82, 'Be': 1.53, 'Na': 2.27, 'Mg': 1.73,
    'K': 2.75, 'Ca': 2.31, 'Rb': 3.03, 'Sr': 2.49,
    'Cs': 3.43, 'Ba': 2.68,
    'Al': 1.84, 'Ga': 1.87, 'Ge': 2.11, 'In': 1.93,
    'Sn': 2.17, 'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07,
    'Sc': 2.18, 'Ti': 2.11, 'V': 2.07, 'Cr': 2.06,
    'Mn': 2.05, 'Fe': 2.04, 'Co': 2.00, 'Ni': 1.97,
    'Cu': 1.96, 'Zn': 2.01,
    'Y': 2.32, 'Zr': 2.23, 'Nb': 2.18, 'Mo': 2.17,
    'Tc': 2.16, 'Ru': 2.13, 'Rh': 2.10, 'Pd': 2.10,
    'Ag': 2.11, 'Cd': 2.18,
    'Hf': 2.23, 'Ta': 2.22, 'W': 2.18, 'Re': 2.16,
    'Os': 2.16, 'Ir': 2.13, 'Pt': 2.09, 'Au': 2.14,
    'Hg': 2.23,
    'La': 2.43, 'Ce': 2.42, 'Pr': 2.40, 'Nd': 2.39,
    'Pm': 2.38, 'Sm': 2.36, 'Eu': 2.35, 'Gd': 2.34,
    'Tb': 2.33, 'Dy': 2.31, 'Ho': 2.30, 'Er': 2.29,
    'Tm': 2.27, 'Yb': 2.26, 'Lu': 2.24,
    'Ac': 2.47, 'Th': 2.45, 'Pa': 2.43, 'U': 2.41,
    'Np': 2.39, 'Pu': 2.37, 'Am': 2.35, 'Cm': 2.35,
}
    # CPK_RADII={
    #     'C':1.50,   'H':1.00,   'S.O':1.70,  'Si':2.10,
    #     'C2':1.60,  'N':1.50,   'S1':1.00,   'Co':2.00,
    #     'C3':1.60,  'C66':1.70, 'F':1.35,    'Ni':2.00,
    #     'C4':1.50,  'N4':1.45,  'Cl':1.75,
    #     'C5/N5':1.70, 'O':1.35, 'S4':1.40,
    #     'C6/N6':1.70, 'O2':1.35, 'Br':1.95,
    #     'C7':1.70,    'P':1.40,  'I':2.15,
    #     'C8':1.50,    'S':1.70,  'B':1.92,
    
    # }

    REGULAR_BOND_TYPE = {

        'O.2': 'O', 'N.2': 'N', 'S.3': 'S',
        'O.3': 'O', 'N.1': 'N', 'S.O2': 'S',
        'O.co2': 'O', 'N.3': 'N', 'P.3': 'P',
        'C.1': 'C', 'N.ar': 'N',
        'C.2': 'C', 'N.am': 'N',
        "C.cat": 'C', 'N.pl3': 'N',
        'C.3': 'C', 'N.4': 'N',
        'C.ar': 'C', 'S.2': 'S',
    }

    BOND_TYPE={
        
        'O.2':'O2', 'N.2':'C6/N6','S.3':'S4',
        'O.3':'O', 'N.1':'N', 'S.O2':'S',
        'O.co2':'O', 'N.3':'C6/N6','P.3':'P',
        'C.1':'C', 'N.ar':'C6/N6',
        'C.2':'C3', 'N.am':'C6/N6',
        "C.cat":'C3', 'N.pl3':'C6/N6',
        'C.3':'C', 'N.4':'N4',
        'C.ar':'C6/N6', 'S.2':'S','H':'H' 
        }
    
    ATOMIC_NUMBERS ={
    '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
             '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}
        

    ATOMIC_WEIGHTS = {
            'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
            'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
            'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
            'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
            'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
            'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
            'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
            'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
            'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
            'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
            'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
            'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
            'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
            'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
            'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
            'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
            'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
            'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
            'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
            'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
            'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
            'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
            'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
            'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
    }

_METAL_ELEMENTS = {
    'Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba', 'Fr', 'Ra',
    'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
}

import numpy.typing as npt


def adjust_indices(indices: npt.ArrayLike, adjustment_num: int=1) -> npt.ArrayLike:
    """
    adjust indices by adjustment_num
    """
    return np.array(indices)-adjustment_num

def adjust_indices_xyz(indices: npt.ArrayLike) -> npt.ArrayLike:
    """
    adjust indices by adjustment_num
    """
    return adjust_indices(indices, adjustment_num=1)

def calc_angle(p1: npt.ArrayLike, p2: npt.ArrayLike, degrees: bool=False) -> float: ###works, name in R: 'angle' , radians
    dot_product=np.dot(p1, p2)
    norm_p1=np.linalg.norm(p1)
    norm_p2=np.linalg.norm(p2)
    thetha=np.arccos(dot_product/(norm_p1*norm_p2))
    if degrees:
        thetha=np.degrees(thetha)   
    return thetha
    
def calc_new_base_atoms(coordinates_array: npt.ArrayLike, atom_indices: npt.ArrayLike):  #help function for calc_coordinates_transformation
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices is 4, the origin will be the middle of the first two atoms.
    """
    new_origin=coordinates_array[atom_indices[0]]
    if (len(atom_indices)==4):
        new_origin=(new_origin+coordinates_array[atom_indices[1]])/2
    new_y=(coordinates_array[atom_indices[-2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-2]]-new_origin))
    coplane=((coordinates_array[atom_indices[-1]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-1]]-new_origin)+0.00000001))
    return (new_origin,new_y,coplane)

def np_cross_and_vstack(plane_1, plane_2):
    cross_plane=np.cross(plane_1, plane_2)
    united_results=np.vstack([plane_1, plane_2, cross_plane])
    return united_results

def calc_basis_vector(origin, y: npt.ArrayLike, coplane: npt.ArrayLike):#help function for calc_coordinates_transformation
    """
    origin: origin of the new basis
    y: y direction of the new basis
    coplane: a vector that is coplanar with the new y direction
    """
    # cross_y_plane=np.cross(coplane,y)
    # coef_mat=np.vstack([y, coplane, cross_y_plane])
    coef_mat=np_cross_and_vstack(coplane, y)
    angle_new_y_coplane=calc_angle(coplane,y)
    cop_ang_x=angle_new_y_coplane-(np.pi/2)
    # result_vector=[0,np.cos(cop_ang_x),0]
    result_vector=[np.cos(cop_ang_x), 0, 0]
    new_x,_,_,_=np.linalg.lstsq(coef_mat,result_vector,rcond=None)
    new_basis=np_cross_and_vstack(new_x, y)
    # new_z=np.cross(new_x,y)
    # new_basis=np.vstack([new_x, y, new_z])
    return new_basis

def transform_row(row_array, new_basis, new_origin, round_digits):
    translocated_row = row_array - new_origin
    return np.dot(new_basis, translocated_row).round(round_digits)



def calc_coordinates_transformation(coordinates_array: npt.ArrayLike, base_atoms_indices: npt.ArrayLike, round_digits:int=4 ,origin:npt.ArrayLike=None) -> npt.ArrayLike:#origin_atom, y_direction_atom, xy_plane_atom
    """
    a function that recives coordinates_array and new base_atoms_indices to transform the coordinates by
    and returns a dataframe with the shifted coordinates
    parameters:
    ----------
    coordinates_array: np.array
        xyz molecule array
    base_atoms_indices: list of nums
        indices of new atoms to shift coordinates by.
    origin: in case you want to change the origin of the new basis, middle of the ring for example. used in npa_df
    returns:
        transformed xyz molecule dataframe
    -------
        
    example:
    -------
    calc_coordinates_transformation(coordinates_array,[2,3,4])
    
    Output:
        atom       x       y       z
      0    H  0.3477 -0.5049 -1.3214
      1    B     0.0     0.0     0.0
      2    B    -0.0  1.5257     0.0
    """
    indices=adjust_indices_xyz(base_atoms_indices)
    new_basis=calc_basis_vector(*calc_new_base_atoms(coordinates_array,indices))    
    if origin is None:
        new_origin=coordinates_array[indices[0]]
    else:
        new_origin=origin

    transformed_coordinates = np.apply_along_axis(lambda x: transform_row(x, new_basis, new_origin, round_digits), 1,
                                                  coordinates_array)
    # transformed_coordinates=np.array([np.dot(new_basis,(row-new_origin)) for row in coordinates_array]).round(round_digits)
    return transformed_coordinates

def preform_coordination_transformation(xyz_df, indices=None):
    xyz_copy=xyz_df.copy()
    coordinates=np.array(xyz_copy[['x','y','z']].values)
    if indices is None:
        xyz_copy[['x','y','z']]=calc_coordinates_transformation(coordinates, [1,2,3])
    else:
      
        xyz_copy[['x','y','z']]=calc_coordinates_transformation(coordinates, indices)
 
    return xyz_copy

def calc_npa_charges(coordinates_array: npt.ArrayLike,charge_array: npt.ArrayLike,  geom_transform_indices=None):##added option for subunits
    """
    a function that recives coordinates and npa charges, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    
    Parameters
    ---------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    charge_array: np.array
        array of npa charges
    base_atoms_indices:list
        3/4 atom indices for coordinates transformation
        
    optional-sub_atoms:list
        calculate npa charges from a set of sub_atoms instead of all atoms.       
    Returns:
    -------
    dipole_df=calc_npa_charges(coordinates_array,charges,base_atoms_indices,sub_atoms)
    Output:
    dipole_df : pd.DataFrame
        output:            dip_x     dip_y     dip_z     total
                       0  0.097437 -0.611775  0.559625  0.834831
    """
    # what is NPA here??
    # indices=adjust_indices(base_atoms_indices)
    # transformed_coordinates=calc_coordinates_transformation(coordinates_array, indices)
    # if sub_atoms:
    #     atom_mask=sub_atoms
    # else:
    #     atom_mask=range(len(charge_array))
    # atom_mask=range(charge_array) if sub_atoms==None else sub_atoms
#TODO: Add option for sub_atoms!
    # Apply geometric transformation if specified
    # print(geom_transform_indices)
    if geom_transform_indices is not None:
        geometric_center = np.mean(coordinates_array[geom_transform_indices], axis=0)
        coordinates_array -= geometric_center

    dipole_xyz = np.vstack([(row[0] * row[1])for row in
                            list(zip(coordinates_array, charge_array))])
    dipole_vector=np.sum(dipole_xyz,axis=0)
    array_dipole=np.hstack([dipole_vector,np.linalg.norm(dipole_vector)])
    dipole_df=pd.DataFrame(array_dipole,index=XYZConstants.DIPOLE_COLUMNS.value).T
 
    return dipole_df

def calc_dipole_gaussian(coordinates_array, gauss_dipole_array, base_atoms_indices ,geometric_transformation_indices=None):
    """
    a function that recives coordinates and gaussian dipole, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    """
    if geometric_transformation_indices:
        # Calculate the geometric center of specified indices
        geometric_center = np.mean(coordinates_array[geometric_transformation_indices], axis=0)
        # Translate all coordinates
        coordinates_array -= geometric_center

    indices=adjust_indices(base_atoms_indices)
    basis_vector=calc_basis_vector(*calc_new_base_atoms(coordinates_array, indices))
    gauss_dipole_array[0,0:3]=np.matmul(basis_vector,gauss_dipole_array[0,0:3])
    dipole_df=pd.DataFrame(gauss_dipole_array,columns=['dipole_x','dipole_y','dipole_z','total'])
    # print(geometric_transformation_indices, dipole_df)
    return dipole_df

def check_imaginary_frequency(info_df):##return True if no complex frequency, called ground.state in R
        bool_imaginary=not any([isinstance(frequency, complex) for frequency in info_df['Frequency']])
        return bool_imaginary



def indices_to_coordinates_vector(coordinates_array,indices):
    """
    a function that recives coordinates_array and indices of two atoms
    and returns the bond vector between them
    """

    if  isinstance(indices[0], tuple):
        bond_vector=[(coordinates_array[index[0]]-coordinates_array[index[1]]) for index in indices]
    else:
        bond_vector= coordinates_array[indices[0]]-coordinates_array[indices[1]]

    return bond_vector

def get_bonds_vector_for_calc_angle(coordinates_array,atoms_indices): ##for calc_angle_between_atoms

    indices=adjust_indices(atoms_indices)#three atoms-angle four atoms-dihedral
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



def direction_atoms_for_sterimol(bonds_df,base_atoms)->list: #help function for sterinol
    """
    a function that return the base atom indices for coordination transformation according to the bonded atoms.
    you can insert two atom indicess-[1,2] output [1,2,8] or the second bonded atom
    if the first one repeats-[1,2,1] output [1,2,3]
    """
    
    base_atoms_copy=base_atoms[0:2]
    origin,direction=base_atoms[0],base_atoms[1]
    bonds_df = bonds_df[~((bonds_df[0] == origin) & (bonds_df[1] == direction)) & 
                              ~((bonds_df[0] == direction) & (bonds_df[1] == origin))]
    
    try :
        base_atoms[2]==origin
        if(any(bonds_df[0]==direction)):
            # take the second atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[0]==direction)][1].iloc[1])
        else:
            # take the first atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[1]==direction)][0].iloc[1])
    except: 
       
        for _, row in bonds_df.iterrows():
            if row[0] == direction:
                base_atoms_copy.append(row[1])
                break
            elif row[1] == direction:
                base_atoms_copy.append(row[0])
                break
    return base_atoms_copy

def get_molecule_connections(bonds_df,source,direction):
    graph=ig.Graph.DataFrame(edges=bonds_df,directed=True)
    paths=graph.get_all_simple_paths(v=source,mode='all')
    with_direction=[path for path in paths if (direction in path)]
    longest_path=np.unique(flatten_list(with_direction))
    return longest_path




def get_specific_bonded_atoms_df(bonds_df,longest_path,coordinates_df):
    """
    a function that returns a dataframe of the atoms that are bonded in the longest path.
    bonded_atoms_df: dataframe
        the atom type and the index of each bond.

       atom_1 atom_2 index_1 index_2
0       C      N       1       4
1       C      N       1       5
2       C      N       2       3
    """
    if longest_path is not None:
        edited_bonds_df=bonds_df[(bonds_df.isin(longest_path))].dropna().reset_index(drop=True)
    else:
        edited_bonds_df=bonds_df
    bonds_array=(np.array(edited_bonds_df)-1).astype(int) # adjust indices? 
    atom_bonds=np.vstack([(coordinates_df.iloc[bond]['atom'].values) for bond in bonds_array]).reshape(-1,2)
    bonded_atoms_df=(pd.concat([pd.DataFrame(atom_bonds),edited_bonds_df],axis=1))
    bonded_atoms_df.columns=[XYZConstants.BONDED_COLUMNS.value]
    return bonded_atoms_df


def remove_atom_bonds(bonded_atoms_df,atom_remove='H'):
    atom_bonds_array=np.array(bonded_atoms_df)
    delete_rows_left=np.where(atom_bonds_array[:,0]==atom_remove)[0] #itterrow [0] is index [1] are the values
    delete_rows_right=np.where(atom_bonds_array[:,1]==atom_remove)[0]
    atoms_to_delete=np.concatenate((delete_rows_left,delete_rows_right))
    new_bonded_atoms_df=bonded_atoms_df.drop((atoms_to_delete),axis=0)
    return new_bonded_atoms_df



def extract_connectivity(xyz_df, threshold_distance=1.82, metals=None,
                         metal_threshold=2.8, max_coordination=6):
    """
    Build a connectivity table from XYZ coordinates.

    Parameters
    ----------
    threshold_distance : float
        Max bond distance (Å) for non-metal pairs.
    metals : None | str | list[str]
        Metal element symbols to treat with relaxed distance rules.
        None → auto-detect from the full periodic-table metal set.
        Pass an empty list [] to disable metal handling entirely.
    metal_threshold : float
        Max bond distance (Å) for metal–ligand pairs.
    max_coordination : int
        Maximum number of bonds kept per metal centre.
    """
    if metals is None:
        active_metals = _METAL_ELEMENTS
    elif isinstance(metals, str):
        active_metals = {metals}
    else:
        active_metals = set(metals)

    coordinates = np.array(xyz_df[['x', 'y', 'z']].values)
    atoms_symbol = np.array(xyz_df['atom'].values)
    distances = pdist(coordinates)
    dist_matrix = squareform(distances)

    dist_df = pd.DataFrame(dist_matrix).stack().reset_index()
    dist_df.columns = ['a1', 'a2', 'value']
    dist_df['first_atom'] = [atoms_symbol[i] for i in dist_df['a1']]
    dist_df['second_atom'] = [atoms_symbol[i] for i in dist_df['a2']]

    remove_list = []
    dist_array = np.array(dist_df)
    special_atoms = {'Cl', 'Br', 'F', 'I'}

    for idx, row in enumerate(dist_array):
        i, j, dist, atom1, atom2 = row
        remove_flag = False

        if i == j:
            remove_flag = True

        if ((atom1 == 'H') and (atom2 not in XYZConstants.NOF_ATOMS.value)) or \
           ((atom1 == 'H') and (atom2 == 'H')) or \
           ((atom1 == 'H' or atom2 == 'H') and float(dist) >= 1.5):
            remove_flag = True

        involves_metal = atom1 in active_metals or atom2 in active_metals

        if not involves_metal:
            if float(dist) >= threshold_distance or float(dist) == 0:
                remove_flag = True
        else:
            if float(dist) > metal_threshold or float(dist) == 0:
                remove_flag = True

        # Halogens bonded between threshold and 2.6 Å are allowed (e.g. C–I, C–Br)
        if (atom1 in special_atoms or atom2 in special_atoms) and \
           (threshold_distance <= float(dist) < 2.6) and not involves_metal:
            remove_flag = False

        if remove_flag:
            remove_list.append(idx)

    dist_df = dist_df.drop(remove_list)
    dist_df[['min_col', 'max_col']] = pd.DataFrame(
        np.sort(dist_df[['a1', 'a2']], axis=1), index=dist_df.index
    )
    dist_df = dist_df.drop(columns=['a1', 'a2']).rename(columns={'min_col': 0, 'max_col': 1})
    dist_df = dist_df.drop_duplicates(subset=[0, 1])

    # Keep only the closest bond per halogen (terminal halogens bond to one atom)
    special_atoms_idxs = {}
    for idx, row in dist_df.iterrows():
        if row['first_atom'] in special_atoms:
            special_atoms_idxs.setdefault(row[0], []).append((idx, row['value']))
        if row['second_atom'] in special_atoms:
            special_atoms_idxs.setdefault(row[1], []).append((idx, row['value']))

    special_atoms_to_remove = []
    for atom_idx, bonds in special_atoms_idxs.items():
        if len(bonds) > 1:
            bonds.sort(key=lambda x: x[1])
            special_atoms_to_remove.extend([i for i, _ in bonds[1:]])
    dist_df = dist_df.drop(special_atoms_to_remove)

    # For each metal centre keep only the max_coordination shortest bonds
    metal_mask = dist_df['first_atom'].isin(active_metals) | dist_df['second_atom'].isin(active_metals)
    metal_bonds = dist_df[metal_mask].copy()
    non_metal = dist_df[~metal_mask].copy()

    kept_metal_idx = []
    if not metal_bonds.empty:
        metal_bonds['_metal_idx'] = metal_bonds.apply(
            lambda r: r[0] if r['first_atom'] in active_metals else r[1], axis=1
        )
        for _, group in metal_bonds.groupby('_metal_idx'):
            kept_metal_idx.extend(group.nsmallest(max_coordination, 'value').index)

    metal_kept = metal_bonds.loc[kept_metal_idx, [0, 1]] if kept_metal_idx else pd.DataFrame(columns=[0, 1])

    final = pd.concat([non_metal[[0, 1]], metal_kept], ignore_index=True)
    return pd.DataFrame(final[[0, 1]] + 1)

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
    # print(center_atom)
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

def nob_atype(xyz_df, bonds_df):
    
    symbols = xyz_df['atom'].values
    
    list_results=[]
    for index,symbol in enumerate(symbols):
        index+=1
        nob = bonds_df[(bonds_df[0] == index) | (bonds_df[1] == index)].shape[0]
        if symbol == 'H':
            result = 'H'
        elif symbol == 'F':
            result = 'F'
        elif symbol == 'P':
            result = 'P'
        elif symbol == 'Cl':
            result = 'Cl'
        elif symbol == 'Br':
            result = 'Br'
        elif symbol == 'I':
            result = 'I'
        elif symbol == 'O':
            if nob < 1.5:
                result = 'O2'
            elif nob > 1.5:
                result = 'O'
        elif symbol == 'S':
            if nob < 2.5:
                result = 'S'
            elif 2.5 < nob < 5.5:
                result = 'S4'
            elif nob > 5.5:
                result = 'S1'
        elif symbol == 'N':
            if nob < 2.5:
                result = 'C6/N6'
            elif nob > 2.5:
                result = 'N'
        elif symbol == 'C':
            if nob < 2.5:
                result = 'C3'
            elif 2.5 < nob < 3.5:
                result = 'C6/N6'
            elif nob > 3.5:
                result = 'C'
        elif symbol in _METAL_ELEMENTS:
            result = symbol
        else:
            result = 'X'

        list_results.append(result)

    return list_results

def get_sterimol_indices(coordinates,bonds_df):
    center=get_center_of_mass(coordinates)
    center_atom=get_closest_atom_to_center(coordinates,center)
    base_atoms=get_sterimol_base_atoms(center_atom,bonds_df)
    return base_atoms

def filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df):
    """
    a function that filter out NOF bonds and H bonds and returns
     a dataframe of the molecule coordinates without them.
    """
    allowed_bonds_indices= pd.concat([bonded_atoms_df['index_1'],bonded_atoms_df['index_2']],axis=1).reset_index(drop=True)
    atom_filter=adjust_indices(np.unique([atom for sublist in allowed_bonds_indices.values.tolist() for atom in sublist]))
    edited_coordinates_df=coordinates_df.loc[atom_filter].reset_index(drop=True)
   
    return edited_coordinates_df



def get_extended_df_for_sterimol(coordinates_df, bonds_df, radii='CPK'):
    """
    A function that adds information to the regular coordinates_df

    Parameters
    ----------
    coordinates_df : dataframe
    bond_type : str
        The bond type of the molecule
    radii : str, optional
        The type of radii to use ('bondi' or 'CPK'), by default 'bondi'

    Returns
    -------
    dataframe
        The extended dataframe with additional columns

    """
    
    bond_type_map_regular = GeneralConstants.REGULAR_BOND_TYPE.value
    bond_type_map=GeneralConstants.BOND_TYPE.value
    ## if radius is cpk mapping should be done on atype, else on atom
    radii_map = GeneralConstants.CPK_RADII.value if radii == 'CPK' else GeneralConstants.BONDI_RADII.value
    
    df = coordinates_df.copy()  # make a copy of the dataframe to avoid modifying the original
    
    
    
    if radii == 'bondi':
        df['atype']=df['atom']
    else:
        df['atype']=nob_atype(coordinates_df, bonds_df)

    df['magnitude'] = calc_magnitude_from_coordinates_array(df[['x', 'z']].astype(float))
    
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
    # Compute number of points per substituent (assumes evenly divided)
    # print('inside calc_B1')
    # n_total = transformed_plane.shape[0]
    # n_subs = edited_coordinates_df.shape[0]
    # n_points = n_total // n_subs if n_subs > 0 else 1

    # # Print debug info
    # print("calc_B1: transformed_plane.shape =", transformed_plane.shape)
    # print("calc_B1: len(extended_df) =", n_subs)
    # print("calc_B1: n_points per substituent =", n_points)
    
    # # Find first index where the absolute value is close to the minimum of avs
    # idx = np.where(np.isclose(np.abs(transformed_plane[:, column_index]), avs.min()))[0][0]
    # # Map the plane index back to the corresponding DataFrame row
    # index_df = idx // n_points
    # print("calc_B1: raw idx =", idx, "mapped index_df =", index_df)
    # idx=index_df
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

def generate_circle(center_x, center_y, radius, n_points=20):
    """
    Generate circle coordinates given a center and radius.
    Returns a DataFrame with columns 'x' and 'y'.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)
    return np.column_stack((x, y))


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
    
    # Generate circles for each substituent and combine them.
    circles = []
    for idx, row in extended_df.iterrows():
        circle_points = generate_circle(row['x'], row['z'], row['radius'], n_points=100)
        circles.append(circle_points)
    plane = np.vstack(circles)  # All circle points combined.
    
    sterimol_df=b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane, b1_planes)

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
    
    b1s, b1s_loc, b1_planes = [], [], []
    sterimol_df=b1s_for_loop_function(extended_df, b1s, b1s_loc, list(new_degree_list), plane, b1_planes)
  
    b1s=sterimol_df['B1']

    b1_b5_angle=sterimol_df['B1_B5_angle']
    plane=sterimol_df['plane']
    

    return [b1s, b1_b5_angle, plane]

import matplotlib.pyplot as plt

def plot_b1_visualization(rotated_plane, extended_df, n_points=100, title="Rotated Plane Visualization"):
    """
    Visualize the rotated plane by plotting:
      - Complete circles (each generated from a substituent),
      - Dashed lines at extreme x and y values,
      - Arrows for the four extreme directions (with the B1 arrow highlighted),
      - A B5 arrow (the farthest point from the origin),
      - An arc indicating the angle between the B1 and B5 arrows.
    
    Parameters
    ----------
    rotated_plane : np.array
        Rotated plane points (stacked complete circles; shape: [n_total_points, 2]).
    extended_df : pd.DataFrame
        DataFrame with columns 'radius' and 'L' (used for annotations).
    n_points : int, optional
        Number of points per circle (default is 20).
    title : str, optional
        Title for the plot.
    """
    # Compute extreme values from all points
    max_x = np.max(rotated_plane[:, 0])
    min_x = np.min(rotated_plane[:, 0])
    max_y = np.max(rotated_plane[:, 1])
    min_y = np.min(rotated_plane[:, 1])
    avs = np.abs([max_x, min_x, max_y, min_y])
    min_val = np.min(avs)
    min_index = np.argmin(avs)
    
    # Determine B1 arrow coordinates based on the minimum extreme
    if min_index == 0:
        b1_coords = np.array([max_x, 0])
    elif min_index == 1:
        b1_coords = np.array([min_x, 0])
    elif min_index == 2:
        b1_coords = np.array([0, max_y])
    else:
        b1_coords = np.array([0, min_y])
    
    # Determine B5 as the farthest point from the origin
    norms_sq = np.sum(rotated_plane**2, axis=1)
    b5_index = np.argmax(norms_sq)
    b5_point = rotated_plane[b5_index]
    b5_value = np.linalg.norm(b5_point)
  
    # Calculate angles for the arrows
    angle_b1 = np.arctan2(b1_coords[1], b1_coords[0]) % (2 * np.pi)
    angle_b5 = np.arctan2(b5_point[1], b5_point[0]) % (2 * np.pi)
    angle_diff = abs(angle_b5 - angle_b1)
    if angle_diff > np.pi:
        angle_diff = 2 * np.pi - angle_diff
    angle_diff_deg = np.degrees(angle_diff)
    
    plt.figure(figsize=(8, 8))
    
    # Plot complete circles.
    n_total = rotated_plane.shape[0]
    n_circles = n_total // n_points
    for i in range(n_circles):
        circle_points = rotated_plane[i * n_points:(i + 1) * n_points, :]
        # Close the circle by appending the first point to the end
        circle_points = np.vstack([circle_points, circle_points[0]])
        plt.plot(circle_points[:, 0], circle_points[:, 1], color='cadetblue', linewidth=1.5)
    
    # Plot dashed extreme lines
    plt.axvline(x=max_x, color='darkred', linestyle='dashed')
    plt.axvline(x=min_x, color='darkred', linestyle='dashed')
    plt.axhline(y=max_y, color='darkgreen', linestyle='dashed')
    plt.axhline(y=min_y, color='darkgreen', linestyle='dashed')
    
    # Draw arrows for each extreme (all black except the B1 arrow highlighted)
    arrow_colors = ['black'] * 4
    arrow_colors[min_index] = '#8FBC8F'
    plt.arrow(0, 0, max_x, 0, head_width=0.1, length_includes_head=True, color=arrow_colors[0])
    plt.arrow(0, 0, min_x, 0, head_width=0.1, length_includes_head=True, color=arrow_colors[1])
    plt.arrow(0, 0, 0, max_y, head_width=0.1, length_includes_head=True, color=arrow_colors[2])
    plt.arrow(0, 0, 0, min_y, head_width=0.1, length_includes_head=True, color=arrow_colors[3])
    
    # Draw the B5 arrow in red
    plt.arrow(0, 0, b5_point[0], b5_point[1], head_width=0.1, length_includes_head=True, color="#CD3333")
    
    # Annotate B1 and B5 values
    plt.text(b1_coords[0] * 0.5, b1_coords[1] * 0.5, f"B1\n{min_val:.2f}", 
             fontsize=12, ha='center', va='bottom', fontweight='bold')
    plt.text(b5_point[0] * 0.66, b5_point[1] * 0.66, f"B5\n{b5_value:.2f}", 
             fontsize=12, ha='center', va='bottom', fontweight='bold')
    
    # Draw an arc between the B1 and B5 arrows to represent the angle difference.
    arc_theta = np.linspace(min(angle_b1, angle_b5), max(angle_b1, angle_b5), 100)
    arc_x = 0.5 * np.cos(arc_theta)
    arc_y = 0.5 * np.sin(arc_theta)
    plt.plot(arc_x, arc_y, color='gray', linewidth=1.5)
    
    # Annotate the angle in degrees at the midpoint of the arc.
    mid_angle = (min(angle_b1, angle_b5) + max(angle_b1, angle_b5)) / 2
    plt.text(0.8 * np.cos(mid_angle), 0.8 * np.sin(mid_angle), f"{angle_diff_deg:.1f}°",
             fontsize=12, ha='center', va='center', fontweight='bold')
    
    plt.title(title)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.show()



# def get_b1s_list(extended_df, scans=90//5):
    
#     b1s,b1s_loc,b1_planes=[],[],[]
#     scans=scans
#     degree_list=list(range(18,108,scans))
#     plane=np.array(extended_df[['x','z']].astype(float))
#     b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane,b1_planes)

#     if b1s:
#         try:
#             back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
#             front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans
#             degree_list=range(back_ang,front_ang+1)
#         except:
            
#             back_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]-scans
#             front_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]+scans
#             degree_list=range(back_ang,front_ang+1)
#     else:
     
#         return [np.array(b1s),np.array(b1s_loc)]
    
#     b1s,b1s_loc,b1_planes=[],[] ,[]
#     b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane,b1_planes)
   
#     return [np.array(b1s),np.array(b1s_loc),b1_planes]

def calc_sterimol(bonded_atoms_df,extended_df,visualize=False):
    edited_coordinates_df=filter_atoms_for_sterimol(bonded_atoms_df,extended_df)
  
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
    # loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
    # get the idx of the row with the  biggest b5 value from edited
    max_row=edited_coordinates_df['B5'].idxmax()
    max_row=int(max_row)
    L=max(edited_coordinates_df['L'].values)
    loc_B5 = edited_coordinates_df['y'].iloc[np.where(edited_coordinates_df['B5']==max(edited_coordinates_df['B5']))[0][0]]
    
    sterimol_df = pd.DataFrame([B1, b5_value, L ,loc_B5,angle], index=['B1', 'B5', 'L','loc_B5','B1_B5_angle'])
    if visualize:
        plot_b1_visualization(best_b1_plane, edited_coordinates_df)

    return sterimol_df.T 


def get_sterimol_df(coordinates_df, bonds_df, base_atoms,connected_from_direction, radii='bondi', sub_structure=True, drop_atoms=None,visualize=False):

    if drop_atoms is not None:
        drop_atoms=adjust_indices(drop_atoms)
        for atom in drop_atoms:
            bonds_df = bonds_df[~((bonds_df[0] == atom) | (bonds_df[1] == atom))]
            ## drop the rows from coordinates_df
            coordinates_df = coordinates_df.drop(atom)

    bonds_direction = direction_atoms_for_sterimol(bonds_df, base_atoms)
    
    new_coordinates_df = preform_coordination_transformation(coordinates_df, bonds_direction)


    if sub_structure:
        if connected_from_direction is None:
            connected_from_direction = get_molecule_connections(bonds_df, base_atoms[0], base_atoms[1])
        else:
            connected_from_direction = connected_from_direction
    else:
        connected_from_direction = None
    
    
    bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction, new_coordinates_df)
    
    extended_df = get_extended_df_for_sterimol(new_coordinates_df, bonds_df, radii)
    
    ###calculations
    
    sterimol_df = calc_sterimol(bonded_atoms_df, extended_df,visualize)
    sterimol_df= sterimol_df.rename(index={0: str(base_atoms[0]) + '-' + str(base_atoms[1])})
   
    sterimol_df = sterimol_df.round(4)
    
    return sterimol_df



def calc_magnitude_from_coordinates_array(coordinates_array: npt.ArrayLike) -> List[float]:
    """
    Calculates the magnitudes of each row in the given coordinates array.

    Parameters
    ----------
    coordinates_array: np.ndarray
        A nx3 array representing the x, y, z coordinates of n atoms.

    Returns
    -------
    magnitudes: List[float]
        A list of magnitudes corresponding to the rows of the input array.

    """
    magnitude = np.linalg.norm(coordinates_array, axis=1)
    return magnitude



class Molecule():

    def __init__(self, molecule_xyz_filename):

   
        self.molecule_name = molecule_xyz_filename.split('.')[0]
        self.molecule_path = os.path.dirname(os.path.abspath(molecule_xyz_filename))
        os.chdir(self.molecule_path)
        
        
        self.xyz_df = get_df_from_file(molecule_xyz_filename)
        
        self.coordinates_array = np.array(self.xyz_df[['x', 'y', 'z']].astype(float))
        self.bonds_df = extract_connectivity(self.xyz_df)
        
        
            



    def process_sterimol_atom_group(self, atoms, radii, sub_structure=True, drop_atoms=None,visualize=False) -> pd.DataFrame:

        connected = get_molecule_connections(self.bonds_df, atoms[0], atoms[1])
        
        return get_sterimol_df(self.xyz_df, self.bonds_df, atoms, connected, radii, sub_structure=sub_structure, drop_atoms=drop_atoms, visualize=visualize)

    def get_sterimol(self, base_atoms: Union[None, Tuple[int, int]] = None, radii: str = 'bondi',sub_structure=True, drop_atoms=None,visualize=False) -> pd.DataFrame:
        """
        Returns a DataFrame with the Sterimol parameters calculated based on the specified base atoms and radii.

        Args:
            base_atoms (Union[None, Tuple[int, int]], optional): The indices of the base atoms to use for the Sterimol calculation. Defaults to None.
            radii (str, optional): The radii to use for the Sterimol calculation. Defaults to 'bondi'.

        Returns:
            pd.DataFrame: A DataFrame with the Sterimol parameters.
            
            to add
            - only_sub- sterimol of only one part.
            - drop some atoms.
        """
        if base_atoms is None:
            base_atoms = get_sterimol_indices(self.xyz_df, self.bonds_df)
        
        if isinstance(base_atoms[0], list):
            # If base_atoms is a list of lists, process each group individually and concatenate the results
            sterimol_list = [self.process_sterimol_atom_group(atoms, radii, sub_structure=sub_structure, drop_atoms=drop_atoms,visualize=visualize) for atoms in base_atoms]
            sterimol_df = pd.concat(sterimol_list, axis=0)

        else:
            # If base_atoms is a single group, just process that group
            sterimol_df = self.process_sterimol_atom_group(base_atoms, radii,sub_structure=sub_structure, drop_atoms=drop_atoms,visualize=visualize)
        return sterimol_df


    def swap_atom_pair(self, pair_indices: Tuple[int, int]) -> pd.DataFrame:
        """
        Swaps the positions of two atoms in the molecule and returns a new DataFrame with the updated coordinates.

        Args:
            pair_indices (Tuple[int, int]): The indices of the atoms to swap.

        Returns:
            pd.DataFrame: A new DataFrame with the updated coordinates.
        """
        pairs = adjust_indices(pair_indices)
        xyz_df = self.xyz_df
        temp = xyz_df.iloc[pairs[0]].copy()
        xyz_df.iloc[pairs[0]] = self.coordinates_array[pairs[1]]
        xyz_df.iloc[pairs[1]] = temp
        return xyz_df
  
 
    def get_coordination_transformation_df(self, base_atoms_indices: List[int]) -> pd.DataFrame:
        """
        Returns a new DataFrame with the coordinates transformed based on the specified base atoms.

        Args:
            base_atoms_indices (List[int]): The indices of the base atoms to use for the transformation.

        Returns:
            pd.DataFrame: A new DataFrame with the transformed coordinates.
        """
        new_coordinates_df = preform_coordination_transformation(self.xyz_df, base_atoms_indices)
        return new_coordinates_df

    ## not working after renumbering for some reason
    
    

    def get_bond_angle(self, atom_indices: List[int]) -> pd.DataFrame:
        """
        Returns a DataFrame with the bond angles calculated based on the specified atom indices.

        Args:
            atom_indices (List[int]): The indices of the atoms to use for the bond angle calculation.

        Returns:
            pd.DataFrame: A DataFrame with the bond angles.
        """
        return get_angle_df(self.coordinates_array, atom_indices)
    
    def get_bond_length_single(self, atom_pair):
        bond_length = calc_single_bond_length(self.coordinates_array, atom_pair)
        bond_length_df = pd.DataFrame([bond_length], index=[f'bond_length_{atom_pair[0]}-{atom_pair[1]}'])
        return bond_length_df

    def get_bond_length(self, atom_pairs):
        """
            Returns a DataFrame with the bond lengths calculated based on the specified atom pairs.

            Args:
                atom_pairs (Union[List[Tuple[int, int]], Tuple[int, int]]): The pairs of atoms to use for the bond length calculation.

            Returns:
                pd.DataFrame: A DataFrame with the bond lengths.
            """
        if isinstance(atom_pairs[0], list):
            # If atom_pairs is a list of lists, process each pair individually and concatenate the results
            bond_length_list = [self.get_bond_length_single(pair) for pair in atom_pairs]
            bond_df = pd.concat(bond_length_list, axis=0)
        else:
            # If atom_pairs is a single pair, just process that pair
            bond_df = self.get_bond_length_single(atom_pairs)
        return bond_df



def dict_to_horizontal_df(data_dict):
    # Initialize an empty DataFrame to store the transformed data
    df_transformed = pd.DataFrame()
    # Loop through each key-value pair in the original dictionary
    for mol, df in data_dict.items():
        transformed_data = {}
        print(df.head())    
        # Loop through each row and column in the DataFrame
        for index, row in df.iterrows():
            
            index_words = set(index.split('_'))
            for col in df.columns:
                # Create a new key using the format: col_index
                try:
                    col_words = set(col.split('_'))
                except:
                    col_words = []
                 # Check if the index and the column have the same words and remove one
                common_words = index_words.intersection(col_words)
                if col != 0 and '0':
                    if common_words:
                        unique_col_words = col_words - common_words
                        unique_index_words = index_words - common_words
                        new_key_parts = ['_'.join(common_words)] if common_words else []
                        new_key_parts.extend([part for part in ['_'.join(unique_col_words), '_'.join(unique_index_words)] if part])
                        new_key = '_'.join(new_key_parts)
                    else:
                        new_key = f"{col}_{index}"
                else:
                    new_key = f"{index}"
                # Store the corresponding value in the transformed_data dictionary
                transformed_data[new_key] = row[col]
        # Convert the dictionary into a DataFrame row with the molecule name as the index
        df_row = pd.DataFrame([transformed_data], index=[mol])
        # Append the row to df_transformed
        df_transformed = pd.concat([df_transformed, df_row], ignore_index=False)
    return df_transformed




class Molecules_xyz():
    
    def __init__(self,molecules_dir_name, renumber=False):
        self.molecules_path=os.path.abspath(molecules_dir_name)
        os.chdir(self.molecules_path) 
        self.molecules=[]
        self.failed_molecules=[]
        for file in os.listdir(): 
            if file.endswith('.xyz'):
                try:
                    self.molecules.append(Molecule(file))
                except:
                    self.failed_molecules.append(file)
                    print(f'Error: {file} could not be processed')
                    failed_file = file.rsplit('.feather', 1)[0] + '.feather_fail'
                # self.molecules.append(Molecule(log_file))
       
    


    def filter_molecules(self, indices):
        self.molecules = [self.molecules[i] for i in indices]
        self.molecules_names = [self.molecules_names[i] for i in indices]

    def get_sterimol_dict(self,atom_indices,radii='CPK'):
        sterimol_dict={}
        for molecule in self.molecules:
            try:
                sterimol_dict[molecule.molecule_name]=molecule.get_sterimol(atom_indices,radii)
                print(f'calculated sterimol for {molecule.molecule_name}')
            except ValueError as e:
                print(f'Error: {e}')
                print(f'failed to calculate for {molecule.molecule_name}')
                sterimol_dict[molecule.molecule_name]=np.nan

        return sterimol_dict
   
    def get_sterimol_df(self,atom_indices,radii='CPK'):
        sterimol_dict=self.get_sterimol_dict(atom_indices,radii)
        sterimol_df=dict_to_horizontal_df(sterimol_dict)
        return sterimol_df

if __name__=='__main__':
    pass

 