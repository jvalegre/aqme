import pandas as pd
import numpy as np
import os
import glob
import re
from enum import Enum
import tkinter as tk
from tkinter import filedialog
from typing import *
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import shutil
import fileinput
import numpy.typing as npt
import sys
import traceback
from scipy.spatial.distance import pdist, squareform

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

class FileExtensions(Enum):
    """
    Hold commonly used file extensions
    """
    SMI='.smi'
    XYZ='xyz'
    CSV='.csv'
    ZIP='.zip'
    PPT='.ppt'
    CIF='.cif'
    MOL='.mol'
    PDB='.pdb'

class XYZConstants(Enum):
    """
    Constants related to XYZ file processing
    """
    DF_COLUMNS=['atom','x','y','z']
    STERIMOL_INDEX = ['B1', 'B5', 'L','loc_B5','B1_B5_angle']
    DIPOLE_COLUMNS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole']
    RING_VIBRATION_COLUMNS = ['cross', 'cross_angle', 'para', 'para_angle']
    RING_VIBRATION_INDEX=['Product','Frequency','Sin_angle']
    VIBRATION_INDEX = ['Frequency', 'Amplitude']
    BONDED_COLUMNS = ['atom_1', 'atom_2', 'index_1', 'index_2']
    NOF_ATOMS = ['N', 'O', 'F']
    STERIC_PARAMETERS = ['B1', 'B5', 'L', 'loc_B1', 'loc_B5','RMSD']
    ELECTROSTATIC_PARAMETERS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole','energy']
    CHARGE_TYPE = ['nbo', 'hirshfeld', 'cm5']

class LinuxCommands(Enum):
    COPY='cp'
    OBABEL_PREFIX='/gpfs0/gaus/users/itamarwa/transfer_to_itamar/build/bin/obabel '
    OBABEL_XYZ_SETTINGS_1=' -O '
    OBABEL_XYZ_SETTINGS_2=' --gen3d'
    XTB_INPUT_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb --input ' # xtb_di.inp 
    XTB_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb '
    XTB_SUFIX=' --ohess --dipole --pop'
    CREST_INPUT_PREFIX='/gpfs0/gaus/projects/crest --input '
    CREST_PREFIX='/gpfs0/gaus/projects/crest '
    GAUSS_SUFIX='/gpfs0/gaus/users/kozuch/home/scripts/gaussian/g16'
    #GAUSS_SUFIX='/gpfs0/gaus/users/edenspec/g16'



def log_exception(location="<unknown>"):
    exc_type, exc_value, tb = sys.exc_info()
    # Go to innermost traceback frame
    while tb.tb_next:
        tb = tb.tb_next
    frame = tb.tb_frame
    filename = frame.f_code.co_filename
    func_name = frame.f_code.co_name
    lineno = tb.tb_lineno
    code_line = ""
    try:
        with open(filename) as f:
            lines = f.readlines()
            code_line = lines[lineno - 1].strip()
    except Exception:
        code_line = "<Could not retrieve code line>"

    print("\n" + "="*60)
    print(f"🔥 Exception in {location}")
    print(f"  File      : {filename}")
    print(f"  Function  : {func_name}")
    print(f"  Line      : {lineno}")
    print(f"  Code      : {code_line}")
    print(f"  Type      : {exc_type.__name__}")
    print(f"  Message   : {exc_value}")
    print("-"*60)


    
def compare_cosine_distance_matrices(matrix1, matrix2):
    def cosine_distance(matrix):
        # Calculate norms
        norms = np.linalg.norm(matrix, axis=1, keepdims=True)
        # Avoid division by zero
        norms[norms == 0] = np.finfo(float).eps
        # Normalize the rows to unit length
        norm_matrix = matrix / norms
        # Use dot product to find cosine similarity and subtract from 1 to get cosine distance
        similarity = np.dot(norm_matrix, norm_matrix.T)
        return 1 - similarity

    # Calculate the cosine distance matrices
    matrix1_distances = cosine_distance(matrix1)
    matrix2_distances = cosine_distance(matrix2)
    
    # Calculate differences between the two distance matrices
    differences = matrix1_distances - matrix2_distances
    average_difference = np.mean(np.abs(differences))

    return average_difference


def get_file_name_list(file_identifier):
    """
    The function gets a file identifier as input and returns a list of all files in the working 
    which contain the identifier in the files name
    ----------
    Parameters
    ----------
    identifier : str.
        The wanted file identifier like 'txt','info','nbo' contained in the filename
    -------
    Returns
    -------
    list
        A list of all files in the working directory with the chosen extension 
    --------
    Examples
    --------
    
    all_files_in_dir=listdir()
    print(all_files_in_dir)
        ['0_1106253-mod-mod.xyz', '0_1106253-mod.xyz', '1106253.cif', '1109098.cif', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz', 'cif_handler.py']
        
    xyz_files_in_dir=get_filename_list('.xyz')
    print(xyz_files_in_dir)
        ['0_1106253-mod-mod.xyz', '0_1106253-mod.xyz', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz']
  
    """
    return [filename for filename in os.listdir() if file_identifier in filename]

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

def adjust_indices(element):
    if isinstance(element, (list, tuple, np.ndarray)):
        return [adjust_indices(sub_element) for sub_element in element]
    elif isinstance(element, (int, np.integer)):
        return int(element) - 1
    elif isinstance(element, float) and element.is_integer():
        return int(element) - 1
    else:
        raise ValueError(f"Unsupported element type: {type(element)} — value: {element}")



def data_to_xyz(dataframe, output_name, comment_line=''):
    """

     a function that recieves a dataframe, output name, and comment line and creates a xyz type file.
     
    parameters

    ---

    dataframe: an array that can contain different classes of data, needs to be 4 colums to run.

    output_name:str, the name for the file created.

    comment_line: str, the headline of the file .
    ---

    examples:
    ---
    """
    if type(dataframe) == pd.DataFrame :
        number_of_atoms=dataframe.shape[0]
        atoms_np_array=dataframe.to_numpy()
    else:
        number_of_atoms=len(dataframe)
        atoms_np_array=dataframe
    with open(output_name, 'w') as xyz_file:
        xyz_file.write("{}\n{}\n".format(number_of_atoms, comment_line))
        for atom_np_array in atoms_np_array:
            try:
                xyz_file.write("{:1} {:11.6} {:11.6} {:11.6} \n".format(*atom_np_array))
            except:
                xyz_file.write("{:1}".format(*atom_np_array))


def change_filename(old_name, new_name):
    # Get the file extension
    file_extension = os.path.splitext(old_name)[1]
    # Combine the new name with the original file extension
    new_filename = new_name + file_extension
    os.rename(old_name, new_filename)

def change_filetype (filename,new_type='xyz'):
    """
    a function that recieves a file name, and a new type, and changes the type-ending of the file's name to the new one.

    parameters
    ---
    filename: str, the file we want to change

    new_type:str, the new ending we want for the file

    returns
    ---
    the same file name with a new type ending

    examples
    ---
    filename='B_THR_127_5Angs_noHOH.pdb'
    new_filename=change_filetype(filename,'xyz')
    OUTPUT:'B_THR_127_5Angs_noHOH.xyz'
    
    """
    split_result=filename.split('.')
    if '.' in new_type:
        new_filename=split_result[0]+new_type
    else:
        new_filename=split_result[0]+'.'+new_type
    return new_filename

def xyz_string_to_df(lines):
    strip_lines=[line.strip().rstrip('\n') for line in lines]
    
def create_molecule_directories():
    list_of_dirs=[name.split('.')[0] for name in os.listdir()]
    for dir_name in list_of_dirs:
        os.mkdir(dir_name)
    return

def delete_file(filename):
    """
    a function that gets a file name and deletes it.
    """
    os.remove(os.path.abspath(filename))
    return

def delete_type_files(file_type='xyz'): ## my help function to delete xyz files
    """
    a function that gets a directory path and file type, and deletes all the files of said type.
    """
    list_of_molecules=[file for file in os.listdir() if file.endswith(file_type)]
    for molecule in list_of_molecules:
        os.remove(os.path.abspath(molecule))
        
        
def move_files_directory(file_type):#need edit
    """
    a function that moves xyz type files from one directory to another.
    help function for xyz_file_generator_library to move files to the new directory created
    A void function
    """
    list_of_dirs=[name for name in os.listdir() if os.path.isdir(os.path.abspath(name))]
    list_of_files=get_file_name_list(file_type)
    print(list_of_files,list_of_dirs)
    for file_name,dir_name in zip(list_of_files,list_of_dirs):
        new_path=os.path.join(os.path.abspath(dir_name),file_name)
        os.replace(os.path.abspath(file_name),new_path)
    return


def edit_filenames_in_directory(directory_path,old: str =None, new: str =None):
    # Check if the specified directory path exists and is a directory
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        # Loop through all files in the directory
        for filename in os.listdir(directory_path):
            file_path = os.path.join(directory_path, filename)
            
            # Check if it's a file (not a subdirectory)
            if os.path.isfile(file_path):
                # Replace spaces, parentheses, brackets, and commas with underscores
                new_filename = filename.replace(' ', '_').replace('(', '').replace(')', '').replace('[', '').replace(']', '').replace(',', '')
                new_filename = filename.replace(old, new)
                new_file_path = os.path.join(directory_path, new_filename)

                # Check if the filename needs to be changed
                if file_path != new_file_path:
                    # Rename the file in place (replace the original)
                    os.rename(file_path, new_file_path)
                    print(f"Renamed '{filename}' to '{new_filename}'")

    else:
        print("The specified directory path does not exist or is not a directory.")


def sort_files_to_directories(file_type='xyz'):
    create_molecule_directories()
    move_files_directory(file_type)
    return


## write help functions for mol2 to dfs
def string_to_int_list(string,splitter=' '):
    string=string.split(splitter)
    return [int(value) for value in string]

def split_to_pairs(string,splitter=' '):
    list_a=string_to_int_list(string)
    chunck=(len(list_a)/2)
    splited=np.array_split(list_a,chunck)
    return [value.tolist() for value in splited]

def split_for_angles(string):
    string_list=string.split('  ')
    return [string_to_int_list(value,' ') for value in string_list]

def choose_filename():
    # Create a root window
    root = tk.Tk()
    # Ask the user to select a file
    file_path = filedialog.askopenfilename()
    root.withdraw()
    # Extract the file name from the file path
    file_name = file_path.split("/")[-1]
    # Print the file name
    return file_name, file_path.split(file_name)[0]

def choose_directory():
    # Create a root window
    root = tk.Tk()
    # Ask the user to select a directory
    directory_path = filedialog.askdirectory()
    root.withdraw()
    # Extract the directory name from the directory path
    directory_name = directory_path.split("/")[-1]
    # Print the directory name
    return directory_name, directory_path

def flatten_list(nested_list_arg: List[list]) -> List:
    """
    Flatten a nested list.
    turn [[1,2],[3,4]] to [1,2,3,4]
    """
    flat_list=[item for sublist in nested_list_arg for item in sublist]
    return flat_list



def submit_to_gaussian_calculation(queue='milo'):
    # Copy file from source to destination
    os.system(LinuxCommands.COPY.value +' '+ LinuxCommands.GAUSS_SUFIX.value+ ' .')
    # Loop over all files in current directory with .com extension
    for file in os.listdir():
        if file.endswith('.com'):
            # Run the command on each file
            output_name=file.split('.')[0]+'.log'
            os.system(f'./g16 {file} {output_name} -q {queue}')
    # Remove the g16 file
    os.remove('g16')
    # Change directory back to the previous directory
    return
                                  
def calculate_distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    """
    Calculate a distance matrix given an array of coordinates.
    coordinates: np.array of x y z coordinates

    """
    num_atoms = len(coordinates)
    distance_matrix = np.zeros((num_atoms, num_atoms))

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            atom_i = coordinates[i]
            atom_j = coordinates[j]

            distance = np.linalg.norm(atom_i - atom_j)
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    return distance_matrix


def validate_geometry(distance_matrix, threshold=0.5):
    num_atoms = distance_matrix.shape[0]
    invalid_atoms = []

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = distance_matrix[i, j]

            if distance < threshold:
                invalid_atoms.append((i, j))

    return invalid_atoms

def fix_coordinates(coordinates, invalid_atoms, displacement=0.3):
    fixed_coordinates = coordinates.copy()
    for atom_i, atom_j in invalid_atoms:
        # Compute the vector between atom_i and atom_j
        vector = fixed_coordinates[atom_j] - fixed_coordinates[atom_i]
        # Normalize the vector
        normalized_vector = vector / np.linalg.norm(vector)
        # Adjust the coordinates of atom_j beyond the threshold
        fixed_coordinates[atom_j] = fixed_coordinates[atom_i] + displacement * normalized_vector

    return fixed_coordinates

def generate_3d_coordinates(smiles):
    # Generate a molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    # Generate 3D coordinates using the ETKDG method
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    # Return the molecule with 3D coordinates and hydrogen atoms
    conformer = mol.GetConformer()
    symbols = ([atom.GetSymbol() for atom in mol.GetAtoms()])
    list_of_lists = np.array([[x] for x in symbols])
    array= np.concatenate((list_of_lists, conformer.GetPositions()), axis=1)
    return array


def smiles_to_xyz_files(smiles_list: List[str], molecule_names: List[str], new_dir: bool = False):
    """
    Convert a list of SMILES strings to XYZ files.

    Args:
        smiles_list (List[str]): List of SMILES strings.
        molecule_names (List[str]): List of corresponding molecule names.
        new_dir (bool, optional): Create a new directory for XYZ files. Defaults to False.
    """
    if new_dir:
        os.mkdir('xyz_files')
        os.chdir('xyz_files')
    for smiles, name in zip(smiles_list, molecule_names):
        coordinates = generate_3d_coordinates(smiles)
        data_to_xyz(coordinates, name + '.xyz')


def replace_path_in_sh_file(filename: str, new_path: str, sh_file_path: str = r'/gpfs0/gaus/users/edenspec/Conformers_code-main/M1_outward_sender'):
    """
    Replace a specific line containing 'path=' in a shell script file with a new path.

    Args:
        filename (str): The name of the shell script file to modify.
        new_path (str): The new path to replace with.
        sh_file_path (str, optional): The directory containing the shell script file. Defaults to the specified path.
    """
    # Change the current working directory to the shell script directory
    os.chdir(sh_file_path)

    # Open the file in place for editing
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            # Check if the line contains 'path='
            if '    path=' in line:
                # Replace the entire line with the new path
                line = f"    path='{new_path}'\n"
            # Print the modified line or the original line if no replacement was done
            print(line, end='')



def find_files_with_extension(
    directory_path: str,
    file_extension: str,
    min_size_bytes: Optional[int] = None,
    max_size_bytes: Optional[int] = None,
    return_non_matching: bool = False
) -> List[str]:
    """
    Find files with a specific extension and optional size range in a directory.

    Args:
        directory_path (str): The path to the directory to search.
        file_extension (str): The file extension to filter (e.g., '.log').
        min_size_bytes (int, optional): The minimum file size in bytes. Defaults to None.
        max_size_bytes (int, optional): The maximum file size in bytes. Defaults to None.
        return_non_matching (bool, optional): Whether to return files that do not meet the criteria. Defaults to False.

    Returns:
        List[str]: A list of file paths that meet the criteria if return_non_matching is False,
                   otherwise, a list of file paths that do not meet the criteria.
    """
    matching_files = []
    non_matching_files = []
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        file_size = os.path.getsize(file_path)

        if filename.endswith(file_extension):
            if (min_size_bytes is None or file_size >= min_size_bytes) and (max_size_bytes is None or file_size <= max_size_bytes):
                matching_files.append(file_path)
            elif return_non_matching:
                non_matching_files.append(file_path)

    return non_matching_files if return_non_matching else matching_files


def move_files_to_directory(file_list: List[str], destination_directory: str, create_directory: bool = False):
    """
    Move a list of files to a destination directory.

    Args:
        file_list (List[str]): List of file paths to move.
        destination_directory (str): The destination directory.
        create_directory (bool, optional): Create the destination directory if it doesn't exist. Defaults to False.
    """
    if create_directory and not os.path.exists(destination_directory):
        os.makedirs(destination_directory)
    for file_path in file_list:
        filename = os.path.basename(file_path)
        destination_path = os.path.join(destination_directory, filename)
        shutil.move(file_path, destination_path)


def delete_files_not_in_list(file_list: List[str]):
    """
    Delete all files that are not in the list.

    Args:
        file_list (List[str]): List of file paths to keep.
    """
    for filename in os.listdir(os.getcwd()):
        file_path = os.path.join(os.getcwd(), filename)
        if file_path not in file_list:
            #and if its not a directory
            if not os.path.isdir(file_path):
                os.remove(file_path)



def move_files_to_directories():
    """
    Move files to directories based on their names.

    This function looks at the current directory and finds files with names in the format "directory_name.log".
    It then moves each file to a subdirectory with the same name as the directory part of the file name.

    If a file with the same name already exists in the destination directory, the function appends a numeric
    suffix to the file name to avoid overwriting.

    Returns:
        None
    """
    # List of directories to create or use existing ones
    list_of_dirs = [name.split('.')[0] for name in os.listdir()]

    for dir_name in list_of_dirs:
        file_name = dir_name + '.log'
        destination_dir = dir_name

        # Check if the source file exists before moving
        if os.path.exists(file_name):
            # Check if the destination file already exists
            if os.path.exists(os.path.join(destination_dir, file_name)):
                # If it exists, rename the file being moved with a numeric suffix
                i = 1
                while True:
                    new_file_name = f"{dir_name}_{i}.log"
                    if not os.path.exists(os.path.join(destination_dir, new_file_name)):
                        file_name = new_file_name
                        break
                    i += 1

            # Move the file to the destination directory with the final file name
            shutil.move(file_name, os.path.join(destination_dir, file_name))

    
# def convert_to_list_or_nested_list(input_str):
#     split_by_space = input_str.split(' ')
    
#     # If there are no spaces, return a flat list
#     if len(split_by_space) == 1:
#         return list(map(int, filter(lambda x: x.strip(), split_by_space[0].split(','))))
    
#     # Otherwise, return a nested list
#     nested_list = []
#     for sublist_str in split_by_space:
#         sublist = list(map(int, sublist_str.split(',')))
#         nested_list.append(sublist)
#     return nested_list


def convert_to_custom_nested_list(input_str):
    """
    Converts a comma-separated string into a nested list based on the format:
    - "1,2,3,4,5,6" -> [[1,2,3,4], 5, 6]
    - "1,2,3,4,5,6 2,3,1,5,6,7" -> [[[1,2,3,4], 5, 6], [[2,3,1,5], 6, 7]]
    """
    print(f"Input string: {input_str}")
    split_by_space = input_str.split(' ')  # Split by space for multiple sections

    def process_sublist(sublist_str):
        """Convert a single comma-separated list to the required nested structure."""
        try:
            elements = list(map(int, sublist_str.split(',')))  # Convert to list of integers
        except:
            elements = list(map(int, sublist_str))
        if len(elements) > 2:  # Ensure we separate all but the last two elements
            return [elements[:-2]] + elements[-2:]
        return elements  # If fewer than 3 elements, return as-is

    # Process each segment and decide if it's a nested list or a single flat list
    if len(split_by_space) == 1:
        # Single segment, no spaces
        return process_sublist(split_by_space[0])
    else:
        # Multiple segments separated by spaces
        nested_list = []
        for sublist_str in split_by_space:
            nested_list.append(process_sublist(sublist_str))
        return nested_list
    

def convert_to_list_or_nested_list(input_str):
    # Remove trailing spaces
    input_str = input_str.rstrip()

    # Use regular expression to split by space, dash, or underscore
    split_by_delimiter = re.split(' |-|_', input_str)
    
    # Filter out empty strings
    split_by_delimiter = list(filter(None, split_by_delimiter))

    # If there's only one element, return a flat list
    if len(split_by_delimiter) == 1:
        return list(map(int, filter(None, re.split(',', split_by_delimiter[0]))))
    
    # Otherwise, return a nested list
    nested_list = []
    for sublist_str in split_by_delimiter:
        sublist = list(map(int, filter(None, re.split(',', sublist_str))))
        nested_list.append(sublist)
    return nested_list

def charge_dict_to_horizontal_df(charge_data: dict) -> pd.DataFrame:
    """
    Transforms a nested dictionary of charge data into a horizontal DataFrame.
    
    The input dictionary should be structured as:
    
        {
          molecule_name: {
              charge_type: DataFrame with a single row and columns representing atoms,
              ...
          },
          ...
        }
    
    For example:
    
        {
          'Br-3': {
              'nbo': pd.DataFrame({'atom_1': [-0.20616], 'atom_2': [-0.19537], 'atom_3': [-0.16316]}, index=['charge']),
              'hirshfeld': pd.DataFrame({'atom_1': [-0.030327], 'atom_2': [-0.026769], 'atom_3': [-0.020907]}, index=['11']),
              'cm5': pd.DataFrame({'atom_1': [-0.084356], 'atom_2': [-0.08244], 'atom_3': [-0.073024]}, index=['12'])
          },
          'Br-35': {
              'nbo': pd.DataFrame({'atom_1': [-0.23143], 'atom_2': [-0.0846], 'atom_3': [-0.19]}, index=['charge']),
              'hirshfeld': pd.DataFrame({'atom_1': [-0.0363], 'atom_2': [0.009787], 'atom_3': [-0.027177]}, index=['11']),
              'cm5': pd.DataFrame({'atom_1': [-0.088308], 'atom_2': [-0.011496], 'atom_3': [-0.077175]}, index=['12'])
          }
        }
    
    The output DataFrame will have molecule names as the index and columns such as:
    
        nbo_atom_1, nbo_atom_2, nbo_atom_3, hirshfeld_atom_1, hirshfeld_atom_2, hirshfeld_atom_3,
        cm5_atom_1, cm5_atom_2, cm5_atom_3
    
    Returns:
        pd.DataFrame: The transformed horizontal DataFrame.
    """
    rows = {}
    # Loop through each molecule and its charge type dictionaries.
    for mol, charge_types in charge_data.items():
        row_data = {}
        # For each charge type (e.g. 'nbo', 'hirshfeld', 'cm5')
        for charge_type, df in charge_types.items():
            # Assume each DataFrame contains a single row.
            for col in df.columns:
                new_key = f"{charge_type}_{col}"
                # Extract the single value (the first and only row)
                value = df.iloc[0][col]
                row_data[new_key] = value
        rows[mol] = row_data
    
    # Create a DataFrame from the dictionary of rows.
    horizontal_df = pd.DataFrame.from_dict(rows, orient='index')
    return horizontal_df

def dict_to_horizontal_df(data_dict):
    """
    Converts a dictionary of DataFrames or Series into a single horizontally concatenated DataFrame.
    Each molecule (key) becomes one row; columns are dynamically generated from index and column names.

    Handles inconsistent input types, numeric indices, and overlapping label names robustly.
    """

    import pandas as pd

    # --- Safety checks ---
    if data_dict is None or not isinstance(data_dict, dict):
        print("[WARN] dict_to_horizontal_df: Expected dict, got None or wrong type.")
        return pd.DataFrame()

    df_transformed = pd.DataFrame()

    for mol, df in data_dict.items():
        if df is None:
            print(f"[WARN] {mol}: received None, skipping.")
            continue

        # Convert Series to DataFrame
        if isinstance(df, pd.Series):
            df = df.to_frame().T

        if not isinstance(df, pd.DataFrame) or df.empty:
            print(f"[WARN] {mol}: empty or invalid DataFrame, skipping.")
            continue

        transformed_data = {}

        # Iterate rows safely
        for index, row in df.iterrows():
            index_str = str(index).strip()  # <- prevents int.split() errors
            try:
                index_words = set(index_str.split('_'))
            except Exception:
                index_words = set()

            for col in df.columns:
                col_str = str(col).strip()
                try:
                    col_words = set(col_str.split('_'))
                except Exception:
                    col_words = set()

                # Identify overlap
                common_words = index_words.intersection(col_words)
                if common_words:
                    # Build minimal unique key
                    unique_col = col_words - common_words
                    unique_idx = index_words - common_words
                    new_key_parts = list(common_words) + list(unique_col) + list(unique_idx)
                    new_key = "_".join(filter(None, new_key_parts))
                else:
                    new_key = f"{col_str}_{index_str}"

                # Store value
                try:
                    transformed_data[new_key] = row[col]
                except Exception as e:
                    print(f"[WARN] {mol}: failed on {new_key} -> {e}")

        # Build DataFrame row for this molecule
        if transformed_data:
            df_row = pd.DataFrame([transformed_data], index=[mol])
            df_transformed = pd.concat([df_transformed, df_row], ignore_index=False)
        else:
            print(f"[INFO] {mol}: no valid transformed data.")

    # Clean up duplicate columns (optional)
    df_transformed = df_transformed.loc[:, ~df_transformed.columns.duplicated()].copy()

    return df_transformed




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

def calc_angle(p1: npt.ArrayLike, p2: npt.ArrayLike, degrees: bool=False) -> float: ###works, name in R: 'angle' , radians
    dot_product=np.dot(p1, p2)
    norm_p1=np.linalg.norm(p1)
    norm_p2=np.linalg.norm(p2)
    thetha=np.arccos(dot_product/(norm_p1*norm_p2))
    if degrees:
        thetha=np.degrees(thetha)   
    return thetha


def check_imaginary_frequency(info_df):##return True if no complex frequency, called ground.state in R
        bool_imaginary=not any([isinstance(frequency, complex) for frequency in info_df['Frequency']])
        return bool_imaginary




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


import plotly.graph_objs as go
from ipywidgets import interact, FloatSlider
import plotly.io as pio

def interactive_corr_heatmap_with_highlights(
    df,
    initial_threshold=0.9,
    highlight_thresh=0.95,
    cell_size=25,
    min_size=300,
    max_size=800,
    *,
    renderer="colab",   # e.g. "colab", "notebook_connected", "vscode", "browser"
    show=True
):
    """
    Plotly correlation heatmap with a built-in slider to control min |r|.
    Works in Colab, Jupyter, and non-notebook contexts.
    """

    # ---- Choose a renderer that works in current environment (Colab-aware) ----
    def _pick_renderer(r):
        if r:  # explicit choice wins
            return r
        if "google.colab" in sys.modules:
            return "colab"
        # Jupyter kernel?
        try:
            from IPython import get_ipython
            ip = get_ipython()
            if ip and 'IPKernelApp' in getattr(ip, 'config', {}):
                return "notebook_connected"
        except Exception:
            pass
        return "browser"

    try:
        pio.renderers.default = _pick_renderer(renderer)
    except Exception:
        pio.renderers.default = "browser"

    # ---- Use only numeric columns ----
    num = df.select_dtypes(include=[np.number])
    if num.shape[1] == 0:
        print("No numeric columns to correlate.")
        return None

    # ---- Build correlation + ordering ----
    corr = num.corr()
    order = corr.abs().sum().sort_values(ascending=False).index
    corr_sorted = corr.loc[order, order]
    n = len(order)
    size = int(min(max(n * cell_size, min_size), max_size))

    # ---- Label features with '*' if they appear in any |r| > highlight_thresh pair ----
    highly_corr_features = set()
    abs_corr = corr_sorted.abs().values
    for i in range(n):
        for j in range(n):
            if i != j and abs_corr[i, j] > highlight_thresh:
                highly_corr_features.add(order[i])
                highly_corr_features.add(order[j])
    labels = [f"{name}{'*' if name in highly_corr_features else ''}" for name in order]

    # ---- Helper to mask correlations below threshold ----
    def masked_z(thr: float):
        m = corr_sorted.where(corr_sorted.abs() >= thr)
        # Keep diagonal at 1.0 even if thr > 1 (rare), avoid showing 0's
        m.values[[np.arange(n)], [np.arange(n)]] = np.diag(corr_sorted.values)
        m = m.where(m != 0)
        return m.values

    # ---- Build figure with slider ----
    thr_values = np.round(np.linspace(0.0, 1.0, 51), 2)  # 0.00 .. 1.00
    init_idx = int(np.argmin(np.abs(thr_values - initial_threshold)))

    fig = go.Figure(
        data=go.Heatmap(
            z=masked_z(initial_threshold),
            x=labels,
            y=labels,
            colorscale="RdBu",
            zmid=0,
            colorbar=dict(title="Corr"),
            hovertemplate=(
                "Feature 1: %{y}<br>"
                "Feature 2: %{x}<br>"
                "r = %{z:.3f}<extra></extra>"
            ),
        )
    )

    steps = []
    for t in thr_values:
        steps.append(dict(
            method="restyle",
            args=[{"z": [masked_z(float(t))]}],
            label=f"{t:.2f}"
        ))

    fig.update_layout(
        title=f"Correlation Heatmap (|r| ≥ {initial_threshold:.2f})",
        xaxis_nticks=n,
        yaxis_nticks=n,
        autosize=False,
        width=size,
        height=size,
        margin=dict(l=30, r=30, t=50, b=30),
        sliders=[dict(
            active=init_idx,
            steps=steps,
            currentvalue=dict(prefix="Min |r| = ", visible=True),
            x=0.5, xanchor="center",
            y=-0.5, yanchor="top",   # << place below the plot + labels
            len=0.8,
            pad=dict(t=0, b=0, l=0, r=0)
        )]
    )

    if show:
        fig.show()
        # plt.show()
        
    return fig

import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

def pick_samples_to_remove_for_distribution(y, n_remove, metric='ks', bins=20, plot=True):
    """
    Picks indices of n_remove samples to remove from y such that the distribution
    of the remaining y is as close as possible to the original.
    Plots before/after histograms, marks removed bins, and uses count (not density).

    Args:
        y (array-like): Target values.
        n_remove (int): Number of samples to remove.
        metric (str): Metric to use ('ks', 'mean', or 'std').
        bins (int or array): Number of bins or explicit bin edges for the histogram plot.
        plot (bool): If True, show the before/after plot.

    Returns:
        List[int]: indices_to_remove in y.
    """
    y = np.asarray(y).ravel()
    if n_remove < 0:
        raise ValueError("n_remove must be >= 0")
    if n_remove > len(y) - 1:
        raise ValueError("n_remove cannot be >= len(y)")

    indices_to_remove = []
    all_indices = set(range(len(y)))
    y_orig = y.copy()

    # --- Greedy selection: at each step remove the point that minimizes the chosen metric ---
    for _ in range(n_remove):
        best_score = None
        best_idx = None
        remaining = list(all_indices - set(indices_to_remove))
        # short-circuit: nothing left to test
        if not remaining:
            break

        for idx in remaining:
            y_test = np.delete(y, indices_to_remove + [idx])
            if metric == 'ks':
                score = ks_2samp(y_orig, y_test).statistic
            elif metric == 'mean':
                score = abs(np.mean(y_orig) - np.mean(y_test))
            elif metric == 'std':
                score = abs(np.std(y_orig, ddof=0) - np.std(y_test, ddof=0))
            else:
                raise ValueError("metric must be one of {'ks','mean','std'}")

            if best_score is None or score < best_score:
                best_score = score
                best_idx = idx

        indices_to_remove.append(best_idx)

    y_remaining = np.delete(y, indices_to_remove)

    if plot:
        # --- Histograms using shared bin edges (so bars line up exactly) ---
        counts, bin_edges = np.histogram(y, bins=bins)
        counts_remain, _ = np.histogram(y_remaining, bins=bin_edges)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        widths = np.diff(bin_edges)

        plt.figure(figsize=(8, 4))
        plt.bar(bin_centers, counts,        width=widths, alpha=0.5, color='tab:blue', edgecolor='k', label='Original')
        plt.bar(bin_centers, counts_remain, width=widths, alpha=0.5, color='tab:red',  edgecolor='k', label='After removal')

        # --- Mark bins from which samples were removed (robust to edge cases) ---
        removed_y = y[indices_to_remove]
        # right=True ensures max values fall into the last bin
        removed_bins = np.digitize(removed_y, bin_edges, right=True) - 1
        # clip to valid range [0, n_bins-1]
        removed_bins = np.clip(removed_bins, 0, len(bin_centers) - 1)

        unique_bins, removed_per_bin = np.unique(removed_bins, return_counts=True)
        # vertical offset for annotations (5% of max height, at least 1)
        y_offset = max(1, int(0.05 * (counts.max() if counts.size else 1)))

        for b, count in zip(unique_bins, removed_per_bin):
            x_ = float(bin_centers[b])
            y_ = int(counts[b]) if counts.size else 0
            plt.annotate(
                f'{count} removed',
                xy=(x_, y_),
                xytext=(x_, y_ + y_offset),
                ha='center', fontsize=10, color='black',
                arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.2, shrinkA=0.5)
            )

        plt.xlabel('y value')
        plt.ylabel('Count')
        plt.title(f'Distribution Before (blue) and After (red) removal of {len(indices_to_remove)} sample(s)')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return indices_to_remove



def add_output_column_csv(csv_filepath, output_column_name, output_values):
    """
    Adds an output column to a CSV file and saves it.

    Parameters:
    - csv_filepath (str): Path to the CSV file.
    - output_column_name (str): Name of the new output column.
    - output_values (list): List of values to add to the new column.
    """
    # Read the existing CSV file
    df = pd.read_csv(csv_filepath, index_col=0)
    
    # Add the new output column
    df[output_column_name] = output_values
    
    # Save the updated DataFrame back to the CSV file
    df.to_csv(csv_filepath)



from typing import Dict, Iterable, Optional, Union, Literal

import numpy as np
import pandas as pd
import json

def save_input_json(entry_widgets, json_path):
    data = {}
    for question, entry in entry_widgets.items():
        # Allow both Tkinter widgets and normal values
        if hasattr(entry, "get"):  # Tkinter Entry, Combobox, etc.
            value = entry.get()
        else:
            value = entry  # plain str, int, float, list, etc.

        # Try to convert to int or float if possible
        try:
            if isinstance(value, str) and value.strip() != "":
                if "." in value:
                    value = float(value)
                else:
                    value = int(value)
        except ValueError:
            pass

        data[question] = value

    with open(json_path, "w") as f:
        json.dump(data, f, indent=4)

    print(f"Input saved to {json_path}")

def _validate_and_build_maps(n: int, renumbering_dict: Dict[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """
    Validate 1-based renumbering and build:
    - old_to_new (0-based): index i (old) -> new index
    - new_to_old (0-based): index j (new) -> old index (a permutation/order vector)
    """
    if not renumbering_dict:
        # Identity
        old_to_new = np.arange(n, dtype=int)
        new_to_old = np.arange(n, dtype=int)
        return old_to_new, new_to_old

    # Basic checks
    for k, v in renumbering_dict.items():
        if not (isinstance(k, int) and isinstance(v, int)):
            raise ValueError("renumbering_dict must map int→int (1-based).")
        if not (1 <= k <= n and 1 <= v <= n):
            raise ValueError(f"Mapping out of range: {k}->{v} for n={n}.")

    # Ensure no duplicate targets
    new_targets = list(renumbering_dict.values())
    if len(set(new_targets)) != len(new_targets):
        raise ValueError("renumbering_dict has duplicate target indices (non-bijective).")

    # Build maps (0-based)
    old_to_new = np.arange(n, dtype=int)          # default identity
    for old1, new1 in renumbering_dict.items():
        old_to_new[old1 - 1] = new1 - 1

    # Invert map to get the permutation "order" (new_to_old)
    # For positions not explicitly moved, identity holds.
    new_to_old = np.empty(n, dtype=int)
    new_to_old.fill(-1)
    for old0, new0 in enumerate(old_to_new):
        if new_to_old[new0] != -1:
            # Should never happen after duplicate-target check
            raise ValueError(f"Non-bijective mapping detected for new index {new0+1}.")
        new_to_old[new0] = old0

    # Safety: all filled?
    if (new_to_old < 0).any():
        raise ValueError("Incomplete inversion of mapping. Check renumbering_dict.")

    return old_to_new, new_to_old


def _reindex_like(obj: Union[pd.Series, pd.DataFrame], order: np.ndarray) -> Union[pd.Series, pd.DataFrame]:
    """Reindex a Series/DataFrame (row-wise) by an integer order vector (0-based)."""
    if hasattr(obj, "iloc"):
        out = obj.iloc[order].copy()
        out.index = range(len(out))  # normalize to 0..n-1
        return out
    raise TypeError(f"Unsupported type for reindexing: {type(obj)}")


def _remap_bonds_df(bonds_df: pd.DataFrame, old_to_new: np.ndarray, assume_0_based: bool = True) -> pd.DataFrame:
    """
    Remap bonds (two integer columns) using old_to_new (0-based map).
    If your bonds are 1-based, set assume_0_based=False.
    """
    out = bonds_df.copy()
    cols = out.columns[:2]
    if not assume_0_based:
        # convert 1-based to 0-based, map, then back to 1-based
        out[cols] = out[cols].astype(int) - 1
        out[cols] = out[cols].applymap(lambda idx: int(old_to_new[idx]))
        out[cols] = out[cols] + 1
    else:
        out[cols] = out[cols].astype(int)
        out[cols] = out[cols].applymap(lambda idx: int(old_to_new[idx]))
    return out




MappingFormat = Literal[
    "old_to_new_1b",  # dict/array maps old (1-based) -> new (1-based)
    "old_to_new_0b",  # dict/array maps old (0-based) -> new (0-based)
    "new_to_old_1b",  # dict/array maps new (1-based) -> old (1-based)
    "new_to_old_0b",  # dict/array maps new (0-based) -> old (0-based)
]

def check_renumbering(
    df_before: pd.DataFrame,
    df_after: pd.DataFrame,
    mapping: Union[Dict[int, int], Iterable[int], np.ndarray],
    *,
    mapping_format: MappingFormat = "old_to_new_1b",
    cols: Optional[Iterable[str]] = None,
    atol: float = 1e-8,
    rtol: float = 1e-5,
    return_details: bool = True,
) -> Union[bool, dict]:
    """
    Check that df_after is df_before permuted according to `mapping`.

    Parameters
    ----------
    df_before, df_after : pd.DataFrame
        Original and renumbered DataFrames (same number of rows).
    mapping : dict or 1D array-like
        A permutation mapping in one of the supported formats (see `mapping_format`).
        Dict keys/values must be ints. Arrays must be length n and contain a permutation.
    mapping_format : {'old_to_new_1b','old_to_new_0b','new_to_old_1b','new_to_old_0b'}
        How to interpret `mapping`.
    cols : iterable of str, optional
        Columns to compare. Defaults to the intersection of columns in both DataFrames.
    atol, rtol : float
        Absolute/relative tolerances for float comparisons (np.allclose).
    return_details : bool
        If True, return a dict report. If False, return a boolean pass/fail.

    Returns
    -------
    bool or dict
        If return_details=False, True only if all chosen columns match under the permutation.
        If return_details=True, returns a report with per-column mismatch stats.
    """
    # --- Basic shape checks ---
    n_before = len(df_before)
    n_after  = len(df_after)
    if n_before != n_after:
        out = {
            "passed": False,
            "reason": f"Row count mismatch: before={n_before}, after={n_after}",
        }
        return out if return_details else False
    n = n_before

    # --- Normalize mapping to arrays: old_to_new and new_to_old (0-based) ---
    # Convert mapping to np.ndarray
    if isinstance(mapping, dict):
        # Build dense array
        arr = np.full(n, -1, dtype=int)
        for k, v in mapping.items():
            if not (isinstance(k, int) and isinstance(v, int)):
                raise ValueError("mapping dict must be int -> int")
            arr_index = k - 1 if mapping_format.endswith("_1b") else k
            arr_value = v - 1 if mapping_format.endswith("_1b") else v
            if not (0 <= arr_index < n and 0 <= arr_value < n):
                raise ValueError(f"Mapping out of range for n={n}: {k}->{v}")
            arr[arr_index] = arr_value
        # Fill unspecified entries with identity (if dict was partial)
        mask = arr < 0
        arr[mask] = np.arange(n, dtype=int)[mask]
    else:
        arr = np.asarray(mapping, dtype=int)
        if arr.ndim != 1 or arr.size != n:
            raise ValueError("mapping array must be 1D with length equal to number of rows")

        # Convert to 0-based if needed
        if mapping_format.endswith("_1b"):
            if arr.min() < 1 or arr.max() > n:
                raise ValueError("1-based mapping has values out of range")
            arr = arr - 1

    # Validate permutation property
    # Determine whether `arr` is old_to_new or new_to_old based on mapping_format
    if mapping_format.startswith("old_to_new"):
        old_to_new = arr.copy()
        # Invert
        if len(np.unique(old_to_new)) != n:
            raise ValueError("Mapping is not a bijection (duplicate targets in old_to_new).")
        new_to_old = np.empty(n, dtype=int)
        new_to_old[old_to_new] = np.arange(n, dtype=int)
    else:
        new_to_old = arr.copy()
        if len(np.unique(new_to_old)) != n:
            raise ValueError("Mapping is not a bijection (duplicate sources in new_to_old).")
        old_to_new = np.empty(n, dtype=int)
        old_to_new[new_to_old] = np.arange(n, dtype=int)

    # Sanity: verify ranges
    if not (np.array_equal(np.sort(old_to_new), np.arange(n)) and
            np.array_equal(np.sort(new_to_old), np.arange(n))):
        raise ValueError("Mapping must be a permutation of [0..n-1].")

    # --- Reorder df_before to the order of df_after using new_to_old ---
    # new_to_old[j] = i  means: row j in df_after equals row i in df_before
    # so df_before.iloc[new_to_old] should line up with df_after (row-wise).
    aligned_before = df_before.iloc[new_to_old].reset_index(drop=True)
    aligned_after  = df_after.reset_index(drop=True)

    # Choose columns to compare
    if cols is None:
        cols = [c for c in aligned_before.columns.intersection(aligned_after.columns)]
    else:
        cols = [c for c in cols if c in aligned_before.columns and c in aligned_after.columns]
    if not cols:
        out = {"passed": False, "reason": "No overlapping columns to compare."}
        return out if return_details else False

    # --- Compare column by column ---
    per_col = {}
    all_ok = True
    for c in cols:
        s1 = aligned_before[c]
        s2 = aligned_after[c]

        # Align dtypes where possible
        # If both numeric → use allclose for floats, equality for ints
        if pd.api.types.is_numeric_dtype(s1) and pd.api.types.is_numeric_dtype(s2):
            # Promote to float for tolerance if any is float-like
            if pd.api.types.is_float_dtype(s1) or pd.api.types.is_float_dtype(s2):
                a = pd.to_numeric(s1, errors="coerce").to_numpy(dtype=float)
                b = pd.to_numeric(s2, errors="coerce").to_numpy(dtype=float)
                mask_nan = np.isnan(a) | np.isnan(b)
                close = np.allclose(a[~mask_nan], b[~mask_nan], rtol=rtol, atol=atol)
                stat = {
                    "type": "float",
                    "equal": bool(close and np.all(mask_nan == np.isnan(b))),
                    "max_abs_diff": float(np.nanmax(np.abs(a - b))) if a.size else 0.0,
                    "n_nan_mismatch": int(np.sum(np.isnan(a) ^ np.isnan(b))),
                }
                ok = stat["equal"]
            else:
                a = s1.to_numpy()
                b = s2.to_numpy()
                eq = (a == b)
                stat = {
                    "type": "int",
                    "equal": bool(np.all(eq)),
                    "n_mismatch": int((~eq).sum()),
                }
                ok = stat["equal"]
        else:
            # Fallback exact match for non-numeric/string/categorical
            a = s1.astype("object").to_numpy()
            b = s2.astype("object").to_numpy()
            eq = (a == b)
            # Treat NaNs as equal
            nanmask = pd.isna(a) & pd.isna(b)
            eq = np.where(nanmask, True, eq)
            stat = {
                "type": "object",
                "equal": bool(np.all(eq)),
                "n_mismatch": int((~eq).sum()),
            }
            ok = stat["equal"]

        per_col[c] = stat
        all_ok &= ok

    if not return_details:
        return bool(all_ok)

    return {
        "passed": bool(all_ok),
        "n_rows": n,
        "mapping_ok": True,
        "mapping_summary": {
            "format_interpreted": mapping_format,
            "old_to_new_example": old_to_new[:min(10, n)].tolist(),
        },
        "columns_checked": list(cols),
        "per_column": per_col,
    }
