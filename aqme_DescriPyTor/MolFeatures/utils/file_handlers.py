import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from . import help_functions as hf
import chardet 

def adjust_working_dir(working_dir=None):
    if working_dir:
        os.chdir(working_dir)


def process_csv_file_pd(csv_filepath):
    raw_data=pd.read_csv(csv_filepath, header=None)
    return raw_data

def process_identifiers_csv_file(csv_filepath):
    # Detect the encoding of the file
    with open(csv_filepath, 'rb') as f:
        result = chardet.detect(f.read())  # or read a larger portion if needed
        encoding = result['encoding']

    # Use the detected encoding to read the CSV
    raw_data = pd.read_csv(csv_filepath, header=None, encoding=encoding)

    names_list=raw_data[0]
    identifier_list=raw_data[1]
    return names_list, identifier_list



def generate_3d_coordinates(smiles):
    # Fix common issues in SMILES strings and test if it's valid
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # Handle invalid SMILES string
        print(f"Invalid SMILES string: {smiles}")
        return None
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates using the ETKDG method
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) == -1:
        # Handle failed 3D coordinate generation
        print(f"Could not generate 3D coordinates for: {smiles}")
        return None
    
    # Get the molecule's conformer
    conformer = mol.GetConformer()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    
    # Combine atomic symbols with 3D positions
    array = np.concatenate((np.array(symbols, dtype=object)[:, np.newaxis], conformer.GetPositions()), axis=1)
    return array



def smiles_to_coordinates(smiles):
    return [generate_3d_coordinates(smile) for smile in smiles]

def save_single_xyz_file(symbols, coordinates, output_filename, comment_line='comment_line'):
    """
    The function gets a the atom symbols and coordinates, and the function creates a xyz file from it.
    ----------
    Parameters
    ----------
    symbols : iterable.
        Array of atom types, each item is a string

    coordinates : iterable.
        Array of cartesian coordinates of all atoms, each item is a float

    output_filename : str.
        The name given to the output xyz file

    comment_line : str. default 'comment_line'
        A line recorded into the second line of the xyz file
    -------
    Returns
    -------
    None
    --------
    Examples
    --------
    symbols=['O', 'H', 'H']
    coordinates=[[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]]
    output_filename='water.xyz'
    comment_line='Water molecule'
    save_single_xyz_file(symbols, coordinates, output_filename, comment_line)
    --water.xyz--
    3
    Water molecule
    O          0.00000        0.00000        0.11779
    H          0.00000        0.75545       -0.47116
    H          0.00000       -0.75545       -0.47116
    --end_of_file--
    """
    with open(output_filename, 'w') as my_file:
        my_file.write('{}\n'.format(len(symbols)))
        my_file.write(comment_line+'\n')
        for loop_index, (x, y, z) in enumerate(coordinates):
            my_file.write('{}\t {:.5f}\t {:.5f}\t {:.5f}\t\n'.format(symbols[loop_index], x, y, z))


def use_sh_script(file_name):
    import subprocess
    os.chdir(r'C:\Users\edens\Documents\GitHub\Automation_code-main\M1_pre_calculations')
    subprocess.call(['sh', file_name])
    # result = subprocess.run(('/'+file_name), shell=True, capture_output=True, text=True)




def write_gaussian_file(filename, functional, basisset, charge, title, task, freeze):
    """
    Writes the Gaussian input file based on the provided parameters.
    """
    name = filename.split('.')[0]
    with open(f'{name}.com', 'w') as my_file:
        my_file.write(f"%mem=100GB\n%nproc=32\n%chk={name}.chk\n")
        if task == 'sp':
            my_file.write(f"#P {functional}/{basisset}  Freq \n\n{title}\n\n{charge}\n")
        elif task == 'opt':
            my_file.write(f"# opt Freq {functional}/{basisset}  pop=(full,nbo) \n\n{title}\n\n{charge}\n") ## opt freq def2tzv pop=npa m062x
        xyz_df = hf.get_df_from_file(filename)
        atoms_np_array = np.array(xyz_df)
        for atom_np_array in atoms_np_array:
            try:
                my_file.write("{:1} {:11.5f} {:11.5f} {:11.5f}\n".format(*atom_np_array))
            except Exception as e:
                print(f'Error writing atom data: {e}')
        if freeze:
            my_file.write(f"\n{freeze}\n")
        my_file.write("\n")

def process_file(file, functional, basisset, charge, title, task, nbo_answer):
    """
    Processes a single file for Gaussian input file creation.
    """
    
    write_gaussian_file(file, functional, basisset, charge, title, task)

def xyz_to_gaussian_file(filename, functional='HF', basisset='6-31G(d)', charge='0 1', title='title', task='sp'):
    process_file(filename, functional, basisset, charge, title, task)

def xyz_files_list_to_gaussian_files(xyz_list, functional='HF', basisset='6-31G(d)', charge='0 1', title='title', task='sp'):
    if not os.path.exists('com'):
        os.mkdir('com')
    for file in xyz_list:
        process_file(file, functional, basisset, charge, title, task)
        name = file.split('.')[0]
        os.rename(f'{name}.com', f'com/{name}.com')
