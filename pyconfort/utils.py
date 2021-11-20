"""
This module contains some classes and functions that are used from other modules
"""
from pathlib import Path
import subprocess
import os

from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem.rdmolops import RemoveHs

import pybel

items = """X
           H                                                                                                  He
          Li Be  B                                                                             C   N   O   F  Ne
          Na Mg Al                                                                            Si   P   S  Cl  Ar
           K Ca Sc                                           Ti  V Cr Mn Fe Co Ni Cu  Zn  Ga  Ge  As  Se  Br  Kr
          Rb Sr  Y                                           Zr Nb Mo Tc Ru Rh Pd Ag  Cd  In  Sn  Sb  Te   I  Xe
          Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta  W Re Os Ir Pt Au  Hg  Tl  Pb  Bi  Po  At  Rn
          Fr Ra Ac Th Pa  U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo
        """
possible_atoms = items.replace('\n',' ').strip().split()
possible_atoms[0] = ''

#class for logging
class Logger:
    """
    Class that wraps a file object to abstract the logging.
    """
    # Class Logger to writargs.input.split('.')[0] output to a file
    def __init__(self, filein, append, suffix='dat'):
        self.log = open(f'{filein}_{append}.{suffix}', 'w')

    def write(self, message):
        """
        Appends a newline character to the message and writes it into the file.
   
        Parameters
        ----------
        message : str
           text to be written in the log file.
        """
        self.log.write(f'{message}\n')

    def fatal(self, message):
        """
        Writes the message to the file. Closes the file and raises an error exit
   
        Parameters
        ----------
        message : str
           text to be written in the log file.
        """
        self.write(message)
        self.finalize()
        raise SystemExit(1)

    def finalize(self):
        """ Closes the file """
        self.log.close()

# OS utils 

def move_file(file, destination):
    """
    Moves a file to a destination folder. If the file exists rewrites it.

    Parameters
    ----------
    file : str
        path towards a file
    destination : str
        path towards a folder
    """
    dest = Path(destination)
    dest.mkdir(exist_ok=True)
    filepath = Path(file)
    filename = filepath.name
    file.rename(dest/filename)

def move_file_from_folder(destination,src,file):
    """
    Moves files from the source folder to the destination folder and creates 
    the destination folders when needed.

    Parameters
    ----------
    destination : str
        path to the destination folder
    src : str
        path to the source folder
    file : str
        full name of the file (file + extension)
    """
    dest = Path(destination)
    source = Path(src)
    dest.mkdir(exist_ok=True,parents=True)
    filepath = source / file
    filepath.rename(dest / file)

# openbabel utils 

#com to xyz to sdf for obabel
def com_2_xyz_2_sdf(input,default_charge,start_point=None):
    """
    com to xyz to sdf for obabel

    Parameters
    ----------
    input : str
        path to the file to convert
    start_point : str, optional
        file(path/name?) to the starting point, by default None

    Returns
    -------
    int?
        charge or None? 
    """
    extension = Path(input).suffix

    if start_point is None:
        if extension in ['com','gjf','xyz']:
            file = Path(input)

    else:
        file = Path(start_point)

    filename = Path.stem
    
    # Create the 'xyz' file and/or get the total charge
    if extension != 'xyz':                                                      #  RAUL: Originally this pointed towards args.input, shouldn't it be to args.file?
        xyz,charge = get_charge_and_xyz_from_com(file)
        xyz_txt = '\n'.join(xyz)
        with open(f'{filename}.xyz','w') as F: 
            F.write(f"{len(xyz)}\n{filename}\n{xyz_txt}\n")
    else:
        charge = default_charge

    xyz_2_sdf(f'{filename}.xyz')

    return charge
def xyz_2_sdf(file,parent_dir=None):
    """
    Creates a .sdf file from a .xyz in the specified directory. If no directory
    is specified then the files are created in the current directory. 

    Parameters
    ----------
    file : str
        filename and extension of an existing .xyz file
    dir : str or pathlib.Path, optional
        a path to the directory where the .xyz file is located
    """
    if parent_dir is None: 
        parent_dir = Path('')
    else:
        parent_dir = Path(parent_dir)
    mol = next(pybel.readfile('xyz',parent_dir/file))
    ofile = Path(file).stem + '.sdf'
    mol.write('sdf', parent_dir/ofile)
def get_charge_and_xyz_from_com(file):
    """
     Takes a .gjf or .com file and retrieves the coordinates of the atoms and the
    total charge. 

    Parameters
    ----------
    file : str or pathlib.Path
        A path pointing to a valid .com or .gjf file

    Returns
    -------
    coordinates : list
        A list of strings (without \\n) that contain the xyz coordinates of the
        .gjf or .com file 
    charge : str
        A str with the number corresponding to the total charge of the .com or 
        .gjf file
    """

    with open(file,"r") as comfile:
        comlines = comfile.readlines()

    _iter = comlines.__iter__()
    
    line = ''
    # Find the command line
    while '#' not in line: 
        line = next(_iter)

    # pass the title lines
    _ = next(_iter)
    while line:
        line = next(_iter).strip()
    
    # Read the charge from the charge | spin line
    charge,spin = next(_iter).strip().split()

    # Store the coordinates until next empty line. 
    coordinates = []
    line = next(_iter).strip()
    while line: 
        coordinates.append(line.strip())
        line = next(_iter).strip()
    
    return coordinates, charge
    
# RDKit Utils 

def set_metal_atomic_number(mol,metal_idx,metal_sym):
    """
    Changes the atomic number of the metal atoms using their indices. 

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        rdkit molecule object
    metal_idx : list
        sorted list that contains the indices of the metal atoms in the molecule
    metal_sym : list
        sorted list (same order as metal_idx) that contains the symbols of the 
        metals in the molecule. 
    """

    for atom in mol.GetAtoms():
        if atom.GetIdx() in metal_idx:
            re_symbol = metal_sym[metal_idx.index(atom.GetIdx())]
            atomic_number = possible_atoms.index(re_symbol)
            atom.SetAtomicNum(atomic_number)

def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_RMSD):
    """
    Takes in two rdkit.Chem.Mol objects and calculates the RMSD between them.
    (As side efect mol1 is left in the aligned state, if heavy is specified 
    the side efect will not happen)

    Parameters
    ----------
    mol1 : rdkit.Chem.Mol
        Probe molecule
    mol2 : rdkit.Chem.Mol
        Target molecule. The probe is aligned to the target to compute the RMSD.
    c1 : int
        conformation of mol1 to use for the RMSD
    c2 : int
        conformation of mol2 to use for the RMSD
    heavy : bool
        If True it will ignore the H atoms when computing the RMSD.
    max_matches_RMSD : int
        the max number of matches found in a SubstructMatch()

    Returns
    -------
    float
        Returns the best RMSD found.
    """
    if heavy:
         mol1 = RemoveHs(mol1)
         mol2 = RemoveHs(mol2)
    return GetBestRMS(mol1,mol2,c1,c2,maxMatches=max_matches_RMSD)

