"""
This module contains some classes and functions that are used from other modules
"""
from pathlib import Path
import subprocess
import os

from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem.rdmolops import RemoveHs

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
def com_2_xyz_2_sdf(args,start_point=None):
    """
    com to xyz to sdf for obabel

    Parameters
    ----------
    args : argparse.args
        [description]
    start_point : str, optional
        file(path/name?) to the starting point, by default None

    Returns
    -------
    int?
        charge or None? 
    """
    input = args.input
    default_charge = args.default_charge

    extension = os.path.splitext(input)[1]

    if start_point is None:
        if extension in ['.com','.gjf','.xyz']:
            file = args.input

    elif start_point is not None:
        file = start_point

    filename = os.path.splitext(file)[0]

    if extension != '.xyz':                                                  #  RAUL: Originally this pointed towards args.input, shouldn't it be to args.file? 
        with open(file,"r") as comfile:
            comlines = comfile.readlines()

        emptylines=[]

        for i, line in enumerate(comlines):
            if len(line.strip()) == 0:
                emptylines.append(i)

        #assigning the charges
        charge_com = comlines[(emptylines[1]+1)].split(' ')[0]

        with open(f'{filename}.xyz','w') as xyzfile:
            xyzfile.write(str(emptylines[2]- (emptylines[1]+2)))
            xyzfile.write('\n')
            xyzfile.write(filename)
            xyzfile.write('\n')
            for i in range((emptylines[1]+2), emptylines[2]):
                xyzfile.write(comlines[i])

    cmd_obabel = ['obabel',                                                 # RAUL: Again this could be done with openbabel's pybel
                  '-ixyz', f'{filename}.xyz', 
                  '-osdf', '-O', f'{filename}.sdf']
    subprocess.run(cmd_obabel)

    if start_point is None:
        if extension in ['.com','.gjf']:
            return charge_com
        else:
            return default_charge

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
