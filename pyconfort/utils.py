"""
This module contains some classes and functions that are used from other modules
"""
import shutil
from pathlib import Path

from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem.rdmolops import RemoveHs
from openbabel import pybel
import os
import yaml
import pandas as pd
from rdkit.Chem import AllChem as Chem

def periodic_table():
    items = """X
            H                                                                                                  He
            Li Be  B                                                                             C   N   O   F  Ne
            Na Mg Al                                                                            Si   P   S  Cl  Ar
            K Ca Sc                                           Ti  V Cr Mn Fe Co Ni Cu  Zn  Ga  Ge  As  Se  Br  Kr
            Rb Sr  Y                                           Zr Nb Mo Tc Ru Rh Pd Ag  Cd  In  Sn  Sb  Te   I  Xe
            Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta  W Re Os Ir Pt Au  Hg  Tl  Pb  Bi  Po  At  Rn
            Fr Ra Ac Th Pa  U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo
            """
    periodic_table = items.replace('\n',' ').strip().split()
    periodic_table[0] = ''
    return periodic_table

## Bondi VDW radii in Angstrom
def bondi():
    bondi = {"Bq": 0.00, "H": 1.09,"He": 1.40,
        "Li":1.81,"Be":1.53,"B":1.92,"C":1.70,"N":1.55,"O":1.52,"F":1.47,"Ne":1.54,
        "Na":2.27,"Mg":1.73,"Al":1.84,"Si":2.10,"P":1.80,"S":1.80,"Cl":1.75,"Ar":1.88,
        "K":2.75,"Ca":2.31,"Ni": 1.63,"Cu":1.40,"Zn":1.39,"Ga":1.87,"Ge":2.11,"As":1.85,"Se":1.90,"Br":1.83,"Kr":2.02,
        "Rb":3.03,"Sr":2.49,"Pd": 1.63,"Ag":1.72,"Cd":1.58,"In":1.93,"Sn":2.17,"Sb":2.06,"Te":2.06,"I":1.98,"Xe":2.16,
        "Cs":3.43,"Ba":2.68,"Pt":1.72,"Au":1.66,"Hg":1.55,"Tl":1.96,"Pb":2.02,"Bi":2.07,"Po":1.97,"At":2.02,"Rn":2.20,
        "Fr":3.48,"Ra":2.83, "U":1.86 }
    return bondi

## covalent radii in Angstrom (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
def rcov():
    rcov = {"H": 0.32,"He": 0.46,
    "Li":1.33,"Be":1.02,"B":0.85,"C":0.75,"N":0.71,"O":0.63,"F":0.64,"Ne":0.67,
    "Na":1.55,"Mg":1.39,"Al":1.26, "Si":1.16,"P":1.11,"S":1.03,"Cl":0.99, "Ar":0.96,
    "K":1.96,"Ca":1.71,"Sc": 1.48, "Ti": 1.36, "V": 1.34, "Cr": 1.22, "Mn":1.19, "Fe":1.16, "Co":1.11, "Ni":1.10,"Zn":1.18, "Ga":1.24, "Ge":1.21, "As":1.21, "Se":1.16, "Br":1.14, "Kr":1.17,
    "Rb":2.10, "Sr":1.85,"Y":1.63, "Zr":1.54, "Nb":1.47, "Mo":1.38, "Tc":1.28, "Ru":1.25,"Rh":1.25,"Pd":1.20,"Ag":1.28,"Cd":1.36, "In":1.42, "Sn":1.40,"Sb":1.40,"Te":1.36,"I":1.33,"Xe":1.31}
    return rcov


#load paramters from yaml file
def load_from_yaml(args_,log):
    """
    Loads the parameters for the calculation from a yaml if specified. Otherwise
    does nothing.

    Parameters
    ----------
    args : argparse.args
        Dataclass
    log : Logger
        Where to log the program progress
    """
    # Variables will be updated from YAML file
    try:
        if args_.varfile is not None:
            if os.path.exists(args_.varfile):
                if os.path.splitext(args_.varfile)[1] == '.yaml':
                    log.write("\no  Importing pyCONFORT parameters from " + args_.varfile)
                    with open(args_.varfile, 'r') as file:
                        try:
                            param_list = yaml.load(file, Loader=yaml.SafeLoader)
                        except yaml.scanner.ScannerError:
                            log.write("\nx  Error while reading " + args_.varfile+ ". Edit the yaml file and try again (i.e. use ':' instead of '=' to specify variables)")
                            sys.exit()
            for param in param_list:
                if hasattr(args_, param):
                    if getattr(args_, param) != param_list[param]:
                        log.write("o  RESET " + param + " from " + str(getattr(args_, param)) + " to " + str(param_list[param]))
                        setattr(args_, param, param_list[param])
                    else:
                        log.write("o  DEFAULT " + param + " : " + str(getattr(args_, param)))
    except UnboundLocalError: # RAUL: Is it just me or this only happens when the file exists, and ens in .yaml and is empty or does not end in .yaml?
        log.write("\no  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code.\n")

    return args_, log

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
def creation_of_dup_csv(csearch, cmin):

    """
    Generates a pandas.DataFrame object with the appropiate columns for the
    conformational search and the minimization.

    Parameters
    ----------
    csearch : str
        Conformational search method. Current valid methods are:
        ['rdkit','fullmonte','summ']
    cmin : str
        Minimization method. Current valid methods are:
        ['xtb','ani']

    Returns
    -------
    pandas.DataFrame
    """
    # Boolean aliases from args
    is_rdkit = csearch=='rdkit'
    is_fullmonte = csearch=='fullmonte'
    is_crest = csearch=='crest'
    is_summ = csearch=='summ'
    is_xtb = cmin == 'xtb'
    is_ani = cmin == 'ani'

    # column blocks definitions
    base_columns = ['Molecule',
                    'RDKit-Initial-samples',
                    'RDKit-energy-window',
                    'RDKit-initial_energy_threshold',
                    'RDKit-RMSD-and-energy-duplicates',
                    'RDKit-Unique-conformers']
    xtb_columns = ['xTB-Initial-samples',
                   'xTB-energy-window',
                   'xTB-initial_energy_threshold',
                   'xTB-RMSD-and-energy-duplicates',
                   'xTB-Unique-conformers']
    ANI_columns = ['ANI-Initial-samples',
                   'ANI-energy-window',
                   'ANI-initial_energy_threshold',
                   'ANI-RMSD-and-energy-duplicates',
                   'ANI-Unique-conformers']
    end_columns_no_min = ['CSEARCH time (seconds)',
                          'Overall charge']
    end_columns_min = ['CSEARCH time (seconds)',
                       'CMIN time (seconds)',
                       'Overall charge']
    fullmonte_columns = ['FullMonte-Unique-conformers',]
                        #'FullMonte-conformers',
                        #'FullMonte-energy-window',
                        #'FullMonte-initial_energy_threshold',
                        #'FullMonte-RMSD-and-energy-duplicates']
    summ_columns = ['summ-conformers',
                    'summ-energy-window',
                    'summ-initial_energy_threshold',
                    'summ-RMSD-and-energy-duplicates',
                    'summ-Unique-conformers']
    crest_columns = ['Molecule','crest-conformers']

    # Check Conformer Search method
    if is_rdkit:

        columns = base_columns
    elif is_fullmonte:
        columns = base_columns + fullmonte_columns
    elif is_summ:
        columns = base_columns + summ_columns
    elif is_crest:
        columns = crest_columns
    else:
        return None

    # Check Minimization Method
    if is_ani:
        columns += ANI_columns
    if is_xtb:  # is_ani and is_xtb will not happen, but this is what was written
        columns += xtb_columns
    if is_ani or is_xtb:
        columns += end_columns_min
    else:
        columns += end_columns_no_min
    return pd.DataFrame(columns=columns)


def move_file(file, source, destination):
    """
    Moves a file to a destination folder. If the file exists rewrites it.

    Parameters
    ----------
    file : str
        path towards a file
    destination : str
        path towards a folder
    """
    # if not os.path.isdir(destination):
    #     os.makedirs(destination)

    # shutil.move(os.path.join(source, file), os.path.join(destination, file))
    pass


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

    # IF YOU RUN THE PROGRAM 2 TIMES WITH CMIN, YOU GET THIS ERROR:
    # FileExistsError: [WinError 183] Cannot create a file when that file already exists
    # FIX IT!
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
        xyz,charge = get_info_com(file)
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

def get_info_com(file):
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

    # in case the keywords are distributed in multiple lines
    while len(line.split()) > 0:
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
# DETECTS DIHEDRALS IN THE MOLECULE
def getDihedralMatches(mol, heavy,log):
    #this is rdkit's "strict" pattern
    pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
    qmol = Chem.MolFromSmarts(pattern)
    matches = mol.GetSubstructMatches(qmol)

    #these are all sets of 4 atoms, uniquify by middle two
    uniqmatches = []
    seen = set()
    for (a,b,c,d) in matches:
        if (b,c) not in seen and (c,b) not in seen:
            if heavy:
                if mol.GetAtomWithIdx(a).GetSymbol() != 'H' and mol.GetAtomWithIdx(d).GetSymbol() != 'H':
                    seen.add((b,c))
                    uniqmatches.append((a,b,c,d))
            if not heavy:
                if mol.GetAtomWithIdx(c).GetSymbol() == 'C' and mol.GetAtomWithIdx(d).GetSymbol() == 'H':
                    pass
                else:
                    seen.add((b,c))
                    uniqmatches.append((a,b,c,d))
    return uniqmatches

def getDihedralMatches_v2(mol, heavy,log): #New version using openbabel
    # If this version is selected, bring the import to the top
    # import pybel #
    from openbabel import pybel #for openbabel>=3.0.0
    AtomInTripleBond = '$(*#*)'
    TerminalAtom = 'D1'
    CF3 = '$(C(F)(F)F)'
    CCl3 = '$(C(Cl)(Cl)Cl)'
    CBr3 = '$(C(Br)(Br)Br)'
    tBut = '$(C([CH3])([CH3])[CH3])'
    # A 3-bonded C with a double bond to (N, O or S)
    # singlgy bonded to not ring bonded to a non-terminal N,O or S.
    CD3_1d = '$([CD3](=[N,O,S])-!@[#7,O,S!D1])'
    CD3_1r = '$([#7,O,S!D1]-!@[CD3]=[N,O,S])' # Backwards version
    # A 3-bonded C with a double bond to (N+)
    # singlgy bonded to not ring bonded to Any non-terminal N
    CD3_2d = '$([CD3](=[N+])-!@[#7!D1])'
    CD3_2r = '$([#7!D1]-!@[CD3]=[N+])' # Backwards version
    Atom1 = '*'
    Atom2 =  f'!{AtomInTripleBond}&!{TerminalAtom}'
    Atom2 += f'&!{CF3}&!{CCl3}&!{CBr3}&!{tBut}'
    Atom2 += f'&!{CD3_1d}&!{CD3_1r}&!{CD3_2d}&!{CD3_2r}'
    Atom3 =  f'!{AtomInTripleBond}&!{TerminalAtom}'
    Atom3 += f'&!{CF3}&!{CCl3}&!{CBr3}&!{tBut}'
    Atom4 = '*'
    pattern = f'{Atom1}~[{Atom2}]-!@[{Atom3}]~{Atom4}'
    smarts = pybel.Smarts(pattern)
    matches = smarts.findall(mol)

    #these are all sets of 4 atoms, uniquify by middle two
    H_atoms = set(pybel.Smarts('[#1]').findall(mol))
    C_atoms = set(pybel.Smarts('[#6]').findall(mol))
    uniqmatches = []
    seen = set()
    for (a,b,c,d) in matches:
        if (b,c) not in seen and (c,b) not in seen:
            if heavy:
                if a not in H_atoms and d not in H_atoms:
                    seen.add((b,c))
                    uniqmatches.append((a,b,c,d))
            if not heavy:
                # So what if a == 'H' and b == 'C'? is that valid ¿?
                if c not in C_atoms or d not in H_atoms:
                    seen.add((b,c))
                    uniqmatches.append((a,b,c,d))



#checks for salts
def check_for_pieces(smi):
    #taking largest component for salts
    pieces = smi.split('.')
    if len(pieces) > 1:
        # take largest component by length
        smi = max(pieces, key=len)
    return smi

def mol_from_sdf_or_mol_or_mol2(input):
    """
    mol from sdf

    Parameters
    ----------
    input : str
        path to a .sdf .mol or .mol2 file

    Returns
    -------
    tuple of lists?
        suppl, IDs, charges
    """
    filename = os.path.splitext(input)[0]
    extension = os.path.splitext(input)[1]

    if extension =='.sdf':
        suppl = Chem.SDMolSupplier(input, removeHs=False)
    elif extension =='.mol':
        suppl = Chem.MolFromMolFile(input, removeHs=False)
    elif extension =='.mol2':
        suppl = Chem.MolFromMol2File(input, removeHs=False)

    IDs,charges = [],[]

    with open(input,"r") as F:
        lines = F.readlines()

    molecule_count = 0
    for i, line in enumerate(lines):
        if line.find('>  <ID>') > -1:
            ID = lines[i+1].split()[0]
            IDs.append(ID)
        if line.find('M  CHG') > -1:
            charge_line =  line.split('  ')
            charge = 0
            for j in range(4,len(charge_line)):
                if (j % 2) == 0:
                    if j == len(charge_line) - 1:
                        charge_line[j] = charge_line[j].split('\n')[0]
                    charge += int(charge_line[j])
            charges.append(charge)
        if line.find('$$$$') > -1:
            molecule_count += 1
            if molecule_count != len(charges):
                charges.append(0)

    if len(IDs) == 0:
        if extension == '.sdf':
            for i in range(len(suppl)):
                IDs.append(f'{filename}_{i}')
        else:
            IDs.append(filename)
    if len(charges) == 0:
        if extension == '.sdf':
            for _ in suppl:
                charges.append(0)
        else:
            charges.append(0)
    return suppl, IDs, charges

#checks the charge on the smi string
def check_charge_smi(smi):
    charge = 0
    for i,smi_letter in enumerate(smi):
        if smi_letter == '+':
            if smi[i+1] == ']':
                charge += 1
            else:
                charge += int(smi[i+1])
        elif smi_letter == '-':
            if smi[i+1] == ']':
                charge -= 1
            else:
                charge -= int(smi[i+1])
    return charge

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
            atomic_number = periodic_table.index(re_symbol)
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
