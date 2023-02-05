######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import os
import subprocess
import sys
import time
import getopt
import numpy as np
import glob
import yaml
import ast
from pathlib import Path
from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem.rdmolops import RemoveHs
from rdkit.Chem import Mol
from rdkit.Chem import AllChem as Chem
from aqme.argument_parser import set_options, var_dict
from rdkit import RDLogger

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15

aqme_version = "1.4.2"
time_run = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
aqme_ref = f"AQME v {aqme_version}, Alegre-Requena, J. V.; Sowndarya, S.; Perez-Soto, R.; Alturaifi, T. M.; Paton, R. S., 2022. https://github.com/jvalegre/aqme"

RDLogger.DisableLog("rdApp.*")


def run_command(command, outfile):
    """
    Runs the subprocess command and outputs to the necessary output file
    """

    p2 = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    txt = [line.decode("utf-8") for line in p2.stdout]
    p2.stdout.close()

    with open(outfile, "w") as f:
        for line in txt:
            f.write(line)
    f.close()


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
    periodic_table = items.replace("\n", " ").strip().split()
    periodic_table[0] = ""

    return periodic_table


# load paramters from yaml file
def load_from_yaml(self):
    """
    Loads the parameters for the calculation from a yaml if specified. Otherwise
    does nothing.
    """

    txt_yaml = f"\no  Importing AQME parameters from {self.varfile}"
    error_yaml = False
    # Variables will be updated from YAML file
    try:
        if os.path.exists(self.varfile):
            if os.path.splitext(self.varfile)[1] in [".yaml", ".yml", ".txt"]:
                with open(self.varfile, "r") as file:
                    try:
                        param_list = yaml.load(file, Loader=yaml.SafeLoader)
                    except yaml.scanner.ScannerError:
                        txt_yaml = f'\nx  Error while reading {self.varfile}. Edit the yaml file and try again (i.e. use ":" instead of "=" to specify variables)'
                        error_yaml = True
        if not error_yaml:
            for param in param_list:
                if hasattr(self, param):
                    if getattr(self, param) != param_list[param]:
                        setattr(self, param, param_list[param])

    except UnboundLocalError:
        txt_yaml = "\nx  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code."

    return self, txt_yaml


# class for logging
class Logger:
    """
    Class that wraps a file object to abstract the logging.
    """

    # Class Logger to writargs.input.split('.')[0] output to a file
    def __init__(self, filein, append, suffix="dat"):
        self.log = open(f"{filein}_{append}.{suffix}", "w")

    def write(self, message):
        """
        Appends a newline character to the message and writes it into the file.

        Parameters
        ----------
        message : str
           Text to be written in the log file.
        """
        self.log.write(f"{message}\n")
        print(f"{message}\n")

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
        """
        Closes the file
        """
        self.log.close()


def move_file(destination, source, file):
    """
    Moves files from the source folder to the destination folder and creates
    the destination folders when needed.

    Parameters
    ----------
    destination : str
        Path to the destination folder
    src : str
        Path to the source folder
    file : str
        Full name of the file (file + extension)
    """

    destination.mkdir(exist_ok=True, parents=True)
    filepath = source / file
    try:
        filepath.rename(destination / file)
    except FileExistsError:
        filepath.replace(destination / file)


def get_info_input(file):
    """
    Takes an input file and retrieves the coordinates of the atoms and the
    total charge.

    Parameters
    ----------
    file : str or pathlib.Path
        A path pointing to a valid .com or .gjf file

    Returns
    -------
    coordinates : list
        A list of strings (without \\n) that contain the xyz coordinates of the .gjf or .com file
    charge : str
        A str with the number corresponding to the total charge of the .com or .gjf file
    """

    with open(file, "r") as input_file:
        input_lines = input_file.readlines()

    _iter = input_lines.__iter__()

    line = ""
    # input for Gaussian calculations
    if str(file).split(".")[1] in ["com", "gjf"]:

        # Find the command line
        while "#" not in line:
            line = next(_iter)

        # in case the keywords are distributed in multiple lines
        while len(line.split()) > 0:
            line = next(_iter)

        # pass the title lines
        _ = next(_iter)
        while line:
            line = next(_iter).strip()

        # Read charge and multiplicity
        charge, mult = next(_iter).strip().split()

        # Store the atom types and coordinates until next empty line.
        atoms_and_coords = []
        line = next(_iter).strip()
        while line:
            atoms_and_coords.append(line.strip())
            line = next(_iter).strip()

    # input for ORCA calculations
    if str(file).split(".")[1] == ".inp":

        # Find the line with charge and multiplicity
        while "* xyz" not in line or "* int" not in line:
            line = next(_iter)

        # Read charge and multiplicity
        charge = line.strip().split()[-2]

        # Store the coordinates until next *
        atoms_and_coords = []
        line = next(_iter).strip()
        while len(line.split()) > 1:
            atoms_and_coords.append(line.strip())
            line = next(_iter).strip()
    return atoms_and_coords, charge, mult


def rules_get_charge(mol, args):
    """
    Automatically sets the charge for metal complexes
    """

    C_group = ["C", "Se", "Ge"]
    N_group = ["N", "P", "As"]
    O_group = ["O", "S", "Se"]
    F_group = ["F", "Cl", "Br", "I"]
    
    
    M_ligands, N_carbenes, bridge_atoms, neighbours = [], [], [], []
    charge_rules = np.zeros(len(mol.GetAtoms()), dtype=int)
    neighbours, metal_found, sanit_step = [], False, True
    for i, atom in enumerate(mol.GetAtoms()):
        # get the neighbours of metal atom and calculate the charge of metal center + ligands
        if atom.GetIdx() in args.metal_idx:
            # a sanitation step is needed to ensure that metals and ligands show correct valences
            if sanit_step:
                Chem.SanitizeMol(mol)
                sanit_step = False
            metal_found = True
            charge_idx = args.metal_idx.index(atom.GetIdx())
            neighbours = atom.GetNeighbors()
            charge_rules[i] = args.metal_oxi[charge_idx]
            for neighbour in neighbours:
                M_ligands.append(neighbour.GetIdx())
                if neighbour.GetTotalValence() == 4:
                    if neighbour.GetSymbol() in C_group:
                        # first, detects carbenes to adjust charge
                        carbene_like = False
                        bridge_ligand = False
                        for inside_neighbour in neighbour.GetNeighbors():
                            if inside_neighbour.GetSymbol() in N_group:
                                if inside_neighbour.GetTotalValence() == 4:
                                    for N_neighbour in inside_neighbour.GetNeighbors():
                                        # this option detects bridge ligands that connect two metals such as M--CN--M
                                        # we use I since the M is still represented as I at this point
                                        if N_neighbour.GetSymbol() == "I":
                                            bridge_ligand = True
                                            bridge_atoms.append(inside_neighbour.GetIdx())
                                    if not bridge_ligand:
                                        carbene_like = True
                                        N_carbenes.append(inside_neighbour.GetIdx())
                        if not carbene_like:
                            charge_rules[i] = charge_rules[i] - 1
                elif neighbour.GetTotalValence() == 3:
                    if neighbour.GetSymbol() in N_group and neighbour.GetFormalCharge() == 0:
                        charge_rules[i] = charge_rules[i] - 1
                elif neighbour.GetTotalValence() == 2:
                    # radical chalcogen atoms (i.e., Cu-OH(rad))
                    if neighbour.GetSymbol() in O_group and neighbour.GetFormalCharge() == 0 and len(neighbour.GetNeighbors()) == 2:
                        charge_rules[i] = charge_rules[i] - 1
                    # double bonded chalcogen atom (i.e., V=O)
                    elif neighbour.GetSymbol() in O_group and neighbour.GetFormalCharge() == 0 and len(neighbour.GetNeighbors()) == 1:
                        charge_rules[i] = charge_rules[i] - 2
                elif neighbour.GetTotalValence() == 1:
                    if neighbour.GetSymbol() in O_group:
                        charge_rules[i] = charge_rules[i] - 2
                    if neighbour.GetSymbol() in F_group:
                        charge_rules[i] = charge_rules[i] - 1

    # for charges not in the metal, neighbours or exceptions (i.e., C=N+ from carbenes or CN from bridge atoms)
    invalid_charged_atoms = M_ligands + N_carbenes + bridge_atoms + args.metal_idx
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetIdx() not in invalid_charged_atoms:
            charge_rules[i] = atom.GetFormalCharge()

    charge = np.sum(charge_rules)
    if not metal_found:
        # for organic molecules when using a list containing organic and organometallics molecules mixed
        charge = Chem.GetFormalCharge(mol)

    return charge, metal_found


def substituted_mol(self, mol, checkI):
    """
    Returns a molecule object in which all metal atoms specified in args.metal_atoms
    are replaced by Iodine and the charge is set depending on the number of
    neighbors.

    """

    self.args.metal_idx = []
    self.args.complex_coord = []
    self.args.metal_sym = []

    for _ in self.args.metal_atoms:
        self.args.metal_idx.append(None)
        self.args.complex_coord.append(None)
        self.args.metal_sym.append(None)

    Neighbors2FormalCharge = dict()
    for i, j in zip(range(2, 9), range(-3, 4)):
        Neighbors2FormalCharge[i] = j

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in self.args.metal_atoms:
            self.args.metal_sym[self.args.metal_atoms.index(symbol)] = symbol
            self.args.metal_idx[self.args.metal_atoms.index(symbol)] = atom.GetIdx()
            self.args.complex_coord[self.args.metal_atoms.index(symbol)] = len(
                atom.GetNeighbors()
            )
            if checkI == "I":
                atom.SetAtomicNum(53)
                n_neighbors = len(atom.GetNeighbors())
                if n_neighbors > 1:
                    formal_charge = Neighbors2FormalCharge[n_neighbors]
                    atom.SetFormalCharge(formal_charge)

    return self.args.metal_idx, self.args.complex_coord, self.args.metal_sym


def set_metal_atomic_number(mol, metal_idx, metal_sym):
    """
    Changes the atomic number of the metal atoms using their indices.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object
    metal_idx : list
        sorted list that contains the indices of the metal atoms in the molecule
    metal_sym : list
        sorted list (same order as metal_idx) that contains the symbols of the metals in the molecule
    """

    for atom in mol.GetAtoms():
        if atom.GetIdx() in metal_idx:
            re_symbol = metal_sym[metal_idx.index(atom.GetIdx())]
            atomic_number = periodic_table().index(re_symbol)
            atom.SetAtomicNum(atomic_number)


def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_rmsd):
    """
    Takes in two rdkit.Chem.Mol objects and calculates the RMSD between them.
    (As side efect mol1 is left in the aligned state, if heavy is specified
    the side efect will not happen)

    Parameters
    ----------
    mol1 : rdkit.Chem.Mol
        Probe molecule
    mol2 : rdkit.Chem.Mol
        Target molecule. The probe is aligned to the target to compute the RMSD
    c1 : int
        Conformation of mol1 to use for the RMSD
    c2 : int
        Conformation of mol2 to use for the RMSD
    heavy : bool
        If True it will ignore the H atoms when computing the RMSD
    max_matches_rmsd : int
        Max number of matches found in a SubstructMatch()

    Returns
    -------
    float
        Returns the best RMSD found
    """

    if heavy:
        mol1 = RemoveHs(mol1)
        mol2 = RemoveHs(mol2)
    return GetBestRMS(mol1, mol2, c1, c2, maxMatches=max_matches_rmsd)


def command_line_args():
    """
    Load default and user-defined arguments specified through command lines. Arrguments are loaded as a dictionary
    """

    # First, create dictionary with user-defined arguments
    kwargs = {}
    available_args = ["help"]
    bool_args = [
        "csearch",
        "cmin",
        "qprep",
        "qcorr",
        "qdescp",
        "vismol",
        "heavyonly",
        "cregen",
        "lowest_only",
        "lowest_n",
        "chk",
    ]

    for arg in var_dict:
        if arg in bool_args:
            available_args.append(f"{arg}")
        else:
            available_args.append(f"{arg} =")

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "h", available_args)
    except getopt.GetoptError as err:
        print(err)
        sys.exit()

    for arg, value in opts:
        if arg.find("--") > -1:
            arg_name = arg.split("--")[1].strip()
        elif arg.find("-") > -1:
            arg_name = arg.split("-")[1].strip()
        if arg_name in bool_args:
            value = True
        if value == "None":
            value = None
        if arg_name in ("h", "help"):
            print(f"o  AQME v {aqme_version} is installed correctly! For more information about the available options, see the documentation in https://github.com/jvalegre/aqme")
            sys.exit()
        else:
            # this "if" allows to use * to select multiple files in multiple OS
            if arg_name.lower() == 'files' and value.find('*') > -1:
                kwargs[arg_name] = glob.glob(value)
            else:
                # this converts the string parameters to lists
                if arg_name.lower() in ["files", "metal_oxi", "metal_atoms", "gen_atoms", "constraints_atoms", "constraints_dist", "constraints_angle", "constraints_dihedral", "atom_types", "cartesians", "nmr_atoms", "nmr_slope", "nmr_intercept"]:
                    if not isinstance(value, list):
                        try:
                            value = ast.literal_eval(value)
                        except (SyntaxError, ValueError):
                            pass
                kwargs[arg_name] = value

    # Second, load all the default variables as an "add_option" object
    args = load_variables(kwargs, "command")

    return args


def load_variables(kwargs, aqme_module, create_dat=True):
    """
    Load default and user-defined variables
    """

    # first, load default values and options manually added to the function
    self = set_options(kwargs)
    # this part loads variables from yaml files (if varfile is used)
    txt_yaml = ""
    if self.varfile is not None:
        self, txt_yaml = load_from_yaml(self)
    if aqme_module != "command":

        self.initial_dir = Path(os.getcwd())

        if not isinstance(self.files, list):
            self.w_dir_main = os.path.dirname(self.files)
            check_files = os.path.basename(self.files)
        elif len(self.files) != 0:
            self.w_dir_main = os.path.dirname(self.files[0])
        else:
            self.w_dir_main = os.getcwd()

        if (
            Path(f"{self.w_dir_main}").exists()
            and os.getcwd() not in f"{self.w_dir_main}"
        ):
            self.w_dir_main = Path(f"{os.getcwd()}/{self.w_dir_main}")
        else:
            self.w_dir_main = Path(self.w_dir_main)

        if self.isom_type is not None:
            if (
                Path(f"{self.isom_inputs}").exists()
                and os.getcwd() not in f"{self.isom_inputs}"
            ):
                self.isom_inputs = Path(f"{os.getcwd()}/{self.isom_inputs}")
            else:
                self.isom_inputs = Path(self.isom_inputs)

        error_setup = False

        if not self.w_dir_main.exists():
            txt_yaml += "\nx  The PATH specified as input in the w_dir_main option might be invalid! Using current working directory"
            error_setup = True

        if error_setup:
            self.w_dir_main = Path(os.getcwd())

        if not isinstance(self.files, list):
            if not isinstance(self.files, Mol):
                self.files = glob.glob(f"{self.w_dir_main}/{check_files}")
            else:
                self.files = [self.files]
        # this function is useful when PATH objects are given
        for i,file in enumerate(self.files):
            if not isinstance(file, Mol):
                self.files[i] = str(file)
            
        # start a log file to track the QCORR module
        if create_dat:
            logger_1, logger_2 = "AQME", "data"
            if aqme_module == "qcorr":
                # detects cycle of analysis (0 represents the starting point)
                self.round_num, self.resume_qcorr = check_run(self.w_dir_main)
                logger_1 = "QCORR-run"
                logger_2 = f"{str(self.round_num)}"

            elif aqme_module == "csearch":
                logger_1 = "CSEARCH"

            elif aqme_module == "cmin":
                logger_1 = "CMIN"

            elif aqme_module == "qprep":
                logger_1 = "QPREP"

            elif aqme_module == "qdescp":
                logger_1 = "QDESCP"

            elif aqme_module == "vismol":
                logger_1 = "VISMOL"

            if txt_yaml not in [
                "",
                f"\no  Importing AQME parameters from {self.varfile}",
                "\nx  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code.\n",
            ]:
                self.log = Logger(self.initial_dir / logger_1, logger_2)
                self.log.write(txt_yaml)
                error_setup = True

            if not error_setup:
                if not self.command_line:
                    self.log = Logger(self.initial_dir / logger_1, logger_2)
                else:
                    # prevents errors when using command lines and running to remote directories
                    path_command = Path(f"{os.getcwd()}")
                    self.log = Logger(path_command / logger_1, logger_2)

                self.log.write(f"AQME v {aqme_version} {time_run} \nCitation: {aqme_ref}\n")

                if self.command_line:
                    self.log.write(f"Command line used in AQME: aqme {' '.join([str(elem) for elem in sys.argv[1:]])}\n")

            if error_setup:
                # this is added to avoid path problems in jupyter notebooks
                self.log.finalize()
                os.chdir(self.initial_dir)
                sys.exit()

    return self


def read_file(initial_dir, w_dir, file):
    """
    Reads through a file and retrieves a list with all the lines.
    """

    os.chdir(w_dir)
    outfile = open(file, "r")
    outlines = outfile.readlines()
    outfile.close()
    os.chdir(initial_dir)

    return outlines


def QM_coords(outlines, min_RMS, n_atoms, program, keywords_line):
    """
    Retrieves atom types and coordinates from QM output files
    """

    atom_types, cartesians, range_lines = [], [], []
    per_tab = periodic_table()
    count_RMS = -1

    if program == "gaussian":
        if "nosymm" in keywords_line.lower():
            target_ori = "Input orientation:"
        else:
            target_ori = "Standard orientation:"

        if min_RMS > -1:
            for i, line in enumerate(outlines):
                if line.find(target_ori) > -1:
                    count_RMS += 1
                if count_RMS == min_RMS:
                    range_lines = [i + 5, i + 5 + n_atoms]
                    break
        else:
            for i in reversed(range(len(outlines))):
                if outlines[i].find(target_ori) > -1:
                    range_lines = [i + 5, i + 5 + n_atoms]
                    break
        if len(range_lines) != 0:
            for i in range(range_lines[0], range_lines[1]):
                massno = int(outlines[i].split()[1])
                if massno < len(per_tab):
                    atom_symbol = per_tab[massno]
                else:
                    atom_symbol = "XX"
                atom_types.append(atom_symbol)
                cartesians.append(
                    [
                        float(outlines[i].split()[3]),
                        float(outlines[i].split()[4]),
                        float(outlines[i].split()[5]),
                    ]
                )

    return atom_types, cartesians


def cclib_atoms_coords(cclib_data):
    """
    Function to convert atomic numbers and coordinate arrays from cclib into
    a format compatible with QPREP.
    """

    atom_numbers = cclib_data["atoms"]["elements"]["number"]
    atom_types = []
    per_tab = periodic_table()
    for atom_n in atom_numbers:
        if atom_n < len(per_tab):
            atom_symbol = per_tab[atom_n]
        else:
            atom_symbol = "XX"
        atom_types.append(atom_symbol)

    cartesians_array = cclib_data["atoms"]["coords"]["3d"]
    cartesians = [
        cartesians_array[i : i + 3] for i in range(0, len(cartesians_array), 3)
    ]

    return atom_types, cartesians


def check_run(w_dir):
    """
    Determines the folder where input files are gonna be generated in QCORR.
    """

    if "failed" in w_dir.as_posix():
        resume_qcorr = True
        for folder in w_dir.as_posix().replace("\\", "/").split("/"):
            if "run_" in folder:
                folder_count = int(folder.split("_")[1]) + 1
    else:
        input_folder = w_dir.joinpath("failed/")
        resume_qcorr = False
        folder_count = 1

        if os.path.exists(input_folder):
            dir_list = os.listdir(input_folder)
            for folder in dir_list:
                if folder.find("run_") > -1:
                    folder_count += 1

    return folder_count, resume_qcorr


def read_xyz_charge_mult(file):
    """
    Reads charge and multiplicity from XYZ files. These parameters should be defined
    in the title lines as charge=X and mult=Y (i.e. FILENAME charge=1 mult=1 Eopt -129384.564)
    """

    charge_xyz, mult_xyz = None, None
    # read charge and mult from xyz files
    with open(file, "r") as F:
        lines = F.readlines()
    for line in lines:
        for keyword in line.strip().split():
            if keyword.lower().find("charge") > -1:
                charge_xyz = int(keyword.split("=")[1])
            elif keyword.lower().find("mult") > -1:
                mult_xyz = int(keyword.split("=")[1])
            elif charge_xyz is not None and mult_xyz is not None:
                break

    if charge_xyz is None:
        charge_xyz = 0
    if mult_xyz is None:
        mult_xyz = 1

    return charge_xyz, mult_xyz


def mol_from_sdf_or_mol_or_mol2(input_file, module):
    """
    mol object from SDF, MOL or MOL2 files
    """

    if module in ["qprep","cmin"]:
        # using sanitize=False to avoid reading problems
        mols = Chem.SDMolSupplier(input_file, removeHs=False, sanitize=False)
        return mols

    elif module == "csearch":

        # using sanitize=True in this case, which is recommended for RDKit calculations
        filename = os.path.splitext(input_file)[0]
        extension = os.path.splitext(input_file)[1]

        if extension.lower() == ".pdb":
            input_file = f'{input_file.split(".")[0]}.sdf'
            extension = ".sdf"

        if extension.lower() == ".sdf":
            mols = Chem.SDMolSupplier(input_file, removeHs=False)
        elif extension.lower() == ".mol":
            mols = [Chem.MolFromMolFile(input_file, removeHs=False)]
        elif extension.lower() == ".mol2":
            mols = [Chem.MolFromMol2File(input_file, removeHs=False)]

        IDs, charges, mults = [], [], []

        with open(input_file, "r") as F:
            lines = F.readlines()

        molecule_count = 0
        for i, line in enumerate(lines):
            if line.find(">  <ID>") > -1:
                ID = lines[i + 1].split()[0]
                IDs.append(ID)
            if line.find(">  <Real charge>") > -1:
                charge = lines[i + 1].split()[0]
                charges.append(charge)
            if line.find(">  <Mult>") > -1:
                mult = lines[i + 1].split()[0]
                mults.append(mult)
            if line.find("$$$$") > -1:
                molecule_count += 1
                if molecule_count != len(charges):
                    charges.append(0)
                if molecule_count != len(mults):
                    mults.append(1)

        suppl = []
        for i, mol in enumerate(mols):
            suppl.append(mol)

        if len(IDs) == 0:
            if len(suppl) > 1:
                for i in range(len(suppl)):
                    IDs.append(f"{filename}_{i+1}")
            else:
                IDs.append(filename)

        if len(charges) == 0:
            for i, mol in enumerate(mols):
                charges.append(Chem.GetFormalCharge(mol))
        if len(mults) == 0:
            for i, mol in enumerate(mols):
                NumRadicalElectrons = 0
                for Atom in mol.GetAtoms():
                    NumRadicalElectrons += Atom.GetNumRadicalElectrons()
                TotalElectronicSpin = NumRadicalElectrons / 2
                mult = int((2 * TotalElectronicSpin) + 1)
                mults.append(mult)

        return suppl, charges, mults, IDs
