######################################################.
#          This file stores functions used           #
#                in multiple modules                 #
######################################################.

import os
import re
import subprocess
import sys
import time
import getopt
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

aqme_version = "1.6.0"
time_run = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
aqme_ref = f"AQME v {aqme_version}, Alegre-Requena, J. V.; Sowndarya, S.; Perez-Soto, R.; Alturaifi, T.; Paton, R. AQME: Automated Quantum Mechanical Environments for Researchers and Educators. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2023, DOI: 10.1002/wcms.1663."

RDLogger.DisableLog("rdApp.*")


def run_command(command, outfile, cwd=None):
    """
    Runs the subprocess command and saves the results in an output file (not shown in the terminal)
    """

    output = open(outfile, "w")
    subprocess.run(command, stdout=output, stderr=subprocess.DEVNULL, cwd=cwd)
    output.close()


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
            if os.path.basename(Path(self.varfile)).split('.')[1] in ["yaml", "yml", "txt"]:
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
        txt_yaml = f"\nx  The specified yaml file containing parameters {self.varfile} was not found or the extension is not compatible ('.yaml', '.yml' or '.txt')! Also, make sure that the params file is in the folder where you are running the code."

    return self, txt_yaml


# class for logging
class Logger:
    """
    Class that wraps a file object to abstract the logging.
    """

    # Class Logger to write output to a file
    def __init__(self, filein, append, suffix="dat", verbose=True):
        if verbose:
            self.log = open(f"{filein}_{append}.{suffix}", "w")
        else:
            self.log = ''

    def write(self, message):
        """
        Appends a newline character to the message and writes it into the file.

        Parameters
        ----------
        message : str
           Text to be written in the log file.
        """
        try:
            self.log.write(f"{message}\n")
        except AttributeError:
            pass
        print(f"{message}\n")

    def finalize(self):
        """
        Closes the file
        """
        try:
            self.log.close()
        except AttributeError:
            pass


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
    if os.path.basename(Path(file)).split(".")[1] in ["com", "gjf"]:

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
    if os.path.basename(Path(file)).split(".")[1] == "inp":

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
            atom.SetFormalCharge(0)


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
    return GetBestRMS(mol1, mol2, c1, c2, maxMatches=max_matches_rmsd) # don't use numThreads=0 or -1 as the documentation says, it fails! (due to multithreading in CSEARCH module?)


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
        "heavyonly",
        "cregen",
        "lowest_only",
        "chk",
        "oldchk",
        "nodup_check",
        "dbstep_calc",
        "robert"
    ]
    list_args = [
        "files",
        "gen_atoms",
        "constraints_atoms",
        "constraints_dist",
        "constraints_angle",
        "constraints_dihedral",
        "atom_types",
        "cartesians",
        "nmr_atoms",
        "nmr_slope",
        "nmr_intercept",
        "qdescp_atoms",
        "geom"
    ]
    int_args = [
        "opt_steps",
        "opt_steps_rdkit",
        "auto_sample",
        "seed",
        "max_matches_rmsd",
        "nsteps_fullmonte",
        "nrot_fullmonte",
        "nprocs",
        "crest_nrun",
    ]
    float_args = [
        "ewin_cmin",
        "ewin_csearch",
        "opt_fmax",
        "degree",
        "rms_threshold",
        "energy_threshold",
        "initial_energy_threshold",
        "max_mol_wt",
        "ewin_sample_fullmonte",
        "ewin_fullmonte",
        "dup_threshold",
        "ro_threshold",
        "amplitude_ifreq",
        "ifreq_cutoff",
        "s2_threshold",
        "vdwfrac",
        "covfrac",
        "bond_thres",
        "angle_thres",
        "dihedral_thres",
        "crest_force",
        "qdescp_temp",
        "qdescp_acc",
        "dbstep_r",
        "crest_nclust",
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

        if arg_name in ("h", "help"):
            print(f"o  AQME v {aqme_version} is installed correctly! For more information about the available options, see the documentation in https://github.com/jvalegre/aqme")
            sys.exit()
        else:
            # this "if" allows to use * to select multiple files in multiple OS
            if arg_name.lower() == 'files':
                value = get_files(value)
                kwargs[arg_name] = value
            else:
                # this converts the string parameters to lists
                if arg_name in bool_args:
                    value = True                    
                elif arg_name.lower() in list_args:
                    value = format_lists(value)
                elif arg_name.lower() in int_args:
                    value = int(value)
                elif arg_name.lower() in float_args:
                    value = float(value)
                elif value == "None":
                    value = None
                elif value == "False":
                    value = False
                elif value == "True":
                    value = True

                kwargs[arg_name] = value

    # Second, load all the default variables as an "add_option" object
    args = load_variables(kwargs, "command")
    
    return args


def format_lists(value):
    '''
    Transforms strings into a list
    '''

    if not isinstance(value, list):
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            # this line fixes issues when using "[X]" or ["X"] instead of "['X']" when using lists
            value = value.replace('[',']').replace(',',']').replace("'",']').split(']')
            while('' in value):
                value.remove('')
    return value


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

        # get PATH for the files option
        self.files = get_files(self.files)

        if not isinstance(self.files, list):
            self.w_dir_main = os.path.dirname(self.files)
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
            txt_yaml += "\nx  The PATH specified as input or files might be invalid!"
            error_setup = True

        if error_setup:
            self.w_dir_main = Path(os.getcwd())
            
        # start a log file to track the AQME modules
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

            if txt_yaml not in [
                "",
                f"\no  Importing AQME parameters from {self.varfile}",
                "\nx  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code.\n",
            ]:
                self.log = Logger(self.initial_dir / logger_1, logger_2, verbose=self.verbose)
                self.log.write(txt_yaml)
                error_setup = True

            if not error_setup:
                if not self.command_line:
                    self.log = Logger(self.initial_dir / logger_1, logger_2, verbose=self.verbose)
                else:
                    # prevents errors when using command lines and running to remote directories
                    path_command = Path(f"{os.getcwd()}")
                    self.log = Logger(path_command / logger_1, logger_2, verbose=self.verbose)

                self.log.write(f"AQME v {aqme_version} {time_run} \nCitation: {aqme_ref}\n")

                if self.command_line:
                    cmd_print = ''
                    cmd_args = sys.argv[1:]
                    for i,elem in enumerate(cmd_args):
                        if elem[0] in ['"',"'"]:
                            elem = elem[1:]
                        if elem[-1] in ['"',"'"]:
                            elem = elem[:-1]
                        if elem != '-h' and elem.split('--')[-1] not in var_dict:
                            # parse single elements of the list as strings (otherwise the commands cannot be reproduced)
                            if cmd_args[i-1] == '--qdescp_atoms':
                                elem = elem[1:-1]
                                elem = elem.replace(', ',',').replace(' ,',',')
                                new_elem = []
                                for smarts_strings in elem.split(','):
                                    new_elem.append(f'{smarts_strings}'.replace("'",''))
                                elem = f'{new_elem}'.replace(" ","")
                            if cmd_args[i-1].split('--')[-1] in var_dict: # check if the previous word is an arg
                                cmd_print += f'"{elem}'
                            if i == len(cmd_args)-1 or cmd_args[i+1].split('--')[-1] in var_dict: # check if the next word is an arg, or last word in command
                                cmd_print += f'"'
                        else:
                            cmd_print += f'{elem}'
                        if i != len(cmd_args)-1:
                            cmd_print += ' '

                    self.log.write(f"Command line used in AQME: python -m aqme {cmd_print}\n")

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

    failed_run_pattern = re.compile(
        r"(^.*)[/\\]failed[/\\]run_(?P<folder_count>\d+)[/\\]"
    )

    run_pattern = re.compile(r"(^.*)[/\\]run_(?P<folder_count>\d+)[/\\]")
    resume_qcorr = False
    folder_count = 1

    failed_run_match = failed_run_pattern.match(w_dir.as_posix())
    run_match = run_pattern.match(w_dir.as_posix())

    if failed_run_match:
        folder_count = int(failed_run_match.group("folder_count")) + 1
        resume_qcorr = True
        return folder_count, resume_qcorr

    if run_match:
        folder_count = int(run_match.group("folder_count"))
        resume_qcorr = False
        return folder_count, resume_qcorr

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
            if keyword.lower().find("charge=") > -1:
                charge_xyz = int(keyword.split("=")[1])
            elif keyword.lower().find("mult=") > -1:
                mult_xyz = int(keyword.split("=")[1])
            elif charge_xyz is not None and mult_xyz is not None:
                break

    if charge_xyz is None:
        charge_xyz = 0
    if mult_xyz is None:
        mult_xyz = 1

    return charge_xyz, mult_xyz


def mol_from_sdf_or_mol_or_mol2(input_file, module, args, low_check=None):

    """
    mol object from SDF, MOL or MOL2 files
    """
    if module in ["qprep","cmin"]:
        # using sanitize=False to avoid reading problems
        mols = Chem.SDMolSupplier(input_file, removeHs=False, sanitize=False)
        if low_check=='lowest_only':
            return [mols[0]]
        elif isinstance(low_check, int):
            check_n = min(len(mols),low_check)
            n_mols = []
            for i in range(check_n):
                n_mols.append(mols[i])
            return n_mols
        elif isinstance(low_check, float):
            n_mols = []
            for i in range(len(mols)):
                if abs(float(mols[i].GetProp('Energy')) - float(mols[0].GetProp('Energy'))) < low_check: # kcal/mol
                    n_mols.append(mols[i])
            return n_mols
        else:
            return mols

    elif module == "csearch":

        # using sanitize=True in this case, which is recommended for RDKit calculations
        filename = os.path.basename(Path(input_file)).split('.')[0]
        extension = os.path.basename(Path(input_file)).split('.')[1]

        if extension.lower() == "pdb":
            input_file = f'{os.path.dirname(Path(input_file))}/{filename}.sdf'
            extension = "sdf"

        if extension.lower() == "sdf":
            mols = Chem.SDMolSupplier(input_file, removeHs=False)
        elif extension.lower() == "mol":
            mols = [Chem.MolFromMolFile(input_file, removeHs=False)]
        elif extension.lower() == "mol2":
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
            path_file = f'{os.path.dirname(input_file)}/{os.path.basename(input_file).split(".")[0]}'
            if len(suppl) > 1:
                for i in range(len(suppl)):
                    IDs.append(f"{path_file}_{i+1}")
            else:
                IDs.append(f'{path_file}')

        if args.charge is not None:
            for i, mol in enumerate(mols):
                charges.append(args.charge)
        elif len(charges) == 0:
            for i, mol in enumerate(mols):
                charges.append(Chem.GetFormalCharge(mol))
        if args.mult is not None:
            mults.append(args.mult)
        elif len(mults) == 0:
            for i, mol in enumerate(mols):
                NumRadicalElectrons = 0
                for Atom in mol.GetAtoms():
                    NumRadicalElectrons += Atom.GetNumRadicalElectrons()
                TotalElectronicSpin = NumRadicalElectrons / 2
                mult = int((2 * TotalElectronicSpin) + 1)
                mults.append(mult)

        return suppl, charges, mults, IDs


def add_prefix_suffix(name, args):
    if args.prefix != "":
        name = f"{args.prefix}_{name}"
    if args.suffix != "":
        name += f'_{args.suffix}'

    return name


def check_files(self,module):
    if module.lower() in ['cmin','qdescp','qprep']:
        format_file = "*.sdf"
    elif module.lower() in ['qcorr']:
        format_file = "*.log"
    no_file_found = False
    if len(self.args.files) == 0:
        self.args.log.write(f'\nx  No files were found! In the "files" option, make sure that 1) the PATH to the files is correct, 2) the PATH doesn\'t start with "/", and 3) you use quotation marks if you are using * (i.e. --files "{format_file}")')
        no_file_found = True
    else:
        for file in self.args.files:
            if not os.path.exists(file):
                no_file_found = True
                self.args.log.write(f'\nx  File {file} was not found! In the "files" option, make sure that 1) the PATH to the files is correct and 2) the PATH doesn\'t start with "/".')
                break
    if no_file_found:
        self.args.log.finalize()
        sys.exit()


def check_xtb(self):
    try:
        subprocess.run(
            ["xtb", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        self.args.log.write("x  xTB is not installed (CSEARCH-CREST and CMIN-xTB cannot be used)! You can install the program with 'conda install -c conda-forge xtb'")
        self.args.log.finalize()
        sys.exit()


def check_crest(self):
    try:
        subprocess.run(
            ["crest", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        self.args.log.write("x  CREST is not installed (CSEARCH-CREST cannot be used)! You can install the program with 'conda install -c conda-forge crest'")
        self.args.log.finalize()
        sys.exit()
 

def get_files(value):
    if not isinstance(value, list):
        value = value.replace('[',']').replace(',',']').split(']')
        while('' in value):
            value.remove('')
    new_value = []
    for val in value:
        if not isinstance(val, Path):
            if (
            Path(f"{val}").exists()
            and os.getcwd() not in f"{val}"
            ):
                new_value.append(f"{os.getcwd()}/{val}")
            elif '*' in val:
                if os.getcwd() not in f"{val}":
                    list_of_val = glob.glob(f"{os.getcwd()}/{val}")
                else:
                    list_of_val = glob.glob(val)
                for ele in list_of_val:
                    new_value.append(ele)
            else:
                new_value.append(val)
        else:
            new_value.append(val.as_posix())
    return new_value
