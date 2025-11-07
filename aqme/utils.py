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
from rdkit.Chem import AllChem as Chem
from aqme.argument_parser import set_options, var_dict
from rdkit import RDLogger

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15

obabel_version = "3.1.1" # this MUST match the meta.yaml
aqme_version = "2.0.0"
time_run = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
aqme_ref = f"AQME v {aqme_version}, Alegre-Requena, J. V.; Sowndarya, S.; Perez-Soto, R.; Alturaifi, T.; Paton, R. AQME: Automated Quantum Mechanical Environments for Researchers and Educators. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2023, 13, e1663 (DOI: 10.1002/wcms.1663)."
xtb_version = '6.7.1'
crest_version = '2.12'

RDLogger.DisableLog("rdApp.*")


def run_command(command, outfile, cwd=None, env=None):
    """Run subprocess command and save output to file.
    
    Executes a command in a subprocess and redirects stdout to a file.
    Stderr is suppressed. The command runs without showing terminal output.
    
    Args:
        command (list): Command and arguments to execute
        outfile (str): Path to output file for stdout
        cwd (str or Path, optional): Working directory for command execution
        env (dict, optional): Environment variables for subprocess
    """
    with open(outfile, "w") as output:
        subprocess.run(
            command, 
            stdout=output, 
            stderr=subprocess.DEVNULL, 
            cwd=cwd, 
            env=env
        )


def periodic_table():
    """Generate periodic table as a list of element symbols.
    
    Creates a list where the index corresponds to atomic number.
    Index 0 is empty string, index 1 is 'H', index 2 is 'He', etc.
    
    Returns:
        list: Element symbols indexed by atomic number (0-118)
    """
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


def load_from_yaml(self):
    """Load calculation parameters from YAML file.
    
    Updates self attributes with values from YAML file if varfile is specified.
    Handles file validation and error reporting.
    
    Args:
        self: Arguments object with varfile attribute
    
    Returns:
        tuple: (self, txt_yaml) - updated args and status message
    """
    txt_yaml = f"\no  Importing AQME parameters from {self.varfile}"
    error_yaml = False
    
    try:
        if os.path.exists(self.varfile):
            if os.path.basename(Path(self.varfile)).split('.')[-1] in ["yaml", "yml", "txt"]:
                with open(self.varfile, "r") as file:
                    try:
                        param_list = yaml.load(file, Loader=yaml.SafeLoader)
                    except yaml.scanner.ScannerError:
                        txt_yaml = (
                            f'\nx  Error while reading {self.varfile}. '
                            f'Edit the yaml file and try again '
                            f'(i.e. use ":" instead of "=" to specify variables)'
                        )
                        error_yaml = True
        
        # Update attributes from YAML
        if not error_yaml:
            for param in param_list:
                if hasattr(self, param):
                    if getattr(self, param) != param_list[param]:
                        setattr(self, param, param_list[param])
    
    except UnboundLocalError:
        txt_yaml = (
            f"\nx  The specified yaml file containing parameters {self.varfile} "
            f"was not found or the extension is not compatible ('.yaml', '.yml' or '.txt')! "
            f"Also, make sure that the params file is in the folder where you are running the code."
        )
    
    return self, txt_yaml


class Logger:
    """Wrapper class for file-based logging with console output.
    
    Provides unified logging to both file and console. Messages are
    automatically appended with newlines and printed to stdout.
    
    Attributes:
        log: File object for writing log messages (or empty string if verbose=False)
    """

    def __init__(self, filein, append, suffix="dat", verbose=True):
        """Initialize logger with output file.
        
        Args:
            filein (str or Path): Base path for log file
            append (str): String to append to filename
            suffix (str): File extension (default: "dat")
            verbose (bool): If True, create log file; if False, disable file logging
        """
        if verbose:
            self.log = open(f"{filein}_{append}.{suffix}", "w")
        else:
            self.log = ''

    def write(self, message):
        """Write message to log file and print to console.
        
        Appends newline to message automatically. Prints to stdout
        regardless of file logging status.
        
        Args:
            message (str): Text to log
        """
        try:
            self.log.write(f"{message}\n")
        except AttributeError:
            pass
        print(f"{message}\n")

    def finalize(self):
        """Close the log file.
        
        Safely closes file handle. Does nothing if file logging is disabled.
        """
        try:
            self.log.close()
        except AttributeError:
            pass


def move_file(destination, source, file):
    """Move file from source to destination folder.
    
    Creates destination directory if it doesn't exist. Replaces file
    if it already exists at destination.
    
    Args:
        destination (Path): Destination folder path
        source (Path): Source folder path
        file (str): Filename (including extension)
    """
    destination.mkdir(exist_ok=True, parents=True)
    filepath = source / file
    try:
        filepath.rename(destination / file)
    except FileExistsError:
        filepath.replace(destination / file)


def set_destination(self, module):
    """Set up the destination folder for output files.
    
    Uses destination argument if provided, otherwise creates folder
    in initial directory named after the module.
    
    Args:
        self: Arguments object with destination settings
        module (str): Module name (CSEARCH, CMIN, etc.)
    
    Returns:
        pathlib.Path: Destination folder path
    """
    if self.args.destination is None:
        destination = self.args.initial_dir.joinpath(module)
    elif self.args.initial_dir.joinpath(self.args.destination).exists():
        destination = Path(self.args.initial_dir.joinpath(self.args.destination))
    else:
        destination = Path(self.args.destination)
    
    return destination


def _parse_gaussian_input(input_lines):
    """Parse Gaussian .com/.gjf file for coordinates and charge/mult.
    
    Args:
        input_lines (list): Lines from Gaussian input file
    
    Returns:
        tuple: (atoms_and_coords, charge, mult)
    """
    _iter = input_lines.__iter__()
    line = ""
    
    # Find the command line
    while "#" not in line:
        line = next(_iter)
    
    # Skip keywords on multiple lines
    while len(line.split()) > 0:
        line = next(_iter)
    
    # Pass the title lines
    _ = next(_iter)
    while line:
        line = next(_iter).strip()
    
    # Read charge and multiplicity
    charge, mult = next(_iter).strip().split()
    
    # Store the atom types and coordinates until next empty line
    atoms_and_coords = []
    line = next(_iter).strip()
    while line:
        atoms_and_coords.append(line.strip())
        line = next(_iter).strip()
    
    return atoms_and_coords, charge, mult


def _parse_orca_input(input_lines):
    """Parse ORCA .inp file for coordinates and charge/mult.
    
    Args:
        input_lines (list): Lines from ORCA input file
    
    Returns:
        tuple: (atoms_and_coords, charge, mult)
    """
    _iter = input_lines.__iter__()
    line = ""
    
    # Find the line with charge and multiplicity
    while "* xyz" not in line and "* int" not in line:
        line = next(_iter)
    
    # Read charge and multiplicity
    charge = line.strip().split()[-2]
    mult = line.strip().split()[-1]
    
    # Store the coordinates until next *
    atoms_and_coords = []
    line = next(_iter).strip()
    while len(line.split()) > 1:
        atoms_and_coords.append(line.strip())
        line = next(_iter).strip()
    
    return atoms_and_coords, charge, mult


def get_info_input(file):
    """Extract coordinates and charge from quantum chemistry input files.
    
    Supports Gaussian (.com, .gjf) and ORCA (.inp) file formats.
    
    Parameters
    ----------
    file : str or pathlib.Path
        Path to a valid .com, .gjf, or .inp file

    Returns
    -------
    coordinates : list
        List of strings (without \\n) containing xyz coordinates
    charge : str
        Total charge of the system
    mult : str
        Spin multiplicity of the system
    """
    with open(file, "r") as input_file:
        input_lines = input_file.readlines()
    
    file_extension = os.path.basename(Path(file)).split(".")[-1]
    
    if file_extension in ["com", "gjf"]:
        return _parse_gaussian_input(input_lines)
    elif file_extension == "inp":
        return _parse_orca_input(input_lines)
    else:
        raise ValueError(f"Unsupported file extension: {file_extension}")

def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_rmsd):
    """Calculate RMSD between two molecule conformations.
    
    Computes best RMSD by aligning mol1 to mol2. Note: mol1 is modified
    (left in aligned state) unless heavy=True.
    
    Args:
        mol1 (rdkit.Chem.Mol): Probe molecule to align
        mol2 (rdkit.Chem.Mol): Target molecule for alignment
        c1 (int): Conformer ID of mol1 (-1 for default)
        c2 (int): Conformer ID of mol2 (-1 for default)
        heavy (bool): If True, ignore hydrogens in RMSD calculation
        max_matches_rmsd (int): Maximum substructure matches to consider
    
    Returns:
        float: Best RMSD value found (Angstroms)
    
    Note:
        numThreads must be 1 to not overload the multithreading
    """
    if heavy:
        mol1 = RemoveHs(mol1)
        mol2 = RemoveHs(mol2)
    
    return GetBestRMS(
        mol1, mol2, c1, c2, 
        maxMatches=max_matches_rmsd, 
        numThreads=1
    )


def _get_argument_categories():
    """Define argument categories for command line parsing.
    
    Returns:
        tuple: (bool_args, list_args, int_args, float_args)
    """
    bool_args = [
        "csearch", "cmin", "qprep", "qcorr", "qdescp",
        "heavyonly", "single_system", "cregen", "lowest_only",
        "chk", "oldchk", "nodup_check", "robert", "debug", "pytest_testing"
    ]
    
    list_args = [
        "files", "gen_atoms", "constraints_atoms", "constraints_dist",
        "constraints_angle", "constraints_dihedral", "atom_types", "cartesians",
        "nmr_atoms", "nmr_slope", "nmr_intercept", "qdescp_atoms", "geom"
    ]
    
    int_args = [
        "opt_steps", "opt_steps_rdkit", "seed", "max_matches_rmsd",
        "nprocs", "crest_runs", "sample"
    ]
    
    float_args = [
        "ewin_cmin", "ewin_csearch", "opt_fmax", "rms_threshold",
        "energy_threshold", "initial_energy_threshold", "max_mol_wt",
        "dup_threshold", "ro_threshold", "amplitude_ifreq", "ifreq_cutoff",
        "s2_threshold", "vdwfrac", "covfrac", "bond_thres", "angle_thres",
        "dihedral_thres", "crest_force", "qdescp_temp", "qdescp_acc",
        "dbstep_r", "crest_nclust", "vbur_radius"
    ]
    
    return bool_args, list_args, int_args, float_args


def _build_available_args(bool_args):
    """Build list of available command line arguments.
    
    Args:
        bool_args (list): List of boolean argument names
    
    Returns:
        list: Available argument names with appropriate formatting
    """
    available_args = ["help"]
    
    for arg in var_dict:
        if arg in bool_args:
            available_args.append(f"{arg}")
        else:
            available_args.append(f"{arg} =")
    
    return available_args


def _parse_argument_value(arg_name, value, bool_args, list_args, int_args, float_args):
    """Convert argument value to appropriate type.
    
    Args:
        arg_name (str): Argument name
        value (str): Raw value from command line
        bool_args (list): Boolean argument names
        list_args (list): List argument names
        int_args (list): Integer argument names
        float_args (list): Float argument names
    
    Returns:
        Converted value in appropriate type
    """
    if arg_name in bool_args:
        return True
    elif arg_name.lower() in list_args:
        return format_lists(value, arg_name)
    elif arg_name.lower() in int_args:
        return int(value)
    elif arg_name.lower() in float_args:
        return float(value)
    elif value == "None":
        return None
    elif value == "False":
        return False
    elif value == "True":
        return True
    else:
        return value


def command_line_args():
    """Load command line arguments for AQME.
    
    Parses arguments from sys.argv, converts them to appropriate types,
    and loads default values for unspecified parameters.
    
    Returns:
        argparse.Namespace: Arguments object with all AQME parameters
    """
    # Get argument categories
    bool_args, list_args, int_args, float_args = _get_argument_categories()
    
    # Build available arguments list
    available_args = _build_available_args(bool_args)
    
    # Parse command line options
    kwargs = {}
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "h", available_args)
    except getopt.GetoptError as err:
        print(err)
        sys.exit()
    
    # Process each argument
    for arg, value in opts:
        # Extract argument name
        if arg.find("--") > -1:
            arg_name = arg.split("--")[1].strip()
        elif arg.find("-") > -1:
            arg_name = arg.split("-")[1].strip()
        
        # Handle help
        if arg_name in ("h", "help"):
            print(
                f"o  AQME v {aqme_version} is installed correctly! "
                f"For more information about the available options, see the documentation in "
                f"https://github.com/jvalegre/aqme"
            )
            sys.exit()
        
        # Process file arguments specially to handle wildcards
        elif arg_name.lower() == 'files':
            kwargs[arg_name] = get_files(value)
        
        # Convert value to appropriate type
        else:
            kwargs[arg_name] = _parse_argument_value(
                arg_name, value, bool_args, list_args, int_args, float_args
            )
    
    # Load all default variables
    args = load_variables(kwargs, "command")
    
    return args


def format_lists(value, arg_name):
    """Transform string arguments into lists.
    
    Handles various input formats and removes extra spaces.
    Special handling for qdescp_atoms which uses different parsing.
    
    Args:
        value (str or list): Input value to convert
        arg_name (str): Argument name for special handling
    
    Returns:
        list: Formatted list of values
    """
    if arg_name.lower() != 'qdescp_atoms':
        if not isinstance(value, list):
            try:
                value = ast.literal_eval(value)
            except (SyntaxError, ValueError):
                # Fix issues with "[X]" or ["X"] instead of "['X']"
                value = value.replace('[', ']').replace(',', ']').replace("'", ']').split(']')
                while '' in value:
                    value.remove('')
    else:
        value = value[1:-1].split(',')
    
    # Remove extra spaces
    value = [str(ele).strip() if isinstance(ele, str) else ele for ele in value]
    
    return value


def _setup_working_directory(self):
    """Setup working directory and file paths.
    
    Args:
        self: Arguments object
    
    Returns:
        tuple: (self, error_setup) - updated args and error flag
    """
    error_setup = False
    
    # Get PATH for the files option
    self.files = get_files(self.files)
    
    # Determine working directory
    if not isinstance(self.files, list):
        self.w_dir_main = os.path.dirname(self.files)
    elif len(self.files) != 0:
        self.w_dir_main = os.path.dirname(self.files[0])
    else:
        self.w_dir_main = os.getcwd()
    
    # Resolve working directory path
    if (Path(f"{self.w_dir_main}").exists() 
        and os.getcwd() not in f"{self.w_dir_main}"):
        self.w_dir_main = Path(f"{os.getcwd()}/{self.w_dir_main}")
    else:
        self.w_dir_main = Path(self.w_dir_main)
    
    # Setup isomerization input path if needed
    if self.isom_type is not None:
        if (Path(f"{self.isom_inputs}").exists() 
            and os.getcwd() not in f"{self.isom_inputs}"):
            self.isom_inputs = Path(f"{os.getcwd()}/{self.isom_inputs}")
        else:
            self.isom_inputs = Path(self.isom_inputs)
    
    # Validate working directory exists
    if not self.w_dir_main.exists():
        error_setup = True
    
    return self, error_setup


def _determine_logger_names(aqme_module, self):
    """Determine logger file names based on module.
    
    Args:
        aqme_module (str): Name of AQME module
        self: Arguments object
    
    Returns:
        tuple: (logger_1, logger_2) - logger name and suffix
    """
    logger_1, logger_2 = "AQME", "data"
    
    if aqme_module == "qcorr":
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
    
    return logger_1, logger_2


def _format_command_line_for_log(cmd_args):
    """Format command line arguments for logging.
    
    Args:
        cmd_args (list): Command line arguments from sys.argv
    
    Returns:
        str: Formatted command line string
    """
    cmd_print = ''
    
    for i, elem in enumerate(cmd_args):
        # Remove quotes from element
        if elem[0] in ['"', "'"]:
            elem = elem[1:]
        if elem[-1] in ['"', "'"]:
            elem = elem[:-1]
        
        # Skip help and non-argument elements
        if elem != '-h' and elem.split('--')[-1] not in var_dict:
            # Special handling for qdescp_atoms
            if cmd_args[i-1] == '--qdescp_atoms':
                elem = elem[1:-1]
                elem = elem.replace(', ', ',').replace(' ,', ',')
                new_elem = []
                for smarts_strings in elem.split(','):
                    new_elem.append(f'{smarts_strings}'.replace("'", ''))
                elem = f'{new_elem}'.replace(" ", "")
            
            # Check if previous word is an argument
            if cmd_args[i-1].split('--')[-1] in var_dict:
                cmd_print += f'"{elem}'
            
            # Check if next word is an arg, or last word in command
            if i == len(cmd_args)-1 or cmd_args[i+1].split('--')[-1] in var_dict:
                cmd_print += f'"'
        else:
            cmd_print += f'{elem}'
        
        if i != len(cmd_args)-1:
            cmd_print += ' '
    
    return cmd_print


def _create_logger(self, aqme_module, logger_1, logger_2, txt_yaml, error_setup):
    """Create and initialize logger for the module.
    
    Args:
        self: Arguments object
        aqme_module (str): Module name
        logger_1 (str): Logger base name
        logger_2 (str): Logger suffix
        txt_yaml (str): YAML loading messages
        error_setup (bool): Whether setup errors occurred
    
    Returns:
        tuple: (self, error_setup) - updated args and error flag
    """
    yaml_error_messages = [
        "",
        f"\no  Importing AQME parameters from {self.varfile}",
        "\nx  The specified yaml file containing parameters was not found! "
        "Make sure that the valid params file is in the folder where you are running the code.\n",
    ]
    
    # Handle YAML errors
    if txt_yaml not in yaml_error_messages:
        self.log = Logger(self.initial_dir / logger_1, logger_2, verbose=self.verbose)
        self.log.write(txt_yaml)
        error_setup = True
        return self, error_setup
    
    # Create logger based on mode
    if not self.command_line:
        self.log = Logger(self.initial_dir / logger_1, logger_2, verbose=self.verbose)
    else:
        # Prevents errors when using command lines and running to remote directories
        path_command = Path(f"{os.getcwd()}")
        self.log = Logger(path_command / logger_1, logger_2, verbose=self.verbose)
    
    # Write header
    self.log.write(f"AQME v {aqme_version} {time_run} \nCitation: {aqme_ref}\n")
    
    # Log command line if used
    if self.command_line:
        cmd_print = _format_command_line_for_log(sys.argv[1:])
        self.log.write(f"Command line used in AQME: python -m aqme {cmd_print}\n")
    else:
        # Format qdescp_atoms for non-command-line usage
        if isinstance(self.qdescp_atoms, list):
            self.qdescp_atoms = [str(pattern) for pattern in self.qdescp_atoms]
    
    return self, error_setup


def load_variables(kwargs, aqme_module, create_dat=True):
    """Load default and user-defined variables for AQME modules.
    
    Combines default values, user-specified values, and YAML file parameters.
    Sets up working directories and creates log files.
    
    Args:
        kwargs (dict): User-defined arguments
        aqme_module (str): Name of AQME module (csearch, cmin, qprep, etc.)
        create_dat (bool): Whether to create log file
    
    Returns:
        argparse.Namespace: Complete arguments object with all parameters
    """
    # Load default values and manually added options
    self = set_options(kwargs)
    
    # Load variables from YAML file if specified
    txt_yaml = ""
    if self.varfile is not None:
        self, txt_yaml = load_from_yaml(self)
    
    # Skip directory setup for command mode
    if aqme_module == "command":
        return self
    
    # Setup initial directory
    self.initial_dir = Path(os.getcwd())
    
    # Setup working directory and validate paths
    self, error_setup = _setup_working_directory(self)
    
    if error_setup:
        txt_yaml += "\nx  The PATH specified as input or files might be invalid!"
        self.w_dir_main = Path(os.getcwd())
    
    # Create log file if requested
    if create_dat:
        logger_1, logger_2 = _determine_logger_names(aqme_module, self)
        self, error_setup = _create_logger(
            self, aqme_module, logger_1, logger_2, txt_yaml, error_setup
        )
        
        # Exit if setup errors occurred
        if error_setup:
            self.log.finalize()
            os.chdir(self.initial_dir)
            sys.exit()
    
    return self


def read_file(initial_dir, w_dir, file):
    """Read file and return all lines.
    
    Changes to working directory, reads file, then returns to initial directory.
    
    Args:
        initial_dir (Path): Directory to return to after reading
        w_dir (Path): Directory containing the file
        file (str): Filename to read
    
    Returns:
        list: Lines from the file (including newline characters)
    """
    os.chdir(w_dir)
    with open(file, "r") as outfile:
        outlines = outfile.readlines()
    os.chdir(initial_dir)
    
    return outlines


def cclib_atoms_coords(cclib_data, geom):
    """Convert cclib atomic data to QPREP-compatible format.
    
    Extracts atom types and coordinates from cclib parsed data.
    
    Args:
        cclib_data (dict): Parsed data from cclib with 'atomnos' and 'atomcoords'
        geom (int): Geometry index to extract (-1 for last geometry)
    
    Returns:
        tuple: (atom_types, cartesians)
            - atom_types (list): Element symbols for each atom
            - cartesians (array): 3D coordinates for specified geometry
    """
    atom_numbers = cclib_data["atomnos"]
    atom_types = []
    per_tab = periodic_table()
    
    for atom_n in atom_numbers:
        if atom_n < len(per_tab):
            atom_symbol = per_tab[atom_n]
        else:
            atom_symbol = "XX"  # Unknown element
        atom_types.append(atom_symbol)
    
    cartesians_array = cclib_data["atomcoords"]
    cartesians = cartesians_array[geom]
    
    return atom_types, cartesians


def check_run(w_dir):
    """Determine QCORR run folder and resume status.
    
    Parses directory path to detect if this is a continuation of a
    failed run or a new run cycle. Used for QCORR workflow management.
    
    Args:
        w_dir (Path): Working directory path to analyze
    
    Returns:
        tuple: (folder_count, resume_qcorr)
            - folder_count (int): Run number (1 for new, N+1 for failed run N)
            - resume_qcorr (bool): True if resuming from failed run
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
    """Extract charge and multiplicity from XYZ file title line.
    
    Searches for charge=X and mult=Y keywords in XYZ title lines.
    Example: "FILENAME charge=1 mult=1 Eopt -129384.564"
    
    Args:
        file (str or Path): Path to XYZ file
    
    Returns:
        tuple: (charge_xyz, mult_xyz)
            - charge_xyz (int): Molecular charge (default: 0)
            - mult_xyz (int): Spin multiplicity (default: 1)
    """
    charge_xyz, mult_xyz = None, None
    
    with open(file, "r") as F:
        lines = F.readlines()
    
    for line in lines:
        for keyword in line.strip().split():
            if keyword.lower().find("charge=") > -1:
                charge_xyz = int(keyword.split("=")[1])
            elif keyword.lower().find("mult=") > -1:
                mult_xyz = int(keyword.split("=")[1])
            
            # Exit early if both found
            if charge_xyz is not None and mult_xyz is not None:
                break
    
    # Set defaults if not found
    if charge_xyz is None:
        charge_xyz = 0
    if mult_xyz is None:
        mult_xyz = 1
    
    return charge_xyz, mult_xyz


def _filter_mols_by_criteria(mols, low_check):
    """Filter molecules based on selection criteria.
    
    Args:
        mols: RDKit molecule supplier
        low_check: Selection criterion (None, 'lowest_only', int, or float)
    
    Returns:
        list: Filtered molecule list
    """
    if low_check == 'lowest_only':
        return [mols[0]]
    
    elif isinstance(low_check, int):
        # Return first N conformers
        check_n = min(len(mols), low_check)
        return [mols[i] for i in range(check_n)]
    
    elif isinstance(low_check, float):
        # Return conformers within energy window (kcal/mol)
        n_mols = []
        for i in range(len(mols)):
            energy_diff = abs(
                float(mols[i].GetProp('Energy')) - 
                float(mols[0].GetProp('Energy'))
            )
            if energy_diff < low_check:
                n_mols.append(mols[i])
        return n_mols
    
    else:
        return mols


def _load_mols_for_qprep_cmin(input_file, low_check):
    """Load molecules for QPREP or CMIN modules.
    
    Args:
        input_file (str): Path to SDF file
        low_check: Selection criterion for filtering
    
    Returns:
        list: Molecule objects
    """
    # Use sanitize=False to avoid reading problems
    mols = Chem.SDMolSupplier(input_file, removeHs=False, sanitize=False)
    
    # Transform invalid SDF files created with GaussView
    if None in mols:
        mols = load_sdf(input_file)
    
    return _filter_mols_by_criteria(mols, low_check)


def _extract_sdf_properties(input_file):
    """Extract ID, charge, and multiplicity from SDF file.
    
    Args:
        input_file (str): Path to SDF file
    
    Returns:
        tuple: (IDs, charges, mults) - lists of properties
    """
    IDs, charges, mults = [], [], []
    
    with open(input_file, "r") as F:
        lines = F.readlines()
    
    for i, line in enumerate(lines):
        if line.find(">  <ID>") > -1:
            ID = ' '.join(lines[i + 1].split()[:-1])
            IDs.append(ID)
        if line.find(">  <Real charge>") > -1:
            charge = lines[i + 1].split()[0]
            charges.append(charge)
        if line.find(">  <Mult>") > -1:
            mult = lines[i + 1].split()[0]
            mults.append(mult)
    
    return IDs, charges, mults


def _generate_default_ids(input_file, num_mols):
    """Generate default IDs when not found in file.
    
    Args:
        input_file (str): Path to input file
        num_mols (int): Number of molecules
    
    Returns:
        list: Generated IDs
    """
    path_file = f'{os.path.dirname(input_file)}/{".".join(os.path.basename(input_file).split(".")[:-1])}'
    
    if num_mols > 1:
        return [f"{path_file}_{i+1}" for i in range(num_mols)]
    else:
        return [f'{path_file}']


def _calculate_charges_from_mols(mols, args):
    """Calculate charges from molecule objects.
    
    Args:
        mols: List of RDKit molecules
        args: Arguments object with charge setting
    
    Returns:
        list: Charges for each molecule
    """
    if args.charge is not None:
        return [args.charge for _ in mols]
    else:
        return [Chem.GetFormalCharge(mol) for mol in mols]


def _calculate_mults_from_mols(mols, args):
    """Calculate multiplicities from molecule objects.
    
    Args:
        mols: List of RDKit molecules
        args: Arguments object with mult setting
    
    Returns:
        list: Multiplicities for each molecule
    """
    if args.mult is not None:
        return [args.mult for _ in mols]
    
    mults = []
    for mol in mols:
        NumRadicalElectrons = sum(
            atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()
        )
        TotalElectronicSpin = NumRadicalElectrons / 2
        mult = int((2 * TotalElectronicSpin) + 1)
        mults.append(mult)
    
    return mults


def _load_mols_for_csearch(input_file, args, keep_xyz):
    """Load molecules for CSEARCH module.
    
    Args:
        input_file (str): Path to input file
        args: Arguments object
        keep_xyz (bool): Keep XYZ input
    
    Returns:
        tuple: (suppl, charges, mults, IDs)
    """
    # Use sanitize=True for RDKit calculations
    filename = '.'.join(os.path.basename(Path(input_file)).split('.')[:-1])
    extension = os.path.basename(Path(input_file)).split('.')[-1]
    
    # Convert PDB to SDF
    if extension.lower() == "pdb":
        input_file = f'{os.path.dirname(Path(input_file))}/{filename}.sdf'
        extension = "sdf"
    
    # Load molecules based on format
    if extension.lower() == "sdf":
        mols = load_sdf(input_file,keep_xyz=keep_xyz)
    elif extension.lower() == "mol":
        mols = [Chem.MolFromMolFile(input_file, removeHs=False)]
    elif extension.lower() == "mol2":
        mols = [Chem.MolFromMol2File(input_file, removeHs=False)]
    
    # Extract properties from file
    IDs, charges, mults = _extract_sdf_properties(input_file)
    
    # Create supplier list
    suppl = [mol for mol in mols]
    
    # Generate IDs if not found
    if len(IDs) == 0:
        IDs = _generate_default_ids(input_file, len(suppl))
    
    # Determine charges
    charges = _calculate_charges_from_mols(mols, args) if len(charges) == 0 else charges
    
    # Determine multiplicities
    mults = _calculate_mults_from_mols(mols, args) if len(mults) == 0 else mults
    
    return suppl, charges, mults, IDs


def mol_from_sdf_or_mol_or_mol2(input_file, module, args, low_check=None, keep_xyz=False):
    """Load molecule objects from SDF, MOL, or MOL2 files.
    
    Handles different file formats and modules (qprep, cmin, csearch).
    For qprep/cmin, can filter conformers by count or energy window.
    
    Parameters
    ----------
    input_file : str
        Path to input file (SDF, MOL, or MOL2)
    module : str
        AQME module name ('qprep', 'cmin', or 'csearch')
    args : argparse.Namespace
        Arguments object with charge/mult settings
    low_check : str, int, float, or None
        Filtering criterion:
        - 'lowest_only': Return only lowest energy conformer
        - int: Return first N conformers
        - float: Return conformers within energy window (kcal/mol)
        - None: Return all conformers
    keep_xyz : bool
        Keep the generated XYZ for CREST runs in CSEARCH when using an XYZ input
    
    Returns
    -------
    For qprep/cmin:
        list: Molecule objects (optionally filtered)
    For csearch:
        tuple: (suppl, charges, mults, IDs)
    """
    if module in ["qprep", "cmin"]:
        return _load_mols_for_qprep_cmin(input_file, low_check)
    
    elif module == "csearch":
        return _load_mols_for_csearch(input_file, args, keep_xyz=keep_xyz)


def load_sdf(input_file, keep_xyz=False):
    """Load molecules from SDF file with error recovery.
    
    Attempts to load SDF file. If loading fails (e.g., from GaussView),
    uses OpenBabel to repair the file. If that fails, converts via XYZ.
    
    Args:
        input_file (str): Path to SDF file
        keep_xyz (bool): Keep XYZ input
    
    Returns:
        list: RDKit molecule objects with hydrogens preserved
    """
    mols = Chem.SDMolSupplier(input_file, removeHs=False)
    
    # Some software don't generate valid mol objects (e.g., GaussView)
    if None in mols:
        # Try repairing with OpenBabel
        command_sdf = [
            "obabel",
            "-isdf",
            f"{input_file}",
            "-osdf",
            f"-O{input_file}",
        ]
        subprocess.run(
            command_sdf,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        
        mols = Chem.SDMolSupplier(input_file, removeHs=False)
        
        # If still failing, convert via XYZ
        if None in mols:
            xyz_file = input_file.replace('.sdf', '.xyz')
            command_xyz = [
                "obabel",
                "-isdf",
                f"{input_file}",
                "-oxyz",
                f"-O{xyz_file}",
            ]
            subprocess.run(
                command_xyz,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            
            mols = [Chem.rdmolfiles.MolFromXYZFile(xyz_file)]
            if not keep_xyz:
                os.remove(xyz_file)
    
    return mols


def add_prefix_suffix(name, args):
    """Add prefix and/or suffix to filename.
    
    Args:
        name (str): Base filename
        args: Arguments object with prefix and suffix attributes
    
    Returns:
        str: Modified filename with prefix/suffix
    """
    if args.prefix != "":
        name = f"{args.prefix}_{name}"
    if args.suffix != "":
        name += f'_{args.suffix}'
    
    return name


def check_files(self, module):
    """Validate that input files exist and are accessible.
    
    Args:
        self: Arguments object with files list and logger
        module (str): Module name for error messages
    """
    if module.lower() in ['cmin', 'qdescp', 'qprep']:
        format_file = "*.sdf"
    elif module.lower() in ['qcorr']:
        format_file = "*.log"
    
    no_file_found = False
    
    if len(self.args.files) == 0:
        self.args.log.write(
            f'\nx  No files were found! In the "files" option, make sure that '
            f'1) the PATH to the files is correct, 2) the PATH doesn\'t start with "/", '
            f'and 3) you use quotation marks if you are using * (i.e. --files "{format_file}")'
        )
        no_file_found = True
    else:
        for file in self.args.files:
            if not os.path.exists(file):
                no_file_found = True
                self.args.log.write(
                    f'\nx  File {file} was not found! In the "files" option, '
                    f'make sure that 1) the PATH to the files is correct and '
                    f'2) the PATH doesn\'t start with "/".'
                )
                break
    
    if no_file_found:
        self.args.log.finalize()
        sys.exit()


def check_xtb(self):
    """Verify xTB is installed and accessible.
    
    Attempts to run xTB help command. Exits with error if not found.
    
    Args:
        self: Object with args.log for error reporting
    """
    try:
        subprocess.run(
            ["xtb", "-h"], 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        self.args.log.write(
            "x  xTB is not installed (CSEARCH-CREST and CMIN-xTB cannot be used)! "
            "You can install the program with 'conda install -y -c conda-forge xtb'"
        )
        self.args.log.finalize()
        sys.exit()


def check_crest(self):
    """Verify CREST is installed and accessible.
    
    Attempts to run CREST help command. Exits with error if not found.
    
    Args:
        self: Object with args.log for error reporting
    """
    try:
        subprocess.run(
            ["crest", "-h"], 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        self.args.log.write(
            "x  CREST is not installed (CSEARCH-CREST cannot be used)! "
            "You can install the program with 'conda install -y -c conda-forge crest'"
        )
        self.args.log.finalize()
        sys.exit()
 

def get_files(value):
    """Process and expand file path specifications.
    
    Handles wildcards, list parsing, and path resolution.
    Converts relative paths to absolute and expands glob patterns.
    
    Args:
        value (str or list): File path(s), potentially with wildcards
    
    Returns:
        list: Expanded list of absolute file paths
    """
    # Parse string to list if needed
    if not isinstance(value, list):
        value = value.replace('[', ']').replace(',', ']').split(']')
        while '' in value:
            value.remove('')
    
    new_value = []
    for val in value:
        if not isinstance(val, Path):
            # Add current directory to relative paths
            if Path(f"{val}").exists() and os.getcwd() not in f"{val}":
                new_value.append(f"{os.getcwd()}/{val}")
            
            # Expand wildcards
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


def _check_openbabel(self):
    """Check if OpenBabel is installed.
    
    Args:
        self: Arguments object with logger
    """
    try:
        command_run_1 = ["obabel", "-H"]
        subprocess.run(
            command_run_1, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        self.args.log.write(
            f"x  Open Babel is not installed! You can install the program with "
            f"'conda install -y -c conda-forge openbabel={obabel_version}'"
        )
        self.args.log.finalize()
        sys.exit()


def _check_rdkit(self):
    """Check if RDKit is installed.
    
    Args:
        self: Arguments object with logger
    """
    try:
        from rdkit.Chem import AllChem as Chem
    except ModuleNotFoundError:
        self.args.log.write(
            "x  RDKit is not installed! You can install the program with "
            "'pip install rdkit' or 'conda install -y -c conda-forge rdkit'"
        )
        self.args.log.finalize()
        sys.exit()


def _check_xtb_crest(self):
    """Check if xTB and CREST are installed for quantum calculations.
    
    Args:
        self: Arguments object with logger and program settings
    """
    # Check xTB installation
    try:
        command_run_1 = ["xtb", "-h"]
        subprocess.run(
            command_run_1, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        self.args.log.write(
            f"x  xTB is not installed! You can install the program with "
            f"'conda install -y -c conda-forge xtb={xtb_version}'"
        )
        self.args.log.finalize()
        sys.exit()
    
    # Check xTB version
    _ = check_version(
        self, 'xTB', 'xtb version', xtb_version, 3,
        f'conda install -y -c conda-forge xtb={xtb_version}'
    )
    
    # Check CREST if needed
    if self.args.program.lower() == 'crest':
        try:
            command_run_1 = ["crest", "-h"]
            subprocess.run(
                command_run_1, 
                stdout=subprocess.DEVNULL, 
                stderr=subprocess.DEVNULL
            )
        except FileNotFoundError:
            self.args.log.write(
                f"x  CREST is not installed! You can install the program with "
                f"'conda install -y -c conda-forge crest={crest_version}'"
            )
            self.args.log.finalize()
            sys.exit()
        
        # Check CREST version
        _ = check_version(
            self, 'CREST', 'Version', crest_version, 1,
            f'conda install -y -c conda-forge crest={crest_version}'
        )


def _check_ani_dependencies(self):
    """Check if ANI-related dependencies are installed.
    
    Args:
        self: Arguments object with logger
    """
    # Check torch
    try:
        import torch
        import warnings
        warnings.filterwarnings('ignore')
    except ModuleNotFoundError:
        self.args.log.write(
            "x  Torch-related modules are not installed! You can install these modules with "
            "'pip install torch torchvision torchani'"
        )
        self.args.log.finalize()
        sys.exit()
    
    # Check ASE
    try:
        import ase
        import ase.optimize
    except ModuleNotFoundError:
        self.args.log.write(
            "x  ASE is not installed! You can install the program with "
            "'pip install ase' or 'conda install -y -c conda-forge ase'"
        )
        self.args.log.finalize()
        sys.exit()
    
    # Check torchani
    try:
        import torchani
    except (ImportError, ModuleNotFoundError):
        self.args.log.write(
            "x  Torchani is not installed! You can install the program with "
            "'pip install torchani'"
        )
        self.args.log.finalize()
        sys.exit()


def check_dependencies(self):
    """Check that all required dependencies are installed.
    
    Validates OpenBabel, RDKit, and program-specific dependencies
    (xTB/CREST for quantum calculations, or torch/ASE/torchani for ANI).
    
    Args:
        self: Arguments object with program settings and logger
    """
    # Check core dependencies
    _check_openbabel(self)
    _check_rdkit(self)
    
    # Check program-specific dependencies
    if self.args.program is not None:
        if self.args.program.lower() in ['xtb', 'crest']:
            _check_xtb_crest(self)
        elif self.args.program.lower() == 'ani':
            _check_ani_dependencies(self)

def check_version(self, program, version_line, target_version, n_split, install_cmd):
    """Verify program version compatibility with AQME.
    
    Checks if installed xTB/CREST version is compatible. Allows versions
    up to 2-4 minor releases ahead of target.
    
    Args:
        self: Object with args.log and args.initial_dir
        program (str): Program name ('xTB' or 'CREST')
        version_line (str): String to search for in version output
        target_version (str): Minimum required version
        n_split (int): Index of version number in split line
        install_cmd (str): Command to install correct version
    """
    file_txt = self.args.initial_dir.joinpath(f'{program.lower()}_internal_test.txt')
    version_found = '0.0.0'
    
    # Run version command
    command_run_1 = [program.lower(), "--version", '>', f'{file_txt}']
    run_command(command_run_1, f'{file_txt}', cwd=self.args.initial_dir)
    subprocess.run(command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Parse version from output
    with open(file_txt, 'r') as datfile:
        lines = datfile.readlines()
        for _, line in enumerate(lines):
            if version_line in line:
                if program.lower() == 'crest':
                    version_found = line.split()[n_split].split(',')[0]
                else:
                    version_found = line.split()[n_split]
                break
    
    if os.path.exists(file_txt):
        os.remove(file_txt)
    
    # Check version compatibility
    lower_version = True
    if int(version_found.split('.')[0]) == int(target_version.split('.')[0]):
        # For X.Y versions: allow up to 2 minor versions ahead
        if (len(target_version.split('.')) == 2 and
            int(version_found.split('.')[1]) >= int(target_version.split('.')[1]) and
            int(version_found.split('.')[1]) < int(target_version.split('.')[1]) + 3):
            lower_version = False
        
        # For X.Y.Z versions: allow up to 4 patch versions ahead
        if (len(target_version.split('.')) == 3 and
            int(version_found.split('.')[1]) == int(target_version.split('.')[1]) and
            int(version_found.split('.')[2]) >= int(target_version.split('.')[2]) and
            int(version_found.split('.')[2]) < int(target_version.split('.')[2]) + 5):
            lower_version = False
    
    if version_found == '0.0.0':
        version_found = 'Unknown'
    
    self.args.log.write(f"{program} version used: {version_found}\n")
    
    if lower_version:
        self.args.log.write(
            f"x  {program} needs to be adjusted to ensure that AQME works as intended! "
            f"You can adjust the version with '{install_cmd}'"
        )
        self.args.log.finalize()
        sys.exit()


def blocking_wrapper(func, *args, **kwargs):
    """Execute function and wait for any subprocess.Popen objects.
    
    Ensures all subprocess.Popen objects returned by the function
    complete before returning. Useful for managing parallel processes.
    
    Args:
        func (callable): Function to execute
        *args: Positional arguments for func
        **kwargs: Keyword arguments for func
    
    Returns:
        Return value from func (after waiting for any Popen objects)
    """
    result = func(*args, **kwargs)
    
    # If result is a list of Popen objects, wait on them
    if isinstance(result, list) and all(hasattr(p, 'wait') for p in result):
        for p in result:
            p.wait()
    
    return result