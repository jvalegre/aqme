#####################################################.
#            This file stores functions             #
#                 used in CSEARCH                   #
#####################################################.

import os
import sys
import subprocess
import pandas as pd
import ast
from pathlib import Path
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolAlign

from aqme.utils import (
    get_info_input,
    mol_from_sdf_or_mol_or_mol2,
    read_xyz_charge_mult,
    add_prefix_suffix,
    periodic_table
)
from aqme.csearch.crest import nci_ts_mol


def csv_2_list(constraints):
    """Convert CSV constraints to a list format.
    
    Args:
        constraints: CSV constraint data, can be string, list or NaN
        
    Returns:
        list: List of constraints, empty list if input was NaN
    """
    try:
        if pd.isnull(constraints):
            constraints = []
    except ValueError:
        pass
    if not isinstance(constraints, list):
        constraints = ast.literal_eval(constraints)
    
    return constraints


def prepare_direct_smi(args):
    """Prepare job inputs from a direct SMILES input.
    
    Args:
        args: Arguments object containing SMILES input configuration
        
    Returns:
        list: List containing single job configuration tuple
        
    Raises:
        SystemExit: If no name is provided for the SMILES input
    """
    if args.name is None:
        args.log.write("\nx  Specify a name ('name' option) when using the 'smi' option!")
        args.log.finalize()
        sys.exit()

    name = add_prefix_suffix(args.name, args)
    
    # Create job configuration tuple
    job_config = (
        args.smi,         # SMILES string
        name,            # Molecule name
        args.charge,     # Charge
        args.mult,       # Multiplicity
        args.constraints_atoms,
        args.constraints_dist,
        args.constraints_angle, 
        args.constraints_dihedral,
        args.complex_type,
        args.geom,
        args.sample
    )
    
    return [job_config]


def prepare_smiles_files(args, csearch_file):
    """Process SMILES data from a file into job configurations.

    Reads SMILES data from a file where each line contains a SMILES string
    followed by a molecule name. Creates job configurations for each molecule.

    Args:
        args: Configuration object containing processing parameters
        csearch_file (str): Path to file containing SMILES data

    Returns:
        list: List of job configuration tuples for each molecule
    """
    # Read and filter out empty lines
    with open(csearch_file) as smifile:
        lines = [line for line in smifile if line.strip()]
        
    job_inputs = []
    for line in lines:
        # Process SMILES and name from each line
        smi, name = prepare_smiles_from_line(line, args)
        
        # Create job configuration
        job_config = (
            smi,                     # SMILES string
            name,                    # Molecule name
            args.charge,             # Charge
            args.mult,               # Multiplicity
            args.constraints_atoms,  # Atom constraints
            args.constraints_dist,   # Distance constraints
            args.constraints_angle,  # Angle constraints
            args.constraints_dihedral, # Dihedral constraints
            args.complex_type,      # Complex type
            args.geom,              # Geometry
            args.sample             # Sample
        )
        job_inputs.append(job_config)

    return job_inputs


def prepare_smiles_from_line(line, args):
    """Extract SMILES and molecule name from a line of text.

    Processes a space-separated line where the first token is a SMILES string
    and the second token is the molecule name. Also handles special SMILES
    syntax for nitrogen chirality.

    Args:
        line (str): Space-separated line containing SMILES and name
        args: Configuration object containing processing parameters

    Returns:
        tuple: (smiles, name) where:
            - smiles (str): Processed SMILES string
            - name (str): Processed molecule name with prefix/suffix
            
    Note:
        N@@ and N@ in SMILES are replaced with N as chiral nitrogen is not supported
    """
    # Split line into tokens
    tokens = line.split()
    if len(tokens) < 2:
        args.log.write("\nx  Error: Line must contain SMILES and name")
        return None, None
        
    smiles = tokens[0]
    name = add_prefix_suffix(tokens[1], args)

    # Handle unsupported chiral nitrogen in SMILES
    if "N@@" in smiles or "N@" in smiles:
        args.log.write(
            f"\nx  WARNING! AQME does not support chiral N atoms in SMILES strings "
            f"(N@@ or N@). These atoms were replaced by N in the SMILES: {smiles}."
        )
        smiles = smiles.replace("N@@", "N").replace("N@", "N")

    return smiles, name


def prepare_csv_files(args, csearch_file):
    """Process molecule data from a CSV file into job configurations.

    Reads molecule data from a CSV file that must contain SMILES strings and molecule names.
    Performs various validations and creates unique job configurations for each molecule.

    Args:
        args: Configuration object containing processing parameters
        csearch_file (str): Path to CSV file containing molecule data

    Returns:
        list: List of job configuration tuples for each unique molecule

    Raises:
        SystemExit: If file is empty, missing required columns, or contains invalid data
    """
    # Validate file exists and is not empty
    if os.path.getsize(csearch_file) == 0:
        args.log.write(f"File {args.input} is empty!")
        args.log.finalize()
        sys.exit()
    
    # Read and validate CSV content
    csv_smiles = pd.read_csv(csearch_file)
    if csv_smiles.empty:
        args.log.write(f"File {args.input} is empty!")
        args.log.finalize()
        sys.exit()
    
    # Validate code_name column
    if "code_name" in csv_smiles.columns and csv_smiles["code_name"].dropna().empty:
        args.log.write(f"File {args.input} has a 'code_name' column with no values.")
        args.log.finalize()
        sys.exit()

    # Validate molecule names
    for name in csv_smiles['code_name']:
        if '*' in str(name):
            args.log.write(
                f"\nx  WARNING! The names provided in the CSV contain * "
                f"(i.e. {name}). Please, remove all the * characters."
            )
            args.log.finalize()
            sys.exit()

    # Process SMILES columns and generate job configurations
    job_inputs = []
    unique_smiles = set()
    has_smiles_column = False
    
    for col_idx, column in enumerate(csv_smiles.columns):
        if column.upper() == "SMILES" or "SMILES_" in column.upper():
            has_smiles_column = True
            
            # Process each row in the SMILES column
            for row_idx in range(len(csv_smiles)):
                mol_config = generate_mol_from_csv(args, csv_smiles, row_idx, col_idx)
                
                if mol_config is not None:
                    smiles, name = mol_config[0], mol_config[1]
                    
                    # Only add unique SMILES
                    if smiles not in unique_smiles:
                        job_inputs.append(mol_config)
                        unique_smiles.add(smiles)
                    else:
                        args.log.write(
                            f'\nx  SMILES "{smiles}" used in {name} is a duplicate, '
                            'it was already used with a different code_name!'
                        )
                       
    if not has_smiles_column:
        args.log.write(
            "\nx  Make sure the CSV file contains a column called 'SMILES', "
            "'smiles' or 'SMILES_' with the SMILES of the molecules!"
        )
        args.log.finalize()
        sys.exit()
        
    return job_inputs


def generate_mol_from_csv(args, csv_smiles, index, column_index):
    """Generate molecule configuration from CSV data.

    Extracts SMILES and molecule properties from a specific row and column in the CSV data,
    performing necessary validations and transformations. Supports optional properties
    like charge, multiplicity, and various constraints.

    Args:
        args: Configuration object containing processing parameters
        csv_smiles (pandas.DataFrame): DataFrame containing molecule data
        index (int): Row index in the DataFrame
        column_index (int): Column index for SMILES data

    Returns:
        tuple: Job configuration tuple if valid, None if invalid data.
            Format: (smiles, name, charge, mult, constraints...)
            
    Raises:
        SystemExit: If required 'code_name' column is missing
    """
    # Get SMILES from specified column
    column_name = csv_smiles.columns[column_index]
    smiles = csv_smiles.loc[index, column_name]
    
    # Skip empty or NaN SMILES
    if pd.isna(smiles) or str(smiles).lower() == 'nan':
        return None
        
    # Handle unsupported chiral nitrogen in SMILES
    if "N@@" in str(smiles) or "N@" in str(smiles):
        args.log.write(
            f"\nx  WARNING! AQME does not support chiral N atoms in SMILES strings "
            f"(N@@ or N@). These atoms were replaced by N in the SMILES: {smiles}."
        )
        smiles = str(smiles).replace("N@@", "N").replace("N@", "N")

    # Process molecule name
    try:
        name = str(csv_smiles.loc[index, "code_name"])
        # Add suffix based on SMILES column name
        if column_name.upper() != "SMILES" and "_" in column_name:
            name += "_" + column_name.split("_")[-1]
    except KeyError:
        args.log.write("\nx  Make sure the CSV file contains a column called 'code_name' with the names of the molecules!")
        args.log.finalize()
        sys.exit()

    # Initialize properties with default values
    charge = args.charge
    mult = args.mult
    constraints_atoms = args.constraints_atoms
    constraints_dist = args.constraints_dist
    constraints_angle = args.constraints_angle
    constraints_dihedral = args.constraints_dihedral
    complex_type = args.complex_type
    geom = args.geom
    sample = args.sample

    # Helper function to get non-NaN values from DataFrame
    def get_csv_value(col_name, default_value):
        if col_name in csv_smiles.columns:
            value = csv_smiles.loc[index, col_name]
            if not pd.isna(value) and str(value).lower() != 'nan':
                return value
        return default_value

    # Get molecule properties with fallbacks to defaults
    charge = get_csv_value("charge", charge)
    mult = get_csv_value("mult", mult)
    complex_type = get_csv_value("complex_type", complex_type)
    sample = get_csv_value("sample", sample)

    # Get and process constraint values
    constraints_atoms = csv_2_list(get_csv_value("constraints_atoms", constraints_atoms))
    constraints_dist = csv_2_list(get_csv_value("constraints_dist", constraints_dist))
    constraints_angle = csv_2_list(get_csv_value("constraints_angle", constraints_angle))
    constraints_dihedral = csv_2_list(get_csv_value("constraints_dihedral", constraints_dihedral))
    geom = csv_2_list(get_csv_value("geom", geom))

    # Create and return job configuration
    return (
        smiles,                 # SMILES string
        name,                   # Molecule name
        charge,                 # Charge
        mult,                   # Multiplicity
        constraints_atoms,      # Atom constraints
        constraints_dist,       # Distance constraints
        constraints_angle,      # Angle constraints
        constraints_dihedral,   # Dihedral constraints
        complex_type,          # Complex type
        geom,                  # Geometry
        sample                 # Sample
    )


def prepare_cdx_files(args, csearch_file):
    """Convert ChemDraw files to SMILES and prepare job configurations.

    Converts molecules from a ChemDraw file to SMILES format and creates job
    configurations for each molecule. Names are generated based on the input
    filename with an index suffix.

    Args:
        args: Configuration object containing processing parameters
        csearch_file (str): Path to ChemDraw file (.cdx)

    Returns:
        list: List of job configuration tuples for each molecule
    """
    # Convert ChemDraw molecules to SMILES format
    molecules = generate_mol_from_cdx(csearch_file)

    job_inputs = []
    # Process each molecule with an index suffix
    for i, (smiles, _) in enumerate(molecules):
        # Generate name from file basename with index
        basename = '.'.join(os.path.basename(Path(csearch_file)).split('.')[:-1])
        name = f"{basename}_{i}"
        name = add_prefix_suffix(name, args)

        # Create job configuration
        job_config = (
            smiles,                 # SMILES string
            name,                   # Molecule name
            args.charge,            # Charge
            args.mult,              # Multiplicity
            args.constraints_atoms,  # Atom constraints
            args.constraints_dist,   # Distance constraints
            args.constraints_angle,  # Angle constraints
            args.constraints_dihedral, # Dihedral constraints
            args.complex_type,      # Complex type
            args.geom,              # Geometry
            args.sample             # Sample
        )
        job_inputs.append(job_config)
    return job_inputs


def generate_mol_from_cdx(csearch_file):
    """Convert ChemDraw file to SMILES format using OpenBabel.

    Uses OpenBabel to convert a ChemDraw file to SMILES format and creates
    RDKit molecules for each SMILES string.

    Args:
        csearch_file (str): Path to ChemDraw file (.cdx)

    Returns:
        list: List of tuples (smiles, rdkit_molecule) for each molecule

    Note:
        Requires OpenBabel to be installed and available in the system path
    """
    # Convert CDX to SMILES using OpenBabel
    cdx_cmd = ["obabel", "-icdx", csearch_file, "-osmi", "-Ocdx.smi"]
    try:
        subprocess.run(
            cdx_cmd, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL,
            check=True
        )

        # Read generated SMILES file
        with open("cdx.smi", "r") as smifile:
            smi_lines = [str(line.strip()) for line in smifile]

        # Clean up temporary file
        os.remove("cdx.smi")

        # Convert SMILES to RDKit molecules
        molecules = []
        for smi in smi_lines:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:  # Check for valid molecules
                molecules.append((smi, mol))

        return molecules

    except subprocess.CalledProcessError:
        raise RuntimeError("Failed to convert ChemDraw file using OpenBabel")
    except FileNotFoundError:
        raise RuntimeError("OpenBabel not found. Please ensure it is installed and in your system path")


def prepare_com_files(args, csearch_file):
    """Process Gaussian input files or XYZ files into job configurations.
    
    Converts Gaussian (.com/.gjf) or XYZ files into a temporary SDF format,
    extracts molecular information, and creates job configurations. Handles
    charge and multiplicity from file or arguments.

    Args:
        args: Configuration object containing processing parameters
        csearch_file (str): Path to Gaussian or XYZ input file

    Returns:
        list: List containing single job configuration tuple
        
    Note:
        Temporary files are created and cleaned up during processing
    """
    # Get filename and extension
    file_path = Path(csearch_file)
    filename = file_path.name
    extension = file_path.suffix.lower()[1:]  # Remove leading dot
    
    # Process Gaussian files
    if extension in ["gjf", "com"]:
        xyz_file, file_charge, file_mult = com_2_xyz(csearch_file)
        charge = args.charge if args.charge is not None else file_charge
        mult = args.mult if args.mult is not None else file_mult
    # Process XYZ files
    else:
        xyz_file = csearch_file
        file_charge, file_mult = read_xyz_charge_mult(csearch_file)
        charge = args.charge if args.charge is not None else file_charge
        mult = args.mult if args.mult is not None else file_mult

    try:
        # Convert to SDF format
        xyz_2_sdf(xyz_file)
        
        # Create SDF path and read molecule
        sdf_path = file_path.with_suffix('.sdf')
        suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(str(sdf_path), "csearch", args)
        
        # Process name and create job configuration
        name = add_prefix_suffix(file_path.stem, args)
        job_config = (
            suppl[0],                # Molecule
            name,                    # Name
            charge,                  # Charge
            mult,                    # Multiplicity
            args.constraints_atoms,  # Atom constraints
            args.constraints_dist,   # Distance constraints
            args.constraints_angle,  # Angle constraints
            args.constraints_dihedral, # Dihedral constraints
            args.complex_type,      # Complex type
            args.geom,              # Geometry
            args.sample             # Sample
        )
        
        return [job_config]
        
    finally:
        # Clean up temporary files
        if extension in ["gjf", "com"]:
            os.remove(xyz_file)
        os.remove(str(sdf_path))


def prepare_pdb_files(args, csearch_file):
    """Process PDB files into job configurations via SDF conversion.
    
    Converts PDB files to SDF format using OpenBabel, then processes them
    into job configurations.

    Args:
        args: Configuration object containing processing parameters
        csearch_file (str): Path to PDB file

    Returns:
        list: List of job configuration tuples for each molecule in PDB file
        
    Note:
        Creates and cleans up temporary SDF files during processing
    """
    file_path = Path(csearch_file)
    sdf_path = file_path.with_suffix('.sdf')
    
    try:
        # Convert PDB to SDF using OpenBabel
        command_pdb = [
            "obabel",
            "-ipdb", str(file_path),
            "-osdf",
            f"-O{sdf_path}"
        ]
        subprocess.run(
            command_pdb, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL,
            check=True
        )
        
        # Process the SDF file
        return prepare_sdf_files(args, csearch_file)
        
    except subprocess.CalledProcessError:
        raise RuntimeError("Failed to convert PDB file using OpenBabel")
    finally:
        # Clean up temporary file
        if sdf_path.exists():
            os.remove(sdf_path)


def prepare_sdf_files(args, csearch_file):
    """Process SDF files into job configurations.
    
    Extracts molecule data from SDF files and creates job configurations
    for each molecule found, preserving charge, multiplicity and identifiers.

    Args:
        args: Configuration object containing processing parameters
        csearch_file (str): Path to SDF file

    Returns:
        list: List of job configuration tuples for each molecule in SDF file
    """
    sdf_path = Path(csearch_file)
    
    # Read molecules from SDF file
    suppl, charges, mults, ids = mol_from_sdf_or_mol_or_mol2(
        str(sdf_path), 
        "csearch", 
        args
    )
    
    # Process base names only
    mol_ids = [Path(id).name for id in ids]
    
    # Create job configurations for each molecule
    job_inputs = []
    for mol, charge, mult, mol_id in zip(suppl, charges, mults, mol_ids):
        name = add_prefix_suffix(mol_id, args)
        job_config = (
            mol,                     # Molecule
            name,                    # Name
            charge,                  # Charge
            mult,                    # Multiplicity
            args.constraints_atoms,  # Atom constraints
            args.constraints_dist,   # Distance constraints
            args.constraints_angle,  # Angle constraints
            args.constraints_dihedral, # Dihedral constraints
            args.complex_type,      # Complex type
            args.geom,              # Geometry
            args.sample             # Sample
        )
        job_inputs.append(job_config)
        
    return job_inputs


def xyz_2_sdf(file):
    """Convert XYZ file to SDF format using OpenBabel.

    Creates a .sdf file from a .xyz file in the same directory as the input.
    Uses OpenBabel for the conversion.

    Args:
        file (str): Path to existing .xyz file
        
    Raises:
        subprocess.CalledProcessError: If OpenBabel conversion fails
        FileNotFoundError: If OpenBabel is not installed
    """
    file_path = Path(file)
    output_path = file_path.with_suffix('.sdf')
    
    try:
        # Convert XYZ to SDF using OpenBabel
        command = [
            "obabel", 
            "-ixyz", str(file_path),
            "-osdf", 
            f"-O{output_path}"
        ]
        subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )
    except subprocess.CalledProcessError:
        raise RuntimeError("Failed to convert XYZ file using OpenBabel")
    except FileNotFoundError:
        raise RuntimeError("OpenBabel not found. Please ensure it is installed")


def check_constraints(constraints_atoms, constraints_dist, 
                     constraints_angle, constraints_dihedral):
    """Check if any constraints are defined.

    Determines if any type of constraint (atoms, distances, angles, dihedrals)
    has been specified.

    Args:
        constraints_atoms (list): Atom constraints
        constraints_dist (list): Distance constraints
        constraints_angle (list): Angle constraints
        constraints_dihedral (list): Dihedral angle constraints

    Returns:
        bool: True if any constraints are defined, False otherwise
    """
    return any([
        len(constraints_atoms) > 0,
        len(constraints_dist) > 0,
        len(constraints_angle) > 0,
        len(constraints_dihedral) > 0
    ])


def com_2_xyz(input_file):
    """Convert Gaussian input file to XYZ format.

    Extracts geometry and electronic structure information from a Gaussian
    input file (.com/.gjf) and writes it to XYZ format.

    Args:
        input_file (str): Path to Gaussian input file

    Returns:
        tuple: (xyz_path, charge, mult) where:
            - xyz_path (str): Path to generated XYZ file
            - charge (int): Total molecular charge
            - mult (int): Spin multiplicity
            
    Note:
        Creates a new XYZ file in the same directory as the input file
    """
    # Setup paths
    file_path = Path(input_file)
    xyz_path = file_path.with_suffix('.xyz')
    
    # Extract geometry and properties
    xyz_lines, charge, mult = get_info_input(input_file)
    
    # Write XYZ file
    with open(xyz_path, "w") as f:
        f.write(f"{len(xyz_lines)}\n")          # Number of atoms
        f.write(f"{file_path.stem}\n")          # Comment line
        f.write("\n".join(xyz_lines) + "\n")    # Coordinates

    return str(xyz_path), charge, mult


def realign_mol(mol, conf, coord_Map, alg_Map, mol_template, maxsteps):
    """Minimize and align a molecule while preserving template atom positions.

    Performs force field minimization and molecular alignment on a molecule,
    keeping certain atoms fixed based on a template structure. The minimization
    uses the UFF force field and the alignment is based on matching atoms.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to be minimized and aligned
        conf (int): Conformation ID for the minimization and alignment
        coord_Map (list): List of atom indices for coordinate constraints
        alg_Map (list): List of atom indices for alignment matching
        mol_template (rdkit.Chem.rdchem.Mol): Template molecule for alignment
        maxsteps (int): Maximum number of force field optimization steps

    Returns:
        tuple: (mol, energy) where:
            - mol (rdkit.Chem.rdchem.Mol): Updated molecule after minimization/alignment
            - energy (float): Final UFF force field energy
            
    Note:
        This function combines minimization and alignment steps and may need
        refactoring to separate these operations in the future.
    """
    # Setup UFF forcefield with conformation
    forcefield = Chem.UFFGetMoleculeForceField(mol, confId=conf)
    
    # Find matching atoms and add distance constraints
    matching_atoms = mol.GetSubstructMatch(mol_template)
    for i, atom_i in enumerate(matching_atoms):
        # Add pairwise distance constraints between matching atoms
        for atom_j in matching_atoms[i + 1:]:
            # Get target distance from coordinate map
            target_dist = coord_Map[atom_i].Distance(coord_Map[atom_j])
            # Add strong distance constraint (force constant = 10000)
            forcefield.AddDistanceConstraint(atom_i, atom_j, target_dist, target_dist, 10000)
    
    # Run energy minimization
    forcefield.Initialize()
    forcefield.Minimize(maxIts=maxsteps)
    
    # Align optimized molecule to template
    rdMolAlign.AlignMol(
        mol,                # Molecule to align
        mol_template,       # Template to align to
        prbCid=conf,        # Probe conformer ID
        refCid=-1,         # Reference conformer ID (-1 = first)
        atomMap=alg_Map,    # Atom mapping for alignment
        reflect=True,       # Try mirror image if needed
        maxIters=100,       # Maximum alignment iterations
    )
    
    # Get final energy after minimization and alignment
    energy = float(forcefield.CalcEnergy())
    return mol, energy


def minimize_rdkit_energy(mol, conf, log, FF, maxsteps):
    """Minimize molecular energy using RDKit force fields.
    
    Attempts to minimize a molecule's energy using either MMFF94 or UFF force
    fields. Falls back to UFF if MMFF fails, and handles minimization failures
    gracefully.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to minimize
        conf (int): Conformer ID to minimize
        log: Logger object for status messages
        FF (str): Force field to use ('MMFF', 'UFF', or 'NO FF')
        maxsteps (int): Maximum number of minimization steps

    Returns:
        float: Final energy of the minimized structure, or 0 if no FF used
        
    Note:
        Falls back to UFF if MMFF fails. Reports non-optimized geometries
        if minimization fails completely.
    """
    if FF.upper() == "NO FF":
        return 0.0

    # Try MMFF94 first if requested
    forcefield = None
    if FF.upper() == "MMFF":
        properties = Chem.MMFFGetMoleculeProperties(mol)
        forcefield = Chem.MMFFGetMoleculeForceField(mol, properties, confId=conf)
        if forcefield is None:
            log.write(f"x  Force field {FF} did not work! Falling back to UFF.")

    # Fall back to UFF if MMFF failed or was not requested
    if FF.upper() == "UFF" or forcefield is None:
        forcefield = Chem.UFFGetMoleculeForceField(mol, confId=conf)

    # Attempt minimization
    try:
        forcefield.Initialize()
        forcefield.Minimize(maxIts=maxsteps)
        return float(forcefield.CalcEnergy())
    except RuntimeError:
        log.write(f"\nx  Geometry minimization failed with {FF}, using non-optimized geometry.")
        return float(forcefield.CalcEnergy())


def getDihedralMatches(mol, heavy):
    """Find unique rotatable bonds and their associated dihedral angles.
    
    Uses RDKit's strict rotatable bond SMARTS pattern to find all possible
    dihedral angles in a molecule, filtering by heavy atoms if requested.
    Excludes certain groups like CF3, CCl3, terminal atoms, etc.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to search for dihedrals
        heavy (bool): If True, only consider dihedrals between heavy atoms

    Returns:
        list: List of 4-tuples of atom indices defining unique dihedral angles
        
    Note:
        Uses RDKit's strict pattern which excludes:
        - Triple bonds
        - Terminal atoms
        - CF3, CCl3, CBr3 groups
        - t-Butyl groups
        - Amide bonds
    """
    # RDKit's strict rotatable bond pattern
    STRICT_PATTERN = (
        r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)"
        r"&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])"
        r"&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])"
        r"&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)"
        r"&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
    )
    
    # Convert SMARTS pattern to molecule query
    query = Chem.MolFromSmarts(STRICT_PATTERN)
    matches = mol.GetSubstructMatches(query)

    # Filter and uniquify matches
    unique_dihedrals = []
    seen_bonds = set()  # Track unique central bonds
    
    for atoms in matches:
        a, b, c, d = atoms
        # Check if central bond is new
        if (b, c) not in seen_bonds and (c, b) not in seen_bonds:
            # Get atom symbols for filtering
            a_symbol = mol.GetAtomWithIdx(a).GetSymbol()
            c_symbol = mol.GetAtomWithIdx(c).GetSymbol()
            d_symbol = mol.GetAtomWithIdx(d).GetSymbol()
            
            # Apply filters based on heavy flag
            if heavy:
                # Only accept if both terminal atoms are non-hydrogen
                if a_symbol != "H" and d_symbol != "H":
                    seen_bonds.add((b, c))
                    unique_dihedrals.append(atoms)
            else:
                # Skip specific C-H bonds but accept others
                if not (c_symbol == "C" and d_symbol == "H"):
                    seen_bonds.add((b, c))
                    unique_dihedrals.append(atoms)
                    
    return unique_dihedrals


def smi_to_mol(
    smi,
    program,
    log,
    seed,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
):
    """Convert SMILES to RDKit molecule with constraints handling.
    
    Creates an RDKit molecule from SMILES string, handling special cases like
    complexes, transition states, and mapped atoms. Supports constraint
    application for conformer generation.

    Args:
        smi (str): SMILES string to convert
        program (str): Program to use for conformer generation
        log: Logger object for status messages
        seed (int): Random seed for reproducibility
        constraints_atoms (list): Atom-based constraints
        constraints_dist (list): Distance constraints
        constraints_angle (list): Angle constraints
        constraints_dihedral (list): Dihedral angle constraints

    Returns:
        tuple: (mol, constraints_atoms, constraints_dist, constraints_angle,
               constraints_dihedral, complex_ts) where:
            - mol (rdkit.Chem.rdchem.Mol): Generated molecule or None if failed
            - constraints_*: Updated constraint lists
            - complex_ts (bool): True if molecule is a complex/TS

    Note:
        For complexes (multi-part SMILES) or molecules with constraints,
        only 'crest' program is supported.
    """
    complex_ts = False
    smi_parts = smi.split(".")

    # Handle complexes and constrained systems
    if (len(smi_parts) > 1 or 
        any(len(c) > 0 for c in [constraints_atoms, constraints_dist,
                                constraints_angle, constraints_dihedral])):
        # Validate program choice for complexes
        if program not in ["crest"]:
            log.write(
                f"\nx  {program} not supported for conformer generation of complexes "
                f"and TSs (your SMILES has {len(smi_parts)} parts, separated by a "
                "period)! Specify: program='crest' for complexes"
            )
            sys.exit()

        # Process complex or TS molecule
        mol_data = nci_ts_mol(
            smi_parts,
            log,
            seed,
            constraints_atoms,
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
        )
        mol, *constraints = mol_data
        complex_ts = True

    else:
        # Process single molecule
        params = Chem.SmilesParserParams()
        params.removeHs = False
        smi = smi_parts[0]

        try:
            # Handle mapped atoms
            if ':' in smi:
                log.write(
                    f"\nx  WARNING! The SMILES string provided ({smi}) contains mapped "
                    "atoms, make sure you include their corresponding H atoms explicitly "
                    "in the SMILES (otherwise they'll be omitted). For example, use "
                    "[C:1]([H])([H])([H])C instead of [C:1]C.\n"
                )

            # Create and process molecule
            mol = Chem.MolFromSmiles(smi, params)
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)

        except Chem.AtomValenceException:
            log.write(
                f"\nx  The SMILES string provided ({smi}) contains errors or the "
                "molecule needs to be drawn differently. For example, N atoms from "
                "ligands of metal complexes should be N+ since they're drawn with "
                "four bonds in ChemDraw, same for O atoms in carbonyl ligands, etc.\n"
            )
            mol = None

        # Keep original constraints for single molecules
        constraints = [constraints_atoms, constraints_dist,
                      constraints_angle, constraints_dihedral]

    return (mol, *constraints, complex_ts)

def substituted_mol(mol, checkI, metal_atoms):
    """Process metal atoms in a molecule, optionally replacing them with iodine.
    
    Identifies metal atoms in a molecule and optionally replaces them with iodine
    atoms, adjusting formal charges to maintain valid valence states. Tracks metal
    atom positions and coordination numbers.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Input molecule
        checkI (str): If "I", replace metals with iodine atoms
        metal_atoms (list): List of metal element symbols to process

    Returns:
        tuple: (metal_idx, metal_sym) where:
            - metal_idx (list): Indices of metal atoms in molecule
            - metal_sym (list): Original symbols of metal atoms
            
    Note:
        When replacing with iodine (checkI="I"), the function attempts to
        find valid formal charges by incrementally adjusting from the base
        coordination number derived charge.
    """
    # Initialize tracking lists
    metal_idx = [None] * len(metal_atoms)
    complex_coord = [None] * len(metal_atoms)
    metal_sym = [None] * len(metal_atoms)

    # Map coordination number to base formal charge
    coord_to_charge = {
        coord: charge 
        for coord, charge in zip(range(2, 9), range(-3, 4))
    }

    # Process each atom in the molecule
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in metal_atoms:
            idx = metal_atoms.index(symbol)
            
            # Record metal atom information
            metal_sym[idx] = symbol
            metal_idx[idx] = atom.GetIdx()
            complex_coord[idx] = len(atom.GetNeighbors())
            
            # Replace with iodine if requested
            if checkI == "I":
                atom.SetAtomicNum(53)  # Atomic number of iodine
                n_neighbors = len(atom.GetNeighbors())
                
                if n_neighbors > 1:
                    # Get base charge from coordination number
                    base_charge = coord_to_charge.get(n_neighbors, 0)
                    
                    # Try different charge states until molecule is valid
                    for charge_adj in range(0, 5):
                        atom.SetFormalCharge(base_charge + charge_adj)
                        try:
                            # Test if molecule is valid with this charge
                            mol_test = Chem.Mol(mol)
                            Chem.SanitizeMol(mol_test)
                            break  # Valid charge found
                        except Chem.AtomValenceException:
                            continue  # Try next charge state

    return metal_idx, metal_sym

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

    return mol