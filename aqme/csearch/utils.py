#####################################################.
#            This file stores functions             #
#                 used in CSEARCH                   #
#####################################################.

import os
import sys
import itertools
import subprocess
import pandas as pd
import ast
import numpy as np
from pathlib import Path
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolAlign
import rdkit.Chem as _rdChem

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
        suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(str(sdf_path), "csearch", args, keep_xyz=True)
        
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
    keep_xyz = True if args.program.lower() == 'crest' else False # keep XYZ as inputs for CREST runs
    suppl, charges, mults, ids = mol_from_sdf_or_mol_or_mol2(
        str(sdf_path), 
        "csearch", 
        args,
        keep_xyz=keep_xyz
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


def realign_mol(mol, conf, coord_Map, alg_Map, mol_template, maxsteps,
                constraints_atoms=None, constraints_dist=None,
                constraints_angle=None, constraints_dihedral=None):
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
    constraints_atoms = constraints_atoms or []
    constraints_dist = constraints_dist or []
    constraints_angle = constraints_angle or []
    constraints_dihedral = constraints_dihedral or []

    # Setup UFF forcefield with conformation
    forcefield = Chem.UFFGetMoleculeForceField(mol, confId=conf)
    
    # Find matching atoms and add distance constraints
    matching_atoms = mol.GetSubstructMatch(mol_template)
    for i, atom_i in enumerate(matching_atoms):
        # Add pairwise distance constraints between matching atoms
        for atom_j in matching_atoms[i + 1:]:
            # Get target distance from coordinate map
            target_dist = coord_Map[atom_i].Distance(coord_Map[atom_j])
            # Add strong distance constraint (force constant = 100000)
            forcefield.AddDistanceConstraint(atom_i, atom_j, target_dist, target_dist, 100000)
    
    # Apply user-defined constraints
    has_constraints = any([constraints_atoms, constraints_dist,
                           constraints_angle, constraints_dihedral])
    if has_constraints:
        apply_rdkit_constraints(
            forcefield,
            constraints_atoms, constraints_dist,
            constraints_angle, constraints_dihedral
        )

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


def minimize_rdkit_energy(mol, conf, log, FF, maxsteps,
                            constraints_atoms=None, constraints_dist=None,
                            constraints_angle=None, constraints_dihedral=None):
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
    constraints_atoms = constraints_atoms or []
    constraints_dist = constraints_dist or []
    constraints_angle = constraints_angle or []
    constraints_dihedral = constraints_dihedral or []

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

    # Apply constraints before minimization
    has_constraints = any([constraints_atoms, constraints_dist,
                           constraints_angle, constraints_dihedral])
    if has_constraints:
        apply_rdkit_constraints(
            forcefield,
            constraints_atoms, constraints_dist,
            constraints_angle, constraints_dihedral
        )

    # Attempt minimization and handle failures
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


# ──────────────────────────────────────────────────────────────────────────────
#  RDKit aggregate (multi-fragment SMILES) helpers (THIS IS THE NEW THING)
# ──────────────────────────────────────────────────────────────────────────────

def get_interaction_points(mol):
    """Find chemical interaction points on a molecule for fragment placement.

    Identifies Hydrogen Bond Donors (HBD), Hydrogen Bond Acceptors (HBA),
    and aromatic ring centroids (Pi) as potential interaction sites.
    Falls back to the molecular centroid if no specific points are found.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule with a 3D conformer

    Returns:
        dict: Keys 'HBD', 'HBA', 'Pi', and optionally 'Fallback';
              values are lists of numpy arrays (3D coordinates).
    """
    positions = mol.GetConformer().GetPositions()  # numpy (N, 3)
    points = {'HBD': [], 'HBA': [], 'Pi': []}

    hbd_pattern = Chem.MolFromSmarts('[N,O,S;!H0]')
    for match in mol.GetSubstructMatches(hbd_pattern):
        points['HBD'].append(positions[match[0]])

    hba_pattern = Chem.MolFromSmarts('[N,O,F]')
    for match in mol.GetSubstructMatches(hba_pattern):
        points['HBA'].append(positions[match[0]])

    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            centroid = np.mean(positions[list(ring)], axis=0)
            points['Pi'].append(centroid)

    if not any(points.values()):
        points['Fallback'] = [np.mean(positions, axis=0)]

    return points


def _get_fibonacci_sphere_point(i, total_samples):
    """Return a unit vector on the Fibonacci sphere for index i.

    Ensures approximately uniform angular coverage when called for
    i = 0, 1, …, total_samples-1.

    Args:
        i (int): 0-based index of the point
        total_samples (int): Total number of points on the sphere

    Returns:
        np.ndarray: Unit vector [x, y, z]
    """
    if total_samples <= 1:
        return np.array([1.0, 0.0, 0.0])
    phi = np.pi * (3.0 - np.sqrt(5.0))  # golden angle in radians
    y = 1.0 - (i / float(total_samples - 1)) * 2.0
    radius = np.sqrt(max(0.0, 1.0 - y * y))
    theta = phi * i
    return np.array([np.cos(theta) * radius, y, np.sin(theta) * radius])


def _random_rotation_matrix(rng):
    """Generate a random 3x3 rotation matrix from three Euler angles.

    Args:
        rng: numpy random Generator (e.g. np.random.default_rng(seed))

    Returns:
        np.ndarray: 3x3 rotation matrix
    """
    angles = rng.uniform(0, 2 * np.pi, 3)
    cx, sx = np.cos(angles[0]), np.sin(angles[0])
    cy, sy = np.cos(angles[1]), np.sin(angles[1])
    cz, sz = np.cos(angles[2]), np.sin(angles[2])
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def _position_fragment_on_base(base_mol, mobile_mol, conformer_idx, total_conformers, seed=0):
    """Position a mobile fragment relative to a base fragment.

    Uses interaction-point pairing (HBD/HBA/Pi/Fallback) and the Fibonacci
    sphere algorithm to deterministically distribute the mobile fragment
    around the base. A reproducible random rotation is applied to the mobile
    fragment *before* translation for additional conformer diversity.

    Args:
        base_mol (rdkit.Chem.rdchem.Mol): Stationary fragment (pre-built 3D)
        mobile_mol (rdkit.Chem.rdchem.Mol): Fragment to be repositioned
            (its conformer is modified in-place before combining)
        conformer_idx (int): 0-based index of the current conformer
        total_conformers (int): Total number of conformers to generate
        seed (int): Base random seed; each conformer uses seed + conformer_idx

    Returns:
        rdkit.Chem.rdchem.Mol: Combined molecule from CombineMols(base_mol, mobile_mol)
    """
    base_pts = get_interaction_points(base_mol)
    mobile_pts = get_interaction_points(mobile_mol)

    # Build HBD/HBA/Pi interaction pairs with ideal target distances
    pairs = []
    for b in base_pts.get('HBD', []):
        for m in mobile_pts.get('HBA', []):
            pairs.append((b, m, 3.5))
    for b in base_pts.get('HBA', []):
        for m in mobile_pts.get('HBD', []):
            pairs.append((b, m, 3.5))
    for b in base_pts.get('Pi', []):
        for m in mobile_pts.get('Pi', []):
            pairs.append((b, m, 5))

    # Fall back to centroid-based placement when no chemical interaction points exist
    if not pairs:
        b_fallback = base_pts.get(
            'Fallback', [np.mean(base_mol.GetConformer().GetPositions(), axis=0)]
        )[0]
        m_fallback = mobile_pts.get(
            'Fallback', [np.mean(mobile_mol.GetConformer().GetPositions(), axis=0)]
        )[0]
        pairs.append((b_fallback, m_fallback, 3.0))

    # Cycle through pairs across conformers for angular diversity
    base_pt, mobile_pt, target_dist = pairs[conformer_idx % len(pairs)]

    # --- Step 1: rotate mobile fragment around its centroid ---
    rng = np.random.default_rng(seed + conformer_idx)
    rot = _random_rotation_matrix(rng)
    conf = mobile_mol.GetConformer()
    positions = conf.GetPositions()
    centroid = np.mean(positions, axis=0)
    rotated_positions = (rot @ (positions - centroid).T).T + centroid
    for atom_idx, pos in enumerate(rotated_positions):
        conf.SetAtomPosition(atom_idx, pos.tolist())

    # Update mobile interaction point after rotation
    mobile_pt_rotated = rot @ (mobile_pt - centroid) + centroid

    # --- Step 2: translate mobile fragment to target interaction distance ---
    direction = _get_fibonacci_sphere_point(conformer_idx, total_conformers)
    target_pt = base_pt + direction * target_dist
    translation = target_pt - mobile_pt_rotated
    for atom_idx in range(mobile_mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(atom_idx))
        conf.SetAtomPosition(atom_idx, (pos + translation).tolist())

    combined = Chem.CombineMols(base_mol, mobile_mol)
    combined = Chem.Mol(combined)
    Chem.SanitizeMol(combined)
    return combined


def rdkit_aggregate_mol(smi_parts, seed, n_conformers):
    """Generate multiple 3D conformers of a molecular aggregate from SMILES fragments.
    
    Embeds each fragment independently with multiple conformers, then builds 
    combinations of per-fragment conformer indices to produce diverse starting 
    geometries. Fragments are positioned using interaction-point pairing 
    and Fibonacci-sphere sampling.

    Args:
        smi_parts (list): List of SMILES strings, one per fragment
        seed (int): Base random seed for reproducibility
        n_conformers (int): Maximum number of aggregate conformers to build

    Returns:
        rdkit.Chem.rdchem.Mol: Combined molecule with multiple conformers
    """
    n_frags = len(smi_parts)
    # Determine the number of conformers per fragment to ensure the 
    # ensemble covers the requested sample size
    confs_per_frag = max(1, int(np.ceil(n_conformers ** (1.0 / n_frags))))

    # Embed each fragment with multiple conformers independently
    fragments = []
    for smi in smi_parts:
        mol = _rdChem.MolFromSmiles(smi)
        if mol is None:
            return None
        mol = _rdChem.AddHs(mol)
        cids = Chem.EmbedMultipleConfs(mol, numConfs=confs_per_frag, randomSeed=seed)
        # Fallback to random coordinates if standard distance geometry fails
        if len(cids) == 0:
            cids = Chem.EmbedMultipleConfs(
                mol, numConfs=confs_per_frag, randomSeed=seed, useRandomCoords=True
            )
        if len(cids) == 0:
            return None
        fragments.append(mol)

    # Sort by heavy-atom count: largest fragment remains stationary as the base
    fragments.sort(key=lambda x: x.GetNumHeavyAtoms(), reverse=True)

    conformer_mols = []
    for conf_idx in range(n_conformers):
        # Create a single-conformer copy for each fragment by cycling 
        # through its available internal conformation pool
        frag_copies = []
        for frag in fragments:
            n_available = frag.GetNumConformers()
            selected_cid = conf_idx % n_available
            positions = frag.GetConformer(selected_cid).GetPositions()

            copy = _rdChem.RWMol(frag)
            copy.RemoveAllConformers()
            new_conf = _rdChem.Conformer(frag.GetNumAtoms())
            for atom_i, pos in enumerate(positions):
                new_conf.SetAtomPosition(atom_i, pos.tolist())
            copy.AddConformer(new_conf, assignId=True)
            frag_copies.append(copy)

        # Assemble the fragments using Fibonacci-sphere distribution
        combined = frag_copies[0]
        for mobile in frag_copies[1:]:
            combined = _position_fragment_on_base(
                combined, mobile, conf_idx, n_conformers, seed=seed
            )
        conformer_mols.append(combined)

    if not conformer_mols:
        return None

    # Merge all single-conformer results into a single multi-conformer object
    result_mol = Chem.RWMol(conformer_mols[0])
    for single_conf_mol in conformer_mols[1:]:
        result_mol.AddConformer(single_conf_mol.GetConformer(0), assignId=True)

    result_mol = result_mol.GetMol()
    try:
        Chem.SanitizeMol(result_mol)
    except Exception:
        pass
        
    result_mol.SetProp("AggregateSmiles", ".".join(smi_parts))
    return result_mol


def smi_to_mol(
    smi,
    program,
    log,
    seed,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    sample=25,
):
    """Convert SMILES to RDKit molecule with constraints handling.

    Creates an RDKit molecule from SMILES string, handling special cases like
    aggregates, complexes, transition states, and mapped atoms. Supports
    constraint application for conformer generation.

    Args:
        smi (str): SMILES string to convert
        program (str): Program to use for conformer generation
        log: Logger object for status messages
        seed (int): Random seed for reproducibility
        constraints_atoms (list): Atom-based constraints
        constraints_dist (list): Distance constraints
        constraints_angle (list): Angle constraints
        constraints_dihedral (list): Dihedral angle constraints
        sample (int, optional): Number of pre-built conformers for rdkit
            aggregates. Defaults to 25.

    Returns:
        tuple: (mol, constraints_atoms, constraints_dist, constraints_angle,
               constraints_dihedral, complex_ts) where:
            - mol (rdkit.Chem.rdchem.Mol): Generated molecule or None if failed
            - constraints_*: Updated constraint lists
            - complex_ts (bool): True if molecule is a complex/TS

    Note:
        For multi-fragment SMILES (containing '.') program='rdkit' embeds
        fragments independently and positions them via interaction-point
        pairing and Fibonacci-sphere sampling (complex_ts=False, so all
        downstream RDKit minimisation/filtering runs normally).
        program='crest' uses nci_ts_mol (complex_ts=True).
    """
    complex_ts = False
    smi_parts = smi.split(".")

    # Handle multi-fragment SMILES (aggregates / NCI complexes / TSs)
    if len(smi_parts) > 1:
        if program == "rdkit":
            log.write(
                f"\no  Building RDKit aggregate conformers for "
                f"{len(smi_parts)}-fragment SMILES"
            )
            mol = rdkit_aggregate_mol(smi_parts, seed, sample)
            if mol is None:
                log.write(
                    f"\nx  RDKit embedding failed for one or more fragments in "
                    f"'{smi}'. Try program='crest' or check your SMILES.\n"
                )
            # complex_ts stays False: the pre-built multi-conformer mol goes
            # straight through the normal rdkit pipeline (embed_conf bypassed,
            # minimize_rdkit_energy / conformer_filters run as usual)
            # Translate atom map numbers to indices in the combined molecule
            # (mirrors the same translation done for single molecules below)
            if mol is not None:
                map_to_idx = {
                    atom.GetAtomMapNum(): atom.GetIdx()
                    for atom in mol.GetAtoms()
                    if atom.GetAtomMapNum() > 0
                }
                if map_to_idx:
                    constraints_atoms = _translate_constraint_indices(
                        constraints_atoms, map_to_idx, n_indices=1, log=log)
                    constraints_dist = _translate_constraint_indices(
                        constraints_dist, map_to_idx, n_indices=2, log=log)
                    constraints_angle = _translate_constraint_indices(
                        constraints_angle, map_to_idx, n_indices=3, log=log)
                    constraints_dihedral = _translate_constraint_indices(
                        constraints_dihedral, map_to_idx, n_indices=4, log=log)
            constraints = [constraints_atoms, constraints_dist,
                           constraints_angle, constraints_dihedral]

        elif program == "crest":
            # Use CREST-specific nci_ts_mol for NCI complexes and TSs
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
            log.write(
                f"\nx  {program} not supported for conformer generation of "
                f"aggregates (your SMILES has {len(smi_parts)} parts separated "
                "by a period). Use program='rdkit' or program='crest'."
            )
            sys.exit()

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

            # Build map_num -> atom_idx dictionary from mapped atoms
            map_to_idx = {}
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num > 0:
                    map_to_idx[map_num] = atom.GetIdx()

            # Translate constraint indices if mapped atoms exist
            has_constraints = any([constraints_atoms, constraints_dist,
                                   constraints_angle, constraints_dihedral])
            if map_to_idx:
                constraints_atoms = _translate_constraint_indices(
                    constraints_atoms, map_to_idx, n_indices=1, log=log
                )
                constraints_dist = _translate_constraint_indices(
                    constraints_dist, map_to_idx, n_indices=2, log=log
                )
                constraints_angle = _translate_constraint_indices(
                    constraints_angle, map_to_idx, n_indices=3, log=log
                )
                constraints_dihedral = _translate_constraint_indices(
                    constraints_dihedral, map_to_idx, n_indices=4, log=log
                )
            elif has_constraints:
                log.write(
                    "\nx  Constraints were specified but no atom map numbers were found "
                    "in the molecule. Constraint indices must correspond to atom map numbers "
                    "(e.g. [C:1][N:2]). Stopping."
                )
                sys.exit()

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


def detect_haptic_rings(mol, metal_idx):
    """Detect 5- and 6-member non-metal rings bound to one metal through 2+ atoms."""
    haptic_rings = []
    metal_idx = {idx for idx in (metal_idx or []) if idx is not None}

    if not metal_idx:
        return haptic_rings

    for ring_order, ring in enumerate(mol.GetRingInfo().AtomRings()):
        raw_ring_atoms = list(ring)
        ring_atoms = [atom_idx for atom_idx in raw_ring_atoms if atom_idx not in metal_idx]
        ring_size = len(ring_atoms)
        if ring_size not in [5, 6]:
            continue

        for idx in metal_idx:
            linked_atoms = [
                atom_idx for atom_idx in ring_atoms
                if mol.GetBondBetweenAtoms(atom_idx, idx) is not None
            ]
            if len(linked_atoms) >= 2:
                haptic_rings.append(
                    {
                        "ring_atoms": ring_atoms,
                        "metal_idx": idx,
                        "ring_size": ring_size,
                        "ring_order": ring_order,
                        "raw_ring_atoms": raw_ring_atoms,
                    }
                )
                break

    return haptic_rings


def select_haptic_rings(haptic_rings):
    """Keep non-overlapping haptic rings, preferring the largest shared ring."""
    selected = []
    used_atoms = set()
    candidates = sorted(
        haptic_rings,
        key=lambda ring: (-ring["ring_size"], ring["ring_order"]),
    )

    for ring in candidates:
        ring_atom_set = set(ring["ring_atoms"])
        if ring_atom_set.isdisjoint(used_atoms):
            selected.append(ring)
            used_atoms.update(ring_atom_set)

    return sorted(selected, key=lambda ring: ring["ring_order"])


def label_haptic_ring_atoms(mol, haptic_rings, metal_idx=None):
    """Label only haptic ring atoms with 1000, 1001... then 2000, 2001..."""
    overwritten_maps = []
    metal_idx = set(metal_idx or [])

    for ring_idx, ring in enumerate(haptic_rings, start=1):
        map_base = ring_idx * 1000
        for atom_offset, atom_idx in enumerate(ring["ring_atoms"]):
            if atom_idx in metal_idx:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            previous_map = atom.GetAtomMapNum()
            if previous_map:
                overwritten_maps.append(previous_map)
            atom.SetAtomMapNum(map_base + atom_offset)

    return overwritten_maps


def generate_haptic_constraints(ring_atoms, ring_size, mapped_idxs=None):
    """
    Generate haptic ring constraints with real RDKit atom indices.

    The 1000+ atom map labels are only user-facing labels and are never used
    as force-field atom indices. The optional mapped_idxs argument is accepted
    for compatibility with older internal calls.
    """
    params = {
        5: {"distance": 1.42, "angle": 108.0, "dihedrals": 3},
        6: {"distance": 1.40, "angle": 120.0, "dihedrals": 4},
    }
    if ring_size not in params:
        return [], [], []

    ring_params = params[ring_size]
    constraints_dist = []
    constraints_angle = []
    constraints_dihedral = []

    for i in range(ring_size):
        constraints_dist.append(
            [ring_atoms[i], ring_atoms[(i + 1) % ring_size], ring_params["distance"]]
        )
        constraints_angle.append(
            [
                ring_atoms[i],
                ring_atoms[(i + 1) % ring_size],
                ring_atoms[(i + 2) % ring_size],
                ring_params["angle"],
            ]
        )

    for i in range(ring_params["dihedrals"]):
        constraints_dihedral.append(
            [
                ring_atoms[i],
                ring_atoms[(i + 1) % ring_size],
                ring_atoms[(i + 2) % ring_size],
                ring_atoms[(i + 3) % ring_size],
                0.0,
            ]
        )

    return constraints_dist, constraints_angle, constraints_dihedral


def _copy_constraints(constraints):
    """Copy constraint lists so automatic constraints do not leak between jobs."""
    if not constraints:
        return []
    copied_constraints = []
    for constraint in constraints:
        if isinstance(constraint, (list, tuple)):
            copied_constraints.append(list(constraint))
        else:
            copied_constraints.append(constraint)
    return copied_constraints


def apply_haptic_constraints(
    mol,
    metal_idx,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    log=None,
):
    """Detect haptic rings, label them, and append real-index constraints."""
    constraints_atoms = _copy_constraints(constraints_atoms)
    constraints_dist = _copy_constraints(constraints_dist)
    constraints_angle = _copy_constraints(constraints_angle)
    constraints_dihedral = _copy_constraints(constraints_dihedral)

    detected_haptic_rings = detect_haptic_rings(mol, metal_idx)
    haptic_rings = select_haptic_rings(detected_haptic_rings)
    if not haptic_rings:
        return constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral

    overwritten_maps = label_haptic_ring_atoms(mol, haptic_rings, metal_idx)

    aromatic_constraints_dist = []
    aromatic_constraints_angle = []
    aromatic_constraints_dihedral = []

    for ring in haptic_rings:
        ring_labels = [
            mol.GetAtomWithIdx(atom_idx).GetAtomMapNum()
            for atom_idx in ring["ring_atoms"]
        ]
        label_dist, label_angle, label_dihedral = generate_haptic_constraints(
            ring_labels, ring["ring_size"]
        )
        aromatic_constraints_dist.extend(label_dist)
        aromatic_constraints_angle.extend(label_angle)
        aromatic_constraints_dihedral.extend(label_dihedral)

        dist, angle, dihedral = generate_haptic_constraints(
            ring["ring_atoms"], ring["ring_size"]
        )
        constraints_dist.extend(dist)
        constraints_angle.extend(angle)
        constraints_dihedral.extend(dihedral)

    if log is not None:
        log.write(
            "\no  AQME detected haptic ring binding to a metal and is applying "
            "automatic constraints for the hapticity link. Hydrogens should be "
            "explicit. The ring atoms do not need to be mapped by the user; "
            "AQME labels only the ring atoms automatically as 1000, 1001, etc."
        )
        if len(haptic_rings) < len(detected_haptic_rings):
            log.write(
                "\no  Overlapping haptic rings were detected; AQME kept the largest "
                "ring for shared atoms."
            )
        log.write(
            f"\no  Automatic hapticity constraints: --constraints_dist "
            f"\"{aromatic_constraints_dist}\" --constraints_angle "
            f"\"{aromatic_constraints_angle}\" --constraints_dihedral "
            f"\"{aromatic_constraints_dihedral}\""
        )
        if overwritten_maps:
            log.write(
                "\no  Existing atom-map numbers on haptic ring atoms were replaced "
                "with 1000+ labels, so the previous labels were lost."
            )

    return constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral


def apply_rdkit_constraints(forcefield, constraints_atoms, constraints_dist, 
                             constraints_angle, constraints_dihedral, force_constant=100000.0):
    """Apply constraints to an RDKit force field object (MMFF by default, fallback to UFF)."""
    # Detect FF type once to avoid repeated attribute lookups in loops
    is_mmff = hasattr(forcefield, "MMFFAddDistanceConstraint")

    for constraint in constraints_atoms:
        atom_idx = constraint[0] if isinstance(constraint, (list, tuple)) else constraint
        forcefield.AddFixedPoint(int(atom_idx))

    for constraint in constraints_dist:
        atom1, atom2, distance = int(constraint[0]), int(constraint[1]), float(constraint[2])
        if is_mmff:
            forcefield.MMFFAddDistanceConstraint(atom1, atom2, False, distance, distance, force_constant)
        else:
            forcefield.UFFAddDistanceConstraint(atom1, atom2, False, distance, distance, force_constant)

    for constraint in constraints_angle:
        atom1, atom2, atom3 = int(constraint[0]), int(constraint[1]), int(constraint[2])
        angle = float(constraint[3])
        if is_mmff:
            forcefield.MMFFAddAngleConstraint(atom1, atom2, atom3, False, angle, angle, force_constant)
        else:
            forcefield.UFFAddAngleConstraint(atom1, atom2, atom3, False, angle, angle, force_constant)

    for constraint in constraints_dihedral:
        atom1, atom2, atom3, atom4 = int(constraint[0]), int(constraint[1]), int(constraint[2]), int(constraint[3])
        dihedral = float(constraint[4])
        if is_mmff:
            forcefield.MMFFAddTorsionConstraint(atom1, atom2, atom3, atom4, False, dihedral, dihedral, force_constant)
        else:
            forcefield.UFFAddTorsionConstraint(atom1, atom2, atom3, atom4, False, dihedral, dihedral, force_constant)

def _translate_constraint_indices(constraints, map_to_idx, n_indices, log=None):
    """Translate atom map numbers to RDKit atom indices in constraints.
    
    Args:
        constraints (list): List of constraints, each a list of values
            where the first n_indices values are atom identifiers
        map_to_idx (dict): Mapping from atom map numbers to atom indices
        n_indices (int): Number of atom indices at the start of each constraint
        log: Logger object for warnings (optional)
        
    Returns:
        list: Constraints with translated atom indices
        
    Notes:
        - If map_to_idx is empty, the code exits with an error because constraint
          indices must correspond to mapped atoms in the SMILES.
        - If any constraint index is not found as a map number, the code exits
          with an error indicating which index is invalid.
    """
    if not constraints:
        return constraints

    if not map_to_idx:
        if log is not None:
            log.write(
                "\nx  WARNING! Constraints were specified but no atom mapping was found "
                "in the SMILES. Constraint indices must correspond to mapped atoms "
                "(e.g. use [C:1], [N:2] in the SMILES and pass [1,2,...] as constraints)."
            )
            if hasattr(log, "finalize"):
                log.finalize()
        sys.exit()

    translated = []
    for constraint in constraints:
        if n_indices == 1 and not isinstance(constraint, (list, tuple)):
            new_constraint = [constraint]
        else:
            new_constraint = list(constraint)
        for i in range(n_indices):
            original = int(new_constraint[i])
            if original not in map_to_idx:
                if log is not None:
                    log.write(
                        f"\nx  WARNING! Constraint index {original} is not correct. It does not "
                        f"correspond to any mapped atom number in the molecule. "
                        f"All constraint indices must match atom map numbers "
                        f"(e.g. [C:1], [N:2]). Available map numbers are: "
                        f"{sorted(map_to_idx.keys())}."
                    )
                    if hasattr(log, "finalize"):
                        log.finalize()
                sys.exit()
            new_constraint[i] = map_to_idx[original]
        translated.append(new_constraint)
    return translated
