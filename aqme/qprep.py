"""
Parameters
----------
   files : mol object, str or list of str, default=None
      This module prepares input QM file(s). Formats accepted: mol object(s), 
      Gaussian or ORCA LOG/OUT output files, JSON, XYZ, SDF, PDB. Also, 
      lists can be used (i.e. [FILE1.log, FILE2.log] or \*.FORMAT such as \*.json).
   atom_types : list of str, default=[]
      (If files is None) List containing the atoms of the system
   cartesians : list of str, default=[]
      (If files is None) Cartesian coordinates used for further processing
   w_dir_main : str, default=os.getcwd()
      Working directory
   destination : str, default=None,
      Directory to create the input file(s)
   varfile : str, default=None
      Option to parse the variables using a yaml file (specify the filename)
   program : str, default=None
      Program required to create the new input files. Current options: 'gaussian', 'orca'
   qm_input : str, default=''
      Keywords line for new input files (i.e. 'B3LYP/6-31G opt freq')
   qm_end : str, default=''
      Final line(s) in the new input files
   charge : int, default=None
      Charge of the calculations used in the following input files. If charge isn't defined, it defaults to 0
   mult : int, default=None
      Multiplicity of the calculations used in the following input files. If mult isn't defined, it defaults to 1
   suffix : str, default=''
      Suffix for the new input files (i.e. FILENAME_SUFFIX.com for FILENAME.log)
   prefix : str, default=''  
      Prefix added to all the names  
   chk : bool, default=False
      Include the chk input line in new input files for Gaussian calculations
   oldchk : bool, default=False
      Include the oldchk input line in new input files for Gaussian calculations
   chk_path : str, default=''
      PATH to store CHK files. For example, if chk_path='root/user/FILENAME.chk, the chk line of the input file would be
      %chk=root/user/FILENAME.chk
   oldchk_path : str, default=''
      PATH to read CHK files with %oldchk. For example, if oldchk_path='root/user/FILENAME.chk, the oldchk line of the input file would be
      %oldchk=root/user/FILENAME.chk
   mem : str, default='4GB'
      Memory for the QM calculations (i) Gaussian: total memory; (ii) ORCA: memory per processor
   nprocs : int, default=None
      Number of processors used in the QM calculations
   gen_atoms : list of str, default=[]
      Atoms included in the gen(ECP) basis set (i.e. ['I','Pd'])
   bs_gen : str, default=''
      Basis set used for gen(ECP) atoms (i.e. 'def2svp')
   bs_nogen : str, default=''
      Basis set used for non gen(ECP) atoms in gen(ECP) calculations (i.e. '6-31G*')
   lowest_only : bool, default=False
      Only create input for the conformer with lowest energy of the SDF file
   lowest_n : int, default=None
      Only create inputs for the n conformers with lowest energy of the SDF file
   e_threshold_qprep : float, default=None
      Only create inputs for conformers below the energy threshold (to the lowest conformer)
      of the SDF file
"""
######################################################.
#        This file stores the QPREP class            #
######################################################.

import os
import subprocess
import sys
import glob
import time
import json
import pandas as pd

from aqme.utils import (
    cclib_atoms_coords,
    read_file,
    move_file,
    load_variables,
    read_xyz_charge_mult,
    mol_from_sdf_or_mol_or_mol2,
    add_prefix_suffix,
    check_files,
    check_dependencies,
    set_destination,
    periodic_table
)

from aqme.csearch.crest import xyzall_2_xyz
from pathlib import Path
from rdkit import Chem

from importlib.resources import files
TEMPLATES_PATH = files("aqme").joinpath("templates")

class qprep:
    """A class for preparing quantum chemistry input files.
    
    This class handles the creation of input files for various quantum chemistry 
    programs (currently Gaussian and ORCA) from multiple input formats. It supports
    file conversion, molecular property extraction, and conformer processing.
    
    Attributes:
        args: Configuration object containing user parameters and settings
        
    Examples:
        >>> prep = QPrep(program='gaussian', qm_input='B3LYP/6-31G opt freq')
        >>> prep.process_files(['molecule.xyz'])
    """
    
    SUPPORTED_FORMATS = ['sdf', 'xyz', 'pdb', 'log', 'out', 'json']
    SUPPORTED_PROGRAMS = ['gaussian', 'orca']
    DEFAULT_NPROCS = 8
    
    def __init__(self, create_dat=True, **kwargs):
        """Initialize QPrep with provided configuration.
        
        Args:
            create_dat (bool, optional): Whether to create data files. Defaults to True.
            **kwargs: Configuration parameters for QM calculations.
            
        Raises:
            SystemExit: If invalid configuration or missing required parameters.
        """
        self.start_time = time.time()
        self.create_dat = create_dat
        
        # Initialize configuration
        self.args = load_variables(kwargs, "qprep", create_dat=create_dat)
        self._validate_initial_setup()
        
        # Set up working environment
        self._setup_environment()
        
        # Process input files
        self._process_input_files()
        
        # Finalize
        if create_dat:
            self._finalize()
            
    def _validate_initial_setup(self):
        """Validate initial configuration and dependencies.
        
        Raises:
            SystemExit: If validation fails
        """
        # Check dependencies and files
        check_dependencies(self)
        check_files(self, 'qprep')
        
        # Validate file format
        file_format = os.path.basename(Path(self.args.files[0])).split('.')[-1].lower()
        if file_format not in self.SUPPORTED_FORMATS:
            self.args.log.write(
                f"\nx  Format {file_format} not supported! Accepted: {', '.join(self.SUPPORTED_FORMATS)}"
            )
            self.args.log.finalize()
            sys.exit()
            
        # Validate QM program
        if not self.args.program or self.args.program.lower() not in self.SUPPORTED_PROGRAMS:
            self.args.log.write(
                '\nx  Specify supported program: "gaussian" or "orca"'
            )
            self.args.log.finalize()
            sys.exit()
            
        # Validate QM input
        if self.args.qm_input == "" and self.create_dat:
            self.args.log.write("x  No keywords specified (qm_input=KEYWORDS_LINE).")
            self.args.log.finalize()
            sys.exit()
            
        # Validate gen/genecp setup
        self._validate_gen_basis()
            
    def _validate_gen_basis(self):
        """Validate basis set specifications for gen/genecp calculations.
        
        Raises:
            SystemExit: If required basis sets are missing
        """
        if not self.args.gen_atoms:
            return
            
        if self.args.bs_nogen == "" and self.create_dat:
            self.args.log.write(
                "x  Gen(ECP) atoms specified but missing bs_nogen=BASIS_SET"
            )
            self.args.log.finalize()
            sys.exit()
            
        if self.args.bs_gen == "" and self.create_dat:
            self.args.log.write(
                "x  Gen(ECP) atoms specified but missing bs_gen=BASIS_SET"
            )
            self.args.log.finalize()
            sys.exit()
            
    def _setup_environment(self):
        """Configure working environment and parameters."""
        # Set processors
        self.args.nprocs = self.args.nprocs or self.DEFAULT_NPROCS
        
        # Set destination
        self.destination = set_destination(self, 'QCALC')
        
        # Validate theory level for Gaussian
        if self.args.program.lower() == 'gaussian' and self.create_dat:
            self.check_level_of_theory()
            
    def _process_input_files(self):
        """Process all input files and generate QM inputs."""
        for file in self.args.files:
            file_format = os.path.basename(Path(file)).split('.')[-1].lower()
            name = os.path.basename(Path(file)).split('.')[0]
            
            if file_format in ['sdf', 'xyz', 'pdb']:
                self._process_molecular_file(file, file_format, name)
            else:
                self._process_qm_file(file, name)
                
    def _process_molecular_file(self, file, file_format, name):
        """Process molecular structure files (SDF/XYZ/PDB).
        
        Args:
            file (str): Input file path
            file_format (str): File format
            name (str): Base filename
        """
        sdf_files = self._convert_to_sdf(file, file_format, name)
        
        for sdf_file in sdf_files:
            try:
                self.sdf_2_com(sdf_file, self.destination, file_format)
                if self.create_dat:
                    self.args.log.write(f"o  {name} processed at {self.destination}")
            except OSError:
                self.args.log.write(f"x  {name} couldn't be processed!")
                continue
                
            # Cleanup temporary files
            if file_format in ['xyz', 'pdb']:
                os.remove(sdf_file)
                
    def _process_qm_file(self, file, name):
        """Process quantum chemistry output files.
        
        Args:
            file (str): Input file path
            name (str): Base filename
        """
        atom_types, cartesians, charge, mult, found_coords = self.qprep_coords(
            file, None, file.split('.')[-1]
        )
        
        if not found_coords:
            return
            
        qprep_data = {
            "atom_types": atom_types,
            "cartesians": cartesians,
            "charge": charge,
            "mult": mult,
            "name": name,
        }
        
        comfile = self.write(qprep_data)
        move_file(self.destination, self.args.w_dir_main, comfile)
        
        if self.create_dat:
            self.args.log.write(f"o  {name} processed at {self.destination}")
            
    def _convert_to_sdf(self, file, file_format, name):
        """Convert input files to SDF format.
        
        Args:
            file (str): Input file path
            file_format (str): File format
            name (str): Base filename
            
        Returns:
            list: Paths to generated SDF files
        """
        sdf_files = []
        
        if file_format == 'xyz':
            sdf_files = self._convert_xyz_to_sdf(file, name)
        elif file_format == 'pdb':
            sdf_files = self._convert_pdb_to_sdf(file)
        else:  # Already SDF
            sdf_files = [file]
            
        return sdf_files
        
    def _convert_xyz_to_sdf(self, file, name):
        """Convert XYZ file to SDF format.
        
        Args:
            file (str): XYZ file path
            name (str): Base filename
            
        Returns:
            list: Paths to generated SDF files
        """
        sdf_files = []
        xyzall_2_xyz(file, f"{self.args.w_dir_main}/{name}")
        
        for conf_file in glob.glob(f"{self.args.w_dir_main}/{name}_conf_*.xyz"):
            charge = self.args.charge if self.args.charge is not None else read_xyz_charge_mult(conf_file)[0]
            mult = self.args.mult if self.args.mult is not None else read_xyz_charge_mult(conf_file)[1]
            
            # Convert to SDF using OpenBabel
            sdf_path = f"{conf_file.split('.xyz')[0]}.sdf"
            command = [
                "obabel", "-ixyz", conf_file, "-osdf", f"-O{sdf_path}",
                "--property", f"Real charge={charge}", ";", f"Mult={mult}"
            ]
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            sdf_files.append(sdf_path)
            os.remove(conf_file)  # Clean up XYZ file
            
        return sdf_files
        
    def _convert_pdb_to_sdf(self, file):
        """Convert PDB file to SDF format.
        
        Args:
            file (str): PDB file path
            
        Returns:
            list: Path to generated SDF file
        """
        sdf_path = f"{file.split('.pdb')[0]}.sdf"
        command = ["obabel", "-ipdb", file, "-osdf", f"-O{sdf_path}"]
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return [sdf_path]
        
    def _finalize(self):
        """Finalize processing and log completion."""
        elapsed_time = round(time.time() - self.start_time, 2)
        self.args.log.write(f"\nTime QPREP: {elapsed_time} seconds\n")
        self.args.log.finalize()

    def sdf_2_com(self, sdf_file, destination, file_format):
        """Process SDF file to generate QM input files.
        
        Handles conversion of SDF structures to Gaussian/ORCA input files,
        with support for conformer selection and energy filtering.
        
        Args:
            sdf_file (str): Path to input SDF file
            destination (Path): Output directory path
            file_format (str): Original file format
            
        Note:
            Uses lowest_only, lowest_n and e_threshold_qprep settings to
            filter conformers if specified.
        """
        sdf_name = os.path.basename(Path(sdf_file)).split(".")[0]
        
        # Determine conformer filtering criteria
        low_check = self._get_conformer_filter()
        
        # Get filtered molecule list
        mols = mol_from_sdf_or_mol_or_mol2(
            sdf_file, "qprep", self.args, low_check=low_check
        )
        
        # Process each molecule/conformer
        for i, mol in enumerate(mols):
            self._process_single_mol(
                mol, sdf_file, file_format, sdf_name, i, destination
            )
            
    def _get_conformer_filter(self):
        """Determine conformer filtering method.
        
        Returns:
            Union[str, int, float, None]: Filter criteria
        """
        if self.args.lowest_only:
            return 'lowest_only'
        if self.args.lowest_n is not None:
            return int(self.args.lowest_n)
        if self.args.e_threshold_qprep is not None:
            return float(self.args.e_threshold_qprep)
        return None
        
    def _process_single_mol(self, mol, sdf_file, file_format, sdf_name, idx, destination):
        """Process single molecule/conformer from SDF.
        
        Args:
            mol (rdkit.Mol): RDKit molecule object
            sdf_file (str): Source SDF file path
            file_format (str): Original file format
            sdf_name (str): Base name from SDF
            idx (int): Molecule/conformer index
            destination (Path): Output directory
        """
        # Extract molecular data
        atom_types, cartesians, charge, mult, _ = self.qprep_coords(
            sdf_file, mol, file_format
        )
        
        # Generate conformer name
        name_conf = (
            f"{sdf_name}_conf_{idx+1}" if "_conf_" not in sdf_name 
            else sdf_name
        )
        
        # Create and write input file
        qprep_data = {
            "atom_types": atom_types,
            "cartesians": cartesians,
            "charge": charge,
            "mult": mult,
            "name": name_conf,
        }
        
        comfile = self.write(qprep_data)
        move_file(destination, self.args.w_dir_main, comfile)

    def get_header(self, qprep_data):
        """Generate input file header section.
        
        Creates program-specific header with resource specs and calculation setup.
        
        Args:
            qprep_data (dict): Molecular data including name and properties
            
        Returns:
            str: Formatted header text
        """
        name_file = add_prefix_suffix(qprep_data["name"], self.args)
        
        if self.args.program.lower() == "gaussian":
            return self._get_gaussian_header(name_file, qprep_data)
        else:
            return self._get_orca_header(name_file, qprep_data)
            
    def _get_gaussian_header(self, name_file, qprep_data):
        """Generate Gaussian-specific header section.
        
        Args:
            name_file (str): Base filename
            qprep_data (dict): Molecular data
            
        Returns:
            str: Formatted Gaussian header
        """
        header = []
        
        # Add checkpoint handling
        if self.args.chk_path:
            header.append(f'%chk={self.args.chk_path}')
        elif self.args.chk:
            header.append(f'%chk={name_file}.chk')
            
        if self.args.oldchk_path:
            header.append(f'%oldchk={self.args.oldchk_path}')
        elif self.args.oldchk:
            header.append(f'%oldchk={name_file}.chk')
            
        # Add resource specs
        header.extend([
            f"%nprocshared={self.args.nprocs}",
            f"%mem={self.args.mem}"
        ])
        
        # Add calculation setup
        if self.args.qm_input[:2] not in ['p ','P ']:
            header.append(f"# {self.args.qm_input}")
        else:  # Handle #p directive
            header.append(f"#{self.args.qm_input}")
            
        # Add title and charge/mult
        header.extend([
            "",
            name_file,
            "",
            f'{qprep_data["charge"]} {qprep_data["mult"]}'
        ])
        
        return "\n".join(header)
        
    def _get_orca_header(self, name_file, qprep_data):
        """Generate ORCA-specific header section.
        
        Args:
            name_file (str): Base filename
            qprep_data (dict): Molecular data
            
        Returns:
            str: Formatted ORCA header
        """
        header = [f'# {name_file}']
        
        # Handle memory specification
        mem_orca = self._convert_memory_to_orca()
        if '%maxcore' not in self.args.qm_input:
            header.append(f"%maxcore {mem_orca}")
            
        # Handle parallel processing
        if not self._has_pal_directive():
            header.append(f"%pal nprocs {self.args.nprocs} end")
            
        # Add calculation setup
        header.extend([
            f"! {self.args.qm_input}",
            f'* xyz {qprep_data["charge"]} {qprep_data["mult"]}'
        ])
        
        return "\n".join(header)
        
    def _convert_memory_to_orca(self):
        """Convert memory specification to ORCA format.
        
        Returns:
            str: Memory value in ORCA format
        """
        if "GB" in self.args.mem:
            return str(int(self.args.mem.split("GB")[0]) * 1000)
        elif "MB" in self.args.mem:
            return self.args.mem.split("MB")[0]
        elif "MW" in self.args.mem:
            return self.args.mem.split("MW")[0]
        return self.args.mem
        
    def _has_pal_directive(self):
        """Check if ORCA parallel directive is present.
        
        Returns:
            bool: True if PAL directive found
        """
        pal_list = ['%pal','pal1','pal3','pal3','pal4','pal5','pal6','pal7','pal8']
        return any(
            kw.rstrip("\n").lower() in pal_list 
            for kw in self.args.qm_input.split()
        )


    def get_tail(self, qprep_data):
        """Generate input file tail section.
        
        Creates program-specific closing section with basis sets and modifiers.
        
        Args:
            qprep_data (dict): Molecular data
            
        Returns:
            str: Formatted tail section
        """
        if self.args.program.lower() != "gaussian":
            return ""
            
        parts = []
        modifysph_line = ""
        
        # Process QM end section
        if self.args.qm_end:
            qm_end_local = self.args.qm_end
            if "modifysph" in qm_end_local.lower():
                qm_end_local, modifysph_line = self._extract_modifysph(qm_end_local)
            parts.append(f"{qm_end_local}\n")
            
        # Handle basis sets
        if self.args.gen_atoms:
            basis_section = self._generate_basis_section(qprep_data["atom_types"])
            if basis_section:
                if self.args.qm_end:
                    parts.append("\n")
                parts.append(basis_section)
                
        # Combine and format
        txt = "\n".join(part.strip() for part in parts if part)
        if txt and modifysph_line:
            txt = f"{txt}\n{modifysph_line}"
            
        return txt
        
    def _extract_modifysph(self, qm_end):
        """Extract modifysph section from QM end text.
        
        Args:
            qm_end (str): Original QM end text
            
        Returns:
            tuple: (remaining_text, modifysph_section)
        """
        lines = qm_end.split("\n")
        for idx, line in enumerate(lines):
            if "modifysph" in line.lower():
                # Get modifysph section
                modifysph = f"{line}\n\n{lines[idx+1]}\n\n"
                # Remove section from original
                del lines[idx:idx+2]
                return "\n".join(lines), modifysph
        return qm_end, ""
        
    def _generate_basis_section(self, atom_types):
        """Generate basis set specification section.
        
        Args:
            atom_types (List[str]): Atomic symbols
            
        Returns:
            str: Formatted basis set section
        """
        ecp_used = []
        ecp_not_used = []
        gen_type = "genecp" if "genecp" in self.args.qm_input.lower() else "gen"
        
        # Categorize atoms
        for element in atom_types:
            if element in self.args.gen_atoms and element not in ecp_used:
                ecp_used.append(element)
            elif element not in self.args.gen_atoms and element not in ecp_not_used:
                ecp_not_used.append(element)
                
        # Build basis set blocks
        blocks = []
        
        # Non-ECP atoms
        if ecp_not_used:
            blocks.append(
                f"{' '.join(ecp_not_used)} 0\n{self.args.bs_nogen}\n****"
            )
            
        # ECP atoms 
        if ecp_used:
            blocks.append(
                f"{' '.join(ecp_used)} 0\n{self.args.bs_gen}\n****"
            )
            
            # Additional block for GENECP
            if gen_type == "genecp":
                blocks.append(
                    f"\n{' '.join(ecp_used)} 0\n{self.args.bs_gen}"
                )
                
        return "\n".join(blocks) if blocks else ""
        
    def write(self, qprep_data):
        """Write quantum chemistry input file.
        
        Creates Gaussian (.com) or ORCA (.inp) input files with molecular
        coordinates and calculation parameters.
        
        Args:
            qprep_data (dict): Molecular data including:
                - name (str): Base filename
                - atom_types (List[str]): Atomic symbols
                - cartesians (np.ndarray): XYZ coordinates
                - charge (int): Molecular charge
                - mult (int): Spin multiplicity
                
        Returns:
            str: Name of created input file
        """
        extension = "com" if self.args.program.lower() == "gaussian" else "inp"
        name_file = add_prefix_suffix(qprep_data["name"], self.args)
        comfile = f'{name_file}.{extension}'
        outpath = self.args.w_dir_main / comfile
        
        # Remove existing file
        if outpath.exists():
            outpath.unlink()
            
        # Generate file content
        header = self.get_header(qprep_data)
        coords = self._format_coordinates(qprep_data)
        tail = self.get_tail(qprep_data)
        
        # Write file
        with open(outpath, "w") as f:
            f.write(header)
            f.write("\n")
            f.write(coords)
            f.write("\n*" if self.args.program.lower() == "orca" else "\n\n")
            f.write(tail)
            f.write("\n\n")
            
        return comfile
        
    def _format_coordinates(self, qprep_data):
        """Format atomic coordinates for input file.
        
        Args:
            qprep_data (dict): Molecular structure data
            
        Returns:
            str: Formatted coordinate block
        """
        coords = []
        for i, (atom, xyz) in enumerate(zip(
            qprep_data["atom_types"], 
            qprep_data["cartesians"]
        )):
            line = "{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}".format(
                atom, xyz[0], xyz[1], xyz[2]
            )
            coords.append(line)
            
        return "\n".join(coords)

    def qprep_coords(self, file, mol, file_format):
        """Extract molecular coordinates and properties from various file formats.
        
        Retrieves atomic coordinates, atom types, charge and multiplicity from:
        - RDKit mol objects 
        - QM output files (LOG/OUT)
        - JSON files
        - Custom atom_types/cartesians lists
        
        Args:
            file (str): Input file path
            mol (rdkit.Mol, optional): RDKit molecule object
            file_format (str): Input file format
            
        Returns:
            tuple: (atom_types, cartesians, charge, mult, found_coords)
            - atom_types (List[str]): Atomic symbols
            - cartesians (np.ndarray): XYZ coordinates
            - charge (int): Molecular charge
            - mult (int): Spin multiplicity
            - found_coords (bool): Whether coordinates were found
        """
        found_coords = False
        charge, mult = None, None
        
        # Use provided atom lists if available
        charge, mult = self._get_charge_mult(charge, mult)
        if self.args.atom_types and self.args.cartesians:
            return (self.args.atom_types, self.args.cartesians, 
                   charge, mult, True)
                   
        # Process based on input type
        if mol is not None:
            atom_types, cartesians, charge, mult = self._extract_mol_data(mol)
        elif file_format in ['log', 'out']:
            atom_types, cartesians, charge, mult = self._extract_qm_data(file)
        elif file_format == 'json':
            atom_types, cartesians, charge, mult = self._extract_json_data(file)
        else:
            atom_types, cartesians = [], []
            
        # Validate and set defaults
        try:
            if not (atom_types and cartesians):
                self.args.log.write(
                    f"x  {file} missing coordinates/atom information"
                )
            else:
                found_coords = True
        except ValueError:  # Array comparison issue
            found_coords = True
            
        charge, mult = self._get_charge_mult(charge, mult)
        return atom_types, cartesians, charge, mult, found_coords
        
    def _extract_mol_data(self, mol):
        """Extract data from RDKit molecule object.
        
        Args:
            mol (rdkit.Mol): RDKit molecule
            
        Returns:
            tuple: (atom_types, cartesians, charge, mult)
        """
        # Get atom types with isotope handling
        atom_types = []
        for atom in mol.GetAtoms():
            isotope = atom.GetIsotope()
            symbol = atom.GetSymbol()
            atom_types.append(
                f"{symbol}(iso={isotope})" if isotope else symbol
            )
            
        # Get coordinates
        cartesians = mol.GetConformers()[0].GetPositions()
        
        # Get charge
        try:
            charge = int(mol.GetProp("Real charge"))
        except KeyError:
            charge = Chem.GetFormalCharge(mol)
            
        # Get multiplicity
        try:
            mult = int(mol.GetProp("Mult"))
        except KeyError:
            n_radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
            mult = int(n_radicals / 2 * 2 + 1)
            
        return atom_types, cartesians, charge, mult
        
    def _extract_qm_data(self, file):
        """Extract data from QM output file.
        
        Args:
            file (str): Path to LOG/OUT file
            
        Returns:
            tuple: (atom_types, cartesians, charge, mult)
        """
        # Read file
        base_dir = os.getcwd() if self.args.command_line else self.args.w_dir_main
        outlines = read_file(os.getcwd(), base_dir, file)
        
        # Detect program and find start
        program = None
        resume_line = 0
        for i, line in enumerate(outlines):
            if 'Gaussian, Inc.' in line:
                program = 'gaussian'
                resume_line = i
                break
            elif 'O   R   C   A' in line:
                program = 'orca'
                resume_line = i
                break
                
        if not program:
            return [], [], None, None
            
        # Get molecule info
        n_atoms = 0
        charge = mult = None
        if program == 'gaussian':
            n_atoms, charge, mult = self._parse_gaussian_output(
                outlines[resume_line:]
            )
            
        atom_types, cartesians = QM_coords(
            outlines, -1, n_atoms, program, ""
        )
        
        return atom_types, cartesians, charge, mult
        
    def _parse_gaussian_output(self, lines):
        """Parse Gaussian output file for molecular info.
        
        Args:
            lines (List[str]): File lines to parse
            
        Returns:
            tuple: (n_atoms, charge, mult)
        """
        n_atoms = 0
        charge = mult = None
        found_n_atoms = False
        resume_line = 0

        for i, line in enumerate(lines):
            if line.find("Gaussian, Inc."):
                program = "gaussian"
                resume_line = i
                break
            elif line[i].find("O   R   C   A"):
                program = "orca"
                resume_line = i
                break

        for i in range(resume_line, len(lines)):
            if program == "gaussian":
                # get charge and mult
                if lines[i].find("Charge = ") > -1:
                    charge = int(lines[i].split()[2])
                    mult = int(lines[i].split()[5].rstrip("\n"))
                # get number of atoms
                elif lines[i].find("Symbolic Z-matrix:") > -1:
                    for j in range(i + 2, len(lines)):
                        if len(lines[j].split()) > 0:
                            n_atoms += 1
                        else:
                            found_n_atoms = True
                            break
                elif found_n_atoms:
                    break
              
        return n_atoms, charge, mult
        
    def _extract_json_data(self, file):
        """Extract data from JSON file.
        
        Args:
            file (str): JSON file path
            
        Returns:
            tuple: (atom_types, cartesians, charge, mult)
        """
        try:
            with open(file) as f:
                data = json.load(f)
            atom_types, cartesians = cclib_atoms_coords(data, -1)
            charge = data['charge']
            mult = data['mult']
        except (AttributeError, KeyError, json.JSONDecodeError):
            return [], [], None, None
            
        return atom_types, cartesians, charge, mult
        
    def _get_charge_mult(self, charge, mult):
        """Get charge and multiplicity, using defaults if needed.
        
        Args:
            charge (int, optional): Molecular charge
            mult (int, optional): Spin multiplicity
            
        Returns:
            tuple: (charge, mult) with defaults applied
        """
        # Set charge
        if self.args.charge is not None:
            charge = self.args.charge
        elif charge is None:
            charge = 0
            
        # Set multiplicity  
        if self.args.mult is not None:
            mult = self.args.mult
        elif mult is None:
            mult = 1
            
        return charge, mult


    def check_level_of_theory(self):
        """Validate quantum chemistry method against known options.
        
        Cross checks the specified functional and basis set against predefined lists.
        Note that these lists are not exhaustive - missing items may still be valid.
        
        Raises:
            Warning: If functional or basis set not found in predefined lists
        """
        # Load reference data
        functional_list = self._load_reference_data('functionals.csv')
        basis_set_list = self._load_reference_data('basis_sets.csv')
        
        # Validate components
        found_func = self._validate_functionals(functional_list)
        found_basis = self._validate_basis_sets(basis_set_list)
        
        # Report warnings
        if not found_func:
            self.args.log.write(
                "x  WARNING! Verify functional correctness. "
                "If valid, please report to add to known functionals."
            )
            
        if not found_basis:
            self.args.log.write(
                "x  WARNING! Verify basis set(s) correctness. "
                "If valid, please report to add to known basis sets."
            )
            
    def _load_reference_data(self, filename):
        """Load and process reference data from CSV.
        
        Args:
            filename (str): Name of CSV file in templates
            
        Returns:
            list: Processed reference data list
        """
        csv_path = TEMPLATES_PATH.joinpath(filename)
        with csv_path.open("rb") as fh:
            data = pd.read_csv(fh)
            
        # Extract and clean data
        values = data[self.args.program].to_numpy().flatten()
        return [x for x in values if str(x) != 'nan']
        
    def _validate_functionals(self, functional_list):
        """Check if input contains known functionals.
        
        Args:
            functional_list (list): Reference list of functionals
            
        Returns:
            bool: True if known functional found
        """
        return any(
            self._is_known_keyword(subkey, functional_list)
            for keyword in self.args.qm_input.split()
            for subkey in keyword.split('/')
            if self._is_method_keyword(subkey)
        )
        
    def _validate_basis_sets(self, basis_set_list):
        """Check if basis sets are known.
        
        Args:
            basis_set_list (list): Reference list of basis sets
            
        Returns:
            bool: True if all basis sets are known
        """
        # Check gen/genecp basis sets if specified
        if self.args.bs_gen:
            if not (
                self._is_known_basis(self.args.bs_gen, basis_set_list) and
                self._is_known_basis(self.args.bs_nogen, basis_set_list)
            ):
                return False
                
        # Check basis sets in input line
        return any(
            self._is_known_basis(subkey, basis_set_list)
            for keyword in self.args.qm_input.split()
            for subkey in keyword.split('/')
            if self._is_method_keyword(subkey)
        )
        
    def _is_method_keyword(self, keyword):
        """Check if keyword is a method specification.
        
        Args:
            keyword (str): Input keyword
            
        Returns:
            bool: True if keyword is method-related
        """
        excluded = ['opt', 'freq', 'scrf', 'pop', 'gen']
        return not any(keyword.count(x) > 0 for x in excluded)
        
    def _is_known_keyword(self, keyword, reference_list):
        """Check if keyword exists in reference list.
        
        Args:
            keyword (str): Keyword to check
            reference_list (list): List of known keywords
            
        Returns:
            bool: True if keyword found in list
        """
        return keyword.upper() in (x.upper() for x in reference_list)
        
    def _is_known_basis(self, basis, basis_list):
        """Check if basis set exists in reference list.
        
        Args:
            basis (str): Basis set name
            basis_list (list): List of known basis sets
            
        Returns:
            bool: True if basis set found in list
        """
        if not basis:
            return False
        return basis.upper() in (x.upper() for x in basis_list)


def QM_coords(outlines, min_RMS, n_atoms, program, keywords_line):
    """Extract atomic coordinates from QM output files.
    
    Parses output files from quantum chemistry programs to extract atomic
    symbols and Cartesian coordinates.
    
    Args:
        outlines (List[str]): Lines from output file
        min_RMS (int): Minimum RMS value to consider (-1 for final structure)
        n_atoms (int): Number of atoms in system
        program (str): QM program name ('gaussian' or 'orca')
        keywords_line (str): QM input keywords
        
    Returns:
        tuple: (atom_types, cartesians)
        - atom_types (List[str]): List of atomic symbols 
        - cartesians (List[List[float]]): XYZ coordinates for each atom
    """
    # Only Gaussian supported for now
    if program != "gaussian":
        return [], []
        
    # Get coordinate section
    target_ori = "Input orientation:" if "nosymm" in keywords_line.lower() else "Standard orientation:"
    range_lines = _find_coordinate_section(outlines, target_ori, min_RMS, n_atoms)
    
    if not range_lines:
        return [], []
        
    # Parse coordinates
    return _parse_gaussian_coordinates(outlines, range_lines[0], range_lines[1])
    
def _find_coordinate_section(outlines, target_ori, min_RMS, n_atoms):
    """Locate coordinate section in output file.
    
    Args:
        outlines (List[str]): File lines
        target_ori (str): Target orientation marker
        min_RMS (int): Minimum RMS value
        n_atoms (int): Number of atoms
        
    Returns:
        List[int]: Start and end line numbers for coordinates
    """
    if min_RMS > -1:
        return _find_rms_coordinates(outlines, target_ori, min_RMS, n_atoms)
    else:
        return _find_final_coordinates(outlines, target_ori, n_atoms)
        
def _find_rms_coordinates(outlines, target_ori, min_RMS, n_atoms):
    """Find coordinates for specific RMS value.
    
    Args:
        outlines (List[str]): File lines
        target_ori (str): Target orientation marker 
        min_RMS (int): Target RMS value
        n_atoms (int): Number of atoms
        
    Returns:
        List[int]: Start and end line numbers
    """
    count_RMS = -1
    for i, line in enumerate(outlines):
        if target_ori in line:
            count_RMS += 1
            if count_RMS == min_RMS:
                return [i + 5, i + 5 + n_atoms]
    return []
    
def _find_final_coordinates(outlines, target_ori, n_atoms): 
    """Find final coordinates in output.
    
    Args:
        outlines (List[str]): File lines
        target_ori (str): Target orientation marker
        n_atoms (int): Number of atoms
        
    Returns:
        List[int]: Start and end line numbers
    """
    for i in reversed(range(len(outlines))):
        if target_ori in outlines[i]:
            return [i + 5, i + 5 + n_atoms]
    return []
    
def _parse_gaussian_coordinates(outlines, start, end):
    """Parse Gaussian coordinate block.
    
    Args:
        outlines (List[str]): File lines
        start (int): Start line number 
        end (int): End line number
        
    Returns:
        tuple: (atom_types, cartesians)
    """
    per_tab = periodic_table()
    atom_types = []
    cartesians = []
    
    for i in range(start, end):
        fields = outlines[i].split()
        
        # Get atom symbol
        massno = int(fields[1])
        atom_symbol = per_tab[massno] if massno < len(per_tab) else "XX"
        atom_types.append(atom_symbol)
        
        # Get coordinates
        coords = [float(x) for x in fields[3:6]]
        cartesians.append(coords)
        
    return atom_types, cartesians