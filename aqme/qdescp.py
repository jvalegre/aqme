"""
Parameters
----------

General
+++++++

   w_dir_main : str, default=os.getcwd()
      Working directory
   destination : str, default=None,
      Directory to create the JSON file(s)
   program : str, default=xtb
      Program required to create the new descriptors. Current options: 'xtb', 'nmr'
   nprocs : int, default=None
      Number of xTB jobs run in parallel with 1 proc each (1 proc for reproducibility 
      in the results). Also, nprocs used in CSEARCH
   qdescp_atoms : list of str, default=[]
      Type of atom or group to calculate atomic properties. This option admits atoms 
      (i.e., qdescp_atoms=['P']) and SMART patterns (i.e., qdescp_atoms=['C=O']) 
   robert : bool, default=True
      Creates a database ready to use in an AQME-ROBERT machine learning workflow,
      combining the input CSV with SMILES/code_name and the calculated xTB/DBSTEP descriptors

xTB and MORFEUS descriptors
+++++++++++++++

   files or input (both options are valid) : list of str, default=''
      Filenames of SDF/PDB/XYZ/CSV files to calculate xTB descriptors. If CSV is selected, a CSV
      with two columns is required (code_name and SMILES), since AQME will generate conformers
      from SMILES with CSEARCH before QDESCP generates descriptors.
   charge : int, default=None
      Charge of the calculations used in the following input files (charges from
      SDF files generated in CSEARCH are read automatically).
   mult : int, default=None
      Multiplicity of the calculations used in the following input files 
      (multiplicities from SDF files generated in CSEARCH are read automatically).
   gfn_version : int, default="2"
      GFN version used in QDESCP to calculate descriptors.
   qdescp_solvent : str, default=None
      Solvent used in the xTB property calculations (ALPB model)
   qdescp_temp : float, default=300
      Temperature required for the xTB property calculations
   qdescp_acc : float, default=0.2
      Accuracy required for the xTB property calculations 
   qdescp_opt : str, default='normal'
      Convergence criteria required for the xTB property calculations 
   boltz : bool, default=True
      Calculation of Boltzmann averaged xTB properties and addition of RDKit 
      molecular descriptors
   xtb_opt : bool, default=True
      Performs an initial xTB geometry optimization before calculating descriptors
   vbur_radius : float, default=3.5
      Adjusts the radius in the buried volume calculations of MORFEUS

NMR simulation
++++++++++++++

   files : list of str, default=''
      Filenames of LOG files to retrieve NMR shifts from Gaussian calculations 
      (\*.log can be used to include all the log files in the working directory)
   boltz : bool, default=True
      Calculation of Boltzmann averaged NMR shifts
   nmr_atoms : list of str, default=[6, 1]
      List containing the atom types (as atomic numbers) to consider. For 
      example, if the user wants to retrieve NMR shifts from C and H atoms 
      nmr_atoms=[6, 1]
   nmr_slope : list of float, default=[-1.0537, -1.0784]
      List containing the slope to apply for the raw NMR shifts calculated with 
      Gaussian. A slope needs to be provided for each atom type in the analysis 
      (i.e., for C and H atoms, the nmr_slope=[-1.0537, -1.0784]). These values 
      can be adjusted using the CHESHIRE repository.
   nmr_intercept : list of float, default=[181.7815, 31.8723]
      List containing the intercept to apply for the raw NMR shifts calculated 
      with Gaussian. An intercept needs to be provided for each atom type in the
      analysis (i.e., for C and H atoms, the nmr_intercept=[-1.0537, -1.0784]). 
      These values can be adjusted using the CHESHIRE repository.
   nmr_experim : str, default=None
      Filename of a CSV containing the experimental NMR shifts. Two columnds are
      needed: A) 'atom_idx' should contain the indexes of the atoms to study as 
      seen in GaussView or other molecular visualizers (i.e., the first atom of 
      the coordinates has index 1); B) 'experimental_ppm' should contain the 
      experimental NMR shifts in ppm observed for the atoms.
"""
######################################################.
#        This file stores the QDESCP class           #
######################################################.

import os
import subprocess
import glob
import sys
import time
import json
import shutil
import concurrent.futures as futures
import numpy as np
from progress.bar import IncrementalBar
import pandas as pd
from pathlib import Path
from aqme.utils import (
    load_variables,
    read_xyz_charge_mult,
    mol_from_sdf_or_mol_or_mol2,
    run_command,
    check_files,
    check_dependencies,
    set_destination,
    load_sdf,
    blocking_wrapper
)
from aqme.qdescp_utils import (
    assign_prefix_atom_props,
    get_rdkit_properties,
    convert_ndarrays,
    read_json,
    calculate_morfeus_descriptors,
    collect_descp_lists,
    get_boltz_props_nmr,
    fix_cols_names,
    dict_to_json,
    full_level_boltz,
    get_mols_qdescp,
    get_matches_idx_n_prefix,
    auto_pattern,
    remove_invalid_smarts,
    update_atom_props_json,
    find_level_names,
    setup_env
)

from aqme.csearch.crest import xyzall_2_xyz


class PropertyCalculator:
    """Handles quantum chemical property calculations.
    
    This class manages the execution of quantum chemistry calculations and property
    collection using xTB and other tools. It handles:
    - XYZ file preparation
    - Charge and multiplicity assignment
    - Property calculation execution
    - Results collection and processing
    
    Attributes:
        args: AQME arguments object with calculation settings
    """
    
    def __init__(self, args):
        """Initialize property calculator.
        
        Args:
            args: AQME arguments object containing calculation parameters
        """
        self.args = args
        
    def calculate_properties(self, xyz_file, charge, mult, name, destination):
        """Run property calculations for a molecule.
        
        Executes full quantum chemistry workflow:
        1. Prepares input files
        2. Runs xTB calculation
        3. Collects properties
        4. Processes results
        
        Args:
            xyz_file (str): Path to XYZ geometry file
            charge (int): Molecular charge
            mult (int): Molecular multiplicity 
            name (str): Base name for output files
            destination (Path): Output directory
            
        Returns:
            tuple: (success status, file paths dict)
        """
        # Prepare working directory
        dat_dir = destination / name
        dat_dir.mkdir(exist_ok=True, parents=True)
        
        # Set up file paths
        files = {
            'xyz': dat_dir / f"{name}.xyz",
            'input': dat_dir / f"{name}_xtb.inp",
            'output': dat_dir / f"{name}_opt.out",
            'json': dat_dir / f"{name}.json"
        }
        
        # Move input file
        shutil.move(xyz_file, str(files['xyz']))
        
        # Create xTB input
        self._create_xtb_input(files['input'])
        
        # Set up environment
        env = setup_env(self)
        
        # Run calculation
        success = self._run_xtb_calculation(
            dat_dir, files, charge, mult, env
        )
        
        return success, {k: str(v) for k,v in files.items()}
        
    def _create_xtb_input(self, input_file):
        """Create xTB input file.
        
        Args:
            input_file (Path): Path to write input file
        """
        with open(input_file, "w", encoding='utf-8') as f:
            f.write("$write\njson=true\n")

    def _run_xtb_calculation(self, dat_dir, files, charge, mult, env):
        """Execute xTB calculation.
        
        Args:
            dat_dir (Path): Working directory
            files (dict): File path dictionary
            charge (int): Molecular charge
            mult (int): Molecular multiplicity
            env (dict): Environment variables
            
        Returns:
            bool: Success status
        """
        if not self.args.xtb_opt:
            return True
            
        command = [
            "xtb",
            str(files['xyz']),
            "--opt", str(self.args.qdescp_opt),
            "--acc", str(self.args.qdescp_acc),
            "--gfn", str(self.args.gfn_version),
            "--chrg", str(int(float(charge))),
            "--uhf", str(int(float(mult)) - 1),
            "--etemp", str(self.args.qdescp_temp),
            "-P", "1"
        ]
        
        if self.args.qdescp_solvent:
            command.extend(["--alpb", self.args.qdescp_solvent])
            
        run_command(command, files['output'], cwd=dat_dir, env=env)
        
        # Handle results
        os.remove(files['xyz'])
        
        if os.path.exists(dat_dir / "xtbopt.xyz"):
            os.rename(dat_dir / "xtbopt.xyz", files['xyz'])
            return True
            
        elif os.path.exists(dat_dir / "xtblast.xyz"):
            os.rename(dat_dir / "xtblast.xyz", files['xyz'])
            return True
            
        return False


class qdescp:
    """Quantum Mechanical Descriptor Calculation and Processing.

    This class handles the generation and processing of quantum mechanical descriptors
    using xTB, as well as analyzes NMR data. It provides functionality for:
    
    - Running xTB calculations for descriptor generation
    - Processing NMR data from Gaussian calculations 
    - Managing Boltzmann averaging of properties
    - Creating ROBERT-compatible descriptor databases
    
    The class supports both single molecule calculations and batch processing of
    multiple structures from various input formats including SDF, PDB, XYZ and CSV files.

    Attributes:
        args: Configuration object holding all calculation parameters
        start_time: Timestamp of calculation start
    """

    def __init__(self, **kwargs):
        """Initialize QDESCP calculator with provided settings.
        
        Args:
            **kwargs: Keyword arguments for configuration including:
            - w_dir_main (str): Working directory path
            - destination (str): Output directory path  
            - program (str): Program to use ('xtb' or 'nmr')
            - nprocs (int): Number of parallel processes
            - qdescp_atoms (list): Atoms/groups for atomic property calculation
            - robert (bool): Create ROBERT-compatible database
            - files (list): Input file paths
            - charge (int): Molecular charge
            - mult (int): Molecular multiplicity
            - gfn_version (int): GFN-xTB version
            - qdescp_solvent (str): Solvent for ALPB model
            - qdescp_temp (float): Temperature for calculations
            - qdescp_acc (float): Calculation accuracy  
            - qdescp_opt (str): Optimization convergence criteria
            - boltz (bool): Calculate Boltzmann averages
            - xtb_opt (bool): Run xTB optimization
        """
        self.start_time = time.time()
        
        # Initialize configuration and components
        self.args = load_variables(kwargs, "qdescp")
        self.property_calc = PropertyCalculator(self.args)
        
        # Set up QDESCP run and check dependencies
        self, destination, smarts_targets, boltz_dir = self.qdescp_set_up()
        _ = check_dependencies(self)

        # full xTB workflow in QDESCP for descriptor generation and collection
        if self.args.program.lower() == "xtb":
            _ = self.qdescp_xtb_workflow(boltz_dir,destination,smarts_targets)

        # full NMR workflow in QDESCP for NMR prediction
        elif self.args.program.lower() == "nmr":
            _ = self.qdescp_nmr_workflow(boltz_dir)
        
        elapsed_time = round(time.time() - self.start_time, 2)
        self.args.log.write(f"\nTime QDESCP: {elapsed_time} seconds\n")
        self.args.log.finalize()
    

    def qdescp_xtb_workflow(self, boltz_dir, destination, smarts_targets):
        """Run the complete xTB workflow for descriptor generation and collection.
        
        This method handles the complete process of generating and collecting quantum
        mechanical descriptors using xTB, including:
        
        1. Initial input validation and setup
        2. Optional conformer generation for CSV inputs
        3. Automatic SMARTS pattern detection
        4. Parallel descriptor calculation
        5. Results collection and processing
        
        Args:
            boltz_dir (Path): Directory for Boltzmann-averaged results
            destination (Path): Main output directory
            smarts_targets (list): List of SMARTS patterns to match
            
        Note:
            - For CSV inputs, conformers are generated before descriptor calculation
            - When no SMARTS patterns are provided, they are auto-detected
            - Invalid SMARTS patterns (< 75% compatibility) are removed
            - Uses parallel processing for efficiency with reproducibility
        """

        # check whether the user have chosen the "input" or "files" option (QDESCP will use "files" from this point on)
        qdescp_files = self.initial_xtb_check()

        self.args.log.write(f"\nStarting QDESCP-{self.args.program} with {len(qdescp_files)} job(s)\n")

        # if the files input is a CSV, first the program generates conformers
        if len(qdescp_files) == 1 and os.path.basename(qdescp_files[0]).split('.')[-1].lower() == 'csv':
            qdescp_files = self.initial_csearch_run(destination, qdescp_files)

        # obtaining mols from input files that will be used to set up atomic descriptors
        mol_list = get_mols_qdescp(qdescp_files)

        # obtaing SMARTS patterns from the input files automatically if no patterns are provided
        if len(smarts_targets) == 0 and len(qdescp_files) > 1:
            smarts_targets = auto_pattern(mol_list,smarts_targets)
            if len(smarts_targets) > 0:
                self.args.log.write(f"\no  Common atoms/groups found in all molecules (repeated once in each): {smarts_targets}")

        # Delete a SMARTS pattern if it is not compatible with more than 75% of the sdf files
        if len(smarts_targets) > 0:
            smarts_targets = remove_invalid_smarts(self,mol_list,smarts_targets)

        # Get descriptors (denovo, interpret, full)
        descp_dict = collect_descp_lists()

        # run all the calculations to generate xTB outputs and JSON files with descriptors

        self.args.invalid_calcs = [] # keep track of unvalid calcs
        bar = IncrementalBar("\no  Number of finished jobs from QDESCP", max=len(qdescp_files))

        # multiprocessing to accelerate and make QDESCP reproducible (since xTB uses 1 processor to be reproducible)
        if not self.args.debug: # errors and try/excepts are not shown in multithreading
            with futures.ThreadPoolExecutor(max_workers=self.args.nprocs) as executor:
                # Submit all tasks at once
                future_tasks = [
                    executor.submit(
                        blocking_wrapper,
                        self.gather_files_and_run,
                        destination,
                        file,
                        descp_dict['atom_props'],
                        smarts_targets,
                        bar
                    )
                    for file in qdescp_files
                ]
                
                # Wait for all tasks to complete
                done, _ = futures.wait(future_tasks, return_when=futures.ALL_COMPLETED)

                # Handle results & exceptions
                for future in done:
                    try:
                        future.result()  # Raises any exception from the thread
                    except Exception as e:
                        self.args.log.write(f"QDESCP raised an exception: {e}\n")

                # When this point is reached, ALL tasks have finished (successfully or not)
        else:
            for file in qdescp_files:
                _ = self.gather_files_and_run(destination, file, descp_dict['atom_props'], smarts_targets, bar)

        bar.finish()

        if self.args.boltz:
            #AQME-ROBERT workflow: Combines the descriptor data from qdescp CSVs with the input CSV and saves the result.
            _ = self.get_boltz_n_save_csv(destination,qdescp_files,descp_dict,boltz_dir,smarts_targets)


    def initial_xtb_check(self):
        """Validate and process input file selection.
        
        This method:
        1. Checks if input is provided via --input or --files
        2. Validates file formats and paths
        3. Converts input to standardized file list
        
        Returns:
            list: List of validated input file paths
            
        Raises:
            SystemExit: If no valid files are found or formats are invalid
        """
        
        valid_input = True
        if self.args.files == [] and self.args.input != '':
            if os.path.basename(self.args.input).split('.')[-1].lower() != "csv":
                self.args.log.write(f"\nx  The format used ({os.path.basename(self.args.input).split('.')[-1]}) is not compatible with the 'input' option! Formats accepted: csv")
                valid_input = False
            if self.args.input[0] == '[' or isinstance(self.args.input, list):
                self.args.log.write(f"\nx  The 'input' option was specified as a list! Please provide only the PATH or name of the CSV (i.e. --input test.csv)")
                valid_input = False
            qdescp_files = [self.args.input]
        elif self.args.files == []:
            self.args.log.write(f'\nx  No files were found! Please provide the correct PATH to your input files (i.e. --files "*.sdf")')
            valid_input = False
        else:
            if os.path.basename(self.args.files[0]).split('.')[-1].lower() != "sdf":
                self.args.log.write(f"\nx  The format used ({os.path.basename(self.args.files[0]).split('.')[-1]}) is not compatible with the 'files' option! Formats accepted: sdf")
                valid_input = False
            qdescp_files = self.args.files

        if not valid_input:
            self.args.log.finalize()
            sys.exit()

        return qdescp_files


    def initial_csearch_run(self, destination, qdescp_files):
        """Generate conformers from SMILES in CSV input.
        
        This method:
        1. Validates CSV input file existence
        2. Sets up conformer search parameters
        3. Runs RDKit conformer generation
        4. Processes and validates generated conformers
        
        Args:
            destination (Path): Output directory path
            qdescp_files (list): List of input files (expecting single CSV)
            
        Returns:
            list: Paths to generated conformer files
            
        Raises:
            SystemExit: If CSV file not found or conformer generation fails
        """
        
        if not os.path.exists(qdescp_files[0]):
            self.args.log.write(f"\nx  The csv_name provided ({qdescp_files[0]}) does not exist! Please specify this name correctly")
            self.args.log.finalize()
            sys.exit()
        # the default number of conformers is reduced to 5 unless overwritten by the user
        if self.args.sample == 25:
            sample_qdescp = 5
        else:
            sample_qdescp = self.args.sample

        # sets the csv_name variable to create the AQME-ROBERT descriptor file
        self.args.csv_name = qdescp_files[0]

        if f'{os.path.basename(destination).upper()}' == 'QDESCP':
            destination_csearch = Path(os.path.dirname(destination)).joinpath('CSEARCH')
        else:
            destination_csearch = destination.joinpath('CSEARCH')

        cmd_csearch = ['python', '-m', 'aqme', '--csearch', '--program', 'rdkit', '--input', 
                    f'{self.args.csv_name}', '--sample', f'{sample_qdescp}', '--destination', f'{destination_csearch}',
                    '--nprocs', f'{self.args.nprocs}','--auto_sample',self.args.auto_sample, '--ff',self.args.ff]

        if self.args.single_system:
            cmd_csearch.append('--single_system')

        # overwrites charge/mult if the user specifies values
        if self.args.charge is not None:
            cmd_csearch = cmd_csearch + ['--charge', f'{self.args.charge}']
        if self.args.mult is not None:
            cmd_csearch = cmd_csearch + ['--mult', f'{self.args.mult}']
        subprocess.run(cmd_csearch)

        # use only the molecules from the input CSV (ignore previous/unrelated CSEARCH runs that generated SDFs)
        qdescp_files = []
        csearch_files = glob.glob(f'{destination_csearch}/*.sdf')
        if len(csearch_files) > 0:
            df_qdescp = pd.read_csv(self.args.csv_name)
            for file in csearch_files:
                # Match if the SDF file starts with the code_name (handles AAA_0_rdkit.sdf for code_name AAA)
                file_base = os.path.basename(file)
                for code_name in df_qdescp['code_name'].astype(str).tolist():
                    if file_base.startswith(f"{code_name}_"):
                        qdescp_files.append(file)
                        break
        if len(qdescp_files) == 0:
            self.args.log.write(f"\nx  WARNING! The CSEARCH conformational search did not produce any results.")
            self.args.log.finalize()
            sys.exit()

        return qdescp_files
    

    def get_boltz_n_save_csv(self, destination, qdescp_files, descp_dict, boltz_dir, smarts_targets):
        """Process molecular properties and generate Boltzmann-averaged data files.
        
        This method:
        1. Processes each valid input file to extract molecular properties
        2. Generates Boltzmann-weighted JSON files for each molecule
        3. Creates final CSV files with combined data
        
        Args:
            destination (Path): Output directory path
            qdescp_files (list): List of input structure files
            descp_dict (dict): Dictionary of descriptors to calculate
            boltz_dir (Path): Directory for Boltzmann-weighted results
            smarts_targets (list): SMARTS patterns for atom matching
        """
        self.args.log.write('\no  Running RDKit and collecting molecular properties (for all inputs)')
        all_prefixes_atoms = []
        
        for _, file in enumerate(qdescp_files):
            if file in self.args.invalid_calcs:
                continue
                
            # Process each valid molecule
            descp_dict_indiv = descp_dict.copy()
            mols = load_sdf(file)
            mol = mols[0]
            
            # Get molecule name from CSV or file
            name = self._get_molecule_name(file)
            
            # Find and process JSON files for each conformer
            json_files = self._find_conformer_jsons(destination, name)
            if json_files:
                _ = self.get_boltz_props(
                    json_files, name, boltz_dir, "xtb", 
                    descp_dict_indiv, smarts_targets, mol, all_prefixes_atoms
                )
            
        # Create final CSV files
        _ = self.write_csv_boltz_data(destination)
        
    def _get_molecule_name(self, file):
        """Extract molecule name from CSV if available, otherwise from filename.
        
        Args:
            file (str): Path to input file
            
        Returns:
            str: Molecule name/identifier
        """
        if self.args.csv_name and os.path.exists(self.args.csv_name):
            df_qdescp = pd.read_csv(self.args.csv_name)
            file_base = os.path.basename(file)
            
            for code_name in df_qdescp['code_name'].astype(str).tolist():
                if file_base.startswith(f"{code_name}_"):
                    return code_name
                    
        return os.path.basename(file).replace('.sdf', '')
        
    def _find_conformer_jsons(self, destination, name):
        """Find JSON files for all conformers of a molecule.
        
        Args:
            destination (Path): Directory containing JSON files
            name (str): Base name of the molecule
            
        Returns:
            list: Paths to conformer JSON files
        """
        return [x for x in glob.glob(f"{destination}/*.json") 
                if os.path.basename(x).startswith(f"{name}_") 
                and "_conf_" in os.path.basename(x)]


    def get_boltz_props(self, json_files, name, boltz_dir, calc_type, descp_dict_indiv, smarts_targets, mol, all_prefixes_atoms):
        """Calculate Boltzmann-weighted properties from conformer results.
        
        This method:
        1. Reads properties from JSON files
        2. Calculates Boltzmann weights based on energies
        3. Averages properties across conformers
        4. Adds RDKit molecular descriptors
        5. Saves results to new JSON file
        
        Args:
            json_files (list): List of JSON files with conformer data
            name (str): Base name for output files
            boltz_dir (Path): Output directory for Boltzmann results
            calc_type (str): Calculation type ('xtb' or 'nmr')
            descp_dict_indiv (dict): Property dictionary for this molecule
            smarts_targets (list): SMARTS patterns to match
            mol: RDKit molecule object for descriptor calculation
            all_prefixes_atoms (list): List of atomic property prefixes
            
        Returns:
            None: Results are saved to JSON file
        """
        # Ensure smarts_targets is a list even if None
        if smarts_targets is None:
            smarts_targets = []

        full_json_data = {}

        energy = []
        for _, json_file in enumerate(json_files):
            json_data = read_json(json_file)
            energy.append(json_data["total energy"] if calc_type.lower() == "xtb" else json_data["optimization"]["scf"]["scf energies"][-1])

        # include the prefixes for names of atomic properties
        if len(json_data['prefixes_atom_prop']) > 0:
            all_prefixes_atoms = all_prefixes_atoms + json_data['prefixes_atom_prop']
            descp_dict_indiv['atom_props'],descp_dict_indiv['interpret_atoms'],descp_dict_indiv['denovo_atoms'] = assign_prefix_atom_props(json_data['prefixes_atom_prop'],descp_dict_indiv['atom_props'],descp_dict_indiv['interpret_atoms'],descp_dict_indiv['denovo_atoms'])

        # get all the properties (full level)
        full_json_data = full_level_boltz(descp_dict_indiv,json_files,energy,smarts_targets,full_json_data)

        # Calculate RDKit descriptors if molecule is provided
        if mol is not None:
            # Calculate all RDKit properties for full_json_data
            full_json_data = get_rdkit_properties(self,full_json_data, mol)

        _ = convert_ndarrays(full_json_data)

        # Save the averaged properties to a file
        _ = dict_to_json(self,os.path.join(boltz_dir, f"{name}_boltz.json"), full_json_data)


    def qdescp_set_up(self):
        """Initialize and validate QDESCP run settings.
        
        This method performs initial setup and validation:
        1. Sets default program to xTB if not specified
        2. Validates program selection (xTB or NMR)
        3. Configures parallel processing settings
        4. Sets sampling parameters
        5. Validates input files and directories
        6. Creates required output directories
        
        Returns:
            tuple: Contains:
                - self: Updated QDESCP instance
                - destination (Path): Configured output directory
                - smarts_targets (list): SMARTS patterns to match
                - boltz_dir (Path): Directory for Boltzmann calculations
  
        Raises:
            SystemExit: If program selection or input files are invalid
        """
        # Program selection validation
        if self.args.program is None:
            self.args.program = "xtb"

        if self.args.program.lower() not in ["xtb", "nmr"]:
            self._error_exit(
                f"Program {self.args.program} not supported for QDESCP. "
                "Use 'xtb' or 'nmr'"
            )

        # Processing configuration
        if self.args.nprocs is None:
            self.args.nprocs = 8

        if self.args.auto_sample == 'auto':
            self.args.auto_sample = 'low'

        # Validate input parameters
        self._validate_inputs()
        
        # Process input files if provided
        if self.args.files and not self.args.input:
            check_files(self, 'qdescp')
            self.args.files = self.get_unique_files()
            
        # Copy SMARTS patterns for atomic descriptors
        smarts_targets = self.args.qdescp_atoms.copy()

        destination = set_destination(self,'QDESCP')

        # print version of xTB
        destination.mkdir(exist_ok=True, parents=True)

        # create folder to store Boltzmann weighted properties
        boltz_dir = Path(f"{destination}/boltz")
        if os.path.exists(f"{boltz_dir}"):
            self.args.log.write(f'\nx  A previous folder of {boltz_dir} already existed, it was removed and replaced with the results of this QDESCP run.')
            shutil.rmtree(f"{boltz_dir}")
        boltz_dir.mkdir(exist_ok=True, parents=True)

        return self,destination,smarts_targets,boltz_dir


    def write_csv_boltz_data(self, destination):
        """Generate CSV files containing Boltzmann-averaged descriptor data.
        
        This method combines descriptor data from JSON files with input CSV data
        and generates multiple output files:
        - Full descriptor database
        - De novo descriptor subset
        - Interpretable descriptor subset
        
        The method handles both AQME-Descriptors and AQME-ROBERT workflows.
        
        Args:
            destination (Path): Directory containing Boltzmann JSON files
        """
        temp_csv = Path(os.getcwd()).joinpath("AQME_run.csv")
        code_names = self._get_code_names(temp_csv)
        
        # Process JSON files
        boltz_data = self._process_boltz_json_files(destination, code_names)
        if boltz_data is None or boltz_data.empty:
            self.args.log.write("x  No descriptors were generated with QDESCP, please check the WARNINGS above.")
            return
            
        # Generate output files
        name_db = 'ROBERT' if self.args.robert and os.path.exists(self.args.csv_name) else 'Descriptors'
        self._generate_output_files(destination, boltz_data, name_db, temp_csv)
        
    def _get_code_names(self, temp_csv):
        """Extract or generate code names for molecules.
        
        Args:
            temp_csv (Path): Path for temporary CSV file
            
        Returns:
            list: List of code names
        """
        if self.args.csv_name and os.path.exists(self.args.csv_name):
            return pd.read_csv(self.args.csv_name)["code_name"].astype(str).tolist()
            
        # Generate code names from file names
        code_names = [os.path.basename(file).replace('.sdf','') for file in self.args.files]
        normalized_names = self._normalize_code_names(code_names)
        
        # Create temporary CSV
        df = pd.DataFrame({
            'code_name': normalized_names,
            'SMILES': [np.nan] * len(code_names)
        })
        df.to_csv(temp_csv, index=False)
        self.args.csv_name = str(temp_csv)
        
        return code_names
        
    def _normalize_code_names(self, code_names):
        """Normalize code names by removing rdkit suffixes.
        
        Args:
            code_names (list): Original code names
            
        Returns:
            list: Normalized code names
        """
        df_name = pd.DataFrame({'code_name': code_names})
        return list(df_name['code_name'].astype(str)
                   .str.replace(r"(_\d+)?_rdkit$", "", regex=True)
                   .str.replace(r"_rdkit$", "", regex=True))
                   
    def _process_boltz_json_files(self, destination, code_names):
        """Process Boltzmann-weighted JSON files into DataFrame.
        
        Args:
            destination (Path): Directory containing JSON files
            code_names (list): List of molecule code names
            
        Returns:
            DataFrame: Combined data from all JSON files
        """
        json_pattern = str(destination) + "/boltz/*_boltz.json"
        boltz_json_files = glob.glob(json_pattern)
        
        dfs = []
        for file in boltz_json_files:
            data = pd.read_json(file, lines=True)
            name_indiv = self._map_json_to_code_name(file, code_names)
            data.insert(loc=0, column='code_name', value=name_indiv)
            dfs.append(data)
            
        return pd.concat(dfs, ignore_index=True) if dfs else None
        
    def _map_json_to_code_name(self, json_file, code_names):
        """Map JSON filename to original code name.
        
        Args:
            json_file (str): Path to JSON file
            code_names (list): List of original code names
            
        Returns:
            str: Mapped code name
        """
        json_basename = os.path.basename(json_file).replace("_boltz.json", "")
        if json_basename.endswith("_rdkit"):
            json_basename = json_basename[:-6]
            
        for code_name in code_names:
            if json_basename == code_name or json_basename.startswith(f"{code_name}_"):
                return code_name
        return json_basename
        
    def _generate_output_files(self, destination, df_full, name_db, temp_csv):
        """Generate final CSV output files.
        
        Args:
            destination (Path): Output directory
            df_full (DataFrame): Combined descriptor data
            name_db (str): Database name (ROBERT or Descriptors)
            temp_csv (Path): Path to temporary CSV
        """
        if not os.path.exists(self.args.csv_name):
            self.args.log.write(
                f"\nx  The input csv_name provided ({self.args.csv_name}) does not exist. "
                f"A combined database for AQME-{name_db} workflows will not be created "
                "(but you can still check the raw descriptors in the QDESCP folder)."
            )
            return
            
        combined_df = self._prepare_combined_dataframe(df_full)
        if combined_df is None:
            return
            
        self._save_descriptor_files(destination, combined_df, name_db)
        _ = self.process_aqme_csv(name_db)
        os.remove(temp_csv) if os.path.exists(temp_csv) else None
        
    def _prepare_combined_dataframe(self, df_full):
        """Prepare combined DataFrame with input and descriptor data.
        
        Args:
            df_full (DataFrame): Descriptor data
            
        Returns:
            DataFrame: Combined data, or None if invalid
        """
        input_df = pd.read_csv(self.args.csv_name)
        input_df = fix_cols_names(input_df)
        
        if 'code_name' not in input_df.columns:
            self.args.log.write(
                f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain "
                "the code_name column. A combined database for AQME-{name_db} workflows "
                "will not be created."
            )
            return None
            
        if 'SMILES' not in input_df.columns:
            return None
            
        df_full["normalized_code_name"] = (
            df_full["code_name"].astype(str)
            .str.replace(r"(_\d+)?_rdkit$", "", regex=True)
            .str.replace(r"_rdkit$", "", regex=True)
        )
        
        input_df["code_name"] = input_df["code_name"].astype(str)
        
        return input_df.merge(
            df_full.drop(columns=["code_name"]),
            left_on="code_name",
            right_on="normalized_code_name",
            how="left"
        ).drop(columns=["normalized_code_name"])
        
    def _save_descriptor_files(self, destination, combined_df, name_db):
        """Save descriptor data to various output files.
        
        Args:
            destination (Path): Output directory
            combined_df (DataFrame): Combined descriptor data
            name_db (str): Database name
        """
        csv_basename = os.path.basename(self.args.csv_name)
        paths = self._setup_output_paths(name_db, csv_basename)
        
        # Save raw data if no specific atoms were selected
        if len(self.args.qdescp_atoms) == 0:
            self._save_raw_data(destination, combined_df, paths)
            
        # Save clean data (without list columns)
        clean_df = combined_df.drop(columns=[
            col for col in combined_df.columns 
            if combined_df[col].apply(lambda x: isinstance(x, list)).any()
        ])
        
        self._save_clean_data(clean_df, paths)
        self._log_success(name_db, csv_basename)
        
    def _setup_output_paths(self, name_db, csv_basename):
        """Set up paths for output files.
        
        Args:
            name_db (str): Database name
            csv_basename (str): Base name for CSV files
            
        Returns:
            dict: Dictionary of output paths
        """
        paths = {
            'full': self.args.initial_dir.joinpath(f'AQME-{name_db}_full_{csv_basename}'),
            'denovo': self.args.initial_dir.joinpath(f'AQME-{name_db}_denovo_{csv_basename}'),
            'interpret': self.args.initial_dir.joinpath(f'AQME-{name_db}_interpret_{csv_basename}')
        }
        
        # Remove existing files
        for path in paths.values():
            if os.path.exists(path):
                os.remove(path)
                
        return paths
        
    def _save_raw_data(self, destination, combined_df, paths):
        """Save raw data including atomic descriptors.
        
        Args:
            destination (Path): Output directory
            combined_df (DataFrame): Combined data
            paths (dict): Output file paths
        """
        dat_dir = Path(destination.joinpath('raw_data'))
        dat_dir.mkdir(exist_ok=True, parents=True)
        
        # Save full data
        combined_df.to_csv(
            Path(dat_dir).joinpath(os.path.basename(paths['full'])),
            index=None, header=True
        )
        
        # Save descriptor subsets
        for level in ['denovo', 'interpret']:
            descriptors = find_level_names(combined_df, level)
            subset_df = combined_df[descriptors]
            subset_df.to_csv(
                Path(dat_dir).joinpath(os.path.basename(paths[level])),
                index=None, header=True
            )
            
    def _save_clean_data(self, clean_df, paths):
        """Save clean data without list columns.
        
        Args:
            clean_df (DataFrame): Cleaned data
            paths (dict): Output file paths
        """
        clean_df.to_csv(paths['full'], index=None, header=True)
        
        for level in ['denovo', 'interpret']:
            descriptors = find_level_names(clean_df, level)
            subset_df = clean_df[descriptors]
            subset_df.to_csv(paths[level], index=None, header=True)
            
    def _log_success(self, name_db, csv_basename):
        """Log successful file creation.
        
        Args:
            name_db (str): Database name
            csv_basename (str): Base name of CSV files
        """
        self.args.log.write(
            f"o  The AQME-{name_db}_full_{csv_basename}, "
            f"AQME-{name_db}_denovo_{csv_basename} and "
            f"AQME-{name_db}_interpret_{csv_basename} databases "
            f"were created in {self.args.initial_dir}"
        )


    def gather_files_and_run(self, destination, file, atom_props, smarts_targets, bar):
        """Process input file(s) through xTB calculation and property collection.
        
        This method handles the complete process for a single input file:
        1. Converts input files to XYZ format if needed
        2. Extracts charge and multiplicity information
        3. Runs xTB calculations for each conformer
        4. Collects and processes properties
        5. Generates JSON output files
        
        Args:
            destination (Path): Output directory path
            file (str): Input file path (XYZ/PDB/SDF format)
            atom_props (list): Atomic properties to collect
            smarts_targets (list): SMARTS patterns to match
            bar (IncrementalBar): Progress bar instance
            
        Note:
            - For XYZ files, conformers are processed directly
            - Other formats are converted to XYZ using OpenBabel
            - Charge/multiplicity are read from files or use defaults
            - Progress is tracked via the progress bar
        """
        # Get base name and extension
        name = '.'.join(os.path.basename(Path(file)).split('.')[:-1])
        ext = os.path.basename(Path(file)).split(".")[-1]
        self.args.log.write(f"\n\n   ----- {name} -----")
        
        # Get conformers and their properties
        xyz_files, xyz_charges, xyz_mults = self._get_conformer_data(file, name, ext)
        
        # Process each conformer
        for xyz_file, charge, mult in zip(xyz_files, xyz_charges, xyz_mults):
            name_xtb = '.'.join(os.path.basename(Path(xyz_file)).split(".")[:-1])
            self._process_single_conformer(
                destination, file, xyz_file, charge, mult, name_xtb, 
                atom_props, smarts_targets
            )
            
        bar.next()
        
    def _get_conformer_data(self, file, name, ext):
        """Extract conformer geometries and properties.
        
        Args:
            file (str): Input file path
            name (str): Base name without extension 
            ext (str): File extension
            
        Returns:
            tuple: Lists of (xyz files, charges, multiplicities)
        """
        xyz_files, xyz_charges, xyz_mults = [], [], []
        
        if ext.lower() == "xyz":
            xyzall_2_xyz(file, name)
            xyz_files, xyz_charges, xyz_mults = self._process_xyz_conformers(file, name)
        else:
            self._convert_to_xyz(file, name, ext)
            xyz_files, xyz_charges, xyz_mults = self._process_other_conformers(file, name)
            
        return xyz_files, xyz_charges, xyz_mults
        
    def _process_xyz_conformers(self, file, name):
        """Process conformers from XYZ input.
        
        Args:
            file (str): Input XYZ file
            name (str): Base name for outputs
            
        Returns:
            tuple: Lists of (xyz files, charges, multiplicities)
        """
        xyz_files, xyz_charges, xyz_mults = [], [], []
        xyz_files_list = [
            x for x in glob.glob(f"{os.path.dirname(Path(file))}/*.xyz") 
            if os.path.basename(x).startswith(f'{name}_conf_')
        ]
        
        for conf_file in xyz_files_list:
            charge = (self.args.charge if self.args.charge is not None 
                     else read_xyz_charge_mult(conf_file)[0])
            mult = (self.args.mult if self.args.mult is not None 
                   else read_xyz_charge_mult(conf_file)[1])
                   
            xyz_files.append(os.path.dirname(os.path.abspath(file)) + "/" + conf_file)
            xyz_charges.append(charge)
            xyz_mults.append(mult)
            
        return xyz_files, xyz_charges, xyz_mults
        
    def _process_other_conformers(self, file, name):
        """Process conformers from non-XYZ input.
        
        Args:
            file (str): Input file path
            name (str): Base name for outputs
            
        Returns:
            tuple: Lists of (xyz files, charges, multiplicities)
        """
        xyz_files, xyz_charges, xyz_mults = [], [], []
        xyz_files_list = [
            x for x in glob.glob(f"{os.path.dirname(Path(file))}/*.xyz") 
            if os.path.basename(x).startswith(f'{name}_conf_')
        ]
        
        # Get charges and multiplicities
        charges = ([self.args.charge] * len(xyz_files_list) if self.args.charge is not None
                  else mol_from_sdf_or_mol_or_mol2(file, "csearch", self.args)[1])
                  
        mults = ([self.args.mult] * len(xyz_files_list) if self.args.mult is not None
                 else mol_from_sdf_or_mol_or_mol2(file, "csearch", self.args)[2])
        
        # Collect conformer data
        for i, conf_file in enumerate(xyz_files_list):
            xyz_files.append(conf_file)
            xyz_charges.append(charges[i])
            xyz_mults.append(mults[i])
            
        return xyz_files, xyz_charges, xyz_mults
        
    def _process_single_conformer(self, destination, file, xyz_file, charge, mult,
                                name_xtb, atom_props, smarts_targets):
        """Process a single conformer through xTB and property collection.
        
        Args:
            destination (Path): Output directory
            file (str): Original input file
            xyz_file (str): Conformer XYZ file
            charge (int): Molecular charge
            mult (int): Multiplicity 
            name_xtb (str): Base name for outputs
            atom_props (list): Properties to collect
            smarts_targets (list): SMARTS patterns
        """
        self.args.log.write(f"\no  Running xTB and collecting properties ({name_xtb})")
        
        xtb_passing, xtb_files_props = self.run_opt_xtb(
            file, xyz_file, charge, mult, name_xtb, destination
        )
        
        path_name = Path(os.path.dirname(file)).joinpath(
            '.'.join(os.path.basename(Path(file)).split(".")[:-1])
        )

        if xtb_passing:
            self.morfeus_properties(
                path_name, atom_props, smarts_targets,
                xtb_files_props, charge, mult, file
            )

        self.cleanup(name_xtb, destination, xtb_passing, xtb_files_props)


    def run_opt_xtb(self, file, xyz_file, charge, mult, name, destination):
        """Run xTB property calculations for a molecule.

        Args:
            file (str): Original input file path
            xyz_file (str): XYZ format geometry file path
            charge (int): Molecular charge
            mult (int): Molecular multiplicity
            name (str): Base name for output files
            destination (Path): Output directory path
            
        Returns:
            tuple: (success status, dict of file paths)
        """
        success, files = self.property_calc.calculate_properties(
            xyz_file, charge, mult, name, destination
        )
        
        if not success and file not in self.args.invalid_calcs:
            self.args.invalid_calcs.append(file)
            self.args.log.write(
                f"x  WARNING! {file} did not finish correctly and no descriptors "
                "will be generated for this system. Common causes: incorrect "
                "CHARGE and/or MULTIPLICITY (adjust with --charge/--mult options)."
            )
            
        return success, files


    def morfeus_properties(self, name_initial, atom_props, smarts_targets, xtb_files_props, charge, mult, file):
        """Calculate and collect MORFEUS molecular descriptors.
        
        This method:
        1. Reads molecular geometry from XYZ file
        2. Calculates MORFEUS descriptors
        3. Assigns atomic properties based on SMARTS patterns
        4. Saves results to JSON file
        
        Args:
            name_initial (str): Base name for input/output files
            atom_props (list): Atomic properties to calculate
            smarts_targets (list): SMARTS patterns for atom matching
            xtb_files_props (dict): Dictionary of file paths 
            charge (int): Molecular charge
            mult (int): Molecular multiplicity
            file (str): Original input file path
            
        Returns:
            None: Results are saved to JSON file
            
        Raises:
            FileNotFoundError: If XYZ file is missing
            Exception: If MORFEUS descriptor calculation fails
        """

        # load initial json and add coordinates
        json_data = {}

        xyz_path = Path(xtb_files_props['xyz']).resolve()
        json_path = Path(xtb_files_props['json']).resolve()

        try:
            with open(xyz_path, "r", encoding='utf-8') as f:
                inputs = f.readlines()
        except FileNotFoundError:
            self.args.log.write(f"x  ERROR! The xyz file for {name_initial} was not found. No descriptors will be generated for this system.")
            return
        
        coordinates = [inputs[i].strip().split()[1:] for i in range(2, int(inputs[0].strip()) + 2)]
        json_data["coordinates"] = coordinates

        # add MORFEUS properties to JSON
        try:
            global_properties_morfeus = calculate_morfeus_descriptors(str(xyz_path),self,charge,mult,smarts_targets,name_initial)
            json_data.update(global_properties_morfeus)
        except Exception as e:
            self.args.log.write(f"x  ERROR! Failed to calculate MORFEUS descriptors for {name_initial}: {e}\n")
            if file not in self.args.invalid_calcs:
                self.args.invalid_calcs.append(file)
            return

        # assign atomic properties to the corresponding atoms
        json_data = self.assign_atomic_properties(json_data,name_initial,atom_props,smarts_targets)

        json_path.parent.mkdir(parents=True, exist_ok=True)
        with json_path.open("w", encoding="utf-8") as outfile:
            json.dump(json_data, outfile)


    def assign_atomic_properties(self, json_data, name_initial, atom_props, smarts_targets):
        """Assign atomic properties based on SMARTS pattern matches.
        
        This method:
        1. Matches SMARTS patterns to atoms
        2. Assigns property prefixes to matched atoms
        3. Updates JSON data with atomic properties
        
        Args:
            json_data (dict): Dictionary of molecular data
            name_initial (str): Base name for molecule
            atom_props (list): Atomic properties to assign
            smarts_targets (list): SMARTS patterns to match
            
        Returns:
            dict: Updated JSON data with atomic properties
            
        Note:
            - Property prefixes are tracked to avoid duplicates
            - Properties are assigned only to matched atoms
            - Pattern matches use RDKit SMARTS matcher
        """

        prefixes_atom_prop = []
        
        pattern_dict = get_matches_idx_n_prefix(self,smarts_targets,name_initial)
        if len(pattern_dict.keys()) > 0:
            for pattern in pattern_dict:
                # Assign atomic descriptors to each identified atom and update database for final JSON file
                prefixes_atom_prop, json_data = update_atom_props_json(pattern_dict[pattern]['sorted_indices'],
                                                    pattern_dict[pattern]['match_names'],
                                                    atom_props,json_data,prefixes_atom_prop,pattern,
                                                    pattern_dict[pattern]['n_types'])

        # updates the prefixes used for atomic props
        json_data['prefixes_atom_prop'] = prefixes_atom_prop

        return json_data
    

    def cleanup(self, name, destination, xtb_passing, xtb_files_props):
        """Clean up calculation files and organize results.
        
        This method:
        1. Moves successful calculation results to destination
        2. Moves failed calculations to failed/ subdirectory
        3. Cleans up temporary files and directories
        
        Args:
            name (str): Base name for files
            destination (Path): Output directory path
            xtb_passing (bool): Whether xTB calculation succeeded
            xtb_files_props (dict): Dictionary of file paths
            
        Note:
            - Successful calculations: JSON and XYZ files preserved
            - Failed calculations: All files moved to failed/ directory
            - Temporary directories are removed
        """

        if xtb_passing: # only move molecules with successful xTB calcs
            # save JSON and XYZ files
            final_json = f"{destination}/{name}.json"
            shutil.move(xtb_files_props['json'], final_json)
            final_xyz = f"{destination}/{name}.xyz"
            shutil.move(xtb_files_props['xyz'], final_xyz)

            # delete xTB raw data
            shutil.rmtree(f"{destination}/{name}")

        else:
            if not os.path.exists(f"{destination}/failed"): 
                Path(f"{destination}/failed").mkdir()
            if os.path.exists(f"{destination}/failed/{name}"): 
                self.args.log.write(f'\nx  A previous folder of {name} already existed in failed, it was removed and replaced with the results of this QDESCP run.')
                shutil.rmtree(f"{destination}/failed/{name}")
            shutil.move(f"{destination}/{name}", f"{destination}/failed/{name}")


    def get_unique_files(self):
        """Filter input files to remove duplicates based on SMILES.
        
        This method:
        1. Reads SMILES strings from SDF files
        2. Identifies duplicate structures
        3. Keeps only unique structures
        4. Warns about duplicates
        
        Returns:
            list: Paths to unique input files
            
        Note:
            - Duplicates are identified by exact SMILES match
            - Files without SMILES are kept
            - Warning is logged for duplicate structures
        """
        unique_files = []
        unique_smiles = []
        for file in self.args.files:
            smi = None
            with open(file, "r", encoding='utf-8') as F:
                lines = F.readlines()
                smi_exist = False
                for i, line in enumerate(lines):
                    if ">  <SMILES>" in line:
                        smi = lines[i + 1].split()[0]
                        if smi not in unique_smiles:
                            unique_smiles.append(smi)
                            unique_files.append(file)
                            smi_exist = True
                if smi_exist:
                    continue
                elif smi is not None:
                    self.args.log.write(f'x  WARNING! "{os.path.basename(file)}" will not be calculated since it has the same SMILES as "{os.path.basename(unique_files[unique_smiles.index(smi)])}"')

        if not unique_smiles:
            unique_files = self.args.files
        return unique_files
        
    def _error_exit(self, message):
        """Log error message and exit the program.
        
        Args:
            message (str): Error message to log before exiting
        """
        self.args.log.write(f"\nx  {message}")
        self.args.log.finalize()
        sys.exit()
        
    def _validate_inputs(self):
        """Validate input parameters for QDESCP calculation.
        
        This method checks:
        - CSV file existence
        - Solvent compatibility
        - Input file validity
        
        Raises:
            SystemExit: If validation fails
        """
        # Check CSV file existence
        if self.args.csv_name is not None:
            csv_path = Path(self.args.csv_name)
            if not csv_path.exists():
                self._error_exit(
                    f"CSV file {csv_path} not found. Please verify the path."
                )
                
        # Check solvent compatibility
        if self.args.qdescp_solvent is not None:
            self._error_exit(
                "PTB calculations do not support solvents. "
                "Please remove the --qdescp_solvent option."
            )
        
    def _convert_to_xyz(self, file, name, ext):
        """Convert input file to XYZ format using OpenBabel.
        
        Args:
            file (str): Input file path
            name (str): Base name for output files
            ext (str): Input file extension
        """
        command_pdb = [
            "obabel",
            f"-i{ext.lower()}",
            file,
            "-oxyz", 
            f"-O{os.path.dirname(os.path.abspath(file))}/{name}_conf_.xyz",
            "-m"
        ]
        subprocess.run(
            command_pdb,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    
    def process_aqme_csv(self, name_db):
        """Process and update AQME CSV result files.
        
        This method processes three types of CSV files:
        - full: All descriptors
        - denovo: De novo descriptors only
        - interpret: Interpretable descriptors only
        
        For each file:
        1. Reads original and generated CSVs
        2. Adds missing entries
        3. Sorts according to input order
        4. Fills missing values
        
        Args:
            name_db (str): Database name prefix for output files
            
        Note:
            - Missing entries are filled with group-wise first values
            - Original file order is preserved
            - Warns if files are not found
        """
        # Process each of the three generated CSVs
        for suffix in ['full', 'denovo', 'interpret']:
            # Try to read the file with the corresponding suffix
            base_filename = os.path.basename(self.args.csv_name)
            csv_file = self.args.initial_dir.joinpath(f'AQME-{name_db}_{suffix}_{base_filename}')
            csv_file = f'{csv_file}'
            
            try:
                # Read the original CSV and the one generated by qdescp.
                csv_temp = pd.read_csv(f'{self.args.csv_name}')
                df_temp = pd.read_csv(csv_file)

                # check that code_name and SMILES are written with the right format
                csv_temp = fix_cols_names(csv_temp)
                df_temp = fix_cols_names(df_temp)

                # Compare and add missing rows.
                if len(df_temp) < len(csv_temp):
                    missing_rows = csv_temp.loc[~csv_temp['code_name'].isin(df_temp['code_name'])]
                    missing_rows[['code_name', 'SMILES']].to_csv(csv_file, mode='a', header=False, index=False)

                    # Sort the data according to the order of 'code_name'
                    order = csv_temp['code_name'].tolist()
                    df_temp = df_temp.sort_values(by='code_name', key=lambda x: x.map({v: i for i, v in enumerate(order)}))
                    df_temp = df_temp.fillna(df_temp.groupby('SMILES').transform('first'))

                    # Overwrite the CSV file with the processed data
                    df_temp.to_csv(csv_file, index=False)

            except FileNotFoundError:
                self.args.log.write(f"Not found {csv_file}. Please check if the file was generated correctly.")


    def qdescp_nmr_workflow(self, boltz_dir):
        """Run NMR workflow for chemical shift prediction.
        
        This method:
        1. Validates input files (must be JSON format)
        2. Processes conformer JSON files
        3. Calculates Boltzmann-weighted NMR properties
        4. Applies empirical corrections if provided
        
        Args:
            boltz_dir (Path): Directory for Boltzmann-averaged results
            
        Raises:
            SystemExit: If input files are not in JSON format
        """
        
        self.args.log.write(f"\nStarting QDESCP-{self.args.program} with {len(self.args.files)} job(s)\n")

        atom_props = ["NMR Chemical Shifts"]
        
        if os.path.basename(Path(self.args.files[0])).split('.')[-1].lower() not in ["json"]:
            self.args.log.write(f"\nx  The format used ({os.path.basename(Path(self.args.files[0])).split('.')[-1]}) is not compatible with QDESCP with NMR! Formats accepted: json")
            self.args.log.finalize()
            sys.exit()

        name = os.path.basename(Path(self.args.files[0])).split("_conf")[0]
        
        json_files = glob.glob(str(os.path.dirname(Path(self.args.files[0]))) + "/" + name + "_conf_*.json")
        get_boltz_props_nmr(json_files, name, boltz_dir, self, atom_props, self.args.nmr_atoms, self.args.nmr_slope, self.args.nmr_intercept, self.args.nmr_experim)
