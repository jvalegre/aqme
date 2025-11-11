"""
Parameters
----------

General
+++++++

   input : str, default=''  
      (If smi is None) Optionally, file containing the SMILES strings and 
      names of the molecules. Current file extensions: .smi, .sdf, .cdx, 
      .csv, .com, .gjf, .mol, .mol2, .xyz, .txt, .yaml, .yml, .rtf  
      For .csv files (i.e. FILENAME.csv), two columns are required, 
      'code_name' with the names and 'SMILES' for the SMILES string  
   program : str, default=None  
      Program required in the conformational sampling. 
      Current options: 'rdkit', 'crest'  
   smi : str, default=None  
      Optionally, define a SMILES string as input  
   name : str, default=None  
      (If smi is defined) optionally, define a name for the system  
   w_dir_main : str, default=os.getcwd()  
      Working directory 
   destination : str, default=None,
     Directory to create the output file(s)   
   varfile : str, default=None  
      Option to parse the variables using a yaml file (specify the filename)  
   charge : int, default=None  
      Charge of the calculations used in the following input files. 
      If charge isn't defined, it automatically reads the charge of the 
      SMILES string  
   mult : int, default=None  
      Multiplicity of the calculations used in the following input files. If 
      mult isn't defined, it automatically reads the multiplicity of the mol 
      object created with the SMILES string. Be careful with the automated 
      calculation of mult from mol objects when using metals!  
   prefix : str, default=''  
      Prefix added to all the names  
   suffix : str, default=''  
      Suffix added to all the names  
   stacksize : str, default='1G'  
      Controls the stack size used (especially relevant for xTB/CREST 
      calculations of large systems, where high stack sizes are needed)  

General RDKit-based
+++++++++++++++++++

   sample : int, default=25
      Number of conformers to keep after the initial RDKit sampling. They are selected using a
      combination of RDKit energies and Butina clustering
   auto_sample : str, default=mid in CSEARCH, low in QDESCP
      Apply automatic calculation of the number of conformers generated initially with RDKit. This number
      of conformers is initially generated and then reduced to the number specified in --sample with
      different filters. There is a sampling factor, which is multiplied by the number of rotatable bonds,
      and a maximum number of conformers allowed to pass to the filters. Options:
      1. Low: good for descriptor generation in machine learning. Base multiplier = 5, max number of confs = 100
      2. Mid: standard, good compromise between number of conformers and computing time. Base multiplier = 10, max number of confs = 250
      3. High: demanding method, more conformers and time. Base multiplier = 20, max number of confs = 500
      4. False: use the number of conformers specified in --sample
   ff : str, default='MMFF'
      Force field used in RDKit optimizations and energy calculations. Current 
      options: MMFF, UFF (if MMFF fails, AQME tries to use UFF automatically), and NO FF (works well with metals
      when UFF doesn't work)
   ewin_csearch : float, default=5.0
      Energy window in kcal/mol to discard conformers (i.e. if a conformer is 
      more than the E window compared to the most stable conformer)
   initial_energy_threshold : float, default=0.0001
      Energy difference in kcal/mol between unique conformers for the first 
      filter of only E
   energy_threshold : float, default=0.25
      Energy difference in kcal/mol between unique conformers for the second 
      filter of E + RMS
   rms_threshold : float, default=0.25
      RMS difference between unique conformers for the second filter of E + RMS
   opt_steps_rdkit : int, default=1000
      Max cycles used in RDKit optimizations
   heavyonly : bool, default=True
      Only consider heavy atoms during RMS calculations for filtering (in the 
      Chem.rdMolAlign.GetBestRMS() RDKit function)
   max_matches_rmsd : int, default=1000
      Max matches during RMS calculations for filtering (maxMatches option in 
      the Chem.rdMolAlign.GetBestRMS() RDKit function)
   max_mol_wt : int, default=0
      Discard systems with molecular weights higher than this parameter 
      (in g/mol). If 0 is set, this filter is off
   max_torsions : int, default=0
      Discard systems with more than this many torsions (relevant to avoid 
      molecules with many rotatable bonds). If 0 is set, this filter is off
   seed : int, default=62609
      Random seed used during RDKit embedding (in the 
      Chem.rdDistGeom.EmbedMultipleConfs() RDKit function)
   geom : list, default=[]
      Geometry rule to pass for the systems. Format: [SMARTS,VALUE]. Geometry rules
      might be atoms, bonds, angles and dihedral. For example, a rule to keep only
      molecules with C-Pd-C atoms at 180 degrees: ['[C][Pd][C]',180]. Multiple rules
      can be used at the same time (['C[Pd]C',180,'C[Pd]N',90]).
      Special rules (--geom ['RULE_NAME']):
         1. ['Ir_squareplanar']
   bond_thres : float, default=0.2
      Threshold used to discard bonds in the geom option (+-0.2 A) 
   angle_thres : float, default=30
      Threshold used to discard angles in the geom option (+-30 degrees) 
   dihedral_thres : float, default=30
      Threshold used to discard dihedral angles in the geom option (+-30 degrees) 

Only organometallic molecules
.............................

   auto_metal_atoms : bool, default=True
     Automatically detect metal atoms for the RDKit conformer generation. Charge 
     and mult should be specified as well since the automatic charge and mult 
     detection might not be precise.
   complex_type : str, default=''
      Forces the metal complexes to adopt a predefined geometry. This option is 
      especially relevant when RDKit predicts wrong complex geometries or gives 
      a mixture of geometries. Current options: squareplanar, squarepyramidal, 
      linear, trigonalplanar
   single_system : bool, default=False
      When using complex_type templates in CSEARCH, keep only one system of all the options. 
      This option is useful to avoid repetition when the complex has two identical
      ligands (i.e. two Cl substituents).

CREST only
++++++++++

   nprocs : int, default=8
      Number of processors used in CREST optimizations
   constraints_atoms : list, default=[]
      Specify constrained atoms as [AT1,AT2,AT3]. An example of multiple constraints with
      atoms 1, 2 and 5 frozen: [1,2,5]
   constraints_dist : list of lists, default=[]
      Specify distance constraints as [AT1,AT2,DIST]. An example of multiple constraints with
      atoms 1 and 2 frozen at a distance of 1.8 Å, and atoms 4 and 5 with distance of 2.0 Å:
      [[1,2,1.8],[4,5,2.0]]
   constraints_angle : list of lists, default=[]
      Specify angle constraints as [AT1,AT2,AT3,ANGLE]. An example of multiple constraints with
      atoms 1, 2 and 3 frozen at an angle of 180 degrees, and atoms 4, 5 and 6 with an angle of 120:
      [[1,2,3,180],[4,5,6,120]]
   constraints_dihedral : list of lists, default=[]
      Specify dihedral constraints as [AT1,AT2,AT3,AT4,DIHEDRAL]. An example of multiple constraints
      with atoms 1, 2, 3 and 4 frozen at a dihedral angle of 180 degrees, and atoms 4, 5, 6 and 7
      with a dihedral angle of 120: [[1,2,3,4,180],[4,5,6,7,120]]
   crest_force : float, default=0.5
      Force constant for constraints in the .xcontrol.sample file for CREST jobs
   crest_keywords : str, default=None
      Define additional keywords to use in CREST that are not included in --chrg, 
      --uhf, -T and -cinp. For example: '--alpb ch2cl2 --nci --cbonds 0.5'
   cregen : bool, default=True
      If True, perform a CREGEN analysis after CREST (filtering options below)
   cregen_keywords : str, default=None
      Additional keywords for CREGEN (i.e. cregen_keywords='--ethr 0.02')
   xtb_keywords : str, default=None
      Define additional keywords to use in the xTB pre-optimization that are not 
      included in -c, --uhf, -P and --input. For example: '--alpb ch2cl2 --gfn 1'
   crest_runs : int, default=1
      Specify as number of runs if multiple starting points from RDKit starting points is required.
"""
# Standard library imports
import concurrent.futures as futures
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Third-party imports
from progress.bar import IncrementalBar
from rdkit.Chem import (
    AllChem as Chem,
    Descriptors,
    Lipinski,
    PropertyMol,
    rdDistGeom,
    rdmolfiles,
)

# Local application imports
from aqme.filter import (
    cluster_conformers,
    conformer_filters,
    filters,
    geom_filter,
)
from aqme.csearch.utils import (
    check_constraints,
    com_2_xyz,
    getDihedralMatches,
    minimize_rdkit_energy,
    prepare_cdx_files,
    prepare_com_files,
    prepare_csv_files,
    prepare_direct_smi,
    prepare_pdb_files,
    prepare_sdf_files,
    prepare_smiles_files,
    realign_mol,
    smi_to_mol,
    substituted_mol,
    set_metal_atomic_number
)
from aqme.csearch.templates import check_metal_neigh, template_embed
from aqme.utils import (
    blocking_wrapper,
    check_crest,
    check_dependencies,
    check_xtb,
    get_files,
    load_sdf,
    load_variables,
    mol_from_sdf_or_mol_or_mol2,
    set_destination,
)
from aqme.csearch.crest import xtb_opt_main


class csearch:
    """Handles geometry generation and conformational search.
    
    This class provides functionality for:
    1. Multiple conformer generation using RDKit or CREST
    2. Structure optimization and refinement
    3. Conformer filtering and clustering
    4. Support for various input formats (SMILES, SDF, MOL2, etc.)
    5. Parallel processing capabilities
    
    For detailed parameter documentation, see module documentation.
    """
    
    SUPPORTED_PROGRAMS = {"rdkit", "crest"}
    SUPPORTED_FORCEFIELDS = {"MMFF", "UFF", "NO FF"}
    DEFAULT_NPROCS = 4
    DEFAULT_AUTO_SAMPLE = "mid"
    
    def __init__(self, **kwargs):
        """Initialize csearch with the given configuration.
        
        Args:
            **kwargs: Configuration parameters, see class docstring
        
        Raises:
            SystemExit: If required parameters are missing or invalid
        """
        self.start_time = time.time()
        self._initialize_args(kwargs)
        self._validate_configuration()
        self._process_input_files()  # Process inputs before setting up output directory
        self._cleanup()
        
    def _initialize_args(self, kwargs):
        """Initialize and validate basic configuration.
        
        Args:
            kwargs: Configuration parameters
        """
        self.args = load_variables(kwargs, "csearch")
        check_dependencies(self)
        
        # Set default values
        if self.args.nprocs is None:
            self.args.nprocs = self.DEFAULT_NPROCS
        if self.args.auto_sample == 'auto':
            self.args.auto_sample = self.DEFAULT_AUTO_SAMPLE
            
    def _validate_configuration(self):
        """Validate program configuration and dependencies."""
        self._validate_program()
        self._validate_forcefield()
        self._validate_input()
        
        if self.args.program.lower() == "crest":
            check_xtb(self)
            check_crest(self)
            
    def _validate_program(self):
        """Verify program selection is valid."""
        program = self.args.program
        if not program or program.lower() not in self.SUPPORTED_PROGRAMS:
            self._error_exit(
                'Program not specified or not supported for CSEARCH! '
                'Specify: program="rdkit" (or "crest")'
            )
            
    def _validate_forcefield(self):
        """Verify force field selection is valid."""
        if self.args.ff.upper() not in self.SUPPORTED_FORCEFIELDS:
            self._error_exit(f"Force field {self.args.ff} not supported!")
            
    def _validate_input(self):
        """Validate input specification."""
        has_smiles = self.args.smi is not None
        has_input = self.args.input != ""
        
        if not has_smiles and not has_input:
            self._error_exit(
                "Program requires either a SMILES or an input file to proceed! "
                "Please look up acceptable file formats. Specify: smi='CCC' "
                "(or input='filename.csv')"
            )
        elif has_smiles and has_input:
            self._error_exit(
                "Program requires either a SMILES or an input file to proceed, "
                "don't use both!"
            )
            
        # Set dummy extension if using SMILES input
        if has_smiles:
            self.args.input = 'no_ext.no_ext'
            
    def _process_input_files(self):
        """Process input files and run conformer search."""
        # Store original directory
        original_dir = os.getcwd()
        
        try:
            if self.args.smi is not None:
                csearch_files = [self.args.name]
            else:
                # Handle input path
                if os.path.isabs(self.args.input):
                    # If absolute path, use it directly and change to its directory
                    input_dir = os.path.dirname(self.args.input)
                    input_file = os.path.basename(self.args.input)
                    os.chdir(input_dir)
                    self.args.input = input_file
                else:
                    # If relative path, resolve from initial directory
                    full_path = os.path.join(self.args.initial_dir, self.args.input)
                    input_dir = os.path.dirname(full_path)
                    input_file = os.path.basename(full_path)
                    os.chdir(input_dir)
                    self.args.input = input_file
                
                # Get list of input files
                csearch_files = self._get_input_files()
                
            # Process each file
            for csearch_file in csearch_files:
                self._process_single_file(csearch_file)
                
        finally:
            # Return to original directory
            os.chdir(original_dir)
            # Setup working directory for outputs
            self._setup_working_directory()
            
    def _setup_working_directory(self):
        """Set up and validate working directory."""
        try:
            if not os.path.isabs(self.args.w_dir_main):
                self.args.w_dir_main = os.path.join(self.args.initial_dir, self.args.w_dir_main)
            if Path(self.args.w_dir_main).exists():
                os.chdir(self.args.w_dir_main)
            else:
                os.makedirs(self.args.w_dir_main, exist_ok=True)
                os.chdir(self.args.w_dir_main)
        except Exception as e:
            self._error_exit(f"Could not set up working directory: {str(e)}")
            
    def _get_input_files(self):
        """Get list of input files to process.
        
        Returns:
            list: Input file paths
            
        Raises:
            SystemExit: If no input files are found
        """
        files = get_files(self.args.input)
        if not files:
            self._error_exit(f"Input file ({self.args.input}) not found!")
        return files
        
    def _process_single_file(self, csearch_file):
        """Process a single input file.
        
        Args:
            csearch_file: Path to input file
        """
        
        # Prepare job inputs
        job_inputs = (prepare_direct_smi(self.args) if self.args.smi is not None 
                     else self.load_jobs(csearch_file))
                     
        self.args.log.write(
            f"\nStarting CSEARCH with {len(job_inputs)} job(s) "
            "(SDF, XYZ, CSV, etc. files might contain multiple jobs/structures inside)\n"
        )
        
        # Run conformer search
        self.run_csearch(job_inputs)
        
        # Clean up empty output files
        self._cleanup_empty_files()
            
    def _cleanup_empty_files(self):
        """Remove any empty output SDF files."""
        for sdf_file in glob.glob(f'{self.csearch_folder}/*.sdf'):
            if os.path.getsize(sdf_file) == 0:
                os.remove(sdf_file)
                
    def _cleanup(self):
        """Perform final cleanup and logging."""
        elapsed_time = round(time.time() - self.start_time, 2)
        self.args.log.write(f"\nTime CSEARCH: {elapsed_time} seconds\n")
        self.args.log.finalize()
        
        # Restore initial directory (for Jupyter notebooks)
        os.chdir(self.args.initial_dir)
        
    def _error_exit(self, message):
        """Log error message and exit.
        
        Args:
            message: Error message to log
        """
        self.args.log.write(f"\nx  {message}")
        self.args.log.finalize()
        sys.exit()

    SUPPORTED_INPUTS = {
        # SMILES-based formats
        "smi", "txt", "yaml", "yml", "rtf",
        # Structure formats
        "sdf", "mol", "mol2", "pdb",
        # Gaussian formats
        "com", "gjf", "xyz",
        # Other formats
        "cdx", "csv"
    }
    
    EXTENSION_TO_HANDLER = {
        # SMILES-based formats
        "smi": prepare_smiles_files,
        "txt": prepare_smiles_files,
        "yaml": prepare_smiles_files,
        "yml": prepare_smiles_files,
        "rtf": prepare_smiles_files,
        # Structure formats
        "sdf": prepare_sdf_files,
        "mol": prepare_sdf_files,
        "mol2": prepare_sdf_files,
        "pdb": prepare_pdb_files,
        # Gaussian formats
        "com": prepare_com_files,
        "gjf": prepare_com_files,
        "xyz": prepare_com_files,
        # Other formats
        "cdx": prepare_cdx_files,
        "csv": prepare_csv_files
    }
    
    def load_jobs(self, csearch_file):
        """Load molecular information for conformer generation.
        
        This method:
        1. Validates input file format
        2. Maps file extension to appropriate handler
        3. Prepares job inputs for conformer generation
        
        Args:
            csearch_file (str or Path): Path to input file
            
        Returns:
            list: List of job inputs for conformer generation
            
        Raises:
            SystemExit: If file format is unsupported or file is not found
        """
        file_format = self._get_file_format(csearch_file)
        self._validate_file_format(file_format)
        
        try:
            prepare_function = self.EXTENSION_TO_HANDLER[file_format]
            return prepare_function(self.args, csearch_file)
        except FileNotFoundError:
            self._error_exit(
                f'File {os.path.basename(csearch_file)} was not found! '
                'In the "input" option, make sure that:\n'
                '1) The PATH to the files is correct\n'
                '2) The PATH doesn\'t start with "/"'
            )
            
    def _get_file_format(self, file_path):
        """Extract file format from file path.
        
        Args:
            file_path (str or Path): Path to input file
            
        Returns:
            str: Lowercase file extension
        """
        return os.path.basename(Path(file_path)).split('.')[-1].lower()
        
    def _validate_file_format(self, file_format):
        """Check if file format is supported.
        
        Args:
            file_format (str): File extension to validate
            
        Raises:
            SystemExit: If format is not supported
        """
        if file_format not in self.SUPPORTED_INPUTS:
            self._error_exit(f"Input filetype (.{file_format}) not currently supported!")

    def run_csearch(self, job_inputs):
        """Run conformer search on all input jobs.
        
        This method handles parallel processing of conformer generation jobs:
        1. Sets up progress tracking
        2. Executes jobs in parallel if possible
        3. Handles job exceptions gracefully
        
        Args:
            job_inputs (list): List of job configurations
        """
        bar = IncrementalBar(
            "o  Number of finished jobs from CSEARCH", 
            max=len(job_inputs)
        )
        
        try:
            self._execute_jobs(job_inputs, bar)
        finally:
            bar.finish()
            
    def _execute_jobs(self, job_inputs, progress_bar):
        """Execute conformer generation jobs.
        
        Args:
            job_inputs (list): List of job configurations
            progress_bar (IncrementalBar): Progress tracking
        """
        use_parallel = (not self.args.debug and 
                       self.args.program.lower() != 'crest')
                       
        if use_parallel:
            self._run_parallel_jobs(job_inputs, progress_bar)
        else:
            self._run_sequential_jobs(job_inputs, progress_bar)
            
    def _run_parallel_jobs(self, job_inputs, progress_bar):
        """Execute jobs in parallel using thread pool.
        
        Args:
            job_inputs (list): List of job configurations
            progress_bar (IncrementalBar): Progress tracking
        """
        with futures.ThreadPoolExecutor(max_workers=self.args.nprocs) as executor:
            csearch_nprocs = 1
            future_tasks = [
                executor.submit(
                    blocking_wrapper,
                    self.compute_confs,
                    job_input,
                    progress_bar,
                    csearch_nprocs
                )
                for job_input in job_inputs
            ]

            # Wait for all tasks to complete
            done, _ = futures.wait(future_tasks, return_when=futures.ALL_COMPLETED)

            # Handle exceptions
            for future in done:
                try:
                    future.result()  # Raises any exceptions caught during execution
                except Exception as e:
                    self.args.log.write(f"CSEARCH raised an exception: {e}")
                    
    def _run_sequential_jobs(self, job_inputs, progress_bar):
        """Execute jobs sequentially.
        
        Args:
            job_inputs (list): List of job configurations
            progress_bar (IncrementalBar): Progress tracking
        """
        for job_input in job_inputs:
            _ = self.compute_confs(job_input, progress_bar, self.args.nprocs)
            
    def compute_confs(self, job_input, progress_bar, nprocs):
        """Generate conformers for a single molecule.
        
        This method:
        1. Processes input molecule data
        2. Converts to RDKit molecule if needed
        3. Handles 3D input formats
        4. Sets up conformer generation
        
        Args:
            job_input (tuple): Job configuration parameters
            progress_bar (IncrementalBar): Progress tracking
            nprocs (int): Number of processors to use
            
        Returns:
            None: Updates are made to the filesystem
        """
        try:
            # Unpack job input and store as instance attributes
            (smi, name, charge, mult, 
             constraints_atoms, constraints_dist,
             constraints_angle, constraints_dihedral, 
             complex_type, geom, sample) = job_input

            csearch_nprocs = nprocs
            valid_template_embed = True
            
            self.args.log.write(
                f"\n   ----- {os.path.basename(Path(name))} -----"
            )
            
            # Process molecule
            is_smiles_input = (
            isinstance(smi,str) or 
            self._is_smiles_format(self.args.input)
                )
            if is_smiles_input:
                (
                    mol,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                    complex_ts
                ) = smi_to_mol(
                    smi, 
                    self.args.program.lower(),
                    self.args.log,
                    self.args.seed,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral
                )
            else:
                mol = smi
                complex_ts = check_constraints(constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral)
                
            # Setup for conformer generation
            self._setup_output_directory()
            self._handle_3d_input(name)

            # Process metals if present
            metal_atoms = self._detect_metal_atoms(mol, charge, mult, name)
            
            if metal_atoms:
                # For metal-containing molecules, substitute metals with iodine first
                metal_idx, metal_sym = substituted_mol(mol, "I", metal_atoms)
                
                # Try metal template first
                if complex_type:
                    valid_template_embed = self._process_metal_complex(
                        mol, name, metal_atoms, metal_idx, complex_type, metal_sym, valid_template_embed,
                        constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral,
                        complex_ts, charge, mult, smi, geom, csearch_nprocs, sample
                    )
                
                # If no template or template failed, use standard conformer generation
                if not complex_type or not valid_template_embed:
                    self.conformer_generation(
                        mol=mol,
                        name=name,
                        constraints_atoms=constraints_atoms,
                        constraints_dist=constraints_dist,
                        constraints_angle=constraints_angle,
                        constraints_dihedral=constraints_dihedral,
                        complex_ts=complex_ts,
                        charge=charge,
                        mult=mult,
                        smi=smi,
                        geom=geom,
                        metal_atoms=metal_atoms,
                        metal_idx=metal_idx,
                        metal_sym=metal_sym,
                        csearch_nprocs=csearch_nprocs,
                        sample=sample
                    )
            else:
                # For regular molecules
                self.conformer_generation(
                    mol=mol,
                    name=name,
                    constraints_atoms=constraints_atoms,
                    constraints_dist=constraints_dist,
                    constraints_angle=constraints_angle,
                    constraints_dihedral=constraints_dihedral,
                    complex_ts=complex_ts,
                    charge=charge,
                    mult=mult,
                    smi=smi,
                    geom=geom,
                    metal_atoms=[],
                    metal_idx=[],
                    metal_sym=[],
                    csearch_nprocs=csearch_nprocs,
                    sample=sample
                )

            return None
            
        finally:
            progress_bar.next()

    def _is_smiles_format(self, input_path):
        """Check if input file is in SMILES-based format.
        
        Args:
            input_path (str): Path to input file
            
        Returns:
            bool: True if file is SMILES-based format
        """
        smiles_formats = {"smi", "csv", "cdx", "txt", "yaml", "yml", "rtf"}
        ext = os.path.basename(Path(input_path)).split(".")[-1].lower()
        return ext in smiles_formats
            
    def _setup_output_directory(self):
        """Create output directory for conformer files."""
        self.csearch_folder = set_destination(self, 'CSEARCH')
        self.csearch_folder.mkdir(exist_ok=True, parents=True)
        
    def _handle_3d_input(self, name):
        """Process 3D input structures for CREST.
        
        Args:
            name (str): Molecule name/identifier
        """
        if (self.args.program.lower() == "crest" and 
            self.args.smi is None):
            
            input_ext = os.path.basename(Path(self.args.input)).split(".")[-1]
            
            if input_ext in ["pdb", "mol2", "mol", "sdf"]:
                self._convert_to_xyz_obabel(name, input_ext)
            elif input_ext in ["gjf", "com"]:
                self._convert_to_xyz_com(name, input_ext)
            elif input_ext == "xyz":
                self._copy_xyz(name)
                
    def _convert_to_xyz_obabel(self, name, ext):
        """Convert structure file to XYZ using OpenBabel.
        
        Args:
            name (str): Molecule name/identifier
            ext (str): Input file extension
        """
        command = [
            "obabel",
            f'-i{ext}',
            f'{name}.{ext}',
            "-oxyz",
            f"-O{name}_{self.args.program.lower()}.xyz"
        ]
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
    def _convert_to_xyz_com(self, name, ext):
        """Convert Gaussian input to XYZ format.
        
        Args:
            name (str): Molecule name/identifier
            ext (str): Input file extension
        """
        xyz_file, _, _ = com_2_xyz(f'{name}.{ext}')
        shutil.move(xyz_file, f"{name}_{self.args.program.lower()}.xyz")
        
    def _copy_xyz(self, name):
        """Copy XYZ file for CREST input.
        
        Args:
            name (str): Molecule name/identifier
        """
        shutil.copy(f"{name}.xyz", f"{name}_{self.args.program.lower()}.xyz")
        
    def _detect_metal_atoms(self, mol, charge, mult, name):
        """Detect metal atoms in molecule.
        
        This method:
        1. Detects metal atoms if auto-detection is enabled
        2. Returns list of detected metals for further processing
        
        Args:
            mol: RDKit molecule object
            charge (int): Charge of the molecule
            mult (int): Multiplicity of the molecule
            name (str): Name of the molecule

        Returns:
            list: Detected metal atoms
        """
        metal_atoms = []
        if self.args.auto_metal_atoms:
            metal_atoms = self.find_metal_atom(mol, charge, mult, name)
        return metal_atoms
        
    def _process_metal_complex(self, mol, name, metal_atoms, metal_idx, complex_type, metal_sym, valid_template_embed,
                                constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral,
                                complex_ts, charge, mult, smi, geom, csearch_nprocs, sample):
        """Process metal complex with templates.
        
        Args:
            mol: RDKit molecule object 
            name (str): Name of the molecule
            metal_atoms (list): Metal atom indices
            metal_idx (list): Metal indices after substitution
            complex_type (str): Complex type
            metal_sym (list): Metal symbols
            valid_template_embed (bool): Whether the embedding worked with templates for metal complexes
        """
        # Validate template type
        if complex_type not in self.ACCEPTED_COMPLEX_TYPES:
            self._log_invalid_template()
            return False
            
        # Check template applicability
        count_metals = 0
        valid_template = True
        
        for idx in metal_idx:
            if idx is not None:
                valid_template = check_metal_neigh(
                    mol, complex_type, idx,
                    self.args.log, valid_template
                )
                count_metals += 1
                
        if count_metals == 1 and valid_template:
            valid_template_embed = self._apply_template(mol, metal_atoms, metal_idx, metal_sym, name, complex_type, valid_template_embed,
                                                        constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral,
                                                        complex_ts, charge, mult, smi, geom, csearch_nprocs, sample)
        
        return valid_template_embed
            
    def _apply_template(self, mol, metal_atoms, metal_idx, metal_sym, name, complex_type, valid_template_embed,
                        constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral,
                        complex_ts, charge, mult, smi, geom, csearch_nprocs, sample):
        """Apply template to metal complex and generate conformers.
        
        Args:
            mol: RDKit molecule object
            metal_atoms (list): Metal atom indices
            metal_idx (list): Metal indices after substitution
            metal_sym (list): Metal symbols
            name (str): Name of the molecule
            complex_type (str): Geometry type of the metal complex
            valid_template_embed (bool): Whether the embedding worked with templates for metal complexes
       """
        template_kwargs = {
            "complex_type": complex_type,
            "metal_idx": metal_idx,
            "maxsteps": self.args.opt_steps_rdkit,
            "heavyonly": self.args.heavyonly,
            "maxmatches": self.args.max_matches_rmsd,
            "mol": mol,
            "name": name
        }
        
        try:
            items = template_embed(self, **template_kwargs)
            
            for mol_obj, name_in, coord_map, alg_map, template, original_atn in zip(*items):
                # Handle single system option
                if (self.args.single_system and 
                    glob.glob(f'{self.csearch_folder}/{name}_*.sdf')):
                    break

                self.conformer_generation(
                    mol=mol_obj,
                    name=name_in,
                    constraints_atoms=constraints_atoms,
                    constraints_dist=constraints_dist,
                    constraints_angle=constraints_angle,
                    constraints_dihedral=constraints_dihedral,
                    complex_ts=complex_ts,
                    charge=charge,
                    mult=mult,
                    smi=smi,
                    geom=geom,
                    metal_atoms=metal_atoms,
                    metal_idx=metal_idx,
                    metal_sym=metal_sym,
                    csearch_nprocs=csearch_nprocs,
                    sample=sample,
                    coord_Map=coord_map,
                    alg_Map=alg_map,
                    mol_template=template,
                    original_atn=original_atn
                )
                
        except RuntimeError:
            valid_template_embed = False
            self._log_template_error(name, complex_type)
        
        return valid_template_embed
            
    def _log_invalid_template(self, complex_type, name):
        """Log message for invalid template type."""
        self.args.log.write(
            f"\nx  The metal template specified in complex_type "
            f"({complex_type}) is not valid! Options: "
            f"{', '.join(self.ACCEPTED_COMPLEX_TYPES)} "
            f"({os.path.basename(Path(name))})"
        )
            
    def _log_template_error(self, name, complex_type):
        """Log message for template application error."""
        self.args.log.write(
            f"\nx  Molecule {name} was not optimized in the "
            f"specified template ({complex_type}). Try using "
            f"charges to make it easier for MM protocols (i.e. using "
            f"[NH3+][Ag][NH3+] instead of [NH3][Ag][NH3], since N+ "
            f"typically has 4 bonds)."
        )
        
    # Pre-defined metal complex geometries
    ACCEPTED_COMPLEX_TYPES = {
        "squareplanar", "squarepyramidal",
        "linear", "trigonalplanar"
    }
            
    # List of transition metals for automatic detection
    TRANSITION_METALS = {
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
        'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
    }

    def find_metal_atom(self, mol, charge, mult, name):
        """Detect transition metal atoms in molecule.
        
        This method scans through all atoms in the molecule to identify
        transition metals and warns about potential charge/multiplicity
        issues for metal complexes.
        
        Args:
            mol: RDKit molecule object to analyze
            charge (int): Molecular charge (can be None)
            mult (int): Multiplicity (can be None) 
            name (str): Molecule name for logging
            
        Returns:
            list: Symbols of detected metal atoms
            
        Note:
            Issues warnings if charge/multiplicity are not explicitly 
            specified for metal-containing molecules
        """
        # Find metal atoms in molecule
        metal_atoms = [
            atom.GetSymbol() 
            for atom in mol.GetAtoms()
            if atom.GetSymbol() in self.TRANSITION_METALS
        ]
        
        if metal_atoms:
            # Log detected metals
            self._log_metal_detection(metal_atoms, name)
            
            # Warn about charge/multiplicity if not specified
            if charge is None or mult is None:
                self._warn_about_metal_params(name, charge_warning=charge is None, 
                                           mult_warning=mult is None)
                
        return metal_atoms
        
    def _log_metal_detection(self, metal_atoms, name):
        """Log detected metal atoms.
        
        Args:
            metal_atoms (list): Detected metal symbols
            name (str): Molecule name
        """
        self.args.log.write(
            f"\no  AQME recognized the following metal atoms: "
            f"{metal_atoms} ({os.path.basename(Path(name))})"
        )
        
    def _warn_about_metal_params(self, name, charge_warning=False, 
                               mult_warning=False):
        """Warn about unspecified charge/multiplicity for metal complexes.
        
        Args:
            name (str): Molecule name
            charge_warning (bool): Whether to warn about charge
            mult_warning (bool): Whether to warn about multiplicity 
        """
        base_name = os.path.basename(Path(name))
        
        if charge_warning:
            self.args.log.write(
                f"\nx  The automated charge calculation might not be precise "
                f"for metal complexes! You should use the charge option "
                f"(or the charge column in CSV inputs) ({base_name})"
            )
            
        if mult_warning:
            self.args.log.write(
                f"\nx  The automated multiplicity calculation might not be "
                f"precise for metal complexes! You should use the mult option "
                f"(or the mult column in CSV inputs) ({base_name})"
            )
    
    def conformer_generation(
        self, mol, name, constraints_atoms, constraints_dist,
        constraints_angle, constraints_dihedral, complex_ts,
        charge, mult, smi, geom, metal_atoms, metal_idx,
        metal_sym, csearch_nprocs, sample, coord_Map=None,
        alg_Map=None, mol_template=None, original_atn=None
    ):
        """Generate 3D conformers for a molecule.
        
        This method handles conformer generation using either CREST or RDKit
        depending on the input type and program selection.
        
        Args:
            mol: RDKit molecule object
            name (str): Molecule identifier
            constraints_atoms (list): Constrained atom indices
            constraints_dist (list): Distance constraints 
            constraints_angle (list): Angle constraints
            constraints_dihedral (list): Dihedral angle constraints
            complex_ts (bool): Whether molecule is a complex/transition state
            charge (int): Molecular charge
            mult (int): Multiplicity
            smi (str): SMILES string or 3D structure
            geom (list): Geometry constraints
            metal_atoms (list): Metal atoms in molecule
            metal_idx (list): Metal atom indices
            metal_sym (list): Metal atom symbols
            csearch_nprocs (int): Number of processors
            sample (int): Number of conformers to generate
            coord_Map: Coordinate mapping (optional)
            alg_Map: Alignment mapping (optional) 
            mol_template: Template molecule (optional)
            original_atn: Original atomic numbers (optional)
            
        Returns:
            None: Results written to files
        """
        # Set default charge/multiplicity if not provided
        charge = charge if charge is not None else Chem.GetFormalCharge(mol)
        mult = mult if mult is not None else (Descriptors.NumRadicalElectrons(mol) + 1)
        # Generate conformers using appropriate method
        if self._should_use_crest():
            valid_structure = True
            status = self._run_crest_sampling(
                name, charge, mult, smi, constraints_atoms,
                constraints_dist, constraints_angle, constraints_dihedral,
                geom, sample, mol
            )
        else:
            valid_structure, status = self._run_rdkit_sampling(
                mol, name, charge, mult, constraints_atoms,
                constraints_dist, constraints_angle, constraints_dihedral,
                complex_ts, coord_Map, alg_Map, mol_template, smi,
                geom, original_atn, metal_atoms, metal_idx,
                metal_sym, csearch_nprocs, sample
            )
            
        # Handle errors and combine results if needed
        self._handle_sampling_result(status, valid_structure, name)
        if self.args.crest_runs > 1:
            self._combine_crest_runs(name)
            
    def _should_use_crest(self):
        """Determine if CREST should be used for conformer sampling directly with no RDKit sampling first.
        
        Returns:
            bool: True if CREST should be used
        """
        is_crest = self.args.program.lower() == "crest"
        is_3d_input = (
            self.args.smi is None and
            os.path.basename(Path(self.args.input)).split(".")[-1] in 
            ["pdb", "mol2", "mol", "sdf", "gjf", "com", "xyz"]
        )
        return is_crest and is_3d_input
        
    def _run_crest_sampling(self, name, charge, mult, smi,
                          constraints_atoms, constraints_dist,
                          constraints_angle, constraints_dihedral,
                          geom, sample, mol):
        """Run conformer sampling using CREST.
        
        Args:
            name (str): Molecule name
            charge (int): Molecular charge
            mult (int): Multiplicity
            smi: SMILES or 3D structure
            constraints_*: Various constraint parameters
            geom (list): Geometry constraints
            sample (int): Number of conformers
            mol: RDKit molecule object
            
        Returns:
            int: Status code from xtb_opt_main
        """
        if self.args.crest_runs == 1:
            return xtb_opt_main(
                f"{name}_{self.args.program.lower()}", self,
                charge, mult, smi, constraints_atoms,
                constraints_dist, constraints_angle,
                constraints_dihedral, 'crest', geom,
                sample, mol=mol
            )
        
        # Multiple CREST runs
        status = None
        for pt in range(1, self.args.crest_runs + 1):
            src = f"{name}_{self.args.program.lower()}.xyz"
            dst = f"{name}_run_{pt}_{self.args.program.lower()}.xyz"
            shutil.copy(src, dst)
            
            status = xtb_opt_main(
                f"{name}_run_{pt}_{self.args.program.lower()}",
                self, charge, mult, smi, constraints_atoms,
                constraints_dist, constraints_angle,
                constraints_dihedral, 'crest', geom,
                sample, mol=mol
            )
        return status
        
    def _run_rdkit_sampling(self, mol, name, charge, mult,
                           *args):
        """Run conformer sampling using RDKit.
        
        Args:
            mol: RDKit molecule object
            name (str): Molecule name
            charge (int): Molecular charge
            mult (int): Multiplicity
            *args: Additional parameters for rdkit_search
            
        Returns:
            tuple: (valid_structure, status)
        """
        csearch_file = self.csearch_folder.joinpath(
            f"{name}_{self.args.program.lower()}{self.args.output}"
        )
        
        valid_structure = filters(mol, self.args.log, self.args.max_mol_wt)
        status = None
        
        if valid_structure:
            try:
                status = self.rdkit_search(
                    mol, name, csearch_file, charge, mult, *args
                )
            except (KeyboardInterrupt, SystemExit):
                raise
                
        return valid_structure, status
        
    def _handle_sampling_result(self, status, valid_structure, name):
        """Handle result of conformer sampling.
        
        Args:
            status (int): Status code from sampling
            valid_structure (bool): Whether structure is valid
            name (str): Molecule name
        """
        if status == -1 or not valid_structure:
            self.args.log.write(
                f"\nx  ERROR: The structure is not valid or no conformers "
                f"were obtained from this SMILES string "
                f"({os.path.basename(Path(name))})"
            )
            
    def _combine_crest_runs(self, name):
        """Combine results from multiple CREST runs.
        
        Args:
            name (str): Molecule name
        """
        csearch_file = self.csearch_folder.joinpath(
            f"{name}_{self.args.program.lower()}{self.args.output}"
        )
        
        with Chem.SDWriter(str(csearch_file)) as sdwriter_rd:
            # Collect all molecules and energies
            molecules = []
            for file in glob.glob(f"{self.csearch_folder}/{name}_run_*{self.args.program.lower()}.sdf"):
                mols = load_sdf(file)
                for mol in mols:
                    molecules.append((float(mol.GetProp('Energy')), mol))
            
            # Sort by energy and write
            for _, mol in sorted(molecules, key=lambda x: x[0]):
                sdwriter_rd.write(mol)

    def rdkit_search(
        self, mol, name, csearch_file, charge, mult,
        constraints_atoms, constraints_dist, constraints_angle,
        constraints_dihedral, complex_ts, coord_Map, alg_Map,
        mol_template, smi, geom, original_atn, metal_atoms,
        metal_idx, metal_sym, csearch_nprocs, sample
    ):
        """Generate and optimize conformers using RDKit and optionally CREST.
        
        This method handles:
        1. Initial RDKit conformer generation
        2. Optional CREST optimization
        3. Conformer filtering and clustering
        
        Args:
            mol: RDKit molecule object
            name (str): Molecule name
            csearch_file (Path): Output file path
            charge (int): Molecular charge
            mult (int): Multiplicity
            constraints_*: Various constraint parameters
            complex_ts (bool): Whether molecule is complex/TS
            coord_Map: Coordinate mapping
            alg_Map: Alignment mapping
            mol_template: Template molecule
            smi (str): SMILES string
            geom (list): Geometry constraints
            original_atn: Original atomic numbers
            metal_*: Metal atom information
            csearch_nprocs (int): Number of processors
            sample (int): Number of conformers
            
        Returns:
            int: Status code (0 for success, -1 for failure)
        """
        # Initial RDKit conformer generation
        status, mol_crest = self._generate_initial_conformers(
            mol, name, csearch_file, charge, mult,
            coord_Map, alg_Map, mol_template, smi, geom,
            original_atn, metal_atoms, metal_idx, metal_sym,
            csearch_nprocs, sample, complex_ts
        )

        # CREST optimization if selected
        if self.args.program.lower() == 'crest':
            if mol_crest is not None:
                mol = mol_crest
            status = self._run_crest_optimization(
                mol, name, csearch_file, charge, mult,
                constraints_atoms, constraints_dist,
                constraints_angle, constraints_dihedral,
                complex_ts, geom, sample, status
            )
            
        return status
        
    def _generate_initial_conformers(
        self, mol, name, csearch_file, charge, mult,
        coord_Map, alg_Map, mol_template, smi, geom,
        original_atn, metal_atoms, metal_idx, metal_sym,
        csearch_nprocs, sample, complex_ts
    ):
        """Generate initial conformers using RDKit.
        
        Args:
            (see rdkit_search for parameter documentation)
            
        Returns:
            Status, mol_crest
        """
        if complex_ts:
            return None, None
            
        # Log conformer generation start
        if self.args.program.lower() == 'rdkit':
            self.args.log.write(
                f"\no  Starting RDKit conformer sampling "
                f"({os.path.basename(Path(name))})"
            )
        elif self._needs_rdkit_init():
            self.args.log.write(
                "\no  Starting initial RDKit-based mol generation from SMILES"
            )
            
        status, mol_crest = self.rdkit_to_sdf(
            mol, name, csearch_file, charge, mult,
            coord_Map, alg_Map, mol_template, smi, geom,
            original_atn, metal_atoms, metal_idx, metal_sym,
            csearch_nprocs, sample
        )
    
        return status, mol_crest
        
    def _needs_rdkit_init(self):
        """Check if RDKit initialization is needed for CREST.
        
        Returns:
            bool: True if RDKit initialization needed
        """
        is_crest = self.args.program.lower() == 'crest'
        ext = os.path.basename(Path(self.args.input)).split(".")[-1]
        return is_crest and ext not in ["pdb", "mol2", "mol", "sdf", "gjf", "com", "xyz"]
        
    def _run_crest_optimization(
        self, mol, name, csearch_file, charge, mult,
        constraints_atoms, constraints_dist, constraints_angle,
        constraints_dihedral, complex_ts, geom, sample, initial_status
    ):
        """Run CREST optimization on conformers.
        
        Args:
            (see rdkit_search for parameter documentation)
            initial_status: Status from RDKit generation
            
        Returns:
            int: Final status code
        """
        stop_xtb_opt = False
        status = initial_status
        
        if not complex_ts:
            stop_xtb_opt, status = self._prepare_crest_input(
                mol, name, csearch_file, sample
            )
        else:
            stop_xtb_opt, status = self._handle_complex_crest(
                mol, name
            )
            
        if not stop_xtb_opt:
            status = self._run_crest_calculations(
                name, charge, mult, mol,
                constraints_atoms, constraints_dist,
                constraints_angle, constraints_dihedral,
                geom, sample, complex_ts
            )
            
        return status
        
    def _prepare_crest_input(self, mol, name, csearch_file, sample):
        """Prepare input for CREST from RDKit conformers.
        
        Args:
            mol: RDKit molecule object
            name (str): Molecule name
            csearch_file (Path): Output file path
            sample (int): Number of conformers
            
        Returns:
            tuple: (stop_flag, status)
        """
        if mol is None:
            return True, -1
            
        if self.args.crest_runs == 1:
            os.remove(str(csearch_file))
            rdmolfiles.MolToXYZFile(mol, f"{name}_crest.xyz")
        else:
            suppl, *_ = mol_from_sdf_or_mol_or_mol2(str(csearch_file), "csearch", self.args)
            os.remove(str(csearch_file))
            cluster_mols = cluster_conformers(self, suppl, "crest", csearch_file, name, sample)
            
            for i, conf_mol in enumerate(cluster_mols):
                rdmolfiles.MolToXYZFile(conf_mol, f"{name}_run_{i}_crest.xyz")
                
        return False, 0
        
    def _handle_complex_crest(self, mol, name):
        """Handle CREST input for complex/TS molecules.
        
        Args:
            mol: RDKit molecule object
            name (str): Molecule name
            
        Returns:
            tuple: (stop_flag, status)
        """
        if mol is None:
            return True, -1
            
        if self.args.crest_runs == 1:
            rdmolfiles.MolToXYZFile(mol, f"{name}_crest.xyz")
        else:
            for pt in range(1, self.args.crest_runs + 1):
                rdmolfiles.MolToXYZFile(
                    mol, f"{name}_run_{pt}_crest.xyz"
                )
                
        return False, 0
        
    def _run_crest_calculations(
        self, name, charge, mult, mol,
        constraints_atoms, constraints_dist,
        constraints_angle, constraints_dihedral,
        geom, sample, complex_ts
    ):
        """Run CREST calculations.
        
        Args:
            (see rdkit_search for parameter documentation)
            
        Returns:
            int: Status code
        """
        if self.args.crest_runs == 1:
            return self._single_crest_run(
                name, charge, mult, mol,
                constraints_atoms, constraints_dist,
                constraints_angle, constraints_dihedral,
                geom, sample, complex_ts
            )
            
        status = None
        for pt in range(1, self.args.crest_runs + 1):
            status = self._single_crest_run(
                f"{name}_run_{pt}",
                charge, mult, mol,
                constraints_atoms, constraints_dist,
                constraints_angle, constraints_dihedral,
                geom, sample, complex_ts
            )
        return status
        
    def _single_crest_run(
        self, name_base, charge, mult, mol,
        constraints_atoms, constraints_dist,
        constraints_angle, constraints_dihedral,
        geom, sample, complex_ts
    ):
        """Run a single CREST calculation.
        
        Args:
            name_base (str): Base name for output
            (remaining args same as rdkit_search)
            
        Returns:
            int: Status code
        """
        return xtb_opt_main(
            f"{name_base}_{self.args.program.lower()}",
            self, charge, mult, None,
            constraints_atoms, constraints_dist,
            constraints_angle, constraints_dihedral,
            'crest', geom, sample,
            complex_ts=complex_ts,
            mol=mol
        )

    # Sampling configuration for different levels
    SAMPLING_CONFIGS = {
        'low': {'factor': 5, 'max_confs': 100},    # For descriptor generation/ML
        'mid': {'factor': 10, 'max_confs': 250},   # Standard compromise
        'high': {'factor': 20, 'max_confs': 500}   # More exhaustive search
    }
    
    def auto_sampling(self, mol, metal_atoms, metal_idx):
        """Automatically determine number of conformers for sampling.
        
        This method calculates the appropriate number of conformers based on:
        1. Molecule complexity (rotatable bonds, rings, etc.)
        2. Presence of metal atoms and their coordination
        3. User-specified sampling level (low/mid/high)
        
        Args:
            mol: RDKit molecule object
            metal_atoms (list): Metal atoms present
            metal_idx (list): Indices of metal atoms
            
        Returns:
            int: Number of conformers to generate
            
        Raises:
            SystemExit: If invalid sampling level specified
        """
        # Validate and get sampling configuration
        config = self._get_sampling_config()
        
        # Calculate base number of conformers
        auto_samples = self._calculate_base_samples(mol, metal_atoms, metal_idx)
        
        # Apply sampling factor and enforce maximum
        if auto_samples == 0:
            auto_samples = config['factor']
        else:
            auto_samples *= config['factor']
            
        return min(auto_samples, config['max_confs'])
        
    def _get_sampling_config(self):
        """Get configuration for specified sampling level.
        
        Returns:
            dict: Sampling configuration parameters
            
        Raises:
            SystemExit: If invalid sampling level
        """
        level = self.args.auto_sample.lower()
        if level not in self.SAMPLING_CONFIGS:
            self._error_exit(
                f"{self.args.auto_sample} is not a valid option for "
                "--auto_sample! Please use 'low', 'mid' or 'high'"
            )
        return self.SAMPLING_CONFIGS[level]
        
    def _calculate_base_samples(self, mol, metal_atoms, metal_idx):
        """Calculate base number of conformers needed.
        
        Accounts for:
        - Metal complex isomers
        - Rotatable bonds
        - OH/NH groups
        - Saturated rings
        
        Args:
            mol: RDKit molecule object
            metal_atoms (list): Metal atoms present
            metal_idx (list): Indices of metal atoms
            
        Returns:
            int: Base number of conformers
        """
        samples = 0
        
        # Metal complex trans/cis isomers
        if metal_atoms and metal_idx:
            samples += 3 * len(metal_idx)
            
        # Add samples for various molecular features
        samples += 3 * Lipinski.NumRotatableBonds(mol)  # C3 rotations
        samples += 3 * Lipinski.NHOHCount(mol)         # OH/NH rotations
        samples += 3 * Lipinski.NumSaturatedRings(mol) # Ring conformations
        
        return samples

    def genConformer_r(
        self, mol, conf, sdwriter,
        update_to_rdkit, coord_Map, alg_Map, mol_template,
        original_atn, metal_atoms, metal_idx,
        metal_sym, ff
    ):
        """Process and write conformer to SDF file.
        
        This method handles:
        1. Metal atom restoration from iodine placeholders
        2. Energy minimization for metal-containing molecules
        3. Conformer writing to SDF file
        
        Args:
            mol: RDKit molecule object
            conf (int): Conformer ID
            sdwriter: RDKit SDWriter object
            update_to_rdkit (bool): Whether to update coordinates
            coord_Map: Coordinate mapping
            alg_Map: Alignment mapping
            mol_template: Template molecule
            original_atn: Original atomic numbers
            geom (list): Geometry constraints
            metal_atoms (list): Metal atoms present
            metal_idx (list): Metal atom indices
            metal_sym (list): Metal atom symbols
            ff (str): Force field to use
            
        Returns:
            int: Status code (1 for success)
        """
        if metal_atoms:
            if (self.args.program.lower() in ["rdkit", "crest"] or 
                update_to_rdkit):
                # Minimize and update energy
                mol, energy = self._minimize_metal_conformer(
                    mol, conf, coord_Map, alg_Map,
                    mol_template, ff
                )
                mol.SetProp("Energy", str(energy))
                
                # Restore metal atoms
                mol = set_metal_atomic_number(mol, metal_idx, metal_sym)
                # Restore special atoms (e.g. for Ir_squareplanar)
                if original_atn is not None:
                    mol.GetAtomWithIdx(original_atn[1]).SetAtomicNum(original_atn[0])
                    
        sdwriter.write(mol, -1)
        return 1
        
    def _minimize_metal_conformer(
        self, mol, conf, coord_Map, alg_Map,
        mol_template, ff
    ):
        """Minimize metal-containing conformer.
        
        Args:
            mol: RDKit molecule object
            conf (int): Conformer ID
            coord_Map: Coordinate mapping
            alg_Map: Alignment mapping
            mol_template: Template molecule
            ff (str): Force field to use
            
        Returns:
            tuple: (molecule, energy)
        """
        if coord_Map is None and alg_Map is None and mol_template is None:
            energy = minimize_rdkit_energy(
                mol, conf, self.args.log, ff,
                self.args.opt_steps_rdkit
            )
            return mol, energy
        else:
            return realign_mol(
                mol, conf, coord_Map, alg_Map,
                mol_template, self.args.opt_steps_rdkit
            )


    def embed_conf(self, mol, initial_confs, coord_Map, alg_Map, mol_template, csearch_nprocs, name):
        """Embed multiple conformers using RDKit's distance geometry.
        
        This method:
        1. Handles 3D input structures appropriately
        2. Sets up conformer embedding parameters
        3. Attempts fallback options if initial embedding fails
        
        Args:
            mol: RDKit molecule object
            initial_confs (int): Number of conformers to generate
            coord_Map: Coordinate mapping for template alignment 
            alg_Map: Atom mapping for template alignment
            mol_template: Template molecule for alignment
            csearch_nprocs (int): Number of processors to use
            name (str): Molecule name for logging
            
        Returns:
            list: Generated conformer IDs
        """
        # Handle special input formats
        is_3d_input = os.path.basename(Path(self.args.input)).split('.')[-1].lower() in {
            "sdf", "mol", "mol2"
        }
        
        if is_3d_input:
            Chem.AssignStereochemistryFrom3D(mol)
            
        # Setup embedding parameters
        embed_kwargs = {
            "ignoreSmoothingFailures": True,
            "randomSeed": self.args.seed,
            "numThreads": csearch_nprocs
        }
        
        # Add coordinate mapping if template provided
        if any(x is not None for x in (coord_Map, alg_Map, mol_template)):
            embed_kwargs["coordMap"] = coord_Map
            
        # Try standard embedding first
        cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)
        
        # Fall back to random coordinate generation if needed
        if len(cids) <= 1 and initial_confs != 1:
            self._log_embed_fallback(initial_confs, name)
            cids = self._try_random_embedding(mol, initial_confs, embed_kwargs)
            
        # Preserve stereochemistry for 3D input
        if is_3d_input:
            for cid in cids:
                Chem.AssignAtomChiralTagsFromStructure(mol, confId=cid)
                
        return cids
        
    def _log_embed_fallback(self, initial_confs, name):
        """Log message about falling back to random coordinates.
        
        Args:
            initial_confs (int): Number of conformers attempted
            name (str): Molecule name
        """
        self.args.log.write(
            f"\nx  Normal RDKit embedding process failed, trying to generate "
            f"conformers with random coordinates (with {initial_confs} "
            f"possibilities) ({os.path.basename(Path(name))})"
        )
        
    def _try_random_embedding(self, mol, initial_confs, base_kwargs):
        """Attempt conformer generation with random coordinates.
        
        Args:
            mol: RDKit molecule object
            initial_confs (int): Number of conformers to generate
            base_kwargs (dict): Base embedding parameters
            
        Returns:
            list: Generated conformer IDs
        """
        embed_kwargs = base_kwargs.copy()
        embed_kwargs.update({
            "useRandomCoords": True,
            "boxSizeMult": 10.0,
            "numZeroFail": 1000
        })
        return rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

    def min_and_E_calc(self, mol, cids, coord_Map, alg_Map, mol_template,
                       ff, geom, metal_atoms, metal_idx, metal_sym):
        """Energy minimization and geometry filtering of conformers.
        
        This method:
        1. Minimizes each conformer using RDKit force fields
        2. Applies geometry filters and constraints
        3. Collects passing conformers and energies
        
        Args:
            mol: RDKit molecule object
            cids (list): Conformer IDs to process
            coord_Map: Coordinate mapping
            alg_Map: Alignment mapping 
            mol_template: Template molecule
            ff (str): Force field to use
            geom (list): Geometry constraints
            metal_atoms (list): Metal atoms
            metal_idx (list): Metal atom indices
            metal_sym (list): Metal atom symbols
            
        Returns:
            tuple: (
                outmols: List of passing molecule objects,
                passing_cids: List of passing conformer IDs,
                cenergy: List of conformer energies
            )
        """
        cenergy, outmols, passing_cids = [], [], []
        
        for _, conf in enumerate(cids):
            # Minimize conformer
            mol, energy = self._minimize_conformer(
                mol, conf, coord_Map, alg_Map,
                mol_template, ff
            )
            
            # Process metal atoms
            mol_ensemb = self._process_metal_atoms(
                mol, metal_atoms, metal_idx, metal_sym
            )
            
            mol_geom = mol.GetConformer(conf)
            
            # Check geometry constraints
            if self._passes_geometry_filters(mol_ensemb, mol_geom, geom):
                self._add_passing_conformer(
                    mol, mol_geom, energy,
                    cenergy, outmols, passing_cids, conf
                )
                
        return outmols, passing_cids, cenergy
        
    def _minimize_conformer(self, mol, conf, coord_Map, alg_Map,
                          mol_template, ff):
        """Minimize a single conformer.
        
        Args:
            (see min_and_E_calc for parameter details)
            
        Returns:
            tuple: (minimized molecule, energy)
        """
        if coord_Map is None and alg_Map is None and mol_template is None:
            energy = minimize_rdkit_energy(
                mol, conf, self.args.log, ff,
                self.args.opt_steps_rdkit
            )
            return mol, energy
        else:
            return realign_mol(
                mol, conf, coord_Map, alg_Map,
                mol_template, self.args.opt_steps_rdkit
            )
            
    def _process_metal_atoms(self, mol, metal_atoms, metal_idx, metal_sym):
        """Process metal atoms in molecule.
        
        Args:
            (see min_and_E_calc for parameter details)
            
        Returns:
            RDKit.Mol: Processed molecule
        """
        mol_ensemb = Chem.Mol(mol)
        if metal_atoms:
            mol_ensemb = set_metal_atomic_number(mol_ensemb, metal_idx, metal_sym)
        return mol_ensemb
        
    def _passes_geometry_filters(self, mol, mol_geom, geom):
        """Check if conformer passes geometry constraints.
        
        Args:
            mol: RDKit molecule object
            mol_geom: Conformer geometry
            geom (list): Geometry constraints
            
        Returns:
            bool: True if passes all constraints
        """
        if not geom:
            return True
            
        # Handle Ir square planar special case
        if geom == ['Ir_squareplanar']:
            return geom_filter(self, mol, mol_geom, geom)
            
        # Handle regular geometry constraints
        passing_geom = False
        for i, ele in enumerate(geom):
            if i % 2 == 0:
                if i == 0:
                    passing_geom = False
                if passing_geom or i == 0:
                    rule = [ele, geom[i+1]]
                    passing_geom = geom_filter(
                        self, mol, mol_geom, rule
                    )
                else:
                    passing_geom = False
                    
        return passing_geom
        
    def _add_passing_conformer(self, mol, mol_geom, energy,
                             cenergy, outmols, passing_cids, conf):
        """Add a passing conformer to results.
        
        Args:
            mol: Source molecule
            mol_geom: Conformer geometry
            energy (float): Conformer energy
            cenergy (list): List of energies to append to
            outmols (list): List of molecules to append to
            passing_cids (list): List of conformer IDs to append to
            conf: Current conformer ID
        """
        cenergy.append(energy)
        
        # Create single conformer molecule
        mol_single_conf = Chem.Mol(mol)
        mol_single_conf.RemoveAllConformers()
        mol_single_conf.AddConformer(mol_geom, assignId=True)
        
        # Convert to PropertyMol
        pmol = PropertyMol.PropertyMol(mol_single_conf)
        outmols.append(pmol)
        passing_cids.append(conf)

    def min_after_embed(
        self, mol, cids, name, csearch_file,
        update_to_rdkit, coord_Map, alg_Map, mol_template,
        charge, mult, ff, smi, geom, original_atn,
        metal_atoms, metal_idx, metal_sym, sample
    ):
        """Process embedded conformers including minimization and filtering.
        
        This method:
        1. Minimizes and filters conformers
        2. Applies geometry constraints
        3. Sorts by energy
        4. Writes selected conformers to file
        5. Optionally clusters similar conformers
        
        Args:
            mol: RDKit molecule object
            cids (list): Conformer IDs
            name (str): Molecule name
            csearch_file (Path): Output file path
            update_to_rdkit (bool): Update coordinates flag
            coord_Map: Coordinate mapping
            alg_Map: Alignment mapping
            mol_template: Template molecule
            charge (int): Molecular charge
            mult (int): Multiplicity
            ff (str): Force field
            smi (str): SMILES string
            geom (list): Geometry constraints
            original_atn: Original atomic numbers
            metal_atoms (list): Metal atoms
            metal_idx (list): Metal indices 
            metal_sym (list): Metal symbols
            sample (int): Number of conformers to keep
            
        Returns:
            tuple: (status, output molecules)
        """
        # Process and filter conformers
        outmols, cenergy = self._process_conformers(
            mol, cids, name, charge, mult, smi, geom,
            coord_Map, alg_Map, mol_template, ff,
            metal_atoms, metal_idx, metal_sym
        )
        
        # Sort and select conformers
        selectedcids_rdkit = self._select_conformers(
            outmols, cenergy, name, ff
        )
        
        # Write selected conformers
        if self.args.program.lower() in ["rdkit", "crest"]:
            status = self._write_conformers(
                outmols, selectedcids_rdkit, csearch_file,
                update_to_rdkit, coord_Map,
                alg_Map, mol_template, original_atn,
                metal_atoms, metal_idx, metal_sym, ff
            )
            
            # Cluster if needed
            outmols = self._cluster_if_needed(
                csearch_file, name, sample
            )
        else:
            status = 1
        return status, outmols
        
    def _process_conformers(
        self, mol, cids, name, charge, mult, smi, geom,
        coord_Map, alg_Map, mol_template, ff,
        metal_atoms, metal_idx, metal_sym
    ):
        """Process and filter conformers.
        
        Args:
            (see min_after_embed for parameter details)
            
        Returns:
            tuple: (outmols, passing_cids, cenergy)
        """
        if geom:
            self.args.log.write(
                f"o  Applying geometry filters ({geom}) "
                f"({os.path.basename(Path(name))})"
            )
            
        outmols, passing_cids, cenergy = self.min_and_E_calc(
            mol, cids, coord_Map, alg_Map, mol_template,
            ff, geom, metal_atoms, metal_idx, metal_sym
        )
        
        # Add properties to passing molecules
        for i, _ in enumerate(passing_cids):
            self._add_mol_properties(
                outmols[i], name, i+1, cenergy[i],
                charge, mult, smi
            )
            
        return outmols, cenergy
        
    def _select_conformers(self, outmols, cenergy, name, ff):
        """Select conformers based on energy and filtering.
        
        Args:
            outmols (list): Molecule objects
            cenergy (list): Conformer energies
            name (str): Molecule name
            ff (str): Force field
            
        Returns:
            list: Selected conformer IDs
        """
        cids = list(range(len(outmols)))
        sorted_all_cids = sorted(cids, key=lambda cid: cenergy[cid])
        
        if ff.upper() == "NO FF":
            return sorted_all_cids
            
        self.args.log.write(
            f"\no  Applying filters to initial conformers "
            f"({os.path.basename(Path(name))})"
        )
        return conformer_filters(self, sorted_all_cids, cenergy, outmols)
        
    def _write_conformers(
        self, outmols, selected_cids, csearch_file,
        update_to_rdkit, coord_Map,
        alg_Map, mol_template, original_atn,
        metal_atoms, metal_idx, metal_sym, ff
    ):
        """Write selected conformers to SDF file.
        
        Args:
            (see min_after_embed for parameter details)
            
        Returns:
            int: Status code (1 for success)
        """
        total = 0
        with Chem.SDWriter(str(csearch_file)) as sdwriter:
            for conf in selected_cids:
                total += self.genConformer_r(
                    outmols[conf], -1,
                    sdwriter, update_to_rdkit, coord_Map,
                    alg_Map, mol_template, original_atn,
                    metal_atoms, metal_idx,
                    metal_sym, ff
                )
        return 1
        
    def _cluster_if_needed(self, csearch_file, name, sample):
        """Cluster conformers if needed.
        
        Args:
            csearch_file (Path): Path to SDF file
            name (str): Molecule name
            sample (int): Target number of conformers
            
        Returns:
            list: Final molecule objects
        """
        suppl, *_ = mol_from_sdf_or_mol_or_mol2(
            str(csearch_file), "csearch", self.args
        )
        if len(suppl) > sample and self.args.auto_cluster:
            return cluster_conformers(
                self, suppl, "rdkit", csearch_file,
                name, sample
            )
            
        return suppl
        
    def _add_mol_properties(self, mol, name, idx, energy,
                          charge, mult, smi):
        """Add properties to molecule object.
        
        Args:
            mol: RDKit molecule object
            name (str): Base name
            idx (int): Conformer index
            energy (float): Conformer energy
            charge (int): Molecular charge
            mult (int): Multiplicity
            smi (str): SMILES string
        """
        mol.SetProp("_Name", f"{name} {idx}")
        mol.SetProp("Energy", str(energy))
        mol.SetProp("Real charge", str(charge))
        mol.SetProp("Mult", str(mult))
        mol.SetProp("SMILES", str(smi))

    def rdkit_to_sdf(
        self,
        mol,
        name,
        csearch_file,
        charge,
        mult,
        coord_Map,
        alg_Map,
        mol_template,
        smi,
        geom,
        original_atn,
        metal_atoms,
        metal_idx,
        metal_sym,
        csearch_nprocs,
        sample
    ):

        """
        Conversion from RDKit to SDF
        """

        mol.SetProp("_Name", name)

        # detects and applies auto-detection of initial number of conformers
        if self.args.auto_sample:
            initial_confs = int(self.auto_sampling(mol,metal_atoms,metal_idx))
        else:
            initial_confs = sample

        update_to_rdkit = False

        rotmatches = getDihedralMatches(mol, self.args.heavyonly)

        if len(rotmatches) > self.args.max_torsions and self.args.max_torsions > 0:
            self.args.log.write(f"\nx  Too many torsions ({len(rotmatches)}). Skipping {name + self.args.output}")

        ff = self.args.ff
        cids = self.embed_conf(mol, initial_confs, coord_Map, alg_Map, mol_template, csearch_nprocs, name)

        # energy minimize all to get more realistic results
        # identify the atoms and decide Force Field
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 36 and self.args.ff == "MMFF":  # up to Kr for MMFF, if not the code will use UFF
                self.args.log.write(f"\nx  {self.args.ff} is not compatible with the molecule, changing to UFF (({os.path.basename(Path(name))}))")
                ff = "UFF"
                break

        try:
            status, mol_crest = self.min_after_embed(
                mol,
                cids,
                name,
                csearch_file,
                update_to_rdkit,
                coord_Map,
                alg_Map,
                mol_template,
                charge,
                mult,
                ff,
                smi,
                geom,
                original_atn,
                metal_atoms,
                metal_idx,
                metal_sym,
                sample
            )
            return status, mol_crest[0]
        except IndexError:
            status = -1
            mol_crest = None
            return status, mol_crest
