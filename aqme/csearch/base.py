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
      Current options: 'rdkit', 'summ', 'fullmonte', 'crest'  
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
      options: MMFF and UFF (if MMFF fails, AQME tries to use UFF automatically)
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
      molecules with C-Pd-C atoms at 180 degrees: ['[C][Pd][C]',180].
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

SUMM only
+++++++++

   degree : float, default=120.0
      Interval of degrees to rotate dihedral angles during SUMM sampling 
      (i.e. 120.0 would create 3 conformers for each dihedral, at 0, 
      120 and 240 degrees)

Fullmonte only
++++++++++++++

   ewin_fullmonte : float, default=5.0
      Energy window in kcal/mol to discard conformers (i.e. if a conformer is 
      more than the E window compared to the most stable conformer)
   ewin_sample_fullmonte : float, default=2.0
      Energy window in kcal/mol to use conformers during the Fullmonte sampling 
      (i.e. conformers inside the E window compared to the most stable conformer 
      are considered as unique in each step of the sampling)
   nsteps_fullmonte : int, default=100
      Number of steps (or conformer batches) to carry during the Fullmonte 
      sampling
   nrot_fullmonte : int, default=3
      Number of dihedrals to rotate simultaneously (picked at random) during 
      each step of the Fullmonte sampling
   ang_fullmonte : float, default=30
      Available angle interval to use in the Fullmonte sampling. For example, if
      the angle is 120.0, the program chooses randomly between 120 and 240 
      degrees (picked at random) during each step of the sampling

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
#####################################################.
#          This file stores the CSEARCH class       #
#             used in conformer generation          #
#####################################################.

import math
import os
import sys
import time
import shutil
import subprocess
import glob
import concurrent.futures as futures
from pathlib import Path
from progress.bar import IncrementalBar

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors as Descriptors
from rdkit.Chem import rdmolfiles, rdMolTransforms, PropertyMol, rdDistGeom, Lipinski

from aqme.filter import (
    filters,
    conformer_filters,
    geom_filter,
    cluster_conformers
    )
from aqme.csearch.utils import (
    prepare_direct_smi,
    prepare_smiles_files,
    prepare_csv_files,
    prepare_cdx_files,
    prepare_com_files,
    prepare_sdf_files,
    prepare_pdb_files,
    minimize_rdkit_energy,
    com_2_xyz,
    check_constraints,
    smi_to_mol,
    getDihedralMatches,
    substituted_mol
    )
from aqme.csearch.templates import template_embed, check_metal_neigh
from aqme.csearch.fullmonte import generating_conformations_fullmonte, realign_mol
from aqme.utils import (
    load_variables,
    set_metal_atomic_number,
    check_xtb,
    check_crest,
    get_files,
    check_dependencies,
    mol_from_sdf_or_mol_or_mol2,
    set_destination,
    load_sdf
    )
from aqme.csearch.crest import xtb_opt_main


class csearch:
    """
    Class absracting the geometry generation and conformational search procedure.
    For further detail on the currently accepted keyword arguments (kwargs) 
    please look at the Parameters section (in the module documentation). 
    """

    def __init__(self, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "csearch")
        
        # check whether dependencies are installed
        _ = check_dependencies(self)

        csearch_program = True
        if self.args.program is None:
            csearch_program = False
        if csearch_program:
            if self.args.program.lower() not in ["rdkit", "summ", "fullmonte", "crest"]:
                csearch_program = False
        if not csearch_program:
            self.args.log.write('\nx  Program not specified or not supported for CSEARCH! Specify: program="rdkit" (or "crest", "summ", "fullmonte")')
            self.args.log.finalize()
            sys.exit()

        if self.args.ff.upper() not in ["MMFF", "UFF"]:
            self.args.log.write(f"x  Force field {self.args.ff} not supported!")
            self.args.log.finalize()
            sys.exit()

        # set number of processors (threads for CSEARCH)
        if self.args.nprocs is None:
            self.args.nprocs = 4

        # default value of auto_sample
        if self.args.auto_sample == 'auto':
            self.args.auto_sample = 'mid'

        if self.args.program.lower() == "crest":
            _ = check_xtb(self)
            _ = check_crest(self)

        if self.args.smi is None and self.args.input == "":
            self.args.log.write("\nx  Program requires either a SMILES or an input file to proceed! Please look up acceptable file formats. Specify: smi='CCC' (or input='filename.csv')")
            self.args.log.finalize()
            sys.exit()
        elif self.args.smi is not None and self.args.input != "":
            self.args.log.write("\nx  Program requires either a SMILES or an input file to proceed, don't use both!")
            self.args.log.finalize()
            sys.exit()
        # specify a dummy extension for inputif smi is not None (avoids errors in LOG printings)
        if self.args.smi is not None:
            self.args.input = 'no_ext.no_ext'

        try:
            if Path(f"{self.args.w_dir_main}").exists():
                os.chdir(self.args.w_dir_main)
        except FileNotFoundError:
            self.args.w_dir_main = Path(f"{os.getcwd()}/{self.args.w_dir_main}")
            os.chdir(self.args.w_dir_main)

        # load files from AQME input
        if self.args.smi is not None:
            csearch_files = [self.args.name]
        else:
            csearch_files = get_files(self.args.input)
            if len(csearch_files) == 0:
                self.args.log.write(f"\nx  Input file ({self.args.input}) not found!")
                self.args.log.finalize()
                sys.exit()

        for csearch_file in csearch_files:
            # load jobs for conformer generation
            if self.args.smi is not None:
                job_inputs = prepare_direct_smi(self.args)

            else:
                job_inputs = self.load_jobs(csearch_file)

            self.args.log.write(f"\nStarting CSEARCH with {len(job_inputs)} job(s) (SDF, XYZ, CSV, etc. files might contain multiple jobs/structures inside)\n")

            # runs the conformer sampling with multiprocessors
            _ = self.run_csearch(job_inputs)

            # removes systems that did not generate any conformers
            for sdf_file in glob.glob(f'{self.csearch_folder}/*.sdf'):
                if os.path.getsize(sdf_file) == 0:
                    os.remove(sdf_file)

        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime CSEARCH: {elapsed_time} seconds\n")
        self.args.log.finalize()

        # this is added to avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

    def load_jobs(self, csearch_file):
        """
        Load information of the different molecules for conformer generation
        """

        SUPPORTED_INPUTS = [
            "smi",
            "sdf",
            "cdx",
            "csv",
            "com",
            "gjf",
            "mol",
            "mol2",
            "xyz",
            "txt",
            "yaml",
            "yml",
            "rtf",
            "pdb",
        ]

        file_format = os.path.basename(Path(csearch_file)).split('.')[-1]
        # Checks
        if file_format.lower() not in SUPPORTED_INPUTS:
            self.args.log.write("\nx  Input filetype not currently supported!")
            self.args.log.finalize()
            sys.exit()

        smi_derivatives = ["smi", "txt", "yaml", "yml", "rtf"]
        Extension2inputgen = dict()
        for key in smi_derivatives:
            Extension2inputgen[key] = prepare_smiles_files
        Extension2inputgen["csv"] = prepare_csv_files
        Extension2inputgen["cdx"] = prepare_cdx_files
        Extension2inputgen["gjf"] = prepare_com_files
        Extension2inputgen["com"] = prepare_com_files
        Extension2inputgen["xyz"] = prepare_com_files
        Extension2inputgen["sdf"] = prepare_sdf_files
        Extension2inputgen["mol"] = prepare_sdf_files
        Extension2inputgen["mol2"] = prepare_sdf_files
        Extension2inputgen["pdb"] = prepare_pdb_files

        # Prepare the jobs
        prepare_function = Extension2inputgen[file_format]
        try:
            job_inputs = prepare_function(self.args, csearch_file)
        except FileNotFoundError:
            self.args.log.write(f'\nx  File {os.path.basename(csearch_file)} was not found! In the "input" option, make sure that 1) the PATH to the files is correct and 2) the PATH doesn\'t start with "/".')
            self.args.log.finalize()
            sys.exit()     

        return job_inputs

    def run_csearch(self, job_inputs):

        bar = IncrementalBar(
            "o  Number of finished jobs from CSEARCH", max=len(job_inputs)
        )

        # multiprocessing to accelerate and make CSEARCH reproducible (since RDKit uses 1 thread to be reproducible)
        if not self.args.debug and self.args.program.lower() != 'crest': # errors and try/excepts are not shown in multithreading
            with futures.ThreadPoolExecutor(
                max_workers=self.args.nprocs,
            ) as executor:
                csearch_nprocs = 1
                for job_input in job_inputs:
                    _ = executor.submit(
                        self.compute_confs, job_input,bar,csearch_nprocs
                        )

        else:
            for job_input in job_inputs:
                _ = self.compute_confs(job_input,bar,self.args.nprocs)

        bar.finish()


    def compute_confs(self,job_input,bar,csearch_nprocs):
        """
        Function to start conformer generation
        """

        # load variables from job_input
        (
            smi,
            name,
            charge,
            mult,
            constraints_atoms,
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
            complex_type,
            geom
        ) = job_input
        
        self.args.log.write(f"\n   ----- {os.path.basename(Path(name))} -----")

        # load mol and other parameters when using SMILES as input
        if self.args.smi is not None or os.path.basename(Path(self.args.input)).split(".")[-1] in ["smi","csv","cdx","txt","yaml","yml","rtf"]:
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
                constraints_dihedral,
            )
            if mol is None:
                self.args.log.write(f"\nx  Failed to convert the provided SMILES ({smi}) to an RDkit Mol object! Please check the starting smiles ({os.path.basename(Path(name))})")
                # if a list of SMILES is provided, the program doesn't stop if one SMILES fails to convert to mol
                if os.path.basename(Path(self.args.input)).split(".")[-1] not in ["csv","cdx","txt","yaml","yml","rtf"]:
                    self.args.log.finalize()
                    sys.exit()
                bar.next()
                return

        else:
            # for 3D input formats, the smi variable represents the mol object
            mol = smi
            if mol is None:
                self.args.log.write(f"\nx  Failed to convert the provided input to an RDkit Mol object! Please check the starting structure ({os.path.basename(Path(name))})")
                if os.path.basename(Path(self.args.input)).split(".")[-1] not in ["csv","cdx","txt","yaml","yml","rtf"]:
                    self.args.log.finalize()
                    sys.exit()
                bar.next()
                return
                
            # check if the optimization is constrained
            complex_ts = check_constraints(self)

        self.csearch_folder = set_destination(self,'CSEARCH')

        self.csearch_folder.mkdir(exist_ok=True, parents=True)

        # for 3D input types
        if self.args.program.lower() in ["crest"] and self.args.smi is None:
            if os.path.basename(Path(self.args.input)).split(".")[-1] in ["pdb", "mol2", "mol", "sdf"]:
                command_pdb = [
                    "obabel",
                    f'-i{os.path.basename(Path(self.args.input)).split(".")[-1]}',
                    f'{name}.{os.path.basename(Path(self.args.input)).split(".")[-1]}',
                    "-oxyz",
                    f"-O{os.path.dirname(Path(name))}/{'.'.join(os.path.basename(Path(name)).split('.')[:-1])}_{self.args.program.lower()}.xyz",
                ]
                subprocess.run(
                    command_pdb,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            elif os.path.basename(Path(self.args.input)).split(".")[-1] in ["gjf", "com"]:
                xyz_file, _, _ = com_2_xyz(f'{name}.{os.path.basename(Path(self.args.input)).split(".")[-1]}')
                shutil.move(xyz_file, f"{name}_{self.args.program.lower()}.xyz")
            elif os.path.basename(Path(self.args.input)).split(".")[-1] == "xyz":
                shutil.copy(f"{name}.xyz", f"{name}_{self.args.program.lower()}.xyz")

        template_opt,valid_template_embed = False,True

        # detects metal atoms
        metal_atoms,metal_idx,complex_coord,metal_sym = [],[],[],[]
        if self.args.auto_metal_atoms:
            metal_atoms = self.find_metal_atom(mol,charge,mult,name)

        # replaces the metal for an I atom
        if len(metal_atoms) >= 1:
            (
                metal_idx,
                complex_coord, # this variable will be used for checking n of ligands when using templates
                metal_sym,
            ) = substituted_mol(mol, "I", metal_atoms)

            # get pre-determined geometries for metal complexes
            accepted_complex_types = [
                "squareplanar",
                "squarepyramidal",
                "linear",
                "trigonalplanar",
            ]
            if complex_type != '' and complex_type not in accepted_complex_types:
                self.args.log.write(f"x  The metal template specified in complex_type ({complex_type}) is not valid! Options: squareplanar, squarepyramidal, linear and trigonalplanar ({os.path.basename(Path(name))})")
                if os.path.basename(Path(self.args.input)).split(".")[-1] not in ["csv","cdx","txt","yaml","yml","rtf"]:
                    self.args.log.finalize()
                    sys.exit()
                bar.next()
                return

            if complex_type in accepted_complex_types:
                count_metals = 0
                valid_template = True
                # check if the specified metal is included in the system
                for metal_idx_ind in metal_idx:
                    if metal_idx_ind is not None:
                        # calculate number of expected neighbours
                        valid_template = check_metal_neigh(mol, complex_type, metal_idx_ind, self.args.log, valid_template)
                        count_metals += 1

                if count_metals == 1 and valid_template:
                    template_opt = True
                    template_kwargs = dict()
                    template_kwargs["complex_type"] = complex_type
                    template_kwargs["metal_idx"] = metal_idx
                    template_kwargs["maxsteps"] = self.args.opt_steps_rdkit
                    template_kwargs["heavyonly"] = self.args.heavyonly
                    template_kwargs["maxmatches"] = self.args.max_matches_rmsd
                    template_kwargs["mol"] = mol
                    template_kwargs["name"] = name
                    template_kwargs["geom"] = geom
                    try:
                        items = template_embed(self, **template_kwargs)
                    except RuntimeError:
                        valid_template_embed = False
                        self.args.log.write(f"\nx  Molecule {name} was not optimized in the specified template ({complex_type}). Try to use charges to make it easier for MM protocols (i.e. using [NH3+][Ag][NH3+] instead of [NH3][Ag][NH3], since N+ typically has 4 bonds).")

                    if valid_template_embed:
                        for mol_obj, name_in, coord_map, alg_map, template, original_atn in zip(*items):
                            _ = self.conformer_generation(
                                mol_obj,
                                name_in,
                                constraints_atoms,
                                constraints_dist,
                                constraints_angle,
                                constraints_dihedral,
                                complex_ts,
                                charge,
                                mult,
                                smi,
                                geom,
                                metal_atoms,
                                metal_idx,
                                metal_sym,
                                csearch_nprocs,
                                coord_map,
                                alg_map,
                                template,
                                original_atn
                            )

                elif count_metals > 1 or count_metals == 0:
                    self.args.log.write(f"\nx  The template specified {complex_type} is not used for systems with more than 1 metal or for organic molecules ({os.path.basename(Path(name))})")

        if not template_opt and valid_template_embed:
            _ = self.conformer_generation(
                mol,
                name,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
                complex_ts,
                charge,
                mult,
                smi,
                geom,
                metal_atoms,
                metal_idx,
                metal_sym,
                csearch_nprocs
            )

        bar.next()

    # automatic detection of metal atoms   
    def find_metal_atom(self,mol,charge,mult,name):
        metal_atoms = [] # for batch jobs such as CSV inputs with many SMILES
        transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo',
                            'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                            'Hg', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in transition_metals:
                metal_atoms.append(atom.GetSymbol())
        if len(metal_atoms) > 0:
            self.args.log.write(f"\no  AQME recognized the following metal atoms: {metal_atoms} ({os.path.basename(Path(name))})")
            if charge is None:
                self.args.log.write(f"\nx  The automated charge calculation might not be precise for metal complexes! You should use the charge option (or the charge column in CSV inputs) ({os.path.basename(Path(name))})")
            if mult is None:
                self.args.log.write(f"\nx  The automated multiplicity calculation might not be precise for metal complexes! You should use the mult option (or the mult column in CSV inputs) ({os.path.basename(Path(name))})")

        return metal_atoms
    
    def conformer_generation(
        self,
        mol,
        name,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
        complex_ts,
        charge,
        mult,
        smi,
        geom,
        metal_atoms,
        metal_idx,
        metal_sym,
        csearch_nprocs,
        coord_Map=None,
        alg_Map=None,
        mol_template=None,
        original_atn=None
    ):
        """
        Function to load mol objects and create 3D conformers

        """

        status = None

        # Set charge and multiplicity
        # user can overwrite charge and mult with the corresponding arguments
        if charge is None:
            charge = Chem.GetFormalCharge(mol)
        if mult is None:
            mult = Descriptors.NumRadicalElectrons(mol) + 1

        # inputs that go through CREST containing 3D coordinates don't require a previous RDKit conformer sampling
        if (
            self.args.program.lower() in ["crest"]
            and self.args.smi is None
            and os.path.basename(Path(self.args.input)).split(".")[-1] in ["pdb","mol2","mol","sdf","gjf","com","xyz"]
        ):

            valid_structure = True
            if self.args.crest_runs == 1:
                status = xtb_opt_main(
                    f"{name}_{self.args.program.lower()}",
                    self,
                    charge,
                    mult,
                    smi,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                    'crest',
                    geom,
                    mol=mol, 
                )
            else:
                for pt in range(1, int(self.args.crest_runs)+1):
                    shutil.copy(f"{name}_{self.args.program.lower()}.xyz", f"{name}_run_{pt}_{self.args.program.lower()}.xyz")
                    status = xtb_opt_main(
                        f"{name}_run_{pt}_{self.args.program.lower()}",
                        self,
                        charge,
                        mult,
                        smi,
                        constraints_atoms,
                        constraints_dist,
                        constraints_angle,
                        constraints_dihedral,
                        'crest',
                        geom,
                        mol=mol, 
                    )

        else:
            csearch_file = self.csearch_folder.joinpath(
                name + "_" + self.args.program.lower() + self.args.output
            )

            valid_structure = filters(
                mol, self.args.log, self.args.max_mol_wt
            )
            if valid_structure:
                try:
                    # the conformational search for RDKit
                    status = self.summ_search(
                        mol,
                        name,
                        csearch_file,
                        charge,
                        mult,
                        constraints_atoms,
                        constraints_dist,
                        constraints_angle,
                        constraints_dihedral,
                        complex_ts,
                        coord_Map,
                        alg_Map,
                        mol_template,
                        smi,
                        geom,
                        original_atn,
                        metal_atoms,
                        metal_idx,
                        metal_sym,
                        csearch_nprocs
                    )
                except (KeyboardInterrupt, SystemExit):
                    raise

        if status == -1 or not valid_structure:
            error_message = f"\nx  ERROR: The structure is not valid or no conformers were obtained from this SMILES string ({os.path.basename(Path(name))})"
            self.args.log.write(error_message)

        #combining all the sdfs from more than one run
        if self.args.crest_runs != 1:
            sdwriter_rd = Chem.SDWriter(f'{csearch_file}')
            file_runs = glob.glob(str(self.csearch_folder)+'/'+ name +'_run_*'+ self.args.program.lower() +'.sdf')
            allenergy, allmols = [], []
            for file in file_runs:
                mols = load_sdf(file)
                for mol in mols:
                    allmols.append(mol)
                    allenergy.append(float(mol.GetProp('Energy')))
            
            allmols_sorted = [mol for _, mol in sorted(zip(allenergy, allmols), key=lambda pair: pair[0])]
            for mol in allmols_sorted:
                sdwriter_rd.write(mol)
            sdwriter_rd.close()

    def summ_search(
        self,
        mol,
        name,
        csearch_file,
        charge,
        mult,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
        complex_ts,
        coord_Map,
        alg_Map,
        mol_template,
        smi,
        geom,
        original_atn,
        metal_atoms,
        metal_idx,
        metal_sym,
        csearch_nprocs
    ):

        """
        Embeds, optimizes and filters RDKit conformers
        """

        # writes sdf for the first RDKit conformer generation
        if not complex_ts:
            if self.args.program.lower() in ['rdkit']:
                self.args.log.write(f"\no  Starting RDKit conformer sampling ({os.path.basename(Path(name))})")
            elif self.args.program.lower() in ['summ','fullmonte']:
                self.args.log.write(f"\no  Starting RDKit-{self.args.program} conformer sampling ({os.path.basename(Path(name))})")
            elif self.args.program.lower() in ['crest'] and os.path.basename(Path(self.args.input)).split(".")[-1] not in ["pdb","mol2","mol","sdf","gjf","com","xyz"]:
                self.args.log.write(f"\no  Starting initial RDKit-based mol generation from SMILES")

            status, rotmatches, ff, mol_crest = self.rdkit_to_sdf(
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
                csearch_nprocs
            )

            # reads the initial SDF files from RDKit and uses dihedral scan if selected
            if status not in [-1, 0]:
                # getting the energy and mols after rotations
                if self.args.program.lower() == "summ" and len(rotmatches) != 0:
                    status = self.dihedral_filter_and_sdf(
                        name, csearch_file, coord_Map, 
                        alg_Map, mol_template, ff, metal_atoms, metal_idx, metal_sym
                    )

        if self.args.program.lower() in ['crest']:
            stop_xtb_opt = False
            if not complex_ts:
                # mol_crest is the RDKit-optimized mol object with all the conformers sorted by E
                if mol_crest is not None:
                    if self.args.crest_runs == 1:
                        os.remove(f'{csearch_file}')
                        rdmolfiles.MolToXYZFile(mol_crest[0], name + "_crest.xyz") # use lowest E conformer
                    else:
                        # clustering to get the most different mol objects
                        suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(f'{csearch_file}', "csearch", self.args)
                        os.remove(f'{csearch_file}')
                        cluster_centroid_mols = cluster_conformers(self,suppl,"crest",csearch_file,name)
                        for i,mol in enumerate(cluster_centroid_mols):
                            rdmolfiles.MolToXYZFile(mol, f'{name}_run_{i}_crest.xyz')                            
                else:
                    stop_xtb_opt = True
                    status = -1
            else:
                # mol is the raw mol object (no optimization with RDKit to avoid problems when using
                # noncovalent complexes and TSs)
                if mol is not None:
                    if self.args.crest_runs == 1:
                        rdmolfiles.MolToXYZFile(mol, name + "_crest.xyz")
                    else:
                        num_start_points = min(int(self.args.crest_runs), len(mol_crest))
                        for pt in range(1, num_start_points+1):
                            rdmolfiles.MolToXYZFile(mol, name + "_run_{0}_crest.xyz".format(pt))
                else:
                    stop_xtb_opt = True
                    status = -1
            if not stop_xtb_opt:
                if self.args.crest_runs == 1:
                    status = xtb_opt_main(
                        f"{name}_{self.args.program.lower()}",
                        self,
                        charge,
                        mult,
                        smi,
                        constraints_atoms,
                        constraints_dist,
                        constraints_angle,
                        constraints_dihedral,
                        'crest',
                        geom,
                        complex_ts=complex_ts,
                        mol=mol, # this is necessary for CREST calculations with constraints 
                        )
                else:
                    num_start_points = min(int(self.args.crest_runs), len(mol_crest))
                    for pt in range(1, num_start_points+1):
                        status = xtb_opt_main(
                            f"{name}_run_{pt}_{self.args.program.lower()}",
                            self,
                            charge,
                            mult,
                            smi,
                            constraints_atoms,
                            constraints_dist,
                            constraints_angle,
                            constraints_dihedral,
                            'crest',
                            geom,
                            complex_ts=complex_ts,
                            mol=mol, # this is necessary for CREST calculations with constraints 
                        )

        return status

    def dihedral_filter_and_sdf(
        self, name, csearch_file, coord_Map, alg_Map, 
        mol_template, ff, metal_atoms, metal_idx, metal_sym
    ):
        """
        Filtering after dihedral scan to sdf
        """

        rotated_energy = []
        # apply filters
        rdmols = load_sdf(csearch_file)

        for i, rd_mol_i in enumerate(rdmols):
            if coord_Map is None and alg_Map is None and mol_template is None:
                energy = minimize_rdkit_energy(
                    rd_mol_i, -1, self.args.log, ff, self.args.opt_steps_rdkit
                )
            else:
                rd_mol_i, energy = realign_mol(
                    rd_mol_i,
                    -1,
                    coord_Map,
                    alg_Map,
                    mol_template,
                    self.args.opt_steps_rdkit,
                )
            rotated_energy.append(energy)

        rotated_cids = list(range(len(rdmols)))
        sorted_rotated_cids = sorted(rotated_cids, key=lambda cid: rotated_energy[cid])

        # filter based on energy window ewin_csearch
        selectedcids_rotated = conformer_filters(self,sorted_rotated_cids,rotated_energy,rdmols)

        mol_select = []
        for i, cid in enumerate(selectedcids_rotated):
            mol_rd = Chem.RWMol(rdmols[cid])
            mol_rd.SetProp("_Name", rdmols[cid].GetProp("_Name") + " " + str(i))
            mol_rd.SetProp("Energy", str(rotated_energy[cid]))
            # setting the metal back instead of I
            if len(metal_atoms) >= 1:
                set_metal_atomic_number(
                    mol_rd, metal_idx, metal_sym
                )
            mol_select.append(mol_rd) 

        # update SDF file
        os.remove(csearch_file)
        sdwriter_upd = Chem.SDWriter(f'{csearch_file}')
        for mol in mol_select:
            sdwriter_upd.write(mol)
        sdwriter_upd.close()
        status = 1

        return status

    def auto_sampling(self, mol, metal_atoms, metal_idx):
        """
        Detects automatically the initial number of conformers for the sampling
        """

        auto_samples = 0
        # low: good for descriptor generation in machine learning
        if self.args.auto_sample.lower() == 'low':
            sampling_factor = 5
            max_confs = 100
        # mid: standard, good compromise between number of conformers and computing time
        elif self.args.auto_sample.lower() == 'mid':
            sampling_factor = 10
            max_confs = 250
        # high: more demanding method, more conformers and time
        elif self.args.auto_sample.lower() == 'high':
            sampling_factor = 20
            max_confs = 500
        else:
            self.args.log.write(f'x  {self.args.auto_sample} is not a valid option for --auto_sample! Please use "low", "mid" or "high"')
            self.args.log.finalize()
            sys.exit()

        if len(metal_atoms) >= 1:
            if len(metal_idx) > 0:
                auto_samples += 3 * len(metal_idx) # this accounts for possible trans/cis isomers in metal complexes

        auto_samples += 3 * (Lipinski.NumRotatableBonds(mol))  # x3, for C3 rotations
        auto_samples += 3 * (Lipinski.NHOHCount(mol))  # x3, for OH/NH rotations
        auto_samples += 3 * (Lipinski.NumSaturatedRings(mol))  # x3, for boat/chair/envelope confs
        if auto_samples == 0:
            auto_samples = sampling_factor
        else:
            auto_samples = sampling_factor * auto_samples

        # set a maximum number of conformers allowed with auto_sampling
        if auto_samples > max_confs:
            auto_samples = max_confs

        return auto_samples

    def genConformer_r(
        self,
        mol,
        conf,
        i,
        matches,
        name,
        sdwriter,
        update_to_rdkit,
        coord_Map,
        alg_Map,
        mol_template,
        original_atn,
        geom,
        metal_atoms,
        metal_idx,
        metal_sym,
        ff
    ):
        """
        If program = RDKit, this replaces iodine back to the metal (if needed) 
        and writes the RDKit SDF files. With program = summ, this function 
        optimizes rotamers
        """

        if i >= len(matches):  # base case, torsions should be set in conf
            if len(metal_atoms) >= 1:
                if self.args.program.lower() in ["rdkit","crest"] or update_to_rdkit:
                    if coord_Map is None and alg_Map is None and mol_template is None:
                        energy = minimize_rdkit_energy(
                            mol,
                            conf,
                            self.args.log,
                            ff,
                            self.args.opt_steps_rdkit,
                        )
                    else:
                        mol, energy = realign_mol(
                            mol,
                            conf,
                            coord_Map,
                            alg_Map,
                            mol_template,
                            self.args.opt_steps_rdkit,
                        )
                    mol.SetProp("Energy", str(energy))

                    # setting the metal back instead of I
                    set_metal_atomic_number(mol, metal_idx, metal_sym)

                    # setting the problematic As atoms back when using the Ir_squareplanar geometry rule
                    if geom == ['Ir_squareplanar']:
                        if original_atn is not None:
                            mol.GetAtomWithIdx(original_atn[1]).SetAtomicNum(original_atn[0])
            
            sdwriter.write(mol, conf)
            return 1

        elif self.args.program.lower() in ["crest"]:
            # setting the metal back instead of I
            set_metal_atomic_number(mol, metal_idx, metal_sym)
            sdwriter.write(mol, conf)

            return 1

        # when SUMM is selected, this cycle generates conformers based on rotation of dihedral angles
        total = 0
        deg = 0
        while deg < 360.0:
            rad = math.pi * deg / 180.0
            rdMolTransforms.SetDihedralRad(
                mol.GetConformer(conf), *matches[i], value=rad
            )
            mol.SetProp("_Name", name)
            total += self.genConformer_r(
                mol,
                conf,
                i + 1,
                matches,
                name,
                sdwriter,
                update_to_rdkit,
                coord_Map,
                alg_Map,
                mol_template,
                original_atn,
                geom,
                metal_atoms,
                metal_idx,
                metal_sym,
                ff
            )
            deg += int(self.args.degree)

        return total

    def embed_conf(self, mol, initial_confs, coord_Map, alg_Map, mol_template, csearch_nprocs, name):
        """
        Function to embed conformers
        """

        is_sdf_mol_or_mol2 = os.path.basename(Path(self.args.input)).split('.')[-1].lower() in [
            "sdf",
            "mol",
            "mol2",
        ]

        if is_sdf_mol_or_mol2:
            Chem.AssignStereochemistryFrom3D(mol)

        embed_kwargs = dict()
        embed_kwargs["ignoreSmoothingFailures"] = True
        embed_kwargs["randomSeed"] = self.args.seed
        embed_kwargs["numThreads"] = csearch_nprocs

        if (coord_Map, alg_Map, mol_template) != (None, None, None):
            embed_kwargs["coordMap"] = coord_Map
        cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

        if len(cids) <= 1 and initial_confs != 1:
            self.args.log.write(f"\nx  Normal RDKit embeding process failed, trying to generate conformers with random coordinates (with {str(initial_confs)} possibilities) ({os.path.basename(Path(name))})")
            embed_kwargs["useRandomCoords"] = True
            embed_kwargs["boxSizeMult"] = 10.0
            embed_kwargs["numZeroFail"] = 1000
            embed_kwargs["numThreads"] = csearch_nprocs
            cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

        if is_sdf_mol_or_mol2:
            # preserving AssignStereochemistryFrom3D
            for cid in cids:
                Chem.AssignAtomChiralTagsFromStructure(mol, confId=cid)

        return cids

    def min_and_E_calc(self, mol, cids, coord_Map, alg_Map, mol_template, 
                       ff, geom, metal_atoms, metal_idx, metal_sym):
        """
        Minimization and E calculation with RDKit after embeding
        """

        cenergy, outmols = [], []

        for _, conf in enumerate(cids):
            if coord_Map is None and alg_Map is None and mol_template is None:
                energy = minimize_rdkit_energy(
                    mol, conf, self.args.log, ff, self.args.opt_steps_rdkit
                )
            else:  # template realign before doing calculations
                mol, energy = realign_mol(
                    mol,
                    conf,
                    coord_Map,
                    alg_Map,
                    mol_template,
                    self.args.opt_steps_rdkit,
                )

            # removes geometries that do not pass the filters (geom option)
            mol_geom = Chem.Mol(mol)
            # setting the metal back instead of I
            if len(metal_atoms) >= 1:
                set_metal_atomic_number(mol_geom, metal_idx, metal_sym)

            passing_geom = geom_filter(self,mol_geom,geom)
            if passing_geom:
                cenergy.append(energy)
                pmol = PropertyMol.PropertyMol(mol)
                outmols.append(pmol)

        return outmols, cenergy

    def min_after_embed(
        self,
        mol,
        cids,
        name,
        csearch_file,
        rotmatches,
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
        metal_sym
    ):
        """
        Minimizes, gets the energy and filters RDKit conformers after embeding
        """

        # gets optimized mol objects and energies
        if geom != []:
            self.args.log.write(f"o  Applying geometry filters ({geom}) ({os.path.basename(Path(name))})")
        outmols, cenergy = self.min_and_E_calc(
            mol, cids, coord_Map, alg_Map, mol_template, ff, geom, metal_atoms, metal_idx, metal_sym
        )

        for i, cid in enumerate(cids):
            outmols[cid].SetProp("_Name", name + " " + str(i + 1))
            outmols[cid].SetProp("Energy", str(cenergy[cid]))
            outmols[cid].SetProp("Real charge", str(charge))
            outmols[cid].SetProp("Mult", str(mult))
            outmols[cid].SetProp("SMILES", str(smi))

        # sorts the energies
        cids = list(range(len(outmols)))
        sorted_all_cids = sorted(cids, key=lambda cid: cenergy[cid])

        self.args.log.write(f"\no  Applying filters to initial conformers ({os.path.basename(Path(name))})")
        selectedcids_rdkit = conformer_filters(self,sorted_all_cids,cenergy,outmols)

        if self.args.program.lower() in ["summ", "rdkit", "crest"]:
            sdwriter_summ = Chem.SDWriter(f'{csearch_file}')            
            # now exhaustively drive torsions of selected conformers
            total = 0
            sdwriter = Chem.SDWriter(f'{csearch_file}')
            for conf in selectedcids_rdkit:
                if self.args.program.lower() == "summ" and not update_to_rdkit:
                    sdwriter_summ.write(outmols[conf], conf)
                    for m in rotmatches:
                        rdMolTransforms.SetDihedralDeg(
                            outmols[conf].GetConformer(conf), *m, 180.0
                        )

                if self.args.program.lower() in ["summ", "rdkit", "crest"]:
                    total += self.genConformer_r(
                        outmols[conf],
                        conf,
                        0,
                        rotmatches,
                        outmols[conf].GetProp("_Name"),
                        sdwriter,
                        update_to_rdkit,
                        coord_Map,
                        alg_Map,
                        mol_template,
                        original_atn,
                        geom,
                        metal_atoms,
                        metal_idx,
                        metal_sym,
                        ff
                    )

            sdwriter.close()
            sdwriter_summ.close()
            status = 1

        # keep only structurally different conformers
        if self.args.program.lower() in ["rdkit","crest"]:
            suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(f'{csearch_file}', "csearch", self.args)
            if len(suppl) > self.args.sample and self.args.auto_cluster:
                cluster_mols_sorted = cluster_conformers(self,suppl,"rdkit",csearch_file,name)
                outmols = cluster_mols_sorted
            else:
                outmols = suppl
        
        if self.args.program.lower() == "fullmonte":
            status = generating_conformations_fullmonte(
                name,
                self.args,
                rotmatches,
                selectedcids_rdkit,
                outmols,
                csearch_file,
                coord_Map,
                alg_Map,
                mol_template,
                ff,
                metal_atoms,
                metal_idx, 
                metal_sym
            )
            # removes the rdkit file
            os.remove(name + "_" + "rdkit" + self.args.output)

        return status, outmols

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
        csearch_nprocs
    ):

        """
        Conversion from RDKit to SDF
        """

        mol.SetProp("_Name", name)

        # detects and applies auto-detection of initial number of conformers
        if self.args.auto_sample:
            initial_confs = int(self.auto_sampling(mol,metal_atoms,metal_idx))
        else:
            initial_confs = self.args.sample

        update_to_rdkit = False

        rotmatches = getDihedralMatches(mol, self.args.heavyonly)

        if len(rotmatches) > self.args.max_torsions and self.args.max_torsions > 0:
            self.args.log.write(f"\nx  Too many torsions ({len(rotmatches)}). Skipping {name + self.args.output}")
        elif self.args.program.lower() == "summ" and len(rotmatches) == 0:
            update_to_rdkit = True
            self.args.log.write(f"\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to SUMM SDF ({os.path.basename(Path(name))})")
        elif self.args.program.lower() == "fullmonte" and len(rotmatches) == 0:
            update_to_rdkit = True
            self.args.log.write(f"\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to FULLMONTE SDF ({os.path.basename(Path(name))})")

        ff = self.args.ff
        if self.args.program.lower() == "rdkit":
            rotmatches = []
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
                rotmatches,
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
                metal_sym
            )
        except IndexError:
            status = -1
            mol_crest = None

        return status, rotmatches, ff, mol_crest
