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
   varfile : str, default=None  
      Option to parse the variables using a yaml file (specify the filename)  
   max_workers : int, default=4  
      Number of simultaneous RDKit jobs run with multiprocessing 
      (WARNING! More than 12 simultaneous jobs might collapse your computer!)  
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

   sample : int, default='auto'
      Number of conformers used initially in the RDKit sampling. If this option 
      isn't specified, AQME automatically calculates (previously benchmarked) an
      approximate number based on number of rotatable bonds, XH (i.e. OH) groups, 
      saturated cycles, etc (see the auto_sampling() function in csearch.py for 
      more information)
   auto_sample : int, default=20
      Base multiplicator number used in the sample option
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

Only organometallic molecules
.............................

   metal_atoms : list of str, default=[]
     Specify metal atom(s) of the system as [ATOM_TYPE]. Multiple metals can 
     be used simultaneously (i.e. ['Pd','Ir']).  This option is important to 
     calculate the charge of metal complexes based on SMILES strings. Requires 
     the use of metal_oxi.
   metal_oxi : list of int, default=[]
     Specify metal oxidation state as [NUMBER]. Multiple metals can be used 
     simultaneously (i.e. [2,3]).
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

Crest only
++++++++++

   nprocs : int, default=2
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
   cregen : bool, default=False
      If True, perform a CREGEN analysis after CREST (filtering options below)
   cregen_keywords : str, default=None
      Additional keywords for CREGEN (i.e. cregen_keywords='--ethr 0.02')
   xtb_keywords : str, default=None
      Define additional keywords to use in the xTB pre-optimization that are not 
      included in -c, --uhf, -P and --input. For example: '--alpb ch2cl2 --gfn 1' 
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
from pathlib import Path
import pandas as pd
import concurrent.futures as futures
import multiprocessing as mp
from progress.bar import IncrementalBar

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors as Descriptors
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, Lipinski

from aqme.filter import filters, ewin_filter, pre_E_filter, RMSD_and_E_filter
from aqme.csearch.utils import (
    prepare_direct_smi,
    prepare_smiles_files,
    prepare_csv_files,
    prepare_cdx_files,
    prepare_com_files,
    prepare_sdf_files,
    prepare_pdb_files,
    creation_of_dup_csv_csearch,
    minimize_rdkit_energy,
    com_2_xyz,
    check_constraints,
    smi_to_mol,
    getDihedralMatches
)
from aqme.csearch.templates import template_embed, check_metal_neigh
from aqme.csearch.fullmonte import generating_conformations_fullmonte, realign_mol
from aqme.utils import (
    rules_get_charge,
    substituted_mol,
    load_variables,
    set_metal_atomic_number
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

        if self.args.program.lower() == "crest":
            try:
                subprocess.run(
                    ["xtb", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
            except FileNotFoundError:
                self.args.log.write("x  xTB is not installed (CREST cannot be used)! You can install the program with 'conda install -c conda-forge xtb'")
                self.args.log.finalize()
                sys.exit()

        csearch_program = True
        if self.args.program is None:
            csearch_program = False
        if csearch_program:
            if self.args.program.lower() not in ["rdkit", "summ", "fullmonte", "crest"]:
                csearch_program = False
        if not csearch_program:
            self.args.log.write("\nx  Program not supported for CSEARCH conformer generation! Specify: program='rdkit' (or summ, fullmonte, crest)")
            self.args.log.finalize()
            sys.exit()

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
            csearch_files = glob.glob(self.args.input)
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
            self.run_csearch(job_inputs)

            # store all the information into a CSV file
            csearch_file_no_path = (
                csearch_file.replace("/", "\\").split("\\")[-1].split(".")[0]
            )
            self.csearch_csv_file = self.args.w_dir_main.joinpath(
                f"CSEARCH-Data-{csearch_file_no_path}.csv"
            )
            self.final_dup_data.to_csv(self.csearch_csv_file, index=False)

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
            ".smi",
            ".sdf",
            ".cdx",
            ".csv",
            ".com",
            ".gjf",
            ".mol",
            ".mol2",
            ".xyz",
            ".txt",
            ".yaml",
            ".yml",
            ".rtf",
            ".pdb",
        ]

        file_format = os.path.splitext(csearch_file)[1]
        # Checks
        if file_format.lower() not in SUPPORTED_INPUTS:
            self.args.log.write("\nx  Input filetype not currently supported!")
            self.args.log.finalize()
            sys.exit()

        smi_derivatives = [".smi", ".txt", ".yaml", ".yml", ".rtf"]
        Extension2inputgen = dict()
        for key in smi_derivatives:
            Extension2inputgen[key] = prepare_smiles_files
        Extension2inputgen[".csv"] = prepare_csv_files
        Extension2inputgen[".cdx"] = prepare_cdx_files
        Extension2inputgen[".gjf"] = prepare_com_files
        Extension2inputgen[".com"] = prepare_com_files
        Extension2inputgen[".xyz"] = prepare_com_files
        Extension2inputgen[".sdf"] = prepare_sdf_files
        Extension2inputgen[".mol"] = prepare_sdf_files
        Extension2inputgen[".mol2"] = prepare_sdf_files
        Extension2inputgen[".pdb"] = prepare_pdb_files

        # Prepare the jobs
        prepare_function = Extension2inputgen[file_format]
        job_inputs = prepare_function(self.args, csearch_file)

        return job_inputs

    def run_csearch(self, job_inputs):
        # create the dataframe to store the data
        self.final_dup_data = creation_of_dup_csv_csearch(self.args.program.lower())

        bar = IncrementalBar(
            "o  Number of finished jobs from CSEARCH", max=len(job_inputs)
        )
        with futures.ProcessPoolExecutor(
            max_workers=self.args.max_workers, mp_context=mp.get_context("spawn")
        ) as executor:
            # Submit a set of asynchronous jobs
            jobs = []
            # Submit the Jobs
            for job_input in job_inputs:
                (
                    smi_,
                    name_,
                    charge_,
                    mult_,
                    constraints_atoms_,
                    constraints_dist_,
                    constraints_angle_,
                    constraints_dihedral_,
                ) = job_input
                job = executor.submit(
                    self.compute_confs(
                        smi_,
                        name_,
                        charge_,
                        mult_,
                        constraints_atoms_,
                        constraints_dist_,
                        constraints_angle_,
                        constraints_dihedral_,
                    )
                )
                jobs.append(job)

                bar.next()

            bar.finish()

    def compute_confs(
        self,
        smi,
        name,
        charge,
        mult,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
    ):
        """
        Function to start conformer generation
        """

        self.args.log.write(f"\n   ----- {name} -----")

        if self.args.smi is not None or self.args.input.split(".")[1] in ["smi","csv","cdx","txt","yaml","yml","rtf"]:
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
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
            )
            if mol is None:
                self.args.log.write(f"\nx  Failed to convert the provided SMILES ({smi}) to an RDkit Mol object! Please check the starting smiles.")
                # if a list of SMILES is provided, the program doesn't stop if one SMILES fails to convert to mol
                if self.args.input.split(".")[1] not in ["csv","cdx","txt","yaml","yml","rtf"]:
                    self.args.log.finalize()
                    sys.exit()
                return

        else:
            # for 3D input formats, the smi variable represents the mol object
            mol = smi
            if mol is None:
                self.args.log.write(f"\nx  Failed to convert the provided input to an RDkit Mol object! Please check the starting structure.")
                self.args.log.finalize()
                sys.exit()
                
            # check if the optimization is constrained
            complex_ts = check_constraints(self)

        # this part allows to use charges and multiplicity from COM/GJF inputs
        if self.args.charge is None:
            self.args.charge = charge
        if self.args.mult is None:
            self.args.mult = mult

        if self.args.destination is None:
            self.csearch_folder = Path(self.args.initial_dir).joinpath(
                f"CSEARCH"
            )
        else:
            if Path(f"{self.args.destination}").exists():
                self.csearch_folder = Path(self.args.destination)
            else:
                self.csearch_folder = Path(self.args.initial_dir).joinpath(
                    self.args.destination
                )

        self.csearch_folder.mkdir(exist_ok=True, parents=True)

        if self.args.program.lower() in ["crest"] and self.args.smi is None:
            if self.args.input.split(".")[1] in ["pdb", "mol2", "mol", "sdf"]:
                command_pdb = [
                    "obabel",
                    f'-i{self.args.input.split(".")[1]}',
                    f'{name}.{self.args.input.split(".")[1]}',
                    "-oxyz",
                    f"-O{name}_{self.args.program.lower()}.xyz",
                ]
                subprocess.run(
                    command_pdb,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            elif self.args.input.split(".")[1] in ["gjf", "com"]:
                xyz_file, _, _ = com_2_xyz(f'{name}.{self.args.input.split(".")[1]}')
                os.move(xyz_file, f"{name}_{self.args.program.lower()}.xyz")
            elif self.args.input.split(".")[1] == "xyz":
                shutil.copy(f"{name}.xyz", f"{name}_{self.args.program.lower()}.xyz")

        template_opt = False
        # replaces the metal for an I atom
        if len(self.args.metal_atoms) >= 1:
            (
                self.args.metal_idx,
                self.args.complex_coord,
                self.args.metal_sym,
            ) = substituted_mol(self, mol, "I")

            # get pre-determined geometries for metal complexes
            accepted_complex_types = [
                "squareplanar",
                "squarepyramidal",
                "linear",
                "trigonalplanar",
            ]
            if self.args.complex_type != '' and self.args.complex_type not in accepted_complex_types:
                self.args.log.write(f"x  The metal template specified in complex_type ({self.args.complex_type}) is not valid! Options: squareplanar, squarepyramidal, linear and trigonalplanar")
                self.args.log.finalize()
                sys.exit()

            if self.args.complex_type in accepted_complex_types:
                count_metals = 0
                valid_template = True
                # check if the specified metal is included in the system
                for metal_idx_ind in self.args.metal_idx:
                    if metal_idx_ind is not None:
                        # calculate number of expected neighbours
                        valid_template = check_metal_neigh(mol, self.args.complex_type, metal_idx_ind, self.args.log, valid_template)
                        count_metals += 1

                if count_metals == 1 and valid_template:
                    template_opt = True
                    template_kwargs = dict()
                    template_kwargs["complex_type"] = self.args.complex_type
                    template_kwargs["metal_idx"] = self.args.metal_idx
                    template_kwargs["maxsteps"] = self.args.opt_steps_rdkit
                    template_kwargs["heavyonly"] = self.args.heavyonly
                    template_kwargs["maxmatches"] = self.args.max_matches_rmsd
                    template_kwargs["mol"] = mol
                    template_kwargs["name"] = name
                    items = template_embed(self, **template_kwargs)

                    total_data = creation_of_dup_csv_csearch(self.args.program.lower())

                    for mol_obj, name_in, coord_map, alg_map, template in zip(*items):
                        data = self.conformer_generation(
                            mol_obj,
                            name_in,
                            constraints_atoms,
                            constraints_dist,
                            constraints_angle,
                            constraints_dihedral,
                            complex_ts,
                            coord_map,
                            alg_map,
                            template,
                        )
                        frames = [total_data, data]
                        total_data = pd.concat(frames, sort=True)

                elif count_metals > 1 or count_metals == 0:
                    self.args.log.write(f"\nx  The template specified {self.args.complex_type} is not used for systems with more than 1 metal or for organic molecueles.")

        if not template_opt:
            total_data = self.conformer_generation(
                mol,
                name,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
                complex_ts
            )

        # Updates the dataframe with infromation about conformer generation
        frames = [self.final_dup_data, total_data]
        self.final_dup_data = pd.concat(frames, ignore_index=True, sort=True)

    def conformer_generation(
        self,
        mol,
        name,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
        complex_ts,
        coord_Map=None,
        alg_Map=None,
        mol_template=None,
    ):
        """
        Function to load mol objects and create 3D conformers

        """

        dup_data = creation_of_dup_csv_csearch(self.args.program.lower())

        dup_data_idx = 0
        start_time = time.time()
        status = None

        # check if metals and oxidation states are both used
        if self.args.metal_atoms != []:
            if self.args.metal_oxi == []:
                self.args.log.write(f"\nx   Metal atoms ({self.args.metal_atoms}) were specified without their corresponding oxidation state (metal_oxi option)")
                self.args.log.finalize()
                sys.exit()

        if self.args.metal_oxi != []:
            if self.args.metal_atoms == []:
                self.args.log.write(f"\nx   Metal oxidation states ({self.args.metal_oxi}) were specified without their corresponding metal atoms (metal_atoms option)")
                self.args.log.finalize()
                sys.exit()

        # Set charge and multiplicity
        metal_found = False
        # user can overwrite charge and mult with the corresponding arguments
        if self.args.charge is not None:
            charge = self.args.charge
        else:
            if not len(self.args.metal_atoms) >= 1:
                charge = Chem.GetFormalCharge(mol)
            else:
                charge, metal_found = rules_get_charge(mol, self.args)

        if self.args.mult is not None:
            mult = self.args.mult
        else:
            mult = Descriptors.NumRadicalElectrons(mol) + 1
            if metal_found:
                # since RDKit gets the multiplicity of the metal with valence 0, the real multiplicity
                # value needs to be adapted with the charge. If multiplicity is different than 1 or 2,
                # the user must specify the value with the mult option
                if (charge % 2) == 1 and charge != 0: # odd charges (i.e. +1, +3, etc)
                    if mult == 1:
                        mult = mult + 1
                    if mult == 2:
                        mult = mult - 1

        dup_data.at[dup_data_idx, "Real charge"] = charge
        dup_data.at[dup_data_idx, "Mult"] = mult

        # inputs that go through CREST containing 3D coordinates don't require a previous RDKit conformer sampling
        if (
            self.args.program.lower() in ["crest"]
            and self.args.smi is None
            and self.args.input.split(".")[1] in ["pdb","mol2","mol","sdf","gjf","com","xyz"]
        ):

            valid_structure = True
            status = xtb_opt_main(
                f"{name}_{self.args.program.lower()}",
                dup_data,
                dup_data_idx,
                self,
                charge,
                mult,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
                'crest',
                mol=mol
            )

        else:
            name = name.replace("/", "\\").split("\\")[-1].split(".")[0]
            self.csearch_file = self.csearch_folder.joinpath(
                name + "_" + self.args.program.lower() + self.args.output
            )
            sdwriter_init = Chem.SDWriter(str(self.csearch_file))

            valid_structure = filters(
                mol, self.args.log, self.args.max_mol_wt
            )
            if valid_structure:
                try:
                    # the conformational search for RDKit
                    status = self.summ_search(
                        mol,
                        name,
                        sdwriter_init,
                        dup_data,
                        dup_data_idx,
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
                    )
                except (KeyboardInterrupt, SystemExit):
                    raise

        if status == -1 or not valid_structure:
            error_message = "\nx  ERROR: The structure is not valid or no conformers were obtained from this SMILES string"
            self.args.log.write(error_message)

        # if self.args.program.lower() == "crest" and valid_structure:
        #     shutil.rmtree(f"{self.csearch_folder}/../crest_xyz")

        n_seconds = round(time.time() - start_time, 2)
        dup_data.at[dup_data_idx, "CSEARCH time (seconds)"] = n_seconds

        return dup_data

    def summ_search(
        self,
        mol,
        name,
        sdwriter,
        dup_data,
        dup_data_idx,
        charge,
        mult,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
        complex_ts,
        coord_Map,
        alg_Map,
        mol_template
    ):

        """
        Embeds, optimizes and filters RDKit conformers
        """

        # writes sdf for the first RDKit conformer generation
        if not complex_ts:
            if self.args.program.lower() in ['rdkit']:
                self.args.log.write(f"\no  Starting RDKit conformer sampling")
            elif self.args.program.lower() in ['summ','fullmonte']:
                self.args.log.write(f"\no  Starting RDKit-{self.args.program} conformer sampling")
            elif self.args.program.lower() in ['crest'] and self.args.input.split(".")[1] not in ["pdb","mol2","mol","sdf","gjf","com","xyz"]:
                self.args.log.write(f"\no  Starting initial RDKit-based mol generation from SMILES")

            status, rotmatches, ff, mol_crest = self.rdkit_to_sdf(
                mol,
                name,
                dup_data,
                dup_data_idx,
                sdwriter,
                charge,
                mult,
                coord_Map,
                alg_Map,
                mol_template
            )
            # this avoids memory issues when using Windows
            try:
                sdwriter.close()
            except RuntimeError:
                pass
            # reads the initial SDF files from RDKit and uses dihedral scan if selected
            if status not in [-1, 0]:
                # getting the energy and mols after rotations
                if self.args.program.lower() == "summ" and len(rotmatches) != 0:
                    status = self.dihedral_filter_and_sdf(
                        name, dup_data, dup_data_idx, coord_Map, alg_Map, mol_template, ff
                    )

        if self.args.program.lower() in ['crest']:
            if not complex_ts:
                # mol_crest is the RDKit-optimized mol object
                rdmolfiles.MolToXYZFile(mol_crest, name + "_crest.xyz")
            else:
                # mol is the raw mol object (no optimization with RDKit to avoid problems when using
                # noncovalent complexes and TSs)
                rdmolfiles.MolToXYZFile(mol, name + "_crest.xyz")

            status = xtb_opt_main(
                    f"{name}_{self.args.program.lower()}",
                    dup_data,
                    dup_data_idx,
                    self,
                    charge,
                    mult,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                    'crest',
                    complex_ts=complex_ts,
                    mol=mol # this is necessary for CREST calculations with constraints
                )

        return status

    def dihedral_filter_and_sdf(
        self, name, dup_data, dup_data_idx, coord_Map, alg_Map, mol_template, ff
    ):
        """
        Filtering after dihedral scan to sdf
        """

        rotated_energy = []
        # apply filters
        with Chem.SDMolSupplier(str(self.csearch_file), removeHs=False) as rdmols:

            if rdmols is None:
                self.args.log.write("\nCould not open " + name + self.args.output)
                self.args.log.finalize()
                sys.exit()

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
            sortedcids_rotated = ewin_filter(
                sorted_rotated_cids,
                rotated_energy,
                dup_data,
                dup_data_idx,
                "summ",
                self.args.ewin_csearch,
            )
            # pre-filter based on energy only
            selectedcids_initial_rotated = pre_E_filter(
                sortedcids_rotated,
                rotated_energy,
                dup_data,
                dup_data_idx,
                "summ",
                self.args.initial_energy_threshold,
            )
            # filter based on energy and RMSD
            selectedcids_rotated = RMSD_and_E_filter(
                rdmols,
                selectedcids_initial_rotated,
                rotated_energy,
                self.args,
                dup_data,
                dup_data_idx,
                "summ",
            )
            mol_select = []
            for i, cid in enumerate(selectedcids_rotated):
                mol_rd = Chem.RWMol(rdmols[cid])
                mol_rd.SetProp("_Name", rdmols[cid].GetProp("_Name") + " " + str(i))
                mol_rd.SetProp("Energy", str(rotated_energy[cid]))
                if len(self.args.metal_atoms) >= 1:
                    set_metal_atomic_number(
                        mol_rd, self.args.metal_idx, self.args.metal_sym
                    )
                mol_select.append(mol_rd) 

        # update SDF file            
        os.remove(self.csearch_file)
        sdwriter_rd = Chem.SDWriter(str(self.csearch_file))
        for mol in mol_select:
            sdwriter_rd.write(mol)
        sdwriter_rd.close()
        status = 1
        return status

    def auto_sampling(self, mol):
        """
        Detects automatically the initial number of conformers for the sampling
        """

        if len(self.args.metal_atoms) >= 1:
            if len(self.args.metal_idx) > 0:
                self.args.auto_sample = (
                    self.args.auto_sample * 3 * len(self.args.metal_idx)
                )  # this accounts for possible trans/cis isomers in metal complexes
        auto_samples = 0
        auto_samples += 3 * (Lipinski.NumRotatableBonds(mol))  # x3, for C3 rotations
        auto_samples += 3 * (Lipinski.NHOHCount(mol))  # x3, for OH/NH rotations
        auto_samples += 3 * (
            Lipinski.NumSaturatedRings(mol)
        )  # x3, for boat/chair/envelope confs
        if auto_samples == 0:
            auto_samples = self.args.auto_sample
        else:
            auto_samples = self.args.auto_sample * auto_samples
        return auto_samples

    def genConformer_r(
        self,
        mol,
        conf,
        i,
        matches,
        sdwriter,
        name,
        update_to_rdkit,
        coord_Map,
        alg_Map,
        mol_template,
    ):
        """
        If program = RDKit, this replaces iodine back to the metal (if needed) 
        and writes the RDKit SDF files. With program = summ, this function 
        optimizes rotamers
        """

        if i >= len(matches):  # base case, torsions should be set in conf
            if len(self.args.metal_atoms) >= 1 and (
                self.args.program.lower() in ["rdkit","crest"] or update_to_rdkit
            ):
                if coord_Map is None and alg_Map is None and mol_template is None:
                    energy = minimize_rdkit_energy(
                        mol,
                        conf,
                        self.args.log,
                        self.args.ff,
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
                set_metal_atomic_number(mol, self.args.metal_idx, self.args.metal_sym)
            try:
                sdwriter.write(mol, conf)
            except (TypeError):
                raise

            # if CREST is used, this RDKit preoptimzed mol object will be employed to initializethe the trajectories
            if self.args.program.lower() in ["crest"]:
                return mol
            else:
                return 1

        if self.args.program.lower() in ["crest"]:
            return mol

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
                sdwriter,
                name,
                update_to_rdkit,
                coord_Map,
                alg_Map,
                mol_template,
            )
            deg += self.args.degree

        return total

    def embed_conf(self, mol, initial_confs, coord_Map, alg_Map, mol_template):
        """
        Function to embed conformers
        """

        is_sdf_mol_or_mol2 = os.path.splitext(self.args.input)[1].lower() in [
            ".sdf",
            ".mol",
            ".mol2",
        ]

        if is_sdf_mol_or_mol2:
            Chem.AssignStereochemistryFrom3D(mol)

        embed_kwargs = dict()
        embed_kwargs["ignoreSmoothingFailures"] = True
        embed_kwargs["randomSeed"] = self.args.seed
        embed_kwargs["numThreads"] = 0

        if (coord_Map, alg_Map, mol_template) != (None, None, None):
            embed_kwargs["coordMap"] = coord_Map
        cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

        if len(cids) <= 1 and initial_confs != 1:
            self.args.log.write(f"\nx  Normal RDKit embeding process failed, trying to generate conformers with random coordinates (with {str(initial_confs)} possibilities)")
            embed_kwargs["useRandomCoords"] = True
            embed_kwargs["boxSizeMult"] = 10.0
            embed_kwargs["numZeroFail"] = 1000
            embed_kwargs["numThreads"] = 1
            cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

        if is_sdf_mol_or_mol2:
            # preserving AssignStereochemistryFrom3D
            for cid in cids:
                Chem.AssignAtomChiralTagsFromStructure(mol, confId=cid)

        return cids

    def min_and_E_calc(self, mol, cids, coord_Map, alg_Map, mol_template, ff):
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
            cenergy.append(energy)
            pmol = PropertyMol.PropertyMol(mol)
            outmols.append(pmol)

        return outmols, cenergy

    def min_after_embed(
        self,
        mol,
        cids,
        name,
        rotmatches,
        dup_data,
        dup_data_idx,
        sdwriter,
        update_to_rdkit,
        coord_Map,
        alg_Map,
        mol_template,
        charge,
        mult,
        ff,
    ):
        """
        Minimizes, gets the energy and filters RDKit conformers after embeding
        """

        # gets optimized mol objects and energies
        outmols, cenergy = self.min_and_E_calc(
            mol, cids, coord_Map, alg_Map, mol_template, ff
        )

        # writing charges and multiplicity after RDKit
        dup_data.at[dup_data_idx, "Mult"] = mult
        dup_data.at[dup_data_idx, "Real charge"] = charge

        for i, cid in enumerate(cids):
            outmols[cid].SetProp("_Name", name + " " + str(i + 1))
            outmols[cid].SetProp("Energy", str(cenergy[cid]))
            outmols[cid].SetProp("Real charge", str(charge))
            outmols[cid].SetProp("Mult", str(mult))

        # sorts the energies
        cids = list(range(len(outmols)))
        sorted_all_cids = sorted(cids, key=lambda cid: cenergy[cid])

        self.args.log.write("\no  Applying filters to initial conformers")

        # filter based on energy window ewin_csearch
        sortedcids_rdkit = ewin_filter(
            sorted_all_cids,
            cenergy,
            dup_data,
            dup_data_idx,
            "rdkit",
            self.args.ewin_csearch,
        )

        # pre-filter based on energy only
        selectedcids_initial_rdkit = pre_E_filter(
            sortedcids_rdkit,
            cenergy,
            dup_data,
            dup_data_idx,
            "rdkit",
            self.args.initial_energy_threshold,
        )

        # filter based on energy and RMSD
        selectedcids_rdkit = RMSD_and_E_filter(
            outmols,
            selectedcids_initial_rdkit,
            cenergy,
            self.args,
            dup_data,
            dup_data_idx,
            "rdkit",
        )

        if self.args.program.lower() in ["summ", "rdkit", "crest"]:
            # now exhaustively drive torsions of selected conformers
            total = 0
            for conf in selectedcids_rdkit:
                if self.args.program.lower() == "summ" and not update_to_rdkit:
                    sdwriter.write(outmols[conf], conf)
                    for m in rotmatches:
                        rdMolTransforms.SetDihedralDeg(
                            outmols[conf].GetConformer(conf), *m, 180.0
                        )
                if self.args.program.lower() in ["summ", "rdkit"]:
                    total += self.genConformer_r(
                        outmols[conf],
                        conf,
                        0,
                        rotmatches,
                        sdwriter,
                        outmols[conf].GetProp("_Name"),
                        update_to_rdkit,
                        coord_Map,
                        alg_Map,
                        mol_template,
                    )
                elif self.args.program.lower() in ["crest"]:
                    mol = self.genConformer_r(
                        outmols[conf],
                        conf,
                        0,
                        rotmatches,
                        sdwriter,
                        outmols[conf].GetProp("_Name"),
                        update_to_rdkit,
                        coord_Map,
                        alg_Map,
                        mol_template,
                    )
                    break

            status = 1

        if self.args.program.lower() == "summ":
            dup_data.at[dup_data_idx, "summ-conformers"] = total

        if self.args.program.lower() == "fullmonte":
            status = generating_conformations_fullmonte(
                name,
                self.args,
                rotmatches,
                selectedcids_rdkit,
                outmols,
                sdwriter,
                dup_data,
                dup_data_idx,
                coord_Map,
                alg_Map,
                mol_template,
                ff,
            )
            # removes the rdkit file
            os.remove(name + "_" + "rdkit" + self.args.output)

        return status, mol

    def rdkit_to_sdf(
        self,
        mol,
        name,
        dup_data,
        dup_data_idx,
        sdwriter,
        charge,
        mult,
        coord_Map,
        alg_Map,
        mol_template
    ):

        """
        Conversion from RDKit to SDF
        """

        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        mol.SetProp("_Name", name)

        # detects and applies auto-detection of initial number of conformers
        if self.args.program.lower() in ['crest']:
            # CREST only uses the most stable conformer from RDKit,
            # and initial_confs can be low to speed up the process
            initial_confs = self.args.auto_sample
        elif self.args.sample == "auto":
            initial_confs = int(self.auto_sampling(mol))
        else:
            initial_confs = int(self.args.sample)

        dup_data.at[dup_data_idx, "Molecule"] = name
        update_to_rdkit = False

        rotmatches = getDihedralMatches(mol, self.args.heavyonly)

        if len(rotmatches) > self.args.max_torsions and self.args.max_torsions > 0:
            self.args.log.write(f"\nx  Too many torsions ({len(rotmatches)}). Skipping {name + self.args.output}")
        elif self.args.program.lower() == "summ" and len(rotmatches) == 0:
            update_to_rdkit = True
            self.args.log.write("\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to SUMM SDF")
        elif self.args.program.lower() == "fullmonte" and len(rotmatches) == 0:
            update_to_rdkit = True
            self.args.log.write("\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to FULLMONTE SDF")

        ff = self.args.ff
        dup_data.at[dup_data_idx, "RDKit-Initial-samples"] = initial_confs
        if self.args.program.lower() == "rdkit":
            rotmatches = []
        cids = self.embed_conf(mol, initial_confs, coord_Map, alg_Map, mol_template)

        # energy minimize all to get more realistic results
        # identify the atoms and decide Force Field
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 36 and self.args.ff == "MMFF":  # up to Kr for MMFF, if not the code will use UFF
                self.args.log.write(f"\nx  {self.args.ff} is not compatible with the molecule, changing to UFF")
                ff = "UFF"

        try:
            status,mol_crest = self.min_after_embed(
                mol,
                cids,
                name,
                rotmatches,
                dup_data,
                dup_data_idx,
                sdwriter,
                update_to_rdkit,
                coord_Map,
                alg_Map,
                mol_template,
                charge,
                mult,
                ff,
            )
        except IndexError:
            status = -1
            mol_crest = None

        sdwriter.close()

        return status, rotmatches, ff, mol_crest
