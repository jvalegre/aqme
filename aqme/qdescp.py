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

xTB descriptors
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
    find_level_names
)

from aqme.csearch.crest import xyzall_2_xyz


class qdescp:
    """
    Class containing all the functions from the QDESCP module
    """

    def __init__(self, **kwargs):

        start_time_overall = time.time()

        # load default and user-specified variables
        self.args = load_variables(kwargs, "qdescp")

        # detects errors and updates variables before the QDESCP run
        self,destination,smarts_targets,boltz_dir = self.qdescp_set_up()

        # check whether dependencies are installed
        _ = check_dependencies(self)

        # full xTB workflow in QDESCP for descriptor generation and collection
        if self.args.program.lower() == "xtb":
            _ = self.qdescp_xtb_workflow(boltz_dir,destination,smarts_targets)

        # full NMR workflow in QDESCP for NMR prediction
        elif self.args.program.lower() == "nmr":
            _ = self.qdescp_nmr_workflow(boltz_dir)
        
        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime QDESCP: {elapsed_time} seconds\n")
        self.args.log.finalize()
    

    def qdescp_xtb_workflow(self,boltz_dir,destination,smarts_targets):
        '''
        Full xTB workflow in QDESCP for descriptor generation and collection
        '''

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
        '''
        Check whether the user have chosen the "input" or "files" option (QDESCP will use "files" from this point on)
        '''
        
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
        '''
        Initial conformer generation for QDESCP runs that start with SMILES strings from CSV inputs
        '''
        
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
                if os.path.basename(file).replace('_rdkit.sdf','') in df_qdescp['code_name'].astype(str).tolist():
                    qdescp_files.append(file)
        if len(qdescp_files) == 0:
            self.args.log.write(f"\nx  WARNING! The CSEARCH conformational search did not produce any results.")
            self.args.log.finalize()
            sys.exit()

        return qdescp_files
    

    def get_boltz_n_save_csv(self,destination,qdescp_files,descp_dict,boltz_dir,smarts_targets):
        self.args.log.write('\no  Running RDKit and collecting molecular properties (for all inputs)')
        all_prefixes_atoms =  []
        for _,file in enumerate(qdescp_files):
            descp_dict_indiv = descp_dict.copy()
            if file not in self.args.invalid_calcs:
                mols = load_sdf(file)
                mol = mols[0]
                name = '.'.join(os.path.basename(Path(file)).split(".")[:-1])
                # to locate difficult names (i.e. with special characters), glob.glob doesn't work, this is needed:
                json_files = [x for x in glob.glob(f"{destination}/*.json") if os.path.basename(x).startswith(f'{name}_conf_')]
                # Generating the JSON files
                _ = self.get_boltz_props(json_files, name, boltz_dir, "xtb", descp_dict_indiv, smarts_targets, mol, all_prefixes_atoms)
            
        # Create the CSV files from the JSON files
        _ = self.write_csv_boltz_data(destination)  # CSV full


    def get_boltz_props(self, json_files, name, boltz_dir, calc_type, descp_dict_indiv, smarts_targets, mol, all_prefixes_atoms):
        """
        Retrieves the properties from json files and gives Boltzmann averaged properties for rdkit, NMR and morfues descriptors.
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
        _ = dict_to_json(os.path.join(boltz_dir, f"{name}_boltz.json"), full_json_data)


    def qdescp_set_up(self):
        '''
        Detects errors and updates variables before the QDESCP run
        '''

        # most users employ QDESCP to generate descriptors with xTB
        if self.args.program is None:
            self.args.program = "xtb"

        if self.args.program.lower() not in ["xtb", "nmr"]:
            self.args.log.write(f"\nx  The program specified ({self.args.program}) is not supported for QDESCP descriptor generation! Specify: program='xtb' (or nmr)")
            self.args.log.finalize()
            sys.exit()

        # set number of processors
        if self.args.nprocs is None:
            self.args.nprocs = 8

        # default value of auto_sample
        if self.args.auto_sample == 'auto':
            self.args.auto_sample = 'low'

        # detect if the csv_name provided exists
        if self.args.csv_name is not None and not os.path.exists(self.args.csv_name):
            self.args.log.write(f"\nx  The csv_name provided ({self.args.csv_name}) does not exist! Please specify this name correctly")
            self.args.log.finalize()
            sys.exit()

        if self.args.qdescp_solvent is not None:
            self.args.log.write(f"\nx  Currently, PTB calculations do not work with solvent! Please, remove the --qdescp_solvent option")
            self.args.log.finalize()
            sys.exit()

        if self.args.files != [] and self.args.input == '':
            # check if the input files are valid
            _ = check_files(self,'qdescp')

            # get unique files to avoid redundancy in calculations
            self.args.files = self.get_unique_files()

        # copy smarts patterns used to generate atomic descriptors
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
        """
        AQME-ROBERT workflow. 
        Combines the descriptor data from qdescp CSVs with the input CSV and saves the result.

        """

        json_pattern = str(destination) + "/boltz/*_boltz.json"

        boltz_json_files = glob.glob(json_pattern)
        dfs = [] 
        for file in boltz_json_files:
            data = pd.read_json(file, lines=True)

            name_indiv = os.path.basename(file).split("_")[:-2]
            if len(name_indiv) > 1:
                name_indiv = ['_'.join(name_indiv)]
            data.insert(loc=0, column='code_name', value=name_indiv)
            dfs.append(data)
        if dfs != []:
            df_full = pd.concat(dfs, ignore_index=True)

            name_db = 'Descriptors'
            if self.args.csv_name is not None:
                if self.args.robert:
                    name_db = 'ROBERT'

                combined_full_df = pd.DataFrame() #full

                # Read the CSV with descriptors
                input_df = pd.read_csv(self.args.csv_name)

                # check that code_name and SMILES are written with the right format
                input_df = fix_cols_names(input_df)

                if 'code_name' not in input_df.columns:
                    self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the code_name column. A combined database for AQME-{name_db} workflows will not be created.")
                elif 'SMILES' in input_df.columns:
                    # make a normalized column in qdescp_df to match input_df['code_name']
                    df_full["normalized_code_name"] = (
                        df_full["code_name"].astype(str)  # ensure all entries are strings
                        .str.replace(r"(_\d+)?_rdkit$", "", regex=True)  # remove _0_rdkit, _1_rdkit, etc.
                        .str.replace(r"_rdkit$", "", regex=True)  # remove plain _rdkit
    )

                    # to avoid problems with names that are numbers
                    input_df["code_name"] = input_df["code_name"].astype(str)

                    # now merge on the normalized column
                    combined_full_df = input_df.merge(
                        df_full.drop(columns=["code_name"]),  # drop original to avoid duplicates
                        left_on="code_name",
                        right_on="normalized_code_name",
                        how="left"
                    ).drop(columns=["normalized_code_name"])

                    csv_basename = os.path.basename(self.args.csv_name)
                    csv_path = self.args.initial_dir.joinpath(f'AQME-{name_db}_full_{csv_basename}')
                    csv_path_denovo = self.args.initial_dir.joinpath(f'AQME-{name_db}_denovo_{csv_basename}')
                    csv_path_interpret = self.args.initial_dir.joinpath(f'AQME-{name_db}_interpret_{csv_basename}')

                    # delete previous results to avoid problems
                    for csv_path_indiv in [csv_path,csv_path_denovo,csv_path_interpret]:
                        if os.path.exists(csv_path_indiv):
                            os.remove(csv_path_indiv)

                    # Save concatenated DataFrames as CSV files
                    # First, if no atoms were specified in qdescp_atoms, save the atom descriptors for all atoms
                    if len(self.args.qdescp_atoms) == 0:
                        dat_dir = Path(destination.joinpath('raw_data'))
                        dat_dir.mkdir(exist_ok=True, parents=True)

                        combined_full_df.to_csv(f'{Path(dat_dir).joinpath(os.path.basename(csv_path))}', index=None, header=True)

                        # de novo and interpret descriptors
                        descriptors_denovo = find_level_names(combined_full_df,'denovo')
                        descriptors_interpret = find_level_names(combined_full_df,'interpret')

                        combined_denovo_df = combined_full_df[descriptors_denovo]
                        combined_denovo_df.to_csv(f'{Path(dat_dir).joinpath(os.path.basename(csv_path_denovo))}', index=None, header=True)

                        combined_interpret_df = combined_full_df[descriptors_interpret]
                        combined_interpret_df.to_csv(f'{Path(dat_dir).joinpath(os.path.basename(csv_path_interpret))}', index=None, header=True)

                    # drop lists if no atomic descriptors were selected
                    clean_full_df = combined_full_df.drop(columns=[col for col in combined_full_df.columns if combined_full_df[col].apply(lambda x: isinstance(x, list)).any()])
                    clean_full_df.to_csv(f'{csv_path}', index=None, header=True)

                    # de novo and interpret descriptors
                    clean_descriptors_denovo = find_level_names(clean_full_df,'denovo')
                    clean_descriptors_interpret = find_level_names(clean_full_df,'interpret')

                    clean_denovo_df = clean_full_df[clean_descriptors_denovo]
                    clean_denovo_df.to_csv(f'{csv_path_denovo}', index=None, header=True)

                    clean_interpret_df = clean_full_df[clean_descriptors_interpret]
                    clean_interpret_df.to_csv(f'{csv_path_interpret}', index=None, header=True)

                    self.args.log.write(f"o  The AQME-{name_db}_full_{csv_basename}, AQME-{name_db}_denovo_{csv_basename} and AQME-{name_db}_interpret_{csv_basename} databases were created in {self.args.initial_dir}")
                    _ = self.process_aqme_csv(name_db)

            else:
                self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the SMILES column. A combined database for AQME-{name_db} workflows will not be created.")

        else:
            self.args.log.write(f"x  No descriptors were generated with QDESCP, please check the WARNINGS above.")



    def gather_files_and_run(self, destination, file, atom_props, smarts_targets, bar):
        """
        Run all the xTB calculations and collect the properties inside JSON files
        """ 

        xyz_files, xyz_charges, xyz_mults = [], [], []
        name = '.'.join(os.path.basename(Path(file)).split('.')[:-1])
        ext = os.path.basename(Path(file)).split(".")[-1]
        self.args.log.write(f"\n\n   ----- {name} -----")
        if ext.lower() == "xyz":
            # separate the parent XYZ file into individual XYZ files
            xyzall_2_xyz(file, name)
            # to locate difficult names (i.e. with special characters), glob.glob doesn't work, this is needed:
            xyz_files_list = [x for x in glob.glob(f"{os.path.dirname(Path(file))}/*.xyz") if os.path.basename(x).startswith(f'{name}_conf_')]

            for conf_file in xyz_files_list:
                if self.args.charge is None:
                    charge_xyz, _ = read_xyz_charge_mult(conf_file)
                else:
                    charge_xyz = self.args.charge
                if self.args.mult is None:
                    _, mult_xyz = read_xyz_charge_mult(conf_file)
                else:
                    mult_xyz = self.args.mult
                xyz_files.append(
                    os.path.dirname(os.path.abspath(file)) + "/" + conf_file
                )
                xyz_charges.append(charge_xyz)
                xyz_mults.append(mult_xyz)

        else:
            command_pdb = [
                "obabel",
                f"-i{ext.lower()}",
                file,
                "-oxyz",
                f"-O{os.path.dirname(os.path.abspath(file))}/{name}_conf_.xyz",
                "-m",
            ]
            subprocess.run(
                command_pdb,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

            # to locate difficult names (i.e. with special characters), glob.glob doesn't work, this is needed:
            xyz_files_list = [x for x in glob.glob(f"{os.path.dirname(Path(file))}/*.xyz") if os.path.basename(x).startswith(f'{name}_conf_')]

            if self.args.charge is None:
                _, charges, _, _ = mol_from_sdf_or_mol_or_mol2(file, "csearch", self.args)
            else:
                charges = [self.args.charge] * len(xyz_files_list)
            if self.args.mult is None:
                _, _, mults, _ = mol_from_sdf_or_mol_or_mol2(file, "csearch", self.args)
            else:
                mults = [self.args.mult] * len(xyz_files_list)

            for count, f in enumerate(xyz_files_list):
                xyz_files.append(f)
                xyz_charges.append(charges[count])
                xyz_mults.append(mults[count])

        for xyz_file, charge, mult in zip(xyz_files, xyz_charges, xyz_mults):
            name_xtb = '.'.join(os.path.basename(Path(xyz_file)).split(".")[:-1])
            self.args.log.write(f"\no  Running xTB and collecting properties ({name_xtb})")

            # if xTB fails during any of the calculations (UnboundLocalError), xTB assigns weird
            # qm5 charges (i.e. > +10 or < -10, ValueError), or the json file is not created 
            # for some weird xTB error (FileNotFoundError), that molecule is not used 
            xtb_passing,xtb_files_props = self.run_opt_xtb(file, xyz_file, charge, mult, name_xtb, destination)
            path_name = Path(os.path.dirname(file)).joinpath('.'.join(os.path.basename(Path(file)).split(".")[:-1]))

            if xtb_passing:
                # collect all the properties from the output files
                _ = self.morfeus_properties(path_name, atom_props, smarts_targets, xtb_files_props, charge, mult)

            _ = self.cleanup(name_xtb, destination, xtb_passing, xtb_files_props)

        bar.next()


    def run_opt_xtb(self, file, xyz_file, charge, mult, name, destination):
        """
        Runs different types of single point xTB calculations
        """

        xtb_passing = True

        dat_dir = destination / name
        dat_dir.mkdir(exist_ok=True, parents=True)

        # dictionary with the PATHs to the different xTB files
        xtb_files_props = {}

        xtb_files_props['xtb_xyz_path'] = str(dat_dir) + "/{0}.xyz".format(name)
        shutil.move(xyz_file, xtb_files_props['xtb_xyz_path'])

        xtb_input_file = str(dat_dir) + "/{0}_xtb.inp".format(name)
        with open(xtb_input_file, "wt", encoding='utf-8') as f:
            f.write("$write\n")
            f.write("json=true\n")

        xtb_files_props['xtb_opt'] = str(dat_dir) + "/{0}.out".format(name+'_opt')
        xtb_files_props['xtb_json'] = str(dat_dir) + "/{0}.json".format(name)

        os.environ["OMP_STACKSIZE"] = self.args.stacksize
        # run xTB/CREST with 1 processor
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"

        # this avoids problems when parsing charges and mults from SDF files
        mult = int(float(mult))
        charge = int(float(charge))

        # initial xTB optimization
        if self.args.xtb_opt:

            command_opt = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--opt",
                str(self.args.qdescp_opt),
                "--acc",
                str(self.args.qdescp_acc),
                "--gfn",
                str(self.args.gfn_version),
                "--chrg",
                str(charge),
                "--uhf",
                str(mult - 1),
                "--etemp",
                str(self.args.qdescp_temp),
                "-P",
                "1",
            ] # optimization
            if self.args.qdescp_solvent is not None:
                command_opt.append("--alpb")
                command_opt.append(f"{self.args.qdescp_solvent}")
            run_command(command_opt, xtb_files_props['xtb_opt'] , cwd=dat_dir)

            # replaces RDKit geometries with xTB geometries
            os.remove(xtb_files_props['xtb_xyz_path'])
            if os.path.exists(str(dat_dir) + "/xtbopt.xyz"): # finished optimizations
                os.rename(str(dat_dir) + "/xtbopt.xyz", xtb_files_props['xtb_xyz_path'])
            elif os.path.exists(str(dat_dir) + "/xtblast.xyz"): # incomplete optimizations
                os.rename(str(dat_dir) + "/xtblast.xyz", xtb_files_props['xtb_xyz_path'])
            else: # failed optimizations
                xtb_passing = False
                if file not in self.args.invalid_calcs:
                    self.args.invalid_calcs.append(file)
                self.args.log.write(f"x  WARNING! {file} did not finish correctly and no descriptors will be generated for this system. Common causes: the CHARGE and/or MULTIPLICITY used are not correct (adjust with the --charge and --mult options).")

            if xtb_passing:
                final_xyz_path = xtb_files_props['xtb_xyz_path']
                with open(final_xyz_path, "r", encoding='utf-8') as xyz_file:
                    self.xyz_coordinates = xyz_file.readlines()

        return xtb_passing,xtb_files_props


    def check_xtb_errors(self,name,file,file_check,xtb_passing):
        '''
        Check if the initial calculation finished OK
        '''
        
        with open(file_check, "r", encoding='utf-8') as opt_file:
            opt_lines = opt_file.readlines()
            for line in opt_lines:
                if '[ERROR] Program stopped' in line:
                    xtb_passing = False
                    if file not in self.args.invalid_calcs:
                        self.args.invalid_calcs.append(file)
                    self.args.log.write(f"x  WARNING! {name} did not finish correctly and no descriptors will be generated for this system. Common causes: the CHARGE and/or MULTIPLICITY used are not correct (adjust with the --charge and --mult options).")

        return xtb_passing


    def morfeus_properties(self,name_initial,atom_props,smarts_targets,xtb_files_props,charge,mult):
        """
        Collects all xTB properties from the files and adds them into a JSON file
        """

        # load initial json and add coordinates
        json_data = {}
        init_path = os.getcwd()
        morfeus_path = os.path.dirname(xtb_files_props['xtb_xyz_path'])
        os.chdir(morfeus_path)

        with open(xtb_files_props['xtb_xyz_path'], "r", encoding='utf-8') as f:
            inputs = f.readlines()
        coordinates = [inputs[i].strip().split()[1:] for i in range(2, int(inputs[0].strip()) + 2)]
        json_data["coordinates"] = coordinates

        # add MORFEUS properties to JSON
        global_properties_morfeus = calculate_morfeus_descriptors(xtb_files_props['xtb_xyz_path'],self,charge,mult,smarts_targets,name_initial)
        json_data.update(global_properties_morfeus)

        # assign atomic properties to the corresponding atoms
        json_data = self.assign_atomic_properties(json_data,name_initial,atom_props,smarts_targets)

        with open(xtb_files_props['xtb_json'], "w", encoding='utf-8') as outfile:
            json.dump(json_data, outfile)

        os.chdir(init_path)


    def assign_atomic_properties(self,json_data,name_initial,atom_props,smarts_targets):
        '''
        Assign atomic properties to the corresponding atoms
        '''

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
        """
        Removes files from the xTB calculations and place json files in the QDESCP folder
        """

        if xtb_passing: # only move molecules with successful xTB calcs
            final_json = f"{destination}/{name}.json"
            shutil.move(xtb_files_props['xtb_json'], final_json)

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
    
    def process_aqme_csv(self, name_db): 
    # Process each of the three generated CSVs.
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


    def qdescp_nmr_workflow(self,boltz_dir):
        '''
        Full NMR workflow in QDESCP for NMR prediction
        '''
        
        self.args.log.write(f"\nStarting QDESCP-{self.args.program} with {len(self.args.files)} job(s)\n")

        atom_props = ["NMR Chemical Shifts"]
        
        if os.path.basename(Path(self.args.files[0])).split('.')[-1].lower() not in ["json"]:
            self.args.log.write(f"\nx  The format used ({os.path.basename(Path(self.args.files[0])).split('.')[-1]}) is not compatible with QDESCP with NMR! Formats accepted: json")
            self.args.log.finalize()
            sys.exit()

        name = os.path.basename(Path(self.args.files[0])).split("_conf")[0]
        
        json_files = glob.glob(str(os.path.dirname(Path(self.args.files[0]))) + "/" + name + "_conf_*.json")
        get_boltz_props_nmr(json_files, name, boltz_dir, self, atom_props, self.args.nmr_atoms, self.args.nmr_slope, self.args.nmr_intercept, self.args.nmr_experim)
