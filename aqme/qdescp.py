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
from rdkit import Chem
from rdkit.Chem import rdFMCS
from pathlib import Path
from aqme.utils import (
    load_variables,
    read_xyz_charge_mult,
    mol_from_sdf_or_mol_or_mol2,
    run_command,
    check_files,
    check_dependencies,
    set_destination,
    periodic_table
)
from aqme.qdescp_utils import (
    get_boltz_props,
    read_fod,
    read_json,
    read_xtb,
    read_wbo,
    read_gfn1,
    calculate_local_CDFT_descriptors,
    calculate_global_CDFT_descriptors,
    calculate_global_morfeus_descriptors,
    calculate_local_morfeus_descriptors,
    get_descriptors,
    get_boltz_props_nmr,
    fix_cols_names,
    remove_atom_descp,
    load_file_formats,
    read_solv,
    read_triplet
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

        # check whether dependencies are installed
        _ = check_dependencies(self)

        # detects errors and updates variables before the QDESCP run
        self,destination,smarts_targets,boltz_dir = self.qdescp_set_up()

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
        valid_input = True
        if self.args.files == [] and self.args.input != '':
            if os.path.basename(self.args.input).split('.')[-1].lower() != "csv":
                self.args.log.write(f"\nx  The format used ({os.path.basename(self.args.input).split('.')[-1]}) is not compatible with the 'input' option! Formats accepted: csv")
                valid_input = False
            if self.args.input[0] == '[' or isinstance(self.args.input, list):
                self.args.log.write(f"\nx  The 'input' option was specified as a list! Please provide only the PATH or name of the CSV (i.e. --input test.csv)")
                valid_input = False
            qdescp_files = [self.args.input]
        else:
            if os.path.basename(self.args.files[0]).split('.')[-1].lower() != "sdf":
                self.args.log.write(f"\nx  The format used ({os.path.basename(self.args.files[0]).split('.')[-1]}) is not compatible with the 'files' option! Formats accepted: sdf")
                valid_input = False
            qdescp_files = self.args.files

        if not valid_input:
            self.args.log.finalize()
            sys.exit()

        self.args.log.write(f"\nStarting QDESCP-{self.args.program} with {len(qdescp_files)} job(s)\n")

        # if the files input is a CSV, first the program generates conformers
        if len(qdescp_files) == 1 and os.path.basename(qdescp_files[0]).split('.')[-1].lower() == 'csv':
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
                        '--nprocs', f'{self.args.nprocs}','--auto_sample',self.args.auto_sample]

            # overwrites charge/mult if the user specifies values
            if self.args.charge is not None:
                cmd_csearch = cmd_csearch + ['--charge', f'{self.args.charge}']
            if self.args.mult is not None:
                cmd_csearch = cmd_csearch + ['--mult', f'{self.args.mult}']

            subprocess.run(cmd_csearch)

            qdescp_files = glob.glob(f'{destination_csearch}/*.sdf')

        # obtaining mols from input files that will be used to set up atomic descriptors
        mol_list = self.get_mols_qdescp(qdescp_files)

        # obtaing SMARTS patterns from the input files automatically if no patterns are provided
        if len(smarts_targets) == 0 and len(qdescp_files) > 1:
            smarts_targets = self.auto_pattern(mol_list,smarts_targets)

        # Delete a SMARTS pattern if it is not compatible with more than 75% of the sdf files
        if len(smarts_targets) > 0:
            smarts_targets = self.remove_invalid_smarts(mol_list,smarts_targets)

        update_atom_props = [] 
        update_denovo_atom_props = [] 
        update_interpret_atom_props = [] 

        # Get descriptors (denovo, interpret, full)
        denovo_descriptors = get_descriptors('denovo')
        interpret_descriptors = get_descriptors('interpret')
        full_descriptors = get_descriptors('full')

        # Descriptors in every level
        denovo_mols = denovo_descriptors['mol']
        denovo_atoms = denovo_descriptors['atoms']

        interpret_mols = denovo_mols + interpret_descriptors['mol']
        interpret_atoms = denovo_atoms + interpret_descriptors['atoms']

        mol_props = interpret_mols + full_descriptors['mol'] 
        atom_props =  interpret_atoms + full_descriptors['atoms']

        update_atom_props, update_denovo_atom_props, update_interpret_atom_props = self.gather_files_and_run(qdescp_files, destination, atom_props, update_atom_props, smarts_targets, denovo_atoms, update_denovo_atom_props, interpret_atoms, update_interpret_atom_props)

        if len(update_atom_props) > 0:
            atom_props = update_atom_props
        if len(update_denovo_atom_props) > 0:
            denovo_atoms = update_denovo_atom_props
        if len(update_interpret_atom_props) > 0:
            interpret_atoms = update_interpret_atom_props

        #Create the CSV files
        qdescp_csv = "QDESCP_full_descriptors.csv"
        qdescp_denovo_csv = "QDESCP_denovo_descriptors.csv"
        qdescp_interpret_csv = "QDESCP_interpret_descriptors.csv"

        if self.args.boltz:
            if self.args.program.lower() == "xtb":
                self.args.log.write('\no  Running RDKit and collecting molecular properties')
                for file in qdescp_files:
                    if file not in self.args.invalid_calcs:
                        mol = Chem.SDMolSupplier(file, removeHs=False)[0]
                        name = '.'.join(os.path.basename(Path(file)).split(".")[:-1])
                        # to locate difficult names (i.e. with special characters), glob.glob doesn't work, this is needed:
                        json_files = [x for x in glob.glob(f"{destination}/*.json") if f'{name}_conf_' in x]

                        # Generating the JSON files
                        _ = get_boltz_props(json_files, name, boltz_dir, "xtb", self, mol_props, atom_props, smarts_targets,
                                            mol=mol, denovo_mols=denovo_mols, denovo_atoms=denovo_atoms, 
                                            interpret_mols=interpret_mols, interpret_atoms=interpret_atoms)
                    
                # Create the CSV files from the JSON files
                folder_raw = Path(destination).joinpath(f'raw_csv_databases')
                _ = self.write_csv_boltz_data(destination, qdescp_csv, folder_raw, atom_props, smarts_targets, json_type="standard")  # CSV full
                _ = self.write_csv_boltz_data(destination, qdescp_denovo_csv, folder_raw, atom_props, smarts_targets, json_type="denovo")  # CSV denovo
                _ = self.write_csv_boltz_data(destination, qdescp_interpret_csv, folder_raw, atom_props, smarts_targets, json_type="interpret")  # CSV interpret
                self.args.log.write(f"o  The {qdescp_denovo_csv}, {qdescp_interpret_csv} and {qdescp_csv} files containing Boltzmann weighted xTB, Morfeus and RDKit descriptors were created in {self.args.initial_dir}")

        #AQME-ROBERT workflow: Combines the descriptor data from qdescp CSVs with the input CSV and saves the result.
        _ = self.combine_and_save_csvs(qdescp_csv, qdescp_denovo_csv, qdescp_interpret_csv, folder_raw)


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

        # default value of auto_sample
        if self.args.auto_sample == 'auto':
            self.args.auto_sample = 'low'

        # detect if the csv_name provided exists
        if self.args.csv_name is not None and not os.path.exists(self.args.csv_name):
            self.args.log.write(f"\nx  The csv_name provided ({self.args.csv_name}) does not exist! Please specify this name correctly")
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
        xtb_version_dat = f'{destination}/xtb_version.dat'
        run_command(['xtb','-version'], xtb_version_dat, cwd=destination)

        xtb_version = None
        with open(xtb_version_dat, 'r') as datfile:
            lines = datfile.readlines()
            for _,line in enumerate(lines):
                if 'xtb version' in line:
                    xtb_version = line.split('version')[1].split()[0]
        os.remove(xtb_version_dat)

        if xtb_version is not None:
            self.args.log.write(f"xTB version used: {xtb_version}\n")
        else:
            self.args.log.write(f"xTB version could not be determined! Please, provide it along the results to allow other researchers to reproduce the results.\n")

        # create folder to store Boltzmann weighted properties
        boltz_dir = Path(f"{destination}/boltz")
        if os.path.exists(f"{boltz_dir}"):
            self.args.log.write(f'\nx  A previous folder of {boltz_dir} already existed, it was removed and replaced with the results of this QDESCP run.')
            shutil.rmtree(f"{boltz_dir}")
        boltz_dir.mkdir(exist_ok=True, parents=True)

        return self,destination,smarts_targets,boltz_dir
    

    def get_mols_qdescp(self,qdescp_files):
        '''
        Obtaining mols from input files
        '''
        
        mol_list = []
        for file in qdescp_files:
            with open(file, "r") as F:
                lines = F.readlines()
                smi_exist = False
                for i, line in enumerate(lines):
                    if ">  <SMILES>" in line:
                        smi = lines[i + 1].split()[0]
                        mol_indiv = Chem.AddHs(Chem.MolFromSmiles(smi))
                        mol_list.append(mol_indiv)
                        smi_exist = True
                        break
                if not smi_exist:
                    mol_indiv = Chem.SDMolSupplier(file, removeHs=False)[0]
                    mol_list.append(mol_indiv)
        
        return mol_list


    def auto_pattern(self,mol_list,smarts_targets):
        '''
        Obtaing SMARTS patterns from the input files automatically if no patterns are provided
        '''

        mcs = rdFMCS.FindMCS(mol_list)
        if mcs is not None:
            common_substructure = Chem.MolFromSmarts(mcs.smartsString)
            # Filter out non-metal atoms
            atom_types = []
            for atom in common_substructure.GetAtoms():
                atom_types.append(atom.GetSymbol())
            for atom_type in set(atom_types):
                if atom_types.count(atom_type) == 1:
                    smarts_targets.append(f'{atom_type}')
        return smarts_targets


    def remove_invalid_smarts(self,mol_list,smarts_targets):
        '''
        Delete a SMARTS pattern if it is not compatible with more than 75% of the sdf files
        '''

        patterns_remove,matches = [],[]
        for pattern in smarts_targets:
            num_matches = len(mol_list)
            for mol_indiv in mol_list:
                try:
                    # we differentiate if is a number for mapped atom or we are looking for smarts pattern in the molecule
                    if not str(pattern).isalpha() and str(pattern).isdigit():
                        for atom in mol_indiv.GetAtoms():
                            if atom.GetAtomMapNum() == int(pattern):
                                pattern_idx = int(atom.GetIdx())
                                matches = ((int(pattern_idx),),)
                    else:
                        matches = mol_indiv.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                except:
                    try: # I tried to make this except more specific for Boost.Python.ArgumentError, but apparently it's not as simple as it looks
                        matches = mol_indiv.GetSubstructMatches(Chem.MolFromSmarts(f'[{pattern}]'))
                    except:
                        if pattern in self.args.qdescp_atoms:
                            self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} was not specified correctly and its corresponding atomic descriptors will not be generated! Make sure the qdescp_atoms option uses this format: \"[C]\" for atoms, \"[C=N]\" for bonds, and so on.")
                        patterns_remove.append(pattern)
                        break
                if len(matches) != 1:
                    num_matches -= 1
                if (num_matches / len(mol_list)) < 0.75:
                    patterns_remove.append(pattern)
                    if pattern in self.args.qdescp_atoms:
                        self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} is not present (or there is more than one match) in 75% or more of the molecules. Atomic descriptors will not be generated for this pattern.")
                    break

        # remove invalid patterns
        for pattern in patterns_remove:
            smarts_targets.remove(pattern)

        # print patterns recognized automatically
        for smarts_target in smarts_targets:
            if smarts_target not in self.args.qdescp_atoms:
                self.args.log.write(f"\no  Pattern {smarts_target} shared in all input files. Using it for atomic descriptor calculations")

        return smarts_targets


    def combine_and_save_csvs(self, qdescp_csv, qdescp_denovo_csv, qdescp_interpret_csv, folder_raw):
        """
        AQME-ROBERT workflow
        Combines the descriptor data from qdescp CSVs with the input CSV and saves the result.
        """

        name_db = 'Descriptors'
        if self.args.csv_name is not None:
            if self.args.robert:
                name_db = 'ROBERT'

            combined_df = pd.DataFrame() #full
            combined_denovo_df = pd.DataFrame()  # denovo
            combined_interpret_df = pd.DataFrame()  # interpret

            # Read the CSV with descriptors
            qdescp_df = pd.read_csv(qdescp_csv)
            qdescp_denovo_df = pd.read_csv(qdescp_denovo_csv)
            qdescp_interpret_df = pd.read_csv(qdescp_interpret_csv)

            input_df = pd.read_csv(self.args.csv_name)

            # check that code_name and SMILES are written with the right format
            input_df = fix_cols_names(input_df)

            if 'code_name' not in input_df.columns:
                self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the code_name column. A combined database for AQME-{name_db} workflows will not be created.")
            elif 'SMILES' in input_df.columns:
                for i, input_name in enumerate(input_df['code_name']):
                    # concatenate with qdescp_df
                    qdescp_col = input_df.loc[i].to_frame().T.reset_index(drop=True)
                    input_col = qdescp_df.loc[(qdescp_df['code_name'] == f'{input_name}_rdkit') | 
                                            (qdescp_df['code_name'] == f'{input_name}') | 
                                            (qdescp_df['code_name'] == f'{input_name}_0_rdkit') | 
                                            (qdescp_df['code_name'] == f'{input_name}_1_rdkit') | 
                                            (qdescp_df['code_name'] == f'{input_name}_2_rdkit')]
                    input_col = input_col.drop(['code_name'], axis=1).reset_index(drop=True)
                    combined_row = pd.concat([qdescp_col, input_col], axis=1)
                    combined_df = pd.concat([combined_df, combined_row], ignore_index=True)

                    # concatenate with  qdescp_denovo_df
                    input_col_denovo = qdescp_denovo_df.loc[(qdescp_denovo_df['code_name'] == f'{input_name}_rdkit') | 
                                                            (qdescp_denovo_df['code_name'] == f'{input_name}') | 
                                                            (qdescp_denovo_df['code_name'] == f'{input_name}_0_rdkit') | 
                                                            (qdescp_denovo_df['code_name'] == f'{input_name}_1_rdkit') | 
                                                            (qdescp_denovo_df['code_name'] == f'{input_name}_2_rdkit')]
                    input_col_denovo = input_col_denovo.drop(['code_name'], axis=1).reset_index(drop=True)
                    combined_row_denovo = pd.concat([qdescp_col, input_col_denovo], axis=1)
                    combined_denovo_df = pd.concat([combined_denovo_df, combined_row_denovo], ignore_index=True)

                    # concatenate with qdescp_interpret_df
                    input_col_interpret = qdescp_interpret_df.loc[(qdescp_interpret_df['code_name'] == f'{input_name}_rdkit') | 
                                                                (qdescp_interpret_df['code_name'] == f'{input_name}') | 
                                                                (qdescp_interpret_df['code_name'] == f'{input_name}_0_rdkit') | 
                                                                (qdescp_interpret_df['code_name'] == f'{input_name}_1_rdkit') | 
                                                                (qdescp_interpret_df['code_name'] == f'{input_name}_2_rdkit')]
                    input_col_interpret = input_col_interpret.drop(['code_name'], axis=1).reset_index(drop=True)
                    combined_row_interpret = pd.concat([qdescp_col, input_col_interpret], axis=1)
                    combined_interpret_df = pd.concat([combined_interpret_df, combined_row_interpret], ignore_index=True)

                csv_basename = os.path.basename(self.args.csv_name)
                csv_path = self.args.initial_dir.joinpath(f'AQME-{name_db}_full_{csv_basename}')
                csv_path_denovo = self.args.initial_dir.joinpath(f'AQME-{name_db}_denovo_{csv_basename}')
                csv_path_interpret = self.args.initial_dir.joinpath(f'AQME-{name_db}_interpret_{csv_basename}')

                # Save concatenated DataFrames as CSV files
                combined_df.to_csv(f'{csv_path}', index=None, header=True)
                combined_denovo_df.to_csv(f'{csv_path_denovo}', index=None, header=True)
                combined_interpret_df.to_csv(f'{csv_path_interpret}', index=None, header=True)

                self.args.log.write(f"o  The AQME-{name_db}_full_{csv_basename}, AQME-{name_db}_denovo_{csv_basename} and AQME-{name_db}_interpret_{csv_basename} databases were created in {self.args.initial_dir} (the original QDESCP databases were moved to {folder_raw})")

                # move the QDESCP CSV files into the raw data folder
                for csv_file in [qdescp_csv, qdescp_denovo_csv, qdescp_interpret_csv]:
                    shutil.move(csv_file, folder_raw.joinpath(csv_file))
            else:
                self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the SMILES column. A combined database for AQME-{name_db} workflows will not be created.")

            _ = self.process_aqme_csv(name_db)

    def write_csv_boltz_data(self, destination, qdescp_csv, folder_raw, atom_props, smarts_targets, json_type="standard"):
        """
        Concatenate the values for all calculations
        """
        
        if json_type == "denovo":
            json_pattern = str(destination) + "/boltz/*_denovo_boltz.json"
        elif json_type == "interpret":
            json_pattern = str(destination) + "/boltz/*_interpret_boltz.json"
        else:
            json_pattern = str(destination) + "/boltz/*_full_boltz.json"

        boltz_json_files = glob.glob(json_pattern)
        dfs = [] 
        for file in boltz_json_files:
            data = pd.read_json(file, lines=True)

            name_indiv = os.path.basename(file).split("_")[:-2]
            if len(name_indiv) > 1:
                name_indiv = ['_'.join(name_indiv)]
            data.insert(loc=0, column='code_name', value=name_indiv)
            dfs.append(data)
        temp = pd.concat(dfs, ignore_index=True)

        # first, create raw files that will store all the information, including atomic descriptors in lists
        if not os.path.exists(folder_raw):
            folder_raw.mkdir(exist_ok=True, parents=True)
        temp.to_csv(folder_raw.joinpath(f'Raw_{qdescp_csv}'), index=False)

        # in the main folder, if there were no SMARTS matches, remove atomic descps since they're lists
        if len(smarts_targets) == 0:
            temp = remove_atom_descp(temp,atom_props)
        temp.to_csv(qdescp_csv, index=False)
    
    
    def gather_files_and_run(self, qdescp_files, destination, atom_props, update_atom_props, smarts_targets, denovo_atoms, update_denovo_atom_props, interpret_atoms, update_interpret_atom_props):
        """
        Load all the input files, execute xTB calculations, gather descriptors and clean up scratch data
        """

        # keep track of unvalid calcs
        self.args.invalid_calcs = []

        bar = IncrementalBar(
            "\no  Number of finished jobs from QDESCP", max=len(qdescp_files)
        )

        # multiprocessing to accelerate QDESCP (since xTB uses 1 processor to be reproducible)
        if not self.args.debug: # errors and try/excepts are not shown in multithreading
            with futures.ThreadPoolExecutor(
                max_workers=self.args.nprocs,
            ) as executor:
                for file in qdescp_files:
                    _ = executor.submit(
                        self.xtb_complete, destination, file, atom_props, smarts_targets, bar, update_atom_props, denovo_atoms, update_denovo_atom_props, interpret_atoms, update_interpret_atom_props
                        )
        else:
            for file in qdescp_files:
                _ = self.xtb_complete(destination, file, atom_props, smarts_targets, bar, update_atom_props, denovo_atoms, update_denovo_atom_props, interpret_atoms, update_interpret_atom_props)

        bar.finish()

        return update_atom_props, update_denovo_atom_props, update_interpret_atom_props


    def xtb_complete(self, destination, file, atom_props, smarts_targets, bar, update_atom_props, denovo_atoms, update_denovo_atom_props, interptret_atoms, update_interpret_atom_props):
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
            xyz_files_list = [x for x in glob.glob(f"{os.path.dirname(Path(file))}/*.xyz") if f'{name}_conf_' in x]

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
            xyz_files_list = [x for x in glob.glob(f"{os.path.dirname(Path(file))}/*.xyz") if f'{name}_conf_' in x]

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
            self.args.log.write(f"\no  Running xTB and collecting properties")

            # if xTB fails during any of the calculations (UnboundLocalError), xTB assigns weird
            # qm5 charges (i.e. > +10 or < -10, ValueError), or the json file is not created 
            # for some weird xTB error (FileNotFoundError), that molecule is not used 
            xtb_passing,xtb_files_props = self.run_sp_xtb(file, xyz_file, charge, mult, name_xtb, destination)
            path_name = Path(os.path.dirname(file)).joinpath('.'.join(os.path.basename(Path(file)).split(".")[:-1]))

            if xtb_passing:
                #standard
                update_atom_props = self.collect_xtb_properties(path_name, atom_props, update_atom_props, smarts_targets, xtb_files_props)

                #denovo
                update_denovo_atom_props = self.collect_xtb_properties(path_name, denovo_atoms, update_denovo_atom_props, smarts_targets, xtb_files_props)

                #interpret
                update_interpret_atom_props = self.collect_xtb_properties(path_name, interptret_atoms, update_interpret_atom_props, smarts_targets, xtb_files_props)

            _ = self.cleanup(name_xtb, destination, xtb_passing, xtb_files_props, move_folder=True)
            _ = self.merge_results(destination,xtb_files_props)
        bar.next()

        return update_atom_props, update_denovo_atom_props, update_interpret_atom_props


    def run_sp_xtb(self, file, xyz_file, charge, mult, name, destination):
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
        with open(xtb_input_file, "wt") as f:
            f.write("$write\n")
            f.write("json=true\n")

        xtb_files_props['xtb_opt'] = str(dat_dir) + "/{0}.out".format(name+'_opt')
        xtb_files_props['xtb_out'] = str(dat_dir) + "/{0}.out".format(name)
        xtb_files_props['xtb_json'] = str(dat_dir) + "/{0}.json".format(name)
        xtb_files_props['xtb_wbo'] = str(dat_dir) + "/{0}.wbo".format(name)
        xtb_files_props['xtb_gfn1'] = str(dat_dir) + "/{0}.gfn1".format(name)
        xtb_files_props['xtb_Nminus1'] = str(dat_dir) + "/{0}.Nminus1".format(name)
        xtb_files_props['xtb_Nminus2'] = str(dat_dir) + "/{0}.Nminus2".format(name)
        xtb_files_props['xtb_Nplus1'] = str(dat_dir) + "/{0}.Nplus1".format(name)
        xtb_files_props['xtb_Nplus2'] = str(dat_dir) + "/{0}.Nplus2".format(name)
        xtb_files_props['xtb_fukui'] = str(dat_dir) + "/{0}.fukui".format(name)
        xtb_files_props['xtb_fod'] = str(dat_dir) + "/{0}.fod".format(name)
        xtb_files_props['xtb_solv'] = str(dat_dir) + "/{0}.solv".format(name)
        xtb_files_props['xtb_triplet'] = str(dat_dir) + "/{0}.triplet".format(name)

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
                self.args.log.write(f"x  WARNING! {file} did not finish correctly and no descriptors will be generated for this system.")

            if xtb_passing:
                final_xyz_path = xtb_files_props['xtb_xyz_path']
                with open(final_xyz_path, "r") as xyz_file:
                    self.xyz_coordinates = xyz_file.readlines()

        if xtb_passing:
            command1 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
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
                "--input",
                str(xtb_input_file),
                "-P",
                "1",
            ] #Single point (file_N)
            if self.args.qdescp_solvent is not None:
                command1.append("--alpb")
                command1.append(f"{self.args.qdescp_solvent}")
            run_command(command1, xtb_files_props['xtb_out'], cwd=dat_dir)

            # check if the initial calculation finished OK
            with open(xtb_files_props['xtb_out'], "r") as opt_file:
                opt_lines = opt_file.readlines()
                for line in opt_lines:
                    if '[ERROR] Program stopped' in line:
                        xtb_passing = False
                        if file not in self.args.invalid_calcs:
                            self.args.invalid_calcs.append(file)
                        self.args.log.write(f"x  WARNING! {file} did not finish correctly and no descriptors will be generated for this system.")

            os.rename(str(dat_dir) + "/xtbout.json", xtb_files_props['xtb_json'])
            os.rename(str(dat_dir) + "/wbo", xtb_files_props['xtb_wbo'])
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

        if xtb_passing:
            command2 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--pop",
                "--gfn",
                "1",
                "--chrg",
                str(charge),
                "--acc",
                str(self.args.qdescp_acc),
                "--uhf",
                str(mult - 1),
                "--etemp",
                str(self.args.qdescp_temp),
                "--vomega",
                "-P",
                "1",
            ] #For cm5 charges and proportions (localgfn1)
            if self.args.qdescp_solvent is not None:
                command2.append("--alpb")
                command2.append(f"{self.args.qdescp_solvent}")
            run_command(command2, xtb_files_props['xtb_gfn1'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)


            command3 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--vfukui",
                "--gfn",
                str(self.args.gfn_version),
                "--chrg",
                str(charge),
                "--acc",
                str(self.args.qdescp_acc),
                "--uhf",
                str(mult - 1),
                "--etemp",
                str(self.args.qdescp_temp),
                "-P",
                "1",
            ] #For fukuis
            if self.args.qdescp_solvent is not None:
                command3.append("--alpb")
                command3.append(f"{self.args.qdescp_solvent}")
            run_command(command3, xtb_files_props['xtb_fukui'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

            command4 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--fod",
                "--gfn",
                str(self.args.gfn_version),
                "--chrg",
                str(charge),
                "--acc",
                str(self.args.qdescp_acc),
                "--uhf",
                str(mult - 1),
                "--etemp",
                str(self.args.qdescp_temp),
                "-P",
                "1",
            ] # for FOD
            if self.args.qdescp_solvent is not None:
                command4.append("--alpb")
                command4.append(f"{self.args.qdescp_solvent}")
            run_command(command4, xtb_files_props['xtb_fod'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

            command5 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--gfn",
                str(self.args.gfn_version),
                "--chrg",
                str(int(charge) +1),
                "--acc",
                str(self.args.qdescp_acc),
                "--uhf",
                str(mult),
                "--etemp",
                str(self.args.qdescp_temp),
                "-P",
                "1",
            ] #file_Nminus1 (N-1)
            if self.args.qdescp_solvent is not None:
                command5.append("--alpb")
                command5.append(f"{self.args.qdescp_solvent}")
            run_command(command5, xtb_files_props['xtb_Nminus1'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

            command6 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--gfn",
                str(self.args.gfn_version),
                "--chrg",
                str(int(charge) +2),
                "--acc",
                str(self.args.qdescp_acc),
                "--uhf",
                str(mult - 1),
                "--etemp",
                str(self.args.qdescp_temp),
                "-P",
                "1",
            ] #file_N-2 (N-2)
            if self.args.qdescp_solvent is not None:
                command6.append("--alpb")
                command6.append(f"{self.args.qdescp_solvent}")
            run_command(command6, xtb_files_props['xtb_Nminus2'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

            command7 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--gfn",
                str(self.args.gfn_version),
                "--chrg",
                str(int(charge) -1),
                "--acc",
                str(self.args.qdescp_acc),
                "--uhf",
                str(int(mult)),
                "--etemp",
                str(self.args.qdescp_temp),
                "-P",
                "1",
            ] #file_Nplus1 (N+1)
            if self.args.qdescp_solvent is not None:
                command7.append("--alpb")
                command7.append(f"{self.args.qdescp_solvent}")
            run_command(command7, xtb_files_props['xtb_Nplus1'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

            command8 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
                "--gfn",
                str(self.args.gfn_version),
                "--chrg",
                str(int(charge)-2),
                "--acc",
                str(self.args.qdescp_acc),
                "--uhf",
                str(int(mult)-1),
                "--etemp",
                str(self.args.qdescp_temp),
                "-P",
                "1",
            ] #file_Nplus2 (N+2)
            if self.args.qdescp_solvent is not None:
                command8.append("--alpb")
                command8.append(f"{self.args.qdescp_solvent}")
            run_command(command8, xtb_files_props['xtb_Nplus2'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

            with open(xtb_input_file, "wt") as f:
                f.write("$write\n")
                f.write('gbsa=true\n')

            command9 = [
                "xtb",
                xtb_files_props['xtb_xyz_path'],
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
                "--input",
                str(xtb_input_file),
                "-P",
                "1",
                "--alpb",
                "h2o"
            ] # file solvation
            run_command(command9, xtb_files_props['xtb_solv'], cwd=dat_dir)
            _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

            if int(mult) == 1:
                command10 = [
                    "xtb",
                    xtb_files_props['xtb_xyz_path'],
                    "--acc",
                    str(self.args.qdescp_acc),
                    "--gfn",
                    str(self.args.gfn_version),
                    "--chrg",
                    str(charge),
                    "--uhf",
                    '2',
                    "--etemp",
                    str(self.args.qdescp_temp),
                    "--input",
                    str(xtb_input_file),
                    "-P",
                    "1",
                ] # file triplet
                run_command(command10, xtb_files_props['xtb_triplet'], cwd=dat_dir)
                _ = self.cleanup(name, destination, xtb_passing, xtb_files_props)

        return xtb_passing,xtb_files_props


    def collect_xtb_properties(self,name_initial,atom_props,update_atom_props,smarts_targets,xtb_files_props):
        """
        Collects all xTB properties from the files and puts them in a JSON file
        """
        properties_dict = read_xtb(xtb_files_props['xtb_out'],self)
        localgfn1 = read_gfn1(xtb_files_props['xtb_gfn1'],self)
        properties_FOD = read_fod(xtb_files_props['xtb_fod'],self)
        bonds, wbos = read_wbo(xtb_files_props['xtb_wbo'],self)
        properties_solv = read_solv(xtb_files_props['xtb_solv'])
        properties_triplet = read_triplet(xtb_files_props['xtb_triplet'],properties_dict['Total energy'])
        cdft_descriptors  = calculate_global_CDFT_descriptors(xtb_files_props['xtb_out'], xtb_files_props['xtb_Nminus1'], xtb_files_props['xtb_Nminus2'], xtb_files_props['xtb_Nplus1'], xtb_files_props['xtb_Nplus2'],self)
        localDescriptors = calculate_local_CDFT_descriptors(xtb_files_props['xtb_fukui'], cdft_descriptors,self)
        # create matrix of Wiberg bond-orders
        atoms = properties_dict.get("atoms")
        nat = len(atoms)
        wbo_matrix = np.zeros((nat, nat))
        for i, bond in enumerate(bonds):
            wbo_matrix[(bond[0] - 1)][(bond[1] - 1)] = wbos[i]
            wbo_matrix[(bond[1] - 1)][(bond[0] - 1)] = wbos[i]

        """
		Now add xTB descriptors to existing json files.
		"""
        json_data = read_json(xtb_files_props['xtb_json']) 
        json_data["Wiberg matrix"] = wbo_matrix.tolist()
        list_properties = [properties_dict,properties_FOD,properties_solv,localgfn1,cdft_descriptors,localDescriptors,properties_triplet]
        for properties in list_properties:
            if properties is not None:
                json_data.update(properties)


        """
		Morfeus Descriptors
		"""
        #Global descriptors
        global_properties_morfeus = calculate_global_morfeus_descriptors(xtb_files_props['xtb_xyz_path'],self)
        json_data.update(global_properties_morfeus)
        #Local descriptors
        local_properties_morfeus = calculate_local_morfeus_descriptors(xtb_files_props['xtb_xyz_path'],self)
        json_data.update(local_properties_morfeus)


        """
        Locate SMARTS, qdescp_atoms
        """

        with open(xtb_files_props['xtb_xyz_path'], "r") as f:
            inputs = f.readlines()

        coordinates = [inputs[i].strip().split()[1:] for i in range(2, int(inputs[0].strip()) + 2)]
        json_data["coordinates"] = coordinates

        if len(smarts_targets) > 0:
            # detect SMILES from SDF files generated by CSEARCH or create mol objects from regular SDF files
            sdf_file = f'{name_initial}.sdf'
            with open(sdf_file, "r") as F:
                lines = F.readlines()

            smi_exist = False
            for i, line in enumerate(lines):
                if ">  <SMILES>" in line:
                    smi = lines[i + 1].split()[0]
                    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
                    smi_exist = True
            if not smi_exist:
                mol = Chem.SDMolSupplier(sdf_file, removeHs=False)[0]

            # find the target atoms or groups
            for pattern in smarts_targets:
                matches = []
                idx_set = None
                
                # we differentiate if is a number for mapped atom or we are looking for smarts pattern in the molecule
                if not str(pattern).isalpha() and str(pattern).isdigit():
                    for atom in mol.GetAtoms():
                        if atom.GetAtomMapNum() == int(pattern):
                            idx_set = pattern
                            pattern_idx = int(atom.GetIdx())
                            matches = ((int(pattern_idx),),)
                else: 
                    try:
                        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                    except:
                        try:
                            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(f'[{pattern}]'))
                        except:
                            self.args.log.write(f"x  WARNING! SMARTS pattern was not specified correctly! Make sure the qdescp_atoms option uses this format: \"[C]\" for atoms, \"[C=N]\" for bonds, and so on.")

                if len(matches) == 0:
                    self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} not found in {os.path.basename(name_initial)}, atomic descriptors will not be generated.")
                elif matches[0] == -1:
                    self.args.log.write(f"x  WARNING! Mapped atom {pattern} not found in {os.path.basename(name_initial)}, atomic descriptors will not be generated.")
                elif len(matches) > 1:
                    self.args.log.write(f"x  WARNING! More than one {pattern} patterns were found in {os.path.basename(name_initial)}, atomic descriptors will not be generated.")
                elif len(matches) == 1:
                    # get atom types and sort them to keep the same atom order among different molecules
                    atom_indices = list(matches[0])
                    atom_types = []
                    for atom_idx in atom_indices:
                        atom_types.append(mol.GetAtoms()[atom_idx].GetSymbol())

                    n_types = len(set(atom_types))
                    if n_types == 1:
                        sorted_indices = sorted(atom_indices, key=lambda idx: len(mol.GetAtoms()[idx].GetNeighbors()))
                    elif n_types > 1:
                        sorted_indices = sorted(atom_indices, key=lambda idx: mol.GetAtoms()[idx].GetAtomicNum())

                    # Generate unique match names for each atom type in the functional group
                    match_names = []
                    atom_counters = {}

                    # Separate atoms when functional groups are used
                    for atom_idx in sorted_indices:
                        atom_type = mol.GetAtoms()[atom_idx].GetSymbol()

                        # Check if there's more than one atom of the same type
                        n_atoms_of_type = sum(1 for idx in sorted_indices if mol.GetAtoms()[idx].GetSymbol() == atom_type)

                        # Initialize counter if it doesn't exist for this atom type
                        if atom_type not in atom_counters:
                            atom_counters[atom_type] = 1
                        # If the pattern is a single atom number, we include that number in the match name
                        if str(pattern).isdigit():  # This means the pattern is an atom number
                            match_name = f'{atom_type}{pattern}'
                        else:
                            # If it's a SMARTS pattern or more than one atom
                            if len(smarts_targets) == 1 and smarts_targets[0] in periodic_table():
                                # Case where it's just an atom type without SMARTS
                                if n_atoms_of_type == 1:
                                    match_name = f'{atom_type}'
                                else:
                                    match_name = f'{atom_type}_{atom_counters[atom_type]}'
                            else:
                                if pattern[0] == '#':
                                    match_name = f'{atom_type}'
                                # Regular SMARTS pattern handling
                                elif n_atoms_of_type == 1:
                                    if idx_set is None:
                                        match_name = f'{pattern}_{atom_type}'
                                    else:
                                        match_name = f'{pattern}_{atom_type}{idx_set}'
                                else:
                                    if idx_set is None:
                                        match_name = f'{pattern}_{atom_type}_{atom_counters[atom_type]}'
                                    else:
                                        match_name = f'{pattern}_{atom_type}{idx_set}_{atom_counters[atom_type]}'
                        
                        atom_counters[atom_type] += 1  # Increment counter for that atom type

                        match_names.append(match_name)

                    # Assign atomic descriptors to each identified atom
                    for atom_idx, match_name in zip(sorted_indices, match_names):
                        idx_xtb = atom_idx
                        for prop in atom_props:
                            try:
                                json_data[f'{match_name}_{prop}'] = json_data[prop][idx_xtb]
                                if f'{match_name}_{prop}' not in update_atom_props:
                                    update_atom_props.append(f'{match_name}_{prop}')
                            except (KeyError,IndexError): # prevents missing values
                                pass

                    # Adding max and min values for functional groups with the same two atoms
                    if len(match_names) > 1 and n_types == 1:
                        for prop in atom_props:
                            prop_values = []
                            for prop_name in match_names:
                                prop_values.append(json_data[f'{prop_name}_{prop}'])
                            json_data[f'{pattern}_max_{prop}'] = max(prop_values)
                            json_data[f'{pattern}_min_{prop}'] = min(prop_values)
                            if f'{pattern}_max_{prop}' not in update_atom_props:
                                update_atom_props.append(f'{pattern}_max_{prop}')
                            if f'{pattern}_min_{prop}' not in update_atom_props:
                                update_atom_props.append(f'{pattern}_min_{prop}')

        with open(xtb_files_props['xtb_json'], "w") as outfile:
            json.dump(json_data, outfile)
        
        return update_atom_props    

    def cleanup(self, name, destination, xtb_passing, xtb_files_props, move_folder=False):
        """
        Removes files from the xTB calculations that are not relevant and place json files in the 
        QDESCP folder
        """

        if xtb_passing and move_folder: # only move molecules with successful xTB calcs
            final_json = f"{destination}/{name}.json"
            shutil.move(xtb_files_props['xtb_json'], final_json)
        
        # delete xTB files that does not contain useful data
        files = glob.glob(f"{destination}/{name}/*")
        for file in files:
            if name not in os.path.basename(file):
                os.remove(file)
            if '.inp' in os.path.basename(file):
                if move_folder:
                    os.remove(file)
        if os.path.exists(f"{destination}/{name}/.xtboptok"):
            os.remove(f"{destination}/{name}/.xtboptok")

        if move_folder:
            if not os.path.exists(f"{destination}/xtb_data"): 
                Path(f"{destination}/xtb_data").mkdir()
            if os.path.exists(f"{destination}/xtb_data/{name}"): 
                self.args.log.write(f'\nx  A previous folder of {name} already existed, it was removed and replaced with the results of this QDESCP run.')
                shutil.rmtree(f"{destination}/xtb_data/{name}")
            shutil.move(f"{destination}/{name}", f"{destination}/xtb_data/{name}")


    def merge_results(self,destination,xtb_files_props):
        # combine all the results in one single file
        file_formats = load_file_formats()

        name = '.'.join(os.path.basename(xtb_files_props['xtb_out']).split('.')[:-1])
        raw_content = ''
        for file_format in file_formats:
            xtb_file = f"{destination}/xtb_data/{name}/{name}{file_format}"
            if os.path.exists(xtb_file):
                if raw_content != '':
                    raw_content += '\n\n'
                raw_content += f'----- {file_formats[file_format]} -----\n'
                with open(xtb_file, 'r', encoding='utf-8') as file:
                    lines = file.readlines()
                    for line in lines:
                        raw_content += line
                os.remove(xtb_file)

        raw_file = open(f"{destination}/xtb_data/{name}/{name}_All_Calcs.out", "w")
        raw_file.write(f"{raw_content}\n")
        raw_file.close()


    def get_unique_files(self):
        unique_files = []
        unique_smiles = []
        for file in self.args.files:
            with open(file, "r") as F:
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
