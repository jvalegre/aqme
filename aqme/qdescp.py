"""
Parameters
----------

General
+++++++

   w_dir_main : str, default=os.getcwd()
      Working directory
   destination : str, default=None,
      Directory to create the JSON file(s)
   program : str, default=None
      Program required to create the new descriptors. Current options: 'xtb', 'nmr'
   qdescp_atoms : list of str, default=[]
      Type of atom or group to calculate atomic properties. This option admits atoms 
      (i.e., qdescp_atoms=['P']) and SMART patterns (i.e., qdescp_atoms=['C=O']) 
   robert : bool, default=True
      Creates a database ready to use in an AQME-ROBERT machine learning workflow,
      combining the input CSV with SMILES/code_name and the calculated xTB/DBSTEP descriptors

xTB descriptors
+++++++++++++++

   files : list of str, default=''
      Filenames of SDF/PDB/XYZ files to calculate xTB descriptors. If \*.sdf 
      (or other strings that are not lists such as \*.pdb) are specified, 
      the program will look for all the SDF files in the working directory 
      through glob.glob(\*.sdf)
   charge : int, default=None
      Charge of the calculations used in the following input files (charges from
      SDF files generated in CSEARCH are read automatically).
   mult : int, default=None
      Multiplicity of the calculations used in the following input files 
      (multiplicities from SDF files generated in CSEARCH are read automatically).
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
    check_files
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
    calculate_global_CDFT_descriptors_part2,
    calculate_global_morfeus_descriptors,
    calculate_local_morfeus_descriptors,
    get_descriptors
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

        # detect errors and incompatibilities before the QDESCP run
        if self.args.program.lower() not in ["xtb", "nmr"]:
            self.args.log.write("\nx  Program not supported for QDESCP descriptor generation! Specify: program='xtb' (or nmr)")
            self.args.log.finalize()
            sys.exit()

        # get unique files to avoid redundancy in calculations
        self.args.files = self.get_unique_files()

        if self.args.destination is None:
            destination = self.args.initial_dir.joinpath("QDESCP")
        else:
            destination = Path(self.args.destination)

        # retrieve the different files to run in QDESCP
        _ = check_files(self,'qdescp')

        self.args.log.write(f"\nStarting QDESCP-{self.args.program} with {len(self.args.files)} job(s)\n")

        # obtain SMARTS patterns from the input files automatically if no patterns are provided
        smarts_targets = self.auto_pattern_recog()

        # delete SMARTS patterns that are not present in more than 30% of the molecules
        smarts_targets = self.remove_patterns_low(smarts_targets)

        """
        Reccolecting descriptors from XTB and Morfeus
        """
        update_atom_props = [] 
        update_denovo_atom_props = [] 
        update_interpret_atom_props = [] 

        if self.args.program.lower() == "xtb":
            # Get descriptors (denovo, interpret, full)
            denovo_descriptors = get_descriptors('denovo')
            interpret_descriptors = get_descriptors('interpret')
            full_descriptors = get_descriptors('full')

            # Descriptors in every level
            denovo_mols = denovo_descriptors['mol']
            denovo_atoms = denovo_descriptors['atoms']

            interpret_mols = interpret_descriptors['mol']
            interpret_atoms = interpret_descriptors['atoms']

            mol_props = full_descriptors['mol'] + interpret_mols
            atom_props = full_descriptors['atoms'] + interpret_atoms

            update_atom_props, update_denovo_atom_props, update_interpret_atom_props = self.gather_files_and_run(
                destination, atom_props, update_atom_props, smarts_targets, 
                denovo_atoms, update_denovo_atom_props, interpret_atoms, update_interpret_atom_props
    )
            
    #Update the list of descriptors
        if len(update_atom_props) > 0:
            atom_props = update_atom_props
        if len(update_denovo_atom_props) > 0:
            denovo_atoms = update_denovo_atom_props
        if len(update_interpret_atom_props) > 0:
            interpret_atoms = update_interpret_atom_props
        
        # Function to create CSV files and Boltzmann averaged properties.
        qdescp_csv, qdescp_denovo_csv, qdescp_interpret_csv = self.create_csv_files(destination, mol_props, atom_props, smarts_targets, denovo_mols, denovo_atoms, interpret_mols, interpret_atoms)

        #AQME-ROBERT workflow: Combines the descriptor data from qdescp CSVs with the input CSV and saves the result.
        self.combine_and_save_csvs(qdescp_csv, qdescp_denovo_csv, qdescp_interpret_csv)

        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime QDESCP: {elapsed_time} seconds\n")
        self.args.log.finalize()

    def auto_pattern_recog(self):
        '''
        Recognizes automatically patterns that are repeated in the molecule database.
        '''

        smarts_targets = self.args.qdescp_atoms.copy()

        if self.args.csv_name is not None:
            input_df = pd.read_csv(self.args.csv_name)
            if len(smarts_targets) == 0:
                smarts_targets = []
                if 'SMILES' in input_df.columns or 'smiles' in input_df.columns or any(col.startswith('smiles_') for col in input_df.columns):
                    possible_smiles_columns = [col for col in input_df.columns if col.lower().startswith('smiles')]
                    if len(possible_smiles_columns) == 0:
                        self.args.log.write("x  WARNING! No column with SMILES information found in the input CSV file.")
                    else:
                        for col in possible_smiles_columns:
                            if col in input_df.columns:
                                smiles_column = col
                                break
                        smiles_list = input_df[smiles_column].tolist()
                else:
                    self.args.log.write("x  WARNING! No column with SMILES information found in the input CSV file.")
                if len(smiles_list) > 0:
                    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]
                    if len(mols) > 0:
                        mcs = rdFMCS.FindMCS(mols)
                        if mcs is not None:
                            common_substructure = Chem.MolFromSmarts(mcs.smartsString)
                            # Filter out non-metal atoms
                            metal_smarts = []
                            for atom in common_substructure.GetAtoms():
                                if atom.GetSymbol() in ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo',
                                                        'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                                                        'Hg', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']:
                                    metal_smarts.append(f'[{atom.GetSymbol()}]')
                            common_substructure = Chem.MolToSmiles(Chem.MolFromSmarts('.'.join(metal_smarts)))
                            if common_substructure is not None and common_substructure != '':
                                smarts_targets.append(common_substructure)
                                self.args.log.write(f"\nSubstructure {(common_substructure)} found in input files. Using it for atomic descriptor calculations.")
        return smarts_targets
    
    def remove_patterns_low(self,smarts_targets):
        '''
        Delete SMARTS patterns that are not present in more than 30% of the molecules.
        '''

        if len(smarts_targets) > 0:
            mol_list = []
            for file in self.args.files:
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
                        mol_indiv = Chem.SDMolSupplier(file, removeHs=False)
                        mol_list.append(mol_indiv)

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
                            self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} was not specified correctly and its corresponding atomic descriptors will not be generated! Make sure the qdescp_atoms option uses this format: \"[C]\" for atoms, \"[C=N]\" for bonds, and so on.")
                            patterns_remove.append(pattern)
                            break
                    if len(matches) != 1:
                        num_matches -= 1
                    if (num_matches / len(mol_list)) < 0.7:
                        patterns_remove.append(pattern)
                        self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} is not specified correctly or is not present in more than 30% of the molecules. Atomic descriptors will not be generated for this pattern.")
                        break

            # remove invalid patterns
            for pattern in patterns_remove:
                smarts_targets.remove(pattern)

            return smarts_targets
    
    def create_csv_files(self, destination, mol_props, atom_props, smarts_targets, denovo_mols, denovo_atoms, interpret_mols, interpret_atoms):
        """
        Function to create CSV files and Boltzmann averaged properties.
        """

        qdescp_csv = "QDESCP_full_boltz_descriptors.csv"
        qdescp_denovo_csv = "QDESCP_denovo_boltz_descriptors.csv"
        qdescp_interpret_csv = "QDESCP_interpret_boltz_descriptors.csv"

        boltz_dir = Path(f"{destination}/boltz")
        if os.path.exists(f"{boltz_dir}"):
            self.args.log.write(f'\nx  A previous folder of {boltz_dir} already existed, it was removed and replaced with the results of this QDESCP run.')
            shutil.rmtree(f"{boltz_dir}")
        boltz_dir.mkdir(exist_ok=True, parents=True)

        if self.args.boltz:
            if self.args.program.lower() == "xtb":
                self.args.log.write('\no  Running RDKit and collecting molecular properties')
                for file in self.args.files:
                    mol = Chem.SDMolSupplier(file, removeHs=False)[0]
                    name = os.path.basename(Path(file)).split(".")[0]
                    json_files = glob.glob(str(destination) + "/" + name + "_conf_*.json")

                    # Generating the JSON files
                    _ = get_boltz_props(json_files, name, boltz_dir, "xtb", self, mol_props, atom_props, smarts_targets,
                                        mol=mol, denovo_mols=denovo_mols, denovo_atoms=denovo_atoms, 
                                        interpret_mols=interpret_mols, interpret_atoms=interpret_atoms)

                # Create the CSV from the JSON files
                _ = self.write_csv_boltz_data(destination, qdescp_csv, json_type="standard")  # CSV full
                _ = self.write_csv_boltz_data(destination, qdescp_denovo_csv, json_type="denovo")  # CSV denovo
                _ = self.write_csv_boltz_data(destination, qdescp_interpret_csv, json_type="interpret")  # CSV interpret

            elif self.args.program.lower() == "nmr":
                mol_props = None
                atom_props = ["NMR Chemical Shifts"]
                if os.path.basename(Path(self.args.files[0])).split('.')[1].lower() not in ["json"]:
                    self.args.log.write(f"\nx  The format used ({os.path.basename(Path(self.args.files[0])).split('.')[1]}) is not compatible with QDESCP with NMR! Formats accepted: json")
                    self.args.log.finalize()
                    sys.exit()

                name = os.path.basename(Path(self.args.files[0])).split("_conf")[0]
                json_files = glob.glob(str(os.path.dirname(os.path.abspath(self.args.files[0]))) + "/" + name + "_conf_*.json")
                get_boltz_props(json_files, name, boltz_dir, "nmr", self, mol_props, atom_props, smarts_targets, 
                                self.args.nmr_atoms, self.args.nmr_slope, self.args.nmr_intercept, self.args.nmr_experim)

        return qdescp_csv, qdescp_denovo_csv, qdescp_interpret_csv
    
    def combine_and_save_csvs(self, qdescp_csv, qdescp_denovo_csv, qdescp_interpret_csv):
        """
        AQME-ROBERT workflow
        Combines the descriptor data from qdescp CSVs with the input CSV and saves the result.
        """
        name_db = 'Descriptors'
        if self.args.program.lower() == "xtb" and self.args.csv_name is not None:
            if self.args.robert:
                name_db = 'ROBERT'

            combined_df = pd.DataFrame()
            combined_denovo_df = pd.DataFrame()  # denovo
            combined_interpret_df = pd.DataFrame()  # interpret

            # Read the CSV with descriptors
            qdescp_df = pd.read_csv(qdescp_csv)
            qdescp_denovo_df = pd.read_csv(qdescp_denovo_csv)
            qdescp_interpret_df = pd.read_csv(qdescp_interpret_csv)

            input_df = pd.read_csv(self.args.csv_name)

            if 'code_name' not in input_df.columns:
                self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the code_name column. A combined database for AQME-{name_db} workflows will not be created.")
            elif 'SMILES' in input_df.columns or 'smiles' in input_df.columns or 'Smiles' in input_df.columns:
                path_json = os.path.dirname(Path(qdescp_df['Name'][0]))
                for i, input_name in enumerate(input_df['code_name']):
                    # concatenate with qdescp_df
                    qdescp_col = input_df.loc[i].to_frame().T.reset_index(drop=True)
                    input_col = qdescp_df.loc[(qdescp_df['Name'] == f'{path_json}/{input_name}_rdkit_boltz') | 
                                            (qdescp_df['Name'] == f'{path_json}/{input_name}_boltz') | 
                                            (qdescp_df['Name'] == f'{path_json}/{input_name}_0_rdkit_boltz') | 
                                            (qdescp_df['Name'] == f'{path_json}/{input_name}_1_rdkit_boltz') | 
                                            (qdescp_df['Name'] == f'{path_json}/{input_name}_2_rdkit_boltz')]
                    input_col = input_col.drop(['Name'], axis=1).reset_index(drop=True)
                    combined_row = pd.concat([qdescp_col, input_col], axis=1)
                    combined_df = pd.concat([combined_df, combined_row], ignore_index=True)

                    # concatenate with  qdescp_denovo_df
                    input_col_denovo = qdescp_denovo_df.loc[(qdescp_denovo_df['Name'] == f'{path_json}/{input_name}_rdkit_boltz') | 
                                                            (qdescp_denovo_df['Name'] == f'{path_json}/{input_name}_boltz') | 
                                                            (qdescp_denovo_df['Name'] == f'{path_json}/{input_name}_0_rdkit_boltz') | 
                                                            (qdescp_denovo_df['Name'] == f'{path_json}/{input_name}_1_rdkit_boltz') | 
                                                            (qdescp_denovo_df['Name'] == f'{path_json}/{input_name}_2_rdkit_boltz')]
                    input_col_denovo = input_col_denovo.drop(['Name'], axis=1).reset_index(drop=True)
                    combined_row_denovo = pd.concat([qdescp_col, input_col_denovo], axis=1)
                    combined_denovo_df = pd.concat([combined_denovo_df, combined_row_denovo], ignore_index=True)

                    # concatenate with  qdescp_interpret_df
                    input_col_interpret = qdescp_interpret_df.loc[(qdescp_interpret_df['Name'] == f'{path_json}/{input_name}_rdkit_boltz') | 
                                                                (qdescp_interpret_df['Name'] == f'{path_json}/{input_name}_boltz') | 
                                                                (qdescp_interpret_df['Name'] == f'{path_json}/{input_name}_0_rdkit_boltz') | 
                                                                (qdescp_interpret_df['Name'] == f'{path_json}/{input_name}_1_rdkit_boltz') | 
                                                                (qdescp_interpret_df['Name'] == f'{path_json}/{input_name}_2_rdkit_boltz')]
                    input_col_interpret = input_col_interpret.drop(['Name'], axis=1).reset_index(drop=True)
                    combined_row_interpret = pd.concat([qdescp_col, input_col_interpret], axis=1)
                    combined_interpret_df = pd.concat([combined_interpret_df, combined_row_interpret], ignore_index=True)

                # clean and save DataFrames
                combined_df = combined_df.dropna(axis=0)
                combined_denovo_df = combined_denovo_df.dropna(axis=0)
                combined_interpret_df = combined_interpret_df.dropna(axis=0)

                csv_basename = os.path.basename(self.args.csv_name)
                csv_path = self.args.initial_dir.joinpath(f'AQME-{name_db}_full_{csv_basename}')
                csv_path_denovo = self.args.initial_dir.joinpath(f'AQME-{name_db}_denovo_{csv_basename}')
                csv_path_interpret = self.args.initial_dir.joinpath(f'AQME-{name_db}_interpret_{csv_basename}')

                # Save concatenated DataFrames as CSV files
                combined_df.to_csv(f'{csv_path}', index=None, header=True)
                combined_denovo_df.to_csv(f'{csv_path_denovo}', index=None, header=True)
                combined_interpret_df.to_csv(f'{csv_path_interpret}', index=None, header=True)

                self.args.log.write(f"o  The AQME-{name_db}_full_{csv_basename} file containing the database ready for the AQME-{name_db} workflow was successfully created in {self.args.initial_dir}")
                self.args.log.write(f"o  The AQME-{name_db}_denovo_{csv_basename} file containing the database ready for the AQME-{name_db} workflow was successfully created in {self.args.initial_dir}")
                self.args.log.write(f"o  The AQME-{name_db}_interpret_{csv_basename} file containing the database ready for the AQME-{name_db} workflow was successfully created in {self.args.initial_dir}")
            else:
                self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the SMILES column. A combined database for AQME-{name_db} workflows will not be created.")

            _ = self.process_aqme_csv(name_db)
                
    
    def write_csv_boltz_data(self, destination, qdescp_csv, json_type="standard"):
        """
        Concatenate the values for all calculations
        """
        
        # JSON files
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

            data["Name"] = file.split(".json")[0]
            dfs.append(data)

        if len(dfs) > 0:
            temp = pd.concat(dfs, ignore_index=True) 
            temp.to_csv(qdescp_csv, index=False)
            self.args.log.write(f"o  The {qdescp_csv} file containing Boltzmann weighted xTB, RDKit and Morfeus descriptors was successfully created in {self.args.initial_dir}")
        else:
            self.args.log.write(f"x  No CSV file containing Boltzmann weighted descriptors was created. This might happen when using the qdescp_atoms option with an atom/group that is not found in any of the calculations")
    
    def gather_files_and_run(self, destination, atom_props, update_atom_props, smarts_targets, denovo_atoms, update_denovo_atom_props, interpret_atoms, update_interpret_atom_props):
        """
        Load all the input files, execute xTB calculations, gather descriptors and clean up scratch data
        """

        # write input files
        if os.path.basename(Path(self.args.files[0])).split('.')[1].lower() not in ["sdf", "xyz", "pdb"]:
            self.args.log.write(f"\nx  The format used ({os.path.basename(Path(self.args.files[0])).split('.')[1]}) is not compatible with QDESCP with xTB! Formats accepted: sdf, xyz, pdb")
            self.args.log.finalize()
            sys.exit()

        for file in self.args.files:
            xyz_files, xyz_charges, xyz_mults = [], [], []
            name = os.path.basename(Path(file)).split('.')[0]
            ext = os.path.basename(Path(file)).split(".")[1]
            self.args.log.write(f"\n\n   ----- {name} -----")
            if ext.lower() == "xyz":
                # separate the parent XYZ file into individual XYZ files
                xyzall_2_xyz(file, name)
                for conf_file in glob.glob(f"{name}_conf_*.xyz"):
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

                if self.args.charge is None:
                    _, charges, _, _ = mol_from_sdf_or_mol_or_mol2(file, "csearch", self.args)
                else:
                    charges = [self.args.charge] * len(
                        glob.glob(
                            f"{os.path.dirname(os.path.abspath(file))}/{name}_conf_*.xyz"
                        )
                    )
                if self.args.mult is None:
                    _, _, mults, _ = mol_from_sdf_or_mol_or_mol2(file, "csearch", self.args)
                else:
                    mults = [self.args.mult] * len(
                        glob.glob(
                            f"{os.path.dirname(os.path.abspath(file))}/{name}_conf_*.xyz"
                        )
                    )

                for count, f in enumerate(
                    glob.glob(
                        f"{os.path.dirname(os.path.abspath(file))}/{name}_conf_*.xyz"
                    )
                ):
                    xyz_files.append(f)
                    xyz_charges.append(charges[count])
                    xyz_mults.append(mults[count])

            bar = IncrementalBar(
                "\no  Number of finished jobs from QDESCP", max=len(xyz_files)
            )
            # multiprocessing to accelerate QDESCP (since xTB uses 1 processor to be reproducible)
            with futures.ThreadPoolExecutor(
                max_workers=self.args.nprocs,
            ) as executor:
                for xyz_file, charge, mult in zip(xyz_files, xyz_charges, xyz_mults):
                    _ = executor.submit(
                        self.xtb_complete, xyz_file, charge, mult, destination, file, atom_props, smarts_targets, bar, update_atom_props, denovo_atoms, update_denovo_atom_props, interpret_atoms, update_interpret_atom_props
                    )
            bar.finish()

        return update_atom_props, update_denovo_atom_props, update_interpret_atom_props


    def xtb_complete(self, xyz_file, charge, mult, destination, file, atom_props, smarts_targets, bar, update_atom_props, denovo_atoms, update_denovo_atom_props, interptret_atoms, update_interpret_atom_props):
        """
        Run all the xTB calculations and collect the properties inside JSON files
        """ 
        
        name_xtb = os.path.basename(Path(xyz_file)).split(".")[0]
        self.args.log.write(f"\no  Running xTB and collecting properties")

        # if xTB fails during any of the calculations (UnboundLocalError), xTB assigns weird
        # qm5 charges (i.e. > +10 or < -10, ValueError), or the json file is not created 
        # for some weird xTB error (FileNotFoundError), that molecule is not used 
        xtb_passing = True
        try:
            xtb_files_props = self.run_sp_xtb(xyz_file, charge, mult, name_xtb, destination)
            path_name = Path(os.path.dirname(file)).joinpath(os.path.basename(Path(file)).split(".")[0])

            #standar
            update_atom_props = self.collect_xtb_properties(path_name, atom_props, update_atom_props, smarts_targets, xtb_files_props)

            #denovo
            update_denovo_atom_props = self.collect_xtb_properties(path_name, denovo_atoms, update_denovo_atom_props, smarts_targets, xtb_files_props)

            #interpret
            update_interpret_atom_props = self.collect_xtb_properties(path_name, interptret_atoms, update_interpret_atom_props, smarts_targets, xtb_files_props)

        except Exception as e:
            xtb_passing = False

        self.cleanup(name_xtb, destination, xtb_passing, xtb_files_props)
        bar.next()

        return update_atom_props, update_denovo_atom_props, update_interpret_atom_props


    def run_sp_xtb(self, xyz_file, charge, mult, name, destination):
        """
        Runs different types of single point xTB calculations
        """

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
        xtb_files_props['xtb_gfn1'] = str(dat_dir) + "/{0}.gfn1".format(name)
        xtb_files_props['xtb_Nminus1'] = str(dat_dir) + "/{0}.Nminus1".format(name)
        xtb_files_props['xtb_Nminus2'] = str(dat_dir) + "/{0}.Nminus2".format(name)
        xtb_files_props['xtb_Nplus1'] = str(dat_dir) + "/{0}.Nplus1".format(name)
        xtb_files_props['xtb_Nplus2'] = str(dat_dir) + "/{0}.Nplus2".format(name)
        xtb_files_props['xtb_fukui'] = str(dat_dir) + "/{0}.fukui".format(name)
        xtb_files_props['xtb_fod'] = str(dat_dir) + "/{0}.fod".format(name)

        os.environ["OMP_STACKSIZE"] = self.args.stacksize
        # run xTB/CREST with 1 processor
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"

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
                "2",
                "--chrg",
                str(charge),
                "--uhf",
                str(int(mult) - 1),
                "-P",
                "1",
            ]
            if self.args.qdescp_solvent is not None:
                command_opt.append("--alpb")
                command_opt.append(f"{self.args.qdescp_solvent}")
            run_command(command_opt, xtb_files_props['xtb_opt'] , cwd=dat_dir)

            # replaces RDKit geometries with xTB geometries
            os.remove(xtb_files_props['xtb_xyz_path'])
            try:
                os.rename(str(dat_dir) + "/xtbopt.xyz", xtb_files_props['xtb_xyz_path'])
            except FileNotFoundError:
                os.rename(str(dat_dir) + "/xtblast.xyz", xtb_files_props['xtb_xyz_path'])

            final_xyz_path = xtb_files_props['xtb_xyz_path']
            self.final_xyz_path = final_xyz_path
            with open(final_xyz_path, "r") as xyz_file:
                self.xyz_coordinates = xyz_file.readlines() 


        command1 = [
            "xtb",
            xtb_files_props['xtb_xyz_path'],
            "--pop",
            "--wbo",
            "--acc",
            str(self.args.qdescp_acc),
            "--gfn",
            "2",
            "--chrg",
            str(charge),
            "--uhf",
            str(int(mult) - 1),
            "--etemp",
            str(self.args.qdescp_temp),
            "--input",
            str(xtb_input_file),
            "-P",
            "1",
        ]
        if self.args.qdescp_solvent is not None:
            command1.append("--alpb")
            command1.append(f"{self.args.qdescp_solvent}")
        run_command(command1, xtb_files_props['xtb_out'], cwd=dat_dir)

        os.rename(str(dat_dir) + "/xtbout.json", xtb_files_props['xtb_json'])
        os.rename(str(dat_dir) + "/wbo", xtb_files_props['xtb_wbo'])

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
            str(int(mult) - 1),
            "--etemp",
            str(self.args.qdescp_temp),
            "--vomega",
            "-P",
            "1",
        ]
        if self.args.qdescp_solvent is not None:
            command2.append("--alpb")
            command2.append(f"{self.args.qdescp_solvent}")
        run_command(command2, xtb_files_props['xtb_gfn1'], cwd=dat_dir)

        command3 = [
            "xtb",
            xtb_files_props['xtb_xyz_path'],
            "--vfukui",
            "--gfn",
            "2",
            "--chrg",
            str(charge),
            "--acc",
            str(self.args.qdescp_acc),
            "--uhf",
            str(int(mult) - 1),
            "--etemp",
            str(self.args.qdescp_temp),
            "-P",
            "1",
        ]
        if self.args.qdescp_solvent is not None:
            command3.append("--alpb")
            command3.append(f"{self.args.qdescp_solvent}")
        run_command(command3, xtb_files_props['xtb_fukui'], cwd=dat_dir)

        command4 = [
            "xtb",
            xtb_files_props['xtb_xyz_path'],
            "--fod",
            "--gfn",
            "2",
            "--chrg",
            str(charge),
            "--acc",
            str(self.args.qdescp_acc),
            "--uhf",
            str(int(mult) - 1),
            "--etemp",
            str(self.args.qdescp_temp),
            "-P",
            "1",
        ]
        if self.args.qdescp_solvent is not None:
            command4.append("--alpb")
            command4.append(f"{self.args.qdescp_solvent}")
        run_command(command4, xtb_files_props['xtb_fod'], cwd=dat_dir)

        command5 = [
            "xtb",
            xtb_files_props['xtb_xyz_path'],
            "--gfn",
            "1",
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
        ]
        if self.args.qdescp_solvent is not None:
            command5.append("--alpb")
            command5.append(f"{self.args.qdescp_solvent}")
        run_command(command5, xtb_files_props['xtb_Nminus1'], cwd=dat_dir)

        command6 = [
            "xtb",
            xtb_files_props['xtb_xyz_path'],
            "--gfn",
            "1",
            "--chrg",
            str(int(charge) +2),
            "--acc",
            str(self.args.qdescp_acc),
            "--uhf",
            str(int(mult) - 1),
            "--etemp",
            str(self.args.qdescp_temp),
            "-P",
            "1",
        ]
        if self.args.qdescp_solvent is not None:
            command6.append("--alpb")
            command6.append(f"{self.args.qdescp_solvent}")
        run_command(command6, xtb_files_props['xtb_Nminus2'], cwd=dat_dir)

        command7 = [
            "xtb",
            xtb_files_props['xtb_xyz_path'],
            "--gfn",
            "1",
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
        ]
        if self.args.qdescp_solvent is not None:
            command7.append("--alpb")
            command7.append(f"{self.args.qdescp_solvent}")
        run_command(command7, xtb_files_props['xtb_Nplus1'], cwd=dat_dir)

        command8 = [
            "xtb",
            xtb_files_props['xtb_xyz_path'],
            "--gfn",
            "1",
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
        ]
        if self.args.qdescp_solvent is not None:
            command8.append("--alpb")
            command8.append(f"{self.args.qdescp_solvent}")
        run_command(command8, xtb_files_props['xtb_Nplus2'], cwd=dat_dir)

        return xtb_files_props


    def collect_xtb_properties(self,name_initial,atom_props,update_atom_props,smarts_targets,xtb_files_props):
        """
        Collects all xTB properties from the files and puts them in a JSON file
        """
        properties_dict = read_xtb(xtb_files_props['xtb_out'])
        localgfn1 = read_gfn1(xtb_files_props['xtb_gfn1'])
        properties_FOD = read_fod(xtb_files_props['xtb_fod'])
        bonds, wbos = read_wbo(xtb_files_props['xtb_wbo'])
        cdft_descriptors = calculate_global_CDFT_descriptors(xtb_files_props['xtb_gfn1'])
        cdft_descriptors2  = calculate_global_CDFT_descriptors_part2( xtb_files_props['xtb_gfn1'], xtb_files_props['xtb_Nminus1'], xtb_files_props['xtb_Nminus2'], xtb_files_props['xtb_Nplus1'], xtb_files_props['xtb_Nplus2'], cdft_descriptors)
        localDescriptors = calculate_local_CDFT_descriptors(xtb_files_props['xtb_fukui'], cdft_descriptors, cdft_descriptors2)
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
        # From read_xtb
        json_data.update(properties_dict)
        # From read_json
        json_data.update(properties_FOD)
        # From read_gfn1
        json_data.update(localgfn1)
        #CDFT part 1
        json_data.update(cdft_descriptors)
        #CDFT part 2
        json_data.update(cdft_descriptors2)
        # Local Descriptors
        json_data.update(localDescriptors)

        """
		Morfeus Descriptors
		"""
        #Global descriptors
        global_properties_morfeus = calculate_global_morfeus_descriptors(self.final_xyz_path)
        json_data.update(global_properties_morfeus)       
        #Local descriptors
        local_properties_morfeus = calculate_local_morfeus_descriptors(self.final_xyz_path)
        json_data.update(local_properties_morfeus)  

        """
		Locate Smarts, qdescp_atoms
		"""
        with open(xtb_files_props['xtb_xyz_path'], "r") as f:
            inputs = f.readlines()

        coordinates = [inputs[i].strip().split()[1:] for i in range(2, int(inputs[0].strip()) + 2)]
        json_data["Coordinates"] = coordinates

        if len(smarts_targets) > 0:
            # Detect SMILES from SDF files generated by CSEARCH or create molecule objects.
            sdf_file = f'{name_initial}.sdf'
            with open(sdf_file, "r") as F:
                lines = F.readlines()

            smi_exist = False
            for i, line in enumerate(lines):
                if ">  <SMILES>" in line:
                    smi = lines[i + 1].split()[0]
                    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
                    smi_exist = True
                    break
            if not smi_exist:
                mol = Chem.SDMolSupplier(sdf_file, removeHs=False)

            # Search for target atoms or groups
            for pattern in smarts_targets:
                matches = []
                idx_set = None
                
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
                    self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} not found in the system, this molecule will not be used.")
                elif matches[0] == -1:
                    self.args.log.write(f"x  WARNING! Mapped atom {pattern} not found in the system, this molecule will not be used.")
                elif len(matches) > 1:
                    self.args.log.write(f"x  WARNING! More than one {pattern} atom was found in the system, this molecule will not be used.")
                elif len(matches) == 1:

                    atom_indices = list(matches[0])
                    atom_types = []
                    for atom_idx in atom_indices:
                        atom_types.append(mol.GetAtoms()[atom_idx].GetSymbol())

                    n_types = len(set(atom_types))
                    if n_types == 1:
                        sorted_indices = sorted(atom_indices, key=lambda idx: len(mol.GetAtoms()[idx].GetNeighbors()))
                    elif n_types > 1:
                        sorted_indices = sorted(atom_indices, key=lambda idx: mol.GetAtoms()[idx].GetAtomicNum())
                    
                    match_idx = 1
                    match_names = []
                    for atom_idx in sorted_indices:
                        atom_type = mol.GetAtoms()[atom_idx].GetSymbol()
                        if len(matches[0]) == 1:
                            if idx_set is None:
                                match_name = f'{atom_type}'
                            else:
                                match_name = f'{atom_type}{idx_set}'
                        else:
                            if n_types == 1:
                                match_name = f'{pattern}_{atom_type}{match_idx}'
                                match_idx += 1
                            elif n_types > 1:
                                match_name = f'{pattern}_{atom_type}'
                        match_names.append(match_name)

                    # Assign atomic descriptors to each identified atom
                    for atom_idx, match_name in zip(sorted_indices, match_names):
                        idx_xtb = atom_idx
                        for prop in atom_props:
                            json_data[f'{match_name}_{prop}'] = json_data[prop][idx_xtb]
                            if f'{match_name}_{prop}' not in update_atom_props:
                                update_atom_props.append(f'{match_name}_{prop}')

                    # Add maximum and minimum values for functional groups with the same two atoms
                    if len(match_names) > 1 and n_types == 1:
                        for prop in atom_props:
                            prop_values = []
                            for prop_name in match_names:
                                prop_values.append(json_data[f'{prop_name}_{prop}'])
                            json_data[f'{pattern}_max_{prop}'] = max(prop_values)
                            json_data[f'{pattern}_min_{prop}'] = min(prop_values)

        with open(xtb_files_props['xtb_json'], "w") as outfile:
            json.dump(json_data, outfile)

        return update_atom_props

    def cleanup(self, name, destination, xtb_passing, xtb_files_props):
        """
        Removes files from the xTB calculations that are not relevant and place json files in the 
        QDESCP folder
        """

        if xtb_passing: # only move molecules with successful xTB calcs
            final_json = str(destination) + "/" + name + ".json"
            shutil.move(xtb_files_props['xtb_json'], final_json)
        
        # delete xTB files that does not contain useful data
        files = glob.glob(f"{destination}/{name}/*")
        for file in files:
            if name not in os.path.basename(file) or '.inp' in os.path.basename(file):
                os.remove(file)
        if not os.path.exists(f"{destination}/xtb_data"): 
            Path(f"{destination}/xtb_data").mkdir()
        if os.path.exists(f"{destination}/xtb_data/{name}"): 
            self.args.log.write(f'\nx  A previous folder of {name} already existed, it was removed and replaced with the results of this QDESCP run.')
            shutil.rmtree(f"{destination}/xtb_data/{name}")
        shutil.move(f"{destination}/{name}", f"{destination}/xtb_data/{name}")

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
    
    def process_aqme_csv(self,name_db):
        csv_temp = pd.read_csv(f'{self.args.csv_name}')
        df_temp = pd.read_csv(f'AQME-{name_db}_{self.args.csv_name}')

        if len(df_temp) < len(csv_temp):
            missing_rows = csv_temp.loc[~csv_temp['code_name'].isin(df_temp['code_name'])]
            missing_rows[['code_name', 'SMILES']].to_csv(f'AQME-{name_db}_{self.args.csv_name}', mode='a', header=False, index=False)

            order = csv_temp['code_name'].tolist()

            df_temp = df_temp.sort_values(by='code_name', key=lambda x: x.map({v: i for i, v in enumerate(order)}))
            df_temp = df_temp.fillna(df_temp.groupby('SMILES').transform('first'))

            df_temp.to_csv(f'AQME-{name_db}_{self.args.csv_name}', index=False)

        return df_temp
