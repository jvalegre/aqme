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
   boltz : bool, default=True
      Calculation of Boltzmann averaged xTB properties and addition of RDKit 
      molecular descriptors

DBSTEP descriptors
++++++++++++++

   dbstep_r : float, default=3.5
      Radius used in the DBSTEP calculations (in A)

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
import numpy as np
from progress.bar import IncrementalBar
import pandas as pd
from rdkit import Chem
from pathlib import Path
import dbstep.Dbstep as db
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
    read_fukui
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

        if self.args.destination is None:
            destination = self.args.initial_dir.joinpath("QDESCP")
        else:
            destination = Path(self.args.destination)

        # retrieve the different files to run in QDESCP
        _ = check_files(self,'qdescp')

        qdescp_program = True
        if self.args.program is None:
            qdescp_program = False
        if qdescp_program:
            if self.args.program.lower() not in ["xtb", "nmr"]:
                qdescp_program = False
        if not qdescp_program:
            self.args.log.write("\nx  Program not supported for QDESCP descriptor generation! Specify: program='xtb' (or nmr)")
            self.args.log.finalize()
            sys.exit()

        update_atom_props = [] # keeps track of the molecules with suitable atomic properties when using qdescp_atoms

        self.args.log.write(f"\nStarting QDESCP-{self.args.program} with {len(self.args.files)} job(s)\n")

        # run the main xTB workflow
        if self.args.program.lower() == "xtb":
            mol_props = ["total energy","HOMO-LUMO gap/eV","electronic energy","Dipole module/D",
                "Total charge","HOMO","LUMO","Fermi-level/eV","Total dispersion C6",
                "Total dispersion C8","Total polarizability alpha","Total FOD"]
            atom_props = ["partial charges","mulliken charges","cm5 charges","FUKUI+","FUKUI-",
                "FUKUIrad","s proportion","p proportion","d proportion","Coordination numbers",
                "Dispersion coefficient C6","Polarizability alpha","FOD","FOD s proportion",
                "FOD p proportion","FOD d proportion",'DBSTEP_Vbur']
            if len(self.args.qdescp_atoms) == 0:
                atom_props.remove('DBSTEP_Vbur')

            update_atom_props = self.gather_files_and_run(destination,atom_props,update_atom_props)

        if len(update_atom_props) > 0:
            atom_props = update_atom_props

        # Boltzmann averaging of xTB values and DFT-NMR workflow
        qdescp_csv = "QDESCP_boltz_descriptors.csv"
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
                    name = file.replace("/", "\\").split("\\")[-1].split(".")[0]
                    json_files = glob.glob(
                        str(destination) + "/" + name + "_conf_*.json"
                    )
                    _ = get_boltz_props(json_files, name, boltz_dir, "xtb", self, mol_props, atom_props, mol=mol)
                _ = self.write_csv_boltz_data(destination,qdescp_csv)

            elif self.args.program.lower() == "nmr":
                mol_props = None
                atom_props = ["NMR Chemical Shifts"]
                if self.args.files[0].split('.')[1].lower() not in ["json"]:
                    self.args.log.write(f"\nx  The format used ({self.args.files[0].split('.')[1].lower()}) is not compatible with QDESCP with NMR! Formats accepted: json")
                    self.args.log.finalize()
                    sys.exit()

                name = self.args.files[0].replace("/", "\\").split("\\")[-1].split("_conf")[0]
                json_files = glob.glob(
                    str(os.path.dirname(os.path.abspath(self.args.files[0])))
                    + "/"
                    + name
                    + "_conf_*.json"
                )
                get_boltz_props(
                    json_files,
                    name,
                    boltz_dir,
                    "nmr",
                    self,
                    mol_props,
                    atom_props,
                    self.args.nmr_atoms,
                    self.args.nmr_slope,
                    self.args.nmr_intercept,
                    self.args.nmr_experim,
                )

        # AQME-ROBERT workflow
        if self.args.robert and self.args.program.lower() == "xtb":
            if self.args.csv_name is None:
                self.args.log.write(f"\nx  The input csv_name with SMILES and code_name columns are missing. A combined database for AQME-ROBERT workflows will not be created.")
            elif not Path(f"{self.args.csv_name}").exists():
                self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) is not valid. A combined database for AQME-ROBERT workflows will not be created.")
            else:
                combined_df = pd.DataFrame()
                qdescp_df = pd.read_csv(qdescp_csv)
                input_df = pd.read_csv(self.args.csv_name)
                if 'code_name' not in input_df.columns:
                    self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the code_name column. A combined database for AQME-ROBERT workflows will not be created.")
                elif 'SMILES' in input_df.columns or 'smiles' in input_df.columns or 'Smiles' in input_df.columns:
                    path_json = os.path.dirname(Path(qdescp_df['Name'][0]))
                    for i,input_name in enumerate(input_df['code_name']):
                        # match the entries of the two databases using the entry name
                        qdescp_col = input_df.loc[i].to_frame().T.reset_index(drop=True) # transposed, reset index
                        input_col = qdescp_df.loc[(qdescp_df['Name'] == f'{path_json}/{input_name}_rdkit_boltz') | (qdescp_df['Name'] == f'{path_json}/{input_name}_boltz') | (qdescp_df['Name'] == f'{path_json}/{input_name}_0_rdkit_boltz') | (qdescp_df['Name'] == f'{path_json}/{input_name}_1_rdkit_boltz') | (qdescp_df['Name'] == f'{path_json}/{input_name}_2_rdkit_boltz')]
                        input_col = input_col.drop(['Name'], axis=1).reset_index(drop=True)
                        combined_row = pd.concat([qdescp_col,input_col], axis=1)
                        combined_df = pd.concat([combined_df, combined_row], ignore_index=True)
                    combined_df = combined_df.dropna(axis=0)
                    csv_basename = os.path.basename(self.args.csv_name)
                    csv_path = self.args.initial_dir.joinpath(f'AQME-ROBERT_{csv_basename}')
                    _ = combined_df.to_csv(f'{csv_path}', index = None, header=True)
                    self.args.log.write(f"o  The AQME-ROBERT_{csv_basename} file containing the database ready for the AQME-ROBERT workflow was successfully created in {self.args.initial_dir}")
                else:
                    self.args.log.write(f"\nx  The input csv_name provided ({self.args.csv_name}) does not contain the SMILES column. A combined database for AQME-ROBERT workflows will not be created.")

        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime QDESCP: {elapsed_time} seconds\n")
        self.args.log.finalize()

    def write_csv_boltz_data(self, destination, qdescp_csv):
        """
        Concatenate the values for all calculations
        """
        
        boltz_json_files = glob.glob(str(destination) + "/boltz/*.json")
        dfs = [] 
        for file in boltz_json_files:
            data = pd.read_json(file, lines=True)
            data["Name"] = file.split(".json")[0]
            dfs.append(data)
        if len(dfs) > 0:
            temp = pd.concat(dfs, ignore_index=True) 
            temp.to_csv(qdescp_csv, index=False)
            self.args.log.write(f"o  The {qdescp_csv} file containing Boltzmann weighted xTB, DBSTEP and RDKit descriptors was successfully created in {self.args.initial_dir}")
        else:
            self.args.log.write(f"x  No CSV file containing Boltzmann weighted descriptors was created. This might happen when using the qdescp_atoms option with an atom/group that is not found in any of the calculations")

    def gather_files_and_run(self, destination, atom_props, update_atom_props):
        """
        Load all the input files, execute xTB calculations, gather descriptors and clean up scratch data
        """

        bar = IncrementalBar(
            "\no  Number of finished jobs from QDESCP", max=len(self.args.files)
        )
        # write input files
        if self.args.files[0].split('.')[1].lower() not in ["sdf", "xyz", "pdb"]:
            self.args.log.write(f"\nx  The format used ({self.args.files[0].split('.')[1].lower()}) is not compatible with QDESCP with xTB! Formats accepted: sdf, xyz, pdb")
            self.args.log.finalize()
            sys.exit()

        for file in self.args.files:
            xyz_files, xyz_charges, xyz_mults = [], [], []
            name = os.path.basename(file).split('.')[0]
            ext = os.path.basename(file).split(".")[1]
            self.args.log.write(f"\n\n   ----- {name} -----")
            if ext.lower() in ["sdf", "xyz", "pdb"]:
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

                elif ext.lower() == "pdb":
                    command_pdb = [
                        "obabel",
                        "-ipdb",
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

                elif ext.lower() == "sdf":
                    command_sdf = [
                        "obabel",
                        "-isdf",
                        file,
                        "-oxyz",
                        f"-O{os.path.dirname(os.path.abspath(file))}/{name}_conf_.xyz",
                        "-m",
                    ]
                    subprocess.run(
                        command_sdf,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )

            if ext.lower() in ["sdf", "pdb"]:
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

            for xyz_file, charge, mult in zip(xyz_files, xyz_charges, xyz_mults):
                name_xtb = os.path.basename(xyz_file.split(".")[0])
                self.args.log.write(f"\no   Running xTB and collecting properties")
                self.run_sp_xtb(xyz_file, charge, mult, name_xtb, destination)
                path_name = Path(os.path.dirname(file)).joinpath(os.path.basename(file).split(".")[0])
                update_atom_props = self.collect_xtb_properties(path_name, atom_props, update_atom_props)
                self.cleanup(name_xtb, destination)
            bar.next()
        bar.finish()

        return update_atom_props

    def run_sp_xtb(self, xyz_file, charge, mult, name, destination):
        """
        Runs single point xTB calculations to collect properties
        """

        dat_dir = destination / name
        dat_dir.mkdir(exist_ok=True, parents=True)

        self.xtb_xyz = str(dat_dir) + "/{0}.xyz".format(name)
        shutil.move(xyz_file, self.xtb_xyz)

        self.inp = str(dat_dir) + "/{0}_xtb.inp".format(name)
        with open(self.inp, "wt") as f:
            f.write("$write\n")
            f.write("json=true\n")

        self.xtb_out = str(dat_dir) + "/{0}.out".format(name)
        self.xtb_json = str(dat_dir) + "/{0}.json".format(name)
        self.xtb_wbo = str(dat_dir) + "/{0}.wbo".format(name)
        self.xtb_gfn1 = str(dat_dir) + "/{0}.gfn1".format(name)
        self.xtb_fukui = str(dat_dir) + "/{0}.fukui".format(name)
        self.xtb_fod = str(dat_dir) + "/{0}.fod".format(name)

        os.chdir(dat_dir)
        command1 = [
            "xtb",
            self.xtb_xyz,
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
            str(self.inp),
        ]
        if self.args.qdescp_solvent is not None:
            command1.append("--alpb")
            command1.append(f"{self.args.qdescp_solvent}")
        run_command(command1, self.xtb_out)

        os.rename("xtbout.json", self.xtb_json)
        os.rename("wbo", self.xtb_wbo)

        command2 = [
            "xtb",
            self.xtb_xyz,
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
        ]
        if self.args.qdescp_solvent is not None:
            command2.append("--alpb")
            command2.append(f"{self.args.qdescp_solvent}")
        run_command(command2, self.xtb_gfn1)

        command3 = [
            "xtb",
            self.xtb_xyz,
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
        ]
        if self.args.qdescp_solvent is not None:
            command3.append("--alpb")
            command3.append(f"{self.args.qdescp_solvent}")
        run_command(command3, self.xtb_fukui)

        command4 = [
            "xtb",
            self.xtb_xyz,
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
        ]
        if self.args.qdescp_solvent is not None:
            command4.append("--alpb")
            command4.append(f"{self.args.qdescp_solvent}")
        run_command(command4, self.xtb_fod)

        os.chdir(self.args.initial_dir)

    def collect_xtb_properties(self,name_initial,atom_props,update_atom_props):
        """
        Collects all xTB properties from the files and puts them in a JSON file
        """

        (
            _,
            total_charge,
            _,
            homo,
            lumo,
            atoms,
            _,
            _,
            dipole_module,
            Fermi_level,
            transition_dipole_moment,
            covCN,
            C6AA,
            alpha,
            homo_occ,
            lumo_occ,
            born_rad,
            SASA,
            h_bond,
            total_SASA,
            total_C6AA,
            total_C8AA,
            total_alpha,
        ) = read_xtb(self.xtb_out)
        FPLUS, FMINUS, FRAD = read_fukui(self.xtb_fukui)
        MULLIKEN, CM5, s_prop, p_prop, d_prop = read_gfn1(self.xtb_gfn1)
        total_fod, fod, s_prop_fod, p_prop_fod, d_prop_fod = read_fod(self.xtb_fod)
        bonds, wbos = read_wbo(self.xtb_wbo)

        # create matrix of Wiberg bond-orders
        nat = len(atoms)
        wbo_matrix = np.zeros((nat, nat))
        for i, bond in enumerate(bonds):
            wbo_matrix[(bond[0] - 1)][(bond[1] - 1)] = wbos[i]
            wbo_matrix[(bond[1] - 1)][(bond[0] - 1)] = wbos[i]

        """
		Now add xTB descriptors to existing json files.
		"""

        json_data = read_json(self.xtb_json)
        json_data["Dipole module/D"] = dipole_module
        json_data["Total charge"] = total_charge
        json_data["Transition dipole module/D"] = transition_dipole_moment
        json_data["HOMO"] = homo
        json_data["LUMO"] = lumo
        json_data["HOMO occupancy"] = homo_occ
        json_data["LUMO occupancy"] = lumo_occ
        json_data["mulliken charges"] = MULLIKEN
        json_data["cm5 charges"] = CM5
        json_data["FUKUI+"] = FPLUS
        json_data["FUKUI-"] = FMINUS
        json_data["FUKUIrad"] = FRAD
        json_data["s proportion"] = s_prop
        json_data["p proportion"] = p_prop
        json_data["d proportion"] = d_prop
        json_data["Fermi-level/eV"] = Fermi_level
        json_data["Coordination numbers"] = covCN
        json_data["Dispersion coefficient C6"] = C6AA
        json_data["Total dispersion C6"] = total_C6AA
        json_data["Total dispersion C8"] = total_C8AA
        json_data["Polarizability alpha"] = alpha
        json_data["Total polarizability alpha"] = total_alpha
        json_data["Wiberg matrix"] = wbo_matrix.tolist()
        json_data["Born radii"] = born_rad
        json_data["Atomic SASAs"] = SASA
        json_data["Solvent H bonds"] = h_bond
        json_data["Total SASA"] = total_SASA
        json_data["Total FOD"] = total_fod
        json_data["FOD"] = fod
        json_data["FOD s proportion"] = s_prop_fod
        json_data["FOD p proportion"] = p_prop_fod
        json_data["FOD d proportion"] = d_prop_fod

        with open(self.xtb_xyz, "r") as f:
            inputs = f.readlines()

        coordinates = [inputs[i].strip().split()[1:] for i in range(2, int(inputs[0].strip()) + 2)]
        json_data["coordinates"] = coordinates

        if len(self.args.qdescp_atoms) > 0:
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
                mol = Chem.SDMolSupplier(sdf_file, removeHs=False)

            # find the target atoms or groups
            for pattern in self.args.qdescp_atoms:
                matches = []
                try:
                    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                except: # I tried to make this except more specific for Boost.Python.ArgumentError, but apparently it's not as simple as it looks
                    try:
                        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(f'[{pattern}]'))
                    except:
                        self.args.log.write(f"x  WARNING! SMARTS pattern was not specified correctly! Make sure the qdescp_atoms option uses this format: \"[C]\" for atoms, \"[C=N]\" for bonds, and so on.")
      
                if len(matches) == 0:
                    self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} not found in the system, this molecule will not be used.")

                elif len(matches) > 1:
                    self.args.log.write(f"x  WARNING! More than one {pattern} atom was found in the system, this molecule will not be used.")

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
                    
                    match_idx = 1
                    match_names = []
                    # separates atoms when functional groups are used
                    for atom_idx in sorted_indices:
                        idx_dbstep, idx_xtb = None, None
                        idx_dbstep = atom_idx+1 # DBSTEP starts from index 1 (i.e. first atom has idx 1)
                        idx_xtb = atom_idx
                        atom_type = mol.GetAtoms()[atom_idx].GetSymbol()
                        if len(matches[0]) == 1:
                            match_name = f'{atom_type}'
                        else:
                            if n_types == 1:
                                match_name = f'{pattern}_{atom_type}{match_idx}'
                                match_idx += 1
                            elif n_types > 1:
                                match_name = f'{pattern}_{atom_type}'
                        match_names.append(match_name)

                        # calculate DBSTEP descriptors
                        self.args.log.write(f"\no   Running DBSTEP and collecting properties")

                        if idx_dbstep is not None:
                            # calculates buried volume to the type of atom selected
                            try:
                                dbstep_obj = db.dbstep(self.xtb_xyz,atom1=idx_dbstep,r=float(self.args.dbstep_r),commandline=True,verbose=False,volume=True)  
                                json_data[f'{match_name}_DBSTEP_Vbur'] = float(dbstep_obj.bur_vol)
                                if f'{match_name}_DBSTEP_Vbur' not in update_atom_props:
                                    update_atom_props.append(f'{match_name}_DBSTEP_Vbur')
                            except TypeError:
                                self.args.log.write(f'x  WARNING! DBSTEP is not working correctly, DBSTEP properties will not be calculated.')

                            # selects xTB atomic properties
                            for prop in atom_props:
                                if prop != 'DBSTEP_Vbur': # set the value of the atom instead of a list of values
                                    json_data[f'{match_name}_{prop}'] = json_data[prop][idx_xtb]
                                    if f'{match_name}_{prop}' not in update_atom_props:
                                        update_atom_props.append(f'{match_name}_{prop}')

                    # adding max and min values for functional groups with the same two atoms
                    if len(match_names) > 1 and n_types == 1:
                        for prop in atom_props:
                            prop_values = []
                            for prop_name in match_names:
                                prop_values.append(json_data[f'{prop_name}_{prop}'])
                            json_data[f'{pattern}_max_{prop}'] = max(prop_values)
                            if f'{pattern}_max_{prop}' not in update_atom_props:
                                update_atom_props.append(f'{pattern}_max_{prop}')
                            json_data[f'{pattern}_min_{prop}'] = min(prop_values)
                            if f'{pattern}_min_{prop}' not in update_atom_props:
                                update_atom_props.append(f'{pattern}_min_{prop}')

        with open(self.xtb_json, "w") as outfile:
            json.dump(json_data, outfile)
        
        return update_atom_props

    def cleanup(self, name, destination):
        final_json = str(destination) + "/" + name + ".json"
        shutil.move(self.xtb_json, final_json)
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

