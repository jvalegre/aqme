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

XTB descriptors
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
from aqme.utils import (
    load_variables,
    read_xyz_charge_mult,
    mol_from_sdf_or_mol_or_mol2,
    run_command
)
from aqme.qdescp_utils import (
    get_boltz_avg_properties_xtb,
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
    Class containing all the functions from the QDESCP module related to xTB properties for SDF files.
    """

    def __init__(self, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "qdescp")

        if self.args.destination is None:
            destination = self.args.initial_dir.joinpath("QDESCP")
        else:
            destination = Path(self.args.destination)

        if len(self.args.files) == 0:
            self.args.log.write('\nx  No files were found! Make sure you use quotation marks if you are using * (i.e. --files "*.sdf")')
            self.args.log.finalize()
            sys.exit()

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

        if self.args.program.lower() == "xtb":
            self.gather_files_and_run(destination)

        if self.args.boltz == "False":
            self.args.boltz = False
        elif self.args.boltz == "True":
            self.args.boltz = True
        if self.args.boltz:
            boltz_dir = Path(f"{destination}/boltz")
            boltz_dir.mkdir(exist_ok=True, parents=True)
            if self.args.program.lower() == "xtb":
                for file in self.args.files:
                    mol = Chem.SDMolSupplier(file, removeHs=False)[0]
                    name = file.replace("/", "\\").split("\\")[-1].split(".")[0]
                    json_files = glob.glob(
                        str(destination) + "/" + name + "_conf_*.json"
                    )
                    get_boltz_avg_properties_xtb(
                        json_files, name, boltz_dir, "xtb", self, None, None, None, None, mol
                    )
                self.write_csv_boltz_data(destination)

            elif self.args.program.lower() == "nmr":
                if self.args.files[0].split('.')[1].lower() not in ["json"]:
                    self.args.log.write(f"\nx  The format used ({self.args.files[0].split('.')[1].lower()}) is not compatible with QDESCP with NMR! Formats accepted: json")
                    self.args.log.finalize()
                    sys.exit()

                for file in self.args.files:
                    name = file.replace("/", "\\").split("\\")[-1].split("_conf")[0]
                    json_files = glob.glob(
                        str(os.path.dirname(os.path.abspath(file)))
                        + "/"
                        + name
                        + "_conf_*.json"
                    )
                    get_boltz_avg_properties_xtb(
                        json_files,
                        name,
                        boltz_dir,
                        "nmr",
                        self,
                        self.args.nmr_atoms,
                        self.args.nmr_slope,
                        self.args.nmr_intercept,
                        self.args.nmr_experim,
                    )

        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime QDESCP: {elapsed_time} seconds\n")
        self.args.log.finalize()

    def write_csv_boltz_data(self, destination):
        boltz_json_files = glob.glob(str(destination) + "/boltz/*.json")
        dfs = []  # an empty list to store the data frames
        for file in boltz_json_files:
            data = pd.read_json(file, lines=True)  # read data frame from json file
            data["Name"] = file.split(".json")[0]
            dfs.append(data)  # append the data frame to the list

        temp = pd.concat(
            dfs, ignore_index=True
        )  # concatenate all the data frames in the list.
        qdescp_csv = "QDESCP_boltz_descriptors.csv"
        temp.to_csv(qdescp_csv, index=False)
        self.args.log.write(f"o  The {qdescp_csv} file containing Boltzmann weighted xTB and RDKit descriptors was successfully created in {self.args.initial_dir}")

    def gather_files_and_run(self, destination):
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
            name = file.replace("/", "\\").split("\\")[-1].split(".")[0]
            self.args.log.write(f"\n\n   ----- {name} -----")
            if file.split(".")[1].lower() in ["sdf", "xyz", "pdb"]:
                if file.split(".")[1].lower() == "xyz":
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

                elif file.split(".")[1].lower() == "pdb":
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

                elif file.split(".")[1].lower() == "sdf":
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

            if file.split(".")[1].lower() in ["sdf", "pdb"]:
                if self.args.charge is None:
                    _, charges, _, _ = mol_from_sdf_or_mol_or_mol2(file, "csearch")
                else:
                    charges = [self.args.charge] * len(
                        glob.glob(
                            f"{os.path.dirname(os.path.abspath(file))}/{name}_conf_*.xyz"
                        )
                    )
                if self.args.mult is None:
                    _, _, mults, _ = mol_from_sdf_or_mol_or_mol2(file, "csearch")
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
                name = os.path.basename(xyz_file.split(".")[0])
                self.run_sp_xtb(xyz_file, charge, mult, name, destination)
                self.collect_xtb_properties()
                self.cleanup(name, destination)
            bar.next()
        bar.finish()

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

    def collect_xtb_properties(self):
        """
        Collects all xTB properties from the files and puts them in a JSON file
        """

        (
            energy,
            total_charge,
            homo_lumo,
            homo,
            lumo,
            atoms,
            numbers,
            chrgs,
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

        coordinates = [
            inputs[i].strip().split()[1:] for i in range(2, int(inputs[0].strip()) + 2)
        ]

        json_data["coordinates"] = coordinates

        with open(self.xtb_json, "w") as outfile:
            json.dump(json_data, outfile)

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

