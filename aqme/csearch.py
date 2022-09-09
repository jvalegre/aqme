#####################################################.
#          This file storesthe CSEARCH class        #
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

try:
    from rdkit.Chem import AllChem as Chem
    from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, Lipinski
except ModuleNotFoundError:
    print(
        "x  RDKit is not installed! You can install the program with 'conda install -c conda-forge rdkit' or 'pip install rdkit-pypi'"
    )
    sys.exit()
# this is a dummy import just to warn the user if Open babel is not installed
try:
    command_run_1 = ["obabel", "-H"]
    subprocess.run(command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
except FileNotFoundError:
    print(
        "x  Open Babel is not installed! You can install the program with 'conda install -c conda-forge openbabel'"
    )
    sys.exit()
from aqme.filter import filters, ewin_filter, pre_E_filter, RMSD_and_E_filter
from aqme.csearch_utils import (
    prepare_direct_smi,
    prepare_smiles_files,
    prepare_csv_files,
    prepare_cdx_files,
    prepare_com_files,
    prepare_sdf_files,
    prepare_pdb_files,
    template_embed,
    creation_of_dup_csv_csearch,
    minimize_rdkit_energy,
    com_2_xyz,
)
from aqme.fullmonte import generating_conformations_fullmonte, realign_mol
from aqme.utils import (
    rules_get_charge,
    substituted_mol,
    load_variables,
    set_metal_atomic_number,
    getDihedralMatches,
    smi_to_mol,
)
from aqme.crest import crest_opt


class csearch:
    """
    Class containing all the functions from the CSEARCH module.

    Parameters
    ----------
    kwargs : argument class

    Specify any arguments from the CSEARCH module (for a complete list of variables, visit the AQME documentation)
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
                self.args.log.write(
                    "x  xTB is not installed! You can install the program with 'conda install -c conda-forge xtb'"
                )
                sys.exit()

        if self.args.program.lower() not in ["rdkit", "summ", "fullmonte", "crest"]:
            self.args.log.write(
                "\nx  Program not supported for CSEARCH conformer generation! Specify: program='rdkit' (or summ, fullmonte, crest)"
            )
            self.args.log.finalize()
            sys.exit()

        if self.args.smi is None and self.args.input == "":
            self.args.log.write(
                "\nx  Program requires either a SMILES or an input file to proceed! Please look up acceptable file formats. Specify: smi='CCC' (or input='filename.csv')"
            )
            self.args.log.finalize()
            sys.exit()

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

            self.args.log.write(
                f"\nStarting CSEARCH with {len(job_inputs)} job(s) (SDF, XYZ, CSV, etc. files might contain multiple jobs/structures inside)\n"
            )

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

        # if large system increase stack size
        if self.args.stacksize != "1G":
            os.environ["OMP_STACKSIZE"] = self.args.stacksize

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
        self.final_dup_data = creation_of_dup_csv_csearch(self.args.program)

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

        # removing temporary files
        temp_files = [
            "gfn2.out",
            "xTB_opt.traj",
            "ANI1_opt.traj",
            "wbo",
            "xtbrestart",
            "ase.opt",
            "xtb.opt",
            "gfnff_topo",
        ]
        for file in temp_files:
            if os.path.exists(file):
                os.remove(file)

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

        if self.args.smi is not None:
            (
                mol,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
            ) = smi_to_mol(
                smi,
                name,
                self.args.program,
                self.args.log,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
            )

        else:
            if self.args.input.split(".")[1] in [
                "smi",
                "csv",
                "cdx",
                "txt",
                "yaml",
                "yml",
                "rtf",
            ]:
                (
                    mol,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                ) = smi_to_mol(
                    smi,
                    name,
                    self.args.program,
                    self.args.log,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                )
            else:
                (
                    mol,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                ) = (
                    smi,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                )

        if mol is None:
            self.args.log.write(
                f"\nx  Failed to convert the provided SMILES ({smi}) to RDkit Mol object! Please check the starting smiles."
            )
            self.args.log.finalize()
            sys.exit()

        if self.args.charge is None:
            self.args.charge = charge
        if self.args.mult is None:
            self.args.mult = mult

        if self.args.destination is None:
            self.csearch_folder = Path(self.args.initial_dir).joinpath(
                f"CSEARCH/{self.args.program}"
            )
        else:
            if Path(f"{self.args.destination}").exists():
                self.csearch_folder = Path(self.args.destination)
            else:
                self.csearch_folder = Path(self.args.initial_dir).joinpath(
                    self.args.destination
                )

        self.csearch_folder.mkdir(exist_ok=True, parents=True)

        if self.args.program in ["crest"] and self.args.smi is None:
            if self.args.input.split(".")[1] in ["pdb", "mol2", "mol", "sdf"]:
                command_pdb = [
                    "obabel",
                    f'-i{self.args.input.split(".")[1]}',
                    f'{name}.{self.args.input.split(".")[1]}',
                    "-oxyz",
                    f"-O{name}_{self.args.program}.xyz",
                ]
                subprocess.run(
                    command_pdb,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            elif self.args.input.split(".")[1] in ["gjf", "com"]:
                xyz_file, _, _ = com_2_xyz(f'{name}.{self.args.input.split(".")[1]}')
                os.move(xyz_file, f"{name}_{self.args.program}.xyz")
            elif self.args.input.split(".")[1] == "xyz":
                shutil.copy(f"{name}.xyz", f"{name}_{self.args.program}.xyz")

        # Converts each line to an RDKit mol object
        if self.args.verbose:
            self.args.log.write(f"\n   -> Input Molecule {Chem.MolToSmiles(mol)}")

        if self.args.metal_complex:
            for _ in self.args.metal_atoms:
                self.args.metal_idx.append(None)
                self.args.complex_coord.append(None)
                self.args.metal_sym.append(None)

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
            if self.args.complex_type in accepted_complex_types:
                count_metals = 0
                for metal_idx_ind in self.args.metal_idx:
                    if metal_idx_ind is not None:
                        count_metals += 1
                if count_metals == 1:
                    template_kwargs = dict()
                    template_kwargs["complex_type"] = self.args.complex_type
                    template_kwargs["metal_idx"] = self.args.metal_idx
                    template_kwargs["maxsteps"] = self.args.opt_steps_rdkit
                    template_kwargs["heavyonly"] = self.args.heavyonly
                    template_kwargs["maxmatches"] = self.args.max_matches_rmsd
                    template_kwargs["mol"] = mol
                    items = template_embed(self, **template_kwargs)

                    total_data = creation_of_dup_csv_csearch(self.args.program)

                    for mol_obj, name_in, coord_map, alg_map, template in zip(*items):
                        data = self.conformer_generation(
                            mol_obj,
                            name_in,
                            constraints_atoms,
                            constraints_dist,
                            constraints_angle,
                            constraints_dihedral,
                            coord_map,
                            alg_map,
                            template,
                        )
                        frames = [total_data, data]
                        total_data = pd.concat(frames, sort=True)
                else:
                    self.args.log.write(
                        "\nx  Cannot use templates for complexes involving more than 1 metal or for organic molecueles."
                    )
                    total_data = None
            else:
                total_data = self.conformer_generation(
                    mol,
                    name,
                    constraints_atoms,
                    constraints_dist,
                    constraints_angle,
                    constraints_dihedral,
                )
        else:
            total_data = self.conformer_generation(
                mol,
                name,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
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
        coord_Map=None,
        alg_Map=None,
        mol_template=None,
    ):
        """
        Function to load mol objects and create 3D conformers

        """

        dup_data = creation_of_dup_csv_csearch(self.args.program)

        dup_data_idx = 0
        start_time = time.time()
        status = None

        self.args.log.write(f"\n   ----- {name} -----")

        ### Fixing all charges here

        # user can overwrite charge and mult with the corresponding arguments
        if self.args.charge is not None:
            charge = self.args.charge

        elif self.args.charge is None:
            if not self.args.metal_complex:
                charge = Chem.GetFormalCharge(mol)
            else:
                charge = rules_get_charge(mol, self.args, "csearch")
        else:
            charge = 0

        if self.args.mult is not None:
            mult = self.args.mult
        elif self.args.mult is None:
            # if not self.args.metal_complex:
            NumRadicalElectrons = 0
            for Atom in mol.GetAtoms():
                NumRadicalElectrons += Atom.GetNumRadicalElectrons()
            TotalElectronicSpin = NumRadicalElectrons / 2
            mult = int((2 * TotalElectronicSpin) + 1)
        else:
            mult = 1

        dup_data.at[dup_data_idx, "Real charge"] = charge
        dup_data.at[dup_data_idx, "Mult"] = mult

        # inputs that go through CREST containing 3D coordinates don't require a previous RDKit conformer sampling
        if (
            self.args.program in ["crest"]
            and self.args.smi is None
            and self.args.input.split(".")[1]
            in [
                "pdb",
                "mol2",
                "mol",
                "sdf",
                "gjf",
                "com",
                "xyz",
            ]
        ):
            valid_structure = True

            status = crest_opt(
                f"{name}_{self.args.program}",
                dup_data,
                dup_data_idx,
                self.args,
                charge,
                mult,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
            )

        else:
            name = name.replace("/", "\\").split("\\")[-1].split(".")[0]
            self.csearch_file = self.csearch_folder.joinpath(
                name + "_" + self.args.program + self.args.output
            )
            self.sdwriter = Chem.SDWriter(str(self.csearch_file))

            valid_structure = filters(
                mol, self.args.log, self.args.max_mol_wt, self.args.verbose
            )
            if valid_structure:
                try:
                    # the conformational search for RDKit
                    status, update_to_rdkit = self.summ_search(
                        mol,
                        name,
                        dup_data,
                        dup_data_idx,
                        charge,
                        mult,
                        constraints_atoms,
                        constraints_dist,
                        constraints_angle,
                        constraints_dihedral,
                        coord_Map,
                        alg_Map,
                        mol_template,
                    )
                    dup_data.at[dup_data_idx, "status"] = status
                    dup_data.at[dup_data_idx, "update_to_rdkit"] = update_to_rdkit
                except (KeyboardInterrupt, SystemExit):
                    raise

        if status == -1 or not valid_structure:
            error_message = "\nx  ERROR: The structure is not valid or no conformers were obtained from this SMILES string"
            self.args.log.write(error_message)

        # if self.args.program == "crest" and valid_structure:
        #     shutil.rmtree(f"{self.csearch_folder}/../crest_xyz")

        n_seconds = round(time.time() - start_time, 2)
        dup_data.at[dup_data_idx, "CSEARCH time (seconds)"] = n_seconds

        return dup_data

    def summ_search(
        self,
        mol,
        name,
        dup_data,
        dup_data_idx,
        charge,
        mult,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
        coord_Map=None,
        alg_Map=None,
        mol_template=None,
    ):

        """
        Embeds, optimizes and filters RDKit conformers
        """

        # writes sdf for the first RDKit conformer generation
        status, rotmatches, update_to_rdkit, ff = self.rdkit_to_sdf(
            mol,
            name,
            dup_data,
            dup_data_idx,
            self.sdwriter,
            charge,
            mult,
            coord_Map,
            alg_Map,
            mol_template,
            constraints_atoms,
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
        )
        # reads the initial SDF files from RDKit and uses dihedral scan if selected
        if status not in [-1, 0]:
            # getting the energy and mols after rotations
            if self.args.program == "summ" and len(rotmatches) != 0:
                status = self.dihedral_filter_and_sdf(
                    name, dup_data, dup_data_idx, coord_Map, alg_Map, mol_template, ff
                )

        return status, update_to_rdkit

    def dihedral_filter_and_sdf(
        self, name, dup_data, dup_data_idx, coord_Map, alg_Map, mol_template, ff
    ):
        """
        Filtering after dihedral scan to sdf
        """

        rotated_energy = []

        rdmols = Chem.SDMolSupplier(str(self.csearch_file), removeHs=False)

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
            self.args,
            dup_data,
            dup_data_idx,
            self.args.log,
            "summ",
            self.args.ewin_csearch,
        )
        # pre-filter based on energy only
        selectedcids_initial_rotated = pre_E_filter(
            sortedcids_rotated,
            rotated_energy,
            dup_data,
            dup_data_idx,
            self.args.log,
            "summ",
            self.args.initial_energy_threshold,
            self.args.verbose,
        )
        # filter based on energy and RMSD
        selectedcids_rotated = RMSD_and_E_filter(
            rdmols,
            selectedcids_initial_rotated,
            rotated_energy,
            self.args,
            dup_data,
            dup_data_idx,
            self.args.log,
            "summ",
        )

        os.remove(self.csearch_file)
        sdwriter_rd = Chem.SDWriter(str(self.csearch_file))
        for i, cid in enumerate(selectedcids_rotated):
            mol_rd = Chem.RWMol(rdmols[cid])
            mol_rd.SetProp("_Name", rdmols[cid].GetProp("_Name") + " " + str(i))
            mol_rd.SetProp("Energy", str(rotated_energy[cid]))
            if self.args.metal_complex:
                set_metal_atomic_number(
                    mol_rd, self.args.metal_idx, self.args.metal_sym
                )
            sdwriter_rd.write(mol_rd)
        sdwriter_rd.close()
        status = 1
        return status

    def auto_sampling(self, mol):
        """
        Detects automatically the initial number of conformers for the sampling
        """

        if self.args.metal_complex:
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
        If program = RDKit, this replaces iodine back to the metal (for metal_complex = True)
        and writes the RDKit SDF files. With program = summ, this function optimizes rotamers
        """

        if i >= len(matches):  # base case, torsions should be set in conf
            # setting the metal back instead of I
            if self.args.metal_complex and (
                self.args.program == "rdkit" or update_to_rdkit
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
                set_metal_atomic_number(mol, self.args.metal_idx, self.args.metal_sym)
            try:
                sdwriter.write(mol, conf)
            except (TypeError):
                raise

            return 1

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
            self.args.log.write(
                "\no  Normal RDKit embeding process failed, trying to "
                "generate conformers with random coordinates "
                f"(with {str(initial_confs)} possibilities)"
            )
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
            else:  # id template realign before doing calculations
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
            self.args,
            dup_data,
            dup_data_idx,
            self.args.log,
            "rdkit",
            self.args.ewin_csearch,
        )

        # pre-filter based on energy only
        selectedcids_initial_rdkit = pre_E_filter(
            sortedcids_rdkit,
            cenergy,
            dup_data,
            dup_data_idx,
            self.args.log,
            "rdkit",
            self.args.initial_energy_threshold,
            self.args.verbose,
        )

        # filter based on energy and RMSD
        selectedcids_rdkit = RMSD_and_E_filter(
            outmols,
            selectedcids_initial_rdkit,
            cenergy,
            self.args,
            dup_data,
            dup_data_idx,
            self.args.log,
            "rdkit",
        )

        if self.args.program == "summ" or self.args.program == "rdkit":
            # now exhaustively drive torsions of selected conformers
            n_confs = int(
                len(selectedcids_rdkit) * (360 / self.args.degree) ** len(rotmatches)
            )
            if self.args.verbose and len(rotmatches) != 0:
                self.args.log.write(
                    "\no  Systematic generation of " + str(n_confs) + " confomers"
                )

            total = 0
            for conf in selectedcids_rdkit:
                if self.args.program == "summ" and not update_to_rdkit:
                    sdwriter.write(outmols[conf], conf)
                    for m in rotmatches:
                        rdMolTransforms.SetDihedralDeg(
                            outmols[conf].GetConformer(conf), *m, 180.0
                        )
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

            if self.args.verbose and len(rotmatches) != 0:
                self.args.log.write("\no  %d total conformations generated" % total)
            status = 1

        if self.args.program == "summ":
            dup_data.at[dup_data_idx, "summ-conformers"] = total

        if self.args.program == "fullmonte":
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

        return status

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
        mol_template,
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral,
    ):

        """
        Conversion from RDKit to SDF
        """

        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        mol.SetProp("_Name", name)

        # detects and applies auto-detection of initial number of conformers
        if self.args.sample == "auto":
            initial_confs = int(self.auto_sampling(mol))
        else:
            initial_confs = int(self.args.sample)

        dup_data.at[dup_data_idx, "Molecule"] = name
        update_to_rdkit = False

        rotmatches = getDihedralMatches(mol, self.args.heavyonly)

        if len(rotmatches) > self.args.max_torsions and self.args.max_torsions > 0:
            self.args.log.write(
                "\nx  Too many torsions (%d). Skipping %s"
                % (len(rotmatches), (name + self.args.output))
            )
        elif self.args.program == "summ" and len(rotmatches) == 0:
            update_to_rdkit = True
            self.args.log.write(
                "\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to SUMM SDF"
            )
        elif self.args.program == "fullmonte" and len(rotmatches) == 0:
            update_to_rdkit = True
            self.args.log.write(
                "\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to FULLMONTE SDF"
            )

        ff = self.args.ff
        if self.args.program != "crest":
            dup_data.at[dup_data_idx, "RDKit-Initial-samples"] = initial_confs
            if self.args.program == "rdkit":
                rotmatches = []
            cids = self.embed_conf(mol, initial_confs, coord_Map, alg_Map, mol_template)

            # energy minimize all to get more realistic results
            # identify the atoms and decide Force Field
            for atom in mol.GetAtoms():
                if (
                    atom.GetAtomicNum() > 36
                ) and self.args.ff == "MMFF":  # up to Kr for MMFF, if not the code will use UFF
                    self.args.log.write(
                        "\nx  "
                        + self.args.ff
                        + " is not compatible with the molecule, changing to UFF"
                    )
                    ff = "UFF"
            if self.args.verbose:
                self.args.log.write(
                    "\no  Optimizing "
                    + str(len(cids))
                    + " initial conformers with "
                    + ff
                )
                if self.args.program == "summ":
                    self.args.log.write(
                        "\no  Found " + str(len(rotmatches)) + " rotatable torsions"
                    )
                elif self.args.program == "fullmonte":
                    self.args.log.write(
                        "\no  Found " + str(len(rotmatches)) + " rotatable torsions"
                    )
                else:
                    self.args.log.write(
                        "\no  Systematic torsion rotation is set to OFF"
                    )

            try:
                status = self.min_after_embed(
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
        else:
            dup_data.at[dup_data_idx, "Real charge"] = charge
            dup_data.at[dup_data_idx, "Mult"] = mult

            status = crest_opt(
                f"{name}_{self.args.program}",
                dup_data,
                dup_data_idx,
                self.args,
                charge,
                mult,
                constraints_atoms,
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
                mol,
            )

        sdwriter.close()

        return status, rotmatches, update_to_rdkit, ff
