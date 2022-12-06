#####################################################.
#          This file stores the CMIN class          #
#             used in conformer refinement          #
#####################################################.

import os
import sys
import glob
import subprocess
import numpy as np
from pathlib import Path
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors as Descriptors
from rdkit.Chem.PropertyMol import PropertyMol
from progress.bar import IncrementalBar
from rdkit.Geometry import Point3D
import pandas as pd
import time
from ase.units import Hartree
from aqme.utils import load_variables, rules_get_charge, substituted_mol
from aqme.filter import ewin_filter, pre_E_filter, RMSD_and_E_filter
from aqme.cmin_utils import creation_of_dup_csv_cmin, rdkit_sdf_read
from aqme.csearch_crest_utils import xtb_opt_main
from aqme.csearch_utils import prepare_com_files

hartree_to_kcal = 627.509


class cmin:
    """
    Class containing all the functions from the CMIN module.

    Parameters
    ----------
    kwargs : argument class
        Specify any arguments from the CMIN module (for a complete list of variables, visit the AQME documentation)
    """

    def __init__(self, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "cmin")
        if self.args.program.lower() not in ["xtb", "ani"]:
            self.args.log.write("\nx  Program not supported for CMIN refinement! Specify: program='xtb' (or 'ani')")
            self.args.log.finalize()
            sys.exit()

        try:
            os.chdir(self.args.w_dir_main)
        except FileNotFoundError:
            self.args.w_dir_main = Path(f"{os.getcwd()}/{self.args.w_dir_main}")
            os.chdir(self.args.w_dir_main)

        # create the dataframe to store the data
        self.final_dup_data = creation_of_dup_csv_cmin(self.args.program.lower())

        bar = IncrementalBar(
            "\no  Number of finished jobs from CMIN", max=len(self.args.files)
        )

        # retrieves the different files to run in CMIN
        file_format = os.path.splitext(self.args.files[0])[1]
        if file_format.lower() in ['.xyz', '.gjf', '.com']:
            for file in self.args.files:
                prepare_com_files(self.args, file)
            if file_format.lower() in ['.gjf', '.com']:
                files_temp_extra = glob.glob('*.xyz')
            files_cmin = glob.glob('*.sdf')
        elif file_format.lower() == '.pdb':
            for file in self.args.files:
                command_pdb = [
                    "obabel",
                    "-ipdb",
                    f"{file}",
                    "-osdf",
                    f"-O{file.split('.pdb')[0]}.sdf",
                ]
                subprocess.run(
                    command_pdb,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            files_cmin = glob.glob('*.sdf')
        elif file_format.lower() == '.sdf':
            files_cmin = self.args.files
        else:
            self.args.log.write(f"\nx  The input format {file_format} is not supported for CMIN refinement! Formats allowed: SDF, XYZ, COM, GJF and PDB")
            self.args.log.finalize()
            sys.exit()

        for file in files_cmin:
            # load jobs for cmin minimization
            self.mols, self.name = self.load_jobs(file)
            self.args.log.write(f"\n\n   ----- {self.name} -----")

            if self.args.destination is None:
                self.cmin_folder = Path(self.args.w_dir_main).joinpath(
                    f"CMIN"
                )
            else:
                if Path(f"{self.args.destination}").exists():
                    self.cmin_folder = Path(self.args.destination)
                else:
                    self.cmin_folder = Path(self.args.initial_dir).joinpath(
                    self.args.destination)

            self.cmin_folder.mkdir(exist_ok=True, parents=True)

            self.cmin_all_file = self.cmin_folder.joinpath(
                f"{self.name}_{self.args.program.lower()}_all_confs{self.args.output}"
            )
            self.sdwriterall = Chem.SDWriter(str(self.cmin_all_file))

            self.cmin_file = self.cmin_folder.joinpath(
                self.name + "_" + self.args.program.lower() + self.args.output
            )
            self.sdwriter = Chem.SDWriter(str(self.cmin_file))

            # runs the conformer sampling with multiprocessors
            total_data = self.compute_cmin(file)

            frames = [self.final_dup_data, total_data]
            self.final_dup_data = pd.concat(frames, ignore_index=True, sort=True)
            bar.next()
        bar.finish()
        
        # store all the information into a CSV file
        cmin_csv_file = self.args.w_dir_main.joinpath("CMIN-Data.csv")
        self.final_dup_data.to_csv(cmin_csv_file, index=False)

        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime CMIN: {elapsed_time} seconds\n")
        self.args.log.finalize()

        # delete extra temporary files created when using XYZ, GJF, COM and PDB files
        if file_format.lower() in ['.xyz', '.gjf', '.com', '.pdb']:
            if file_format.lower() in ['.gjf', '.com']:
                files_cmin = files_cmin + files_temp_extra
            for temp_file in files_cmin:
                os.remove(temp_file)

        # this is added to avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

    def load_jobs(self, file):

        self.args.log.write(f"\n\no  Multiple minimization of {file} with {self.args.program}")

        # read SDF files from RDKit optimization
        inmols = rdkit_sdf_read(file, self.args)
        name_mol = os.path.basename(file).split(".sdf")[0]

        return inmols, name_mol

    def compute_cmin(self, file):

        dup_data = creation_of_dup_csv_cmin(self.args.program.lower())
        dup_data_idx = 0
        dup_data.at[dup_data_idx, "Molecule"] = self.name
        cenergy, outmols = [], []
        start_time = time.time()

        if self.args.program.lower() == "ani":
            if self.args.charge is not None:
                self.args.log.write("\nx  Charge is automatically calculated for ANI methods, do not use the charge option!")
                self.args.log.finalize()
                sys.exit()
            elif self.args.mult is not None:
                self.args.log.write("\nx  Multiplicity is automatically calculated for ANI methods, do not use the mult option!")
                self.args.log.finalize()
                sys.exit()
            charge,mult,final_mult,dup_data = self.charge_mult_cmin(dup_data, dup_data_idx, 'ani')

        elif self.args.program.lower() == "xtb":
            # sets charge and mult
            file_format = os.path.splitext(file)[1]
            charge_input, mult_input, final_mult = None, None, None
            if file_format.lower() == '.sdf':
                if self.args.charge is None or self.args.mult is None:
                    # read charge and mult from SDF if possible (i.e. charge/mult of SDFs created with CSEARCH)
                    with open(file, "r") as F:
                        lines = F.readlines()
                    charge_found, mult_found = False, False
                    for i, line in enumerate(lines):
                        if line.find(">  <Real charge>") > -1:
                            charge_input = lines[i + 1].split()[0]
                            charge_found = True
                        if line.find(">  <Mult>") > -1:
                            mult_input = lines[i + 1].split()[0]
                            mult_found = True
                        if charge_found and mult_found:
                            break
            if self.args.charge is None and charge_input is None:
                # if no charge/mult was specified or found, the charge is calculated using the mol object
                charge,_,final_mult,_ = self.charge_mult_cmin(dup_data, dup_data_idx, 'xtb')
            elif self.args.charge is None:
                charge = charge_input
            else:
                charge = self.args.charge
            if self.args.mult is None and mult_input is None:
                if final_mult is None:
                    _,_,final_mult,_ = self.charge_mult_cmin(dup_data, dup_data_idx, 'xtb')
                mult = final_mult
            elif self.args.mult is None:
                mult = mult_input
            else:
                mult = self.args.mult

        for i, mol in enumerate(self.mols):
            if mol is not None:
                # ANI calculations use ASE to run
                if self.args.program.lower() == "ani":
                    mol, energy, cmin_valid = self.ani_optimize(mol,charge,mult)
                # xTB calculations use the xTB program directly
                elif self.args.program.lower() == "xtb":
                    # for contrained optimizations
                    complex_ts = False
                    if len(self.args.constraints_atoms) >= 1 or len(self.args.constraints_dist) >= 1 or len(self.args.constraints_angle) >= 1 or len(self.args.constraints_dihedral) >= 1:
                        complex_ts = True
                    name_init = str(open(file, "r").readlines()[0].strip())
                    mol, energy, cmin_valid = xtb_opt_main(
                        f'{self.name}_conf_{i}',
                        dup_data,
                        dup_data_idx,
                        self,
                        charge,
                        mult,
                        self.args.constraints_atoms,
                        self.args.constraints_dist,
                        self.args.constraints_angle,
                        self.args.constraints_dihedral,
                        'xtb',
                        complex_ts=complex_ts,
                        mol=mol,
                        name_init=name_init
                    )
                if cmin_valid:
                    pmol = PropertyMol(mol)
                    outmols.append(pmol)
                    cenergy.append(energy)

        if len(cenergy) >= 1:
            # if SQM energy exists, overwrite RDKit energies and geometries
            cids = list(range(len(outmols)))
            sorted_all_cids = sorted(cids, key=lambda cid: cenergy[cid])

            name_mol = self.name

            for cid in sorted_all_cids:
                outmols[cid].SetProp(
                    "_Name", outmols[cid].GetProp("_Name") + " " + self.args.program.lower()
                )
                outmols[cid].SetProp("Energy", cenergy[cid])
                if self.args.program.lower() == "ani":
                    outmols[cid].SetProp("Real charge", str(np.sum(charge)))
                    outmols[cid].SetProp("Mult", str(final_mult))
                elif self.args.program.lower() == "xtb":
                    outmols[cid].SetProp("Real charge", str(charge))
                    outmols[cid].SetProp("Mult", str(mult))

            write_all_confs = 0
            for cid in sorted_all_cids:
                self.sdwriterall.write(outmols[cid])
                write_all_confs += 1
            self.sdwriterall.close()

            self.args.log.write(f"\no  Applying filters to intial conformers after {self.args.program.lower()} minimization")

            # filter based on energy window ewin_cmin
            sortedcids = ewin_filter(
                sorted_all_cids,
                cenergy,
                self.args,
                dup_data,
                dup_data_idx,
                self.args.log,
                self.args.program.lower(),
                self.args.ewin_cmin,
            )
            # pre-filter based on energy only
            selectedcids_initial = pre_E_filter(
                sortedcids,
                cenergy,
                dup_data,
                dup_data_idx,
                self.args.log,
                self.args.program.lower(),
                self.args.initial_energy_threshold,
            )
            # filter based on energy and RMSD
            selectedcids = RMSD_and_E_filter(
                outmols,
                selectedcids_initial,
                cenergy,
                self.args,
                dup_data,
                dup_data_idx,
                self.args.log,
                self.args.program.lower(),
            )

            if self.args.program.lower() == "xtb":
                dup_data.at[dup_data_idx, "xTB-Initial-samples"] = len(self.mols)
            elif self.args.program.lower() == "ani":
                dup_data.at[dup_data_idx, "ANI-Initial-samples"] = len(self.mols)

            # write the filtered, ordered conformers to external file
            self.write_confs(
                outmols, selectedcids, self.args.log
            )

        dup_data.at[dup_data_idx, "CMIN time (seconds)"] = round(
            time.time() - start_time, 2
        )

        # removing temporary files
        temp_files = [
            "gfn2.out",
            "cmin_opt.traj",
            "wbo",
            "xtbrestart",
            "ase.opt",
            "cmin.opt",
            "gfnff_topo"
        ]
        for file in temp_files:
            if os.path.exists(file):
                os.remove(file)

        return dup_data

    # xTB AND ANI MAIN OPTIMIZATION PROCESS
    def ani_optimize(self, mol, charge, mult):

        # Attempts ANI/xTB imports and exits if the programs are not installed
        try:
            import torch
            import warnings
            warnings.filterwarnings('ignore')

        except ModuleNotFoundError:
            self.args.log.write("x  Torch-related modules are not installed! You can install these modules with 'pip install torch torchvision torchani'")
            self.args.log.finalize()
            sys.exit()
        try:
            import ase
            import ase.optimize
        except ModuleNotFoundError:
            self.args.log.write("x  ASE is not installed! You can install the program with 'conda install -c conda-forge ase' or 'pip install ase'")
            self.args.log.finalize()
            sys.exit()

        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        DEVICE = torch.device("cpu")

        # if a large system is used, you might need to increase the stack size
        os.environ["OMP_STACKSIZE"] = self.args.stacksize

        elements = ""
        for _, atom in enumerate(mol.GetAtoms()):
            elements += atom.GetSymbol()

        cartesians = mol.GetConformers()[0].GetPositions()
        coordinates = torch.tensor(
            [cartesians.tolist()], requires_grad=True, device=DEVICE
        )
        cmin_valid = True
        model = self.get_cmin_model()

        # define ase molecule using ANI calculator
        ase_molecule = ase.Atoms(
            elements, positions=coordinates.tolist()[0], calculator=model.ase()
        )

        # Adjust charge and multiplicity from the input SDFs
        for i, atom in enumerate(ase_molecule):
            # will update only for cdx, smi, and csv formats.
            atom.charge = charge[i]
            atom.magmom = mult[i]

        optimizer = ase.optimize.BFGS(
            ase_molecule, trajectory="cmin_opt.traj", logfile="cmin.opt"
        )
        try:
            optimizer.run(fmax=self.args.opt_fmax, steps=self.args.opt_steps)

        except KeyError:
            self.args.log.write(f"\nx  {self.args.ani_method} could not optimize this molecule (i.e. check if all the atoms used are compatible with ANI)")
            cmin_valid = False
            energy = 0

        if cmin_valid:
            if len(ase.io.Trajectory("cmin_opt.traj", mode="r")) != (self.args.opt_steps + 1):
                species_coords = ase_molecule.get_positions().tolist()
                coordinates = torch.tensor(
                    [species_coords], requires_grad=True, device=DEVICE
                )

            # compute energy:
            species = model.species_to_tensor(elements).to(DEVICE).unsqueeze(0)
            _, ani_energy = model((species, coordinates))
            energy = ani_energy.item() * hartree_to_kcal  # Hartree to kcal/mol

            # update coordinates of mol object
            cartesians = np.array(coordinates.tolist()[0])
            for j in range(mol.GetNumAtoms()):
                [x, y, z] = cartesians[j]
                mol.GetConformer().SetAtomPosition(j, Point3D(x, y, z))

        return mol, energy, cmin_valid

    # generate the CMIN optimization model
    def get_cmin_model(self):
        """
        Function to generate the optimization model for CMIN (using xTB or ANI methods)
        """

        if self.args.program.lower() == "ani":
            try:
                import torchani
            except (ImportError,ModuleNotFoundError):
                self.args.log.write("x  Torchani is not installed! You can install the program with 'pip install torchani'")
                self.args.log.finalize()
                sys.exit()

            model = getattr(torchani.models,self.args.ani_method)()

        elif self.args.program.lower() == "xtb":
            try:
                subprocess.run(
                    ["xtb", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
            except FileNotFoundError:
                self.args.log.write("x  xTB is not installed (CREST cannot be used)! You can install the program with 'conda install -c conda-forge xtb'")
                self.args.log.finalize()
                sys.exit()
    
            model = None

        return model

    # write SDF files for xTB and ANI
    def write_confs(self, conformers, selectedcids, log):
        if len(conformers) > 0:
            write_confs = 0
            for cid in selectedcids:
                self.sdwriter.write(conformers[cid])
                write_confs += 1

            self.sdwriter.close()
        else:
            log.write("x  No conformers found!")


    def charge_mult_cmin(self, dup_data, dup_data_idx, type_cmin):
        """
        Retrieves charge and multiplicity arrays for ANI optimizations.

        Parameters
        ----------
        mol_cmin : Mol object
            Mol used to analyze charges and multiplicity of each atom
        dup_data : pandas.DataFrame
            Dataframe with the information of the molecules studied (i.e. number of conformers, charge, mult, etc.)
        dup_data_idx : int
            Index of the molecule studied in the dup_data dataframe

        Returns
        -------
        charge : list
            List of atomic charges
        mult : list
            List of atomic unpaired electrons
        dup_data : pandas.DataFrame
            Updated dup_data dataframe
        """

        mol_cmin = self.mols[0]

        if type_cmin == 'xtb':
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

            # assigns idx to metal atoms
            self.args.metal_idx, self.args.complex_coord, self.args.metal_sym = substituted_mol(self, mol_cmin, "noI")
            charge, metal_found = rules_get_charge(mol_cmin, self.args)
            mult = None
            final_mult = Descriptors.NumRadicalElectrons(mol_cmin) + 1
            if metal_found:
                # since RDKit gets the multiplicity of the metal with valence 0, the real multiplicity
                # value needs to be adapted with the charge. If multiplicity is different than 1 or 2,
                # the user must specify the value with the mult option
                if (charge % 2) == 1 and charge != 0: # odd charges (i.e. +1, +3, etc)
                    if final_mult == 1:
                        final_mult = final_mult + 1
                    if final_mult == 2:
                        final_mult = final_mult - 1

        elif type_cmin == 'ani':
            charge = []
            mult = []
            for _, atom in enumerate(mol_cmin.GetAtoms()):
                charge.append(atom.GetFormalCharge())
                mult.append(atom.GetNumRadicalElectrons())
            TotalElectronicSpin = np.sum(mult) / 2
            final_mult = int((2 * TotalElectronicSpin) + 1)

        dup_data.at[dup_data_idx, "Overall charge"] = np.sum(charge)
        dup_data.at[dup_data_idx, "Mult"] = final_mult

        return charge,mult,final_mult,dup_data