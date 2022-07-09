#####################################################.
#          This file stores the CMIN class          #
#             used in conformer refinement          #
#####################################################.

import os
import sys
import numpy as np
from pathlib import Path
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.PropertyMol import PropertyMol
from progress.bar import IncrementalBar
from rdkit.Geometry import Point3D
import pandas as pd
import time
try:
    import ase
except ModuleNotFoundError:
    print("x  ASE is not installed! You can install the program with 'conda install -c conda-forge ase' or 'pip install ase'")
    sys.exit()
import ase.optimize
from ase.units import Hartree
from aqme.utils import (
    rules_get_charge,
    load_variables,
    substituted_mol,
    creation_of_dup_csv_cmin,
)
from aqme.filter import ewin_filter, pre_E_filter, RMSD_and_E_filter

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
        if self.args.program not in ["xtb", "ani"]:
            self.args.log.write(
                "\nx  Program not supported for conformer generation! Specify: program='xtb' (or ani)"
            )
            self.args.log.finalize()
            sys.exit()

        try:
            os.chdir(self.args.w_dir_main)
        except FileNotFoundError:
            self.args.w_dir_main = Path(f"{os.getcwd()}/{self.args.w_dir_main}")
            os.chdir(self.args.w_dir_main)

        # create the dataframe to store the data
        self.final_dup_data = creation_of_dup_csv_cmin(self.args.program)

        bar = IncrementalBar(
            "o  Number of finished jobs from CMIN", max=len(self.args.files)
        )
        for file in self.args.files:
            # load jobs for cmin minimization
            self.mols, self.name = self.load_jobs(file)

            if self.args.destination is None:
                self.cmin_folder = Path(self.args.w_dir_main).joinpath(
                    f"CMIN/{self.args.program}"
                )
            else:
                if Path(f"{self.args.destination}").exists():
                    self.cmin_folder = Path(self.args.destination)
                else:
                    self.cmin_folder = Path(f"{os.getcwd()}/{self.args.destination}")

            self.cmin_folder.mkdir(exist_ok=True, parents=True)

            self.cmin_all_file = self.cmin_folder.joinpath(
                f"{self.name}_{self.args.program}_all_confs{self.args.output}"
            )
            self.sdwriterall = Chem.SDWriter(str(self.cmin_all_file))

            self.cmin_file = self.cmin_folder.joinpath(
                self.name + "_" + self.args.program + self.args.output
            )
            self.sdwriter = Chem.SDWriter(str(self.cmin_file))

            # runs the conformer sampling with multiprocessors
            total_data = self.compute_cmin()

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

        # this is added to avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

    def load_jobs(self, file):
        if self.args.verbose:
            if self.args.program == "xtb":
                if self.args.xtb_solvent == "none":
                    method = f"xTB ({self.args.xtb_method}"
                else:
                    method = f"xTB ({self.args.xtb_method} in {self.args.xtb_solvent})"
            if self.args.program == "ani":
                method = f"ANI ({self.args.ani_method})"
            if self.args.program in ["xtb", "ani"]:
                self.args.log.write(
                    f"\no  Multiple minimization of {file} with {method}"
                )

        # read SDF files from RDKit optimization
        inmols = rdkit_sdf_read(file, self.args)
        name_mol = os.path.basename(file).split(".sdf")[0]

        return inmols, name_mol

    def compute_cmin(self):

        dup_data = creation_of_dup_csv_cmin(self.args.program)
        dup_data_idx = 0
        dup_data.at[dup_data_idx, "Molecule"] = self.name
        cenergy, outmols = [], []
        start_time = time.time()

        charge, metal_found = rules_get_charge(self.mols[0], self.args, "cmin")
        mult = []
        for Atom in self.mols[0].GetAtoms():
            mult.append(Atom.GetNumRadicalElectrons())

        TotalElectronicSpin = np.sum(mult) / 2
        final_mult = int((2 * TotalElectronicSpin) + 1)

        dup_data.at[dup_data_idx, "Overall charge"] = np.sum(charge)
        dup_data.at[dup_data_idx, "Mult"] = final_mult

        for _, mol in enumerate(self.mols):
            if mol is not None:
                # optimize this structure and record the energy
                if self.args.metal_complex:
                    # fill the lists with None for every metal in the option
                    for _ in self.args.metal_atoms:
                        self.args.metal_idx.append(None)
                        self.args.complex_coord.append(None)
                        self.args.metal_sym.append(None)

                    (
                        self.args.metal_idx,
                        self.args.complex_coord,
                        self.args.metal_sym,
                    ) = substituted_mol(self, mol, "noI")

                self.mol, energy, ani_incompatible = self.optimize(
                    mol,
                    self.args,
                    self.args.program,
                    self.args.log,
                    dup_data,
                    dup_data_idx,
                    charge,
                    mult,
                )

                if not ani_incompatible:
                    pmol = PropertyMol(mol)
                    outmols.append(pmol)
                    cenergy.append(energy)

        # if SQM energy exists, overwrite RDKit energies and geometries
        cids = list(range(len(outmols)))
        sorted_all_cids = sorted(cids, key=lambda cid: cenergy[cid])

        name_mol = self.name

        for cid in sorted_all_cids:
            outmols[cid].SetProp(
                "_Name", outmols[cid].GetProp("_Name") + " " + self.args.program
            )
            outmols[cid].SetProp("Energy", cenergy[cid])
            if self.args.charge is None:
                outmols[cid].SetProp("Real charge", str(np.sum(charge)))
            else:
                outmols[cid].SetProp("Real charge", str(self.args.charge))
            if self.args.mult is None:
                outmols[cid].SetProp("Mult", str(mult))
            else:
                outmols[cid].SetProp("Mult", str(self.args.mult))

        write_all_confs = 0
        for cid in sorted_all_cids:
            self.sdwriterall.write(outmols[cid])
            write_all_confs += 1
        self.sdwriterall.close()

        if self.args.verbose:
            self.args.log.write(
                f"\no  Writing {str(write_all_confs)} conformers to file {name_mol}_{self.args.program}_all_confs{self.args.output}"
            )

        self.args.log.write(
            f"\no  Applying filters to intial conformers after {self.args.program} minimization"
        )

        # filter based on energy window ewin_cmin
        sortedcids = ewin_filter(
            sorted_all_cids,
            cenergy,
            self.args,
            dup_data,
            dup_data_idx,
            self.args.log,
            self.args.program,
            self.args.ewin_cmin,
        )
        # pre-filter based on energy only
        selectedcids_initial = pre_E_filter(
            sortedcids,
            cenergy,
            dup_data,
            dup_data_idx,
            self.args.log,
            self.args.program,
            self.args.initial_energy_threshold,
            self.args.verbose,
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
            self.args.program,
        )

        if self.args.program == "xtb":
            dup_data.at[dup_data_idx, "xTB-Initial-samples"] = len(self.mols)
        elif self.args.program == "ani":
            dup_data.at[dup_data_idx, "ANI-Initial-samples"] = len(self.mols)

        # write the filtered, ordered conformers to external file
        self.write_confs(
            outmols, selectedcids, name_mol, self.args, self.args.program, self.args.log
        )
        dup_data.at[dup_data_idx, "CMIN time (seconds)"] = round(
            time.time() - start_time, 2
        )
        return dup_data

    def xtb_calc(self, elements, coordinates, args, log, charge, mult):
        """
        Run an xtb optimization and return the optimized coordinates and the energy
        of the molecule.

        Parameters
        ----------
        elements : [type]
                                                                                                                                                                                                   [description]
        coordinates : [type]
                                                                                                                                                                                                   [description]
        args : argparse.args
                                                                                                                                                                                                   [description]
        log : Logger
                                                                                                                                                                                                   [description]

        Returns
        -------
        tuple
                                                                                                                                                                                                   sqm_energy, coordinates
        """

        from xtb.ase.calculator import XTB
        try:
            import torch
        except ModuleNotFoundError:
            print("x  Torch is not installed! You can install the program with 'pip install torch'")
            sys.exit()
        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        DEVICE = torch.device("cpu")

        xtb_calculator = XTB(
            method=args.xtb_method,
            accuracy=args.xtb_accuracy,
            electronic_temperature=args.xtb_electronic_temperature,
            max_iterations=args.xtb_max_iterations,
            solvent=args.xtb_solvent,
        )

        # define ase molecule using GFN2 Calculator
        ase_molecule = ase.Atoms(
            elements, positions=coordinates.tolist()[0], calculator=xtb_calculator
        )

        # Adjust the charge of the metal atoms
        for i, atom in enumerate(ase_molecule):
            # will update only for cdx, smi, and csv formats.
            atom.charge = charge[i]
            atom.magmom = mult[i]

        optimizer = ase.optimize.BFGS(
            ase_molecule, trajectory="xTB_opt.traj", logfile="xtb.opt"
        )

        optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)

        if len(ase.io.Trajectory("xTB_opt.traj", mode="r")) != (args.opt_steps + 1):
            species_coords = ase_molecule.get_positions().tolist()
            coordinates = torch.tensor(
                [species_coords], requires_grad=True, device=DEVICE
            )

        # Now let's compute energy:
        xtb_energy = ase_molecule.get_potential_energy()
        sqm_energy = (xtb_energy / Hartree) * hartree_to_kcal

        return sqm_energy, coordinates

    def ani_calc(self, elements, coordinates, args):
        """
        Run an ANI optimization and return the energy and optimized coordinates.

        Parameters
        ----------
        ase : [type]
                                                                                                                                                                                                   [description]
        torch : [type]
                                                                                                                                                                                                   [description]
        model : [type]
                                                                                                                                                                                                   [description]
        device : [type]
                                                                                                                                                                                                   [description]
        elements : [type]
                                                                                                                                                                                                   [description]
        coordinates : [type]
                                                                                                                                                                                                   [description]
        args : [type]
                                                                                                                                                                                                   [description]

        Returns
        -------
        tuple
                                                                                                                                                                                                   sqm_energy, coordinates
        """

        try:
            import torch
        except ModuleNotFoundError:
            print("x  Torch is not installed! You can install the program with 'pip install torch'")
            sys.exit()
        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        DEVICE = torch.device("cpu")
        # the first try/except avoids bugs related to having installed pytorch but not torch
        try:
            import torchani
        except ImportError:
            print("x  Torch is not installed! You can install the program with 'pip install torch'")
            sys.exit()
        try:
            import torchani
        except ModuleNotFoundError:
            print("x  TorchANI is not installed! You can install the program with 'conda install -c conda-forge torchani' or 'pip install torchani'")
            sys.exit()

        # Select the model
        ANI_method = args.ani_method
        if ANI_method == "ANI1x":
            model = torchani.models.ANI1x()
        if ANI_method == "ANI1ccx":
            model = torchani.models.ANI1ccx()
        if ANI_method == "ANI2x":
            model = torchani.models.ANI2x()
        if ANI_method == "ANI2ccx":
            model = torchani.models.ANI2ccx()
        if ANI_method == "ANI3x":
            model = torchani.models.ANI3x()
        if ANI_method == "ANI3ccx":
            model = torchani.models.ANI3ccx()
        # A bit Fancier
        # model = getattr(torchani.models,ANI_method)()

        species = model.species_to_tensor(elements).to(DEVICE).unsqueeze(0)
        _, ani_energy = model((species, coordinates))

        ase_molecule = ase.Atoms(
            elements, positions=coordinates.tolist()[0], calculator=model.ase()
        )

        optimizer = ase.optimize.BFGS(
            ase_molecule, trajectory="ANI1_opt.traj", logfile="ase.opt"
        )
        optimizer.run(fmax=args.opt_fmax, steps=args.opt_steps)
        if len(ase.io.Trajectory("ANI1_opt.traj", mode="r")) != (args.opt_steps + 1):
            species_coords = ase_molecule.get_positions().tolist()
            coordinates = torch.tensor(
                [species_coords], requires_grad=True, device=DEVICE
            )

        # Now let's compute energy:
        _, ani_energy = model((species, coordinates))
        sqm_energy = ani_energy.item() * hartree_to_kcal  # Hartree to kcal/mol

        return sqm_energy, coordinates

    # xTB AND ANI MAIN OPTIMIZATION PROCESS
    def optimize(self, mol, args, program, log, dup_data, dup_data_idx, charge, mult):

        # Attempt an XTB import, if it fails log it and mock xtb_calc to delay the
        # system exit until it is used.
        try:
            from xtb.ase.calculator import XTB
        except (ModuleNotFoundError, AttributeError):
            log.write("\nx  xTB is not installed correctly - xTB is not available")
            log.finalize()
            sys.exit()
        try:
            import torch
        except ModuleNotFoundError:
            print("x  Torch is not installed! You can install the program with 'pip install torch'")
            sys.exit()
        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        DEVICE = torch.device("cpu")

        # if large system increase stack size
        if args.stacksize != "1G":
            os.environ["OMP_STACKSIZE"] = args.stacksize

        elements = ""
        for _, atom in enumerate(mol.GetAtoms()):
            elements += atom.GetSymbol()

        cartesians = mol.GetConformers()[0].GetPositions()
        coordinates = torch.tensor(
            [cartesians.tolist()], requires_grad=True, device=DEVICE
        )
        ani_incompatible = False
        if program == "ani":
            try:
                energy, coordinates = self.ani_calc(elements, coordinates, args)

            except KeyError:
                log.write(
                    f"\nx  {args.ani_method} could not optimize this molecule (i.e. check of atoms that are not compatible)"
                )
                ani_incompatible = True
                energy = 0
                coordinates = np.zeros((1, 3))

        elif program == "xtb":
            energy, coordinates = self.xtb_calc(
                elements, coordinates, args, log, charge, mult
            )
        else:
            log.write(
                "x  Option not compatible with CMIN (check the available options)!"
            )

        cartesians = np.array(coordinates.tolist()[0])

        # update coordinates of mol object
        for j in range(mol.GetNumAtoms()):
            [x, y, z] = cartesians[j]
            mol.GetConformer().SetAtomPosition(j, Point3D(x, y, z))

        return mol, energy, ani_incompatible

    # WRITE SDF FILES FOR xTB AND ANI1
    def write_confs(self, conformers, selectedcids, name, args, program, log):
        if len(conformers) > 0:
            write_confs = 0
            for cid in selectedcids:
                self.sdwriter.write(conformers[cid])
                write_confs += 1

            if args.verbose:
                log.write(
                    "o  Writing "
                    + str(write_confs)
                    + " conformers to file "
                    + name
                    + "_"
                    + program
                    + args.output
                )
            self.sdwriter.close()
        else:
            log.write("x  No conformers found!")


def rdkit_sdf_read(file, args):
    """
    Reads sdf files and stops the execution if the file was not accesible.                                                                                                                                                                                      rdkit.Chem.Mol objects
    """
    inmols = Chem.SDMolSupplier(file, removeHs=False)

    if inmols is None:
        args.log.write(f"Could not open {file}")
        args.log.finalize()
        sys.exit()
    return inmols
