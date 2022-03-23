#####################################################.
#          This file storesthe CSEARCH class        #
#             used in conformer refinement          #
#####################################################.

import os
import sys
import numpy as np
from pathlib import Path
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.PropertyMol import PropertyMol
from rdkit.Geometry import Point3D
import time
import ase
import ase.optimize
from ase.units import Hartree
import torch
import torchani
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
DEVICE = torch.device("cpu")
from aqme.argument_parser import set_options
from aqme.utils import (
    set_metal_atomic_number,
    Logger,
    rules_get_charge,
    substituted_mol,
    creation_of_dup_csv_cmin,
    load_from_yaml,
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

    def __init__(
        self,
        mols=None,
        name=None,
        w_dir_main=os.getcwd(),
        varfile=None,
        charge_default=0,
        program=None,
        **kwargs,
    ):

        self.mols = mols
        self.name = name
        self.w_dir_main = w_dir_main
        self.charge_default = charge_default
        self.program = program

        if "options" in kwargs:
            self.args = kwargs["options"]
        else:
            self.args = set_options(kwargs)
        self.args.varfile = varfile

        dat_dir = Path(self.w_dir_main + "/CMIN/dat_files")
        dat_dir.mkdir(exist_ok=True, parents=True)
        self.log = Logger(dat_dir / self.name, self.args.output_name)

        if varfile is not None:
            self.args, self.log = load_from_yaml(self.args, self.log)

        self.args.charge = charge_default
        self.args.charge = rules_get_charge(self.mols[0], self.args)

        self.program = self.args.cmin

        cmin_folder = Path(self.w_dir_main).joinpath(f"CMIN/{self.args.cmin}")
        cmin_folder.mkdir(exist_ok=True, parents=True)
        self.cmin_all_file = cmin_folder.joinpath(
            f"{self.name}_{self.program}_all_confs{self.args.output}"
        )
        self.sdwriterall = Chem.SDWriter(str(self.cmin_all_file))

        self.cmin_file = cmin_folder.joinpath(self.name + "_" + self.program + self.args.output)
        self.sdwriter = Chem.SDWriter(str(self.cmin_file))

    def compute_cmin(self):

        dup_data = creation_of_dup_csv_cmin(self.args.cmin)
        dup_data_idx = 0
        dup_data.at[dup_data_idx, "Molecule"] = self.name
        cenergy, outmols = [], []
        start_time = time.time()
        for _,mol in enumerate(self.mols):
            if mol is not None:
                # optimize this structure and record the energy
                if self.args.metal_complex:
                    # fill the lists with None for every metal in the option
                    for _ in self.args.metal:
                        self.args.metal_idx.append(None)
                        self.args.complex_coord.append(None)
                        self.args.metal_sym.append(None)

                    (
                        self.args.metal_idx,
                        self.args.complex_coord,
                        self.args.metal_sym,
                    ) = substituted_mol(self,mol)

                self.mol, energy, ani_incompatible = self.optimize(
                    mol, self.args, self.program, self.log, dup_data, dup_data_idx
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
                "_Name", outmols[cid].GetProp("_Name") + " " + self.program
            )
            outmols[cid].SetProp("Energy", cenergy[cid])

        write_all_confs = 0
        for cid in sorted_all_cids:
            self.sdwriterall.write(outmols[cid])
            write_all_confs += 1
        self.sdwriterall.close()

        if self.args.verbose:
            self.log.write(
                f"\no  Writing {str(write_all_confs)} conformers to file {name_mol}_{self.program}_all_confs{self.args.output}"
            )

        self.log.write(
            f"\no  Applying filters to intial conformers after {self.program} minimization"
        )

        # filter based on energy window ewin_cmin
        sortedcids = ewin_filter(
            sorted_all_cids,
            cenergy,
            self.args,
            dup_data,
            dup_data_idx,
            self.log,
            self.program,
            self.args.ewin_cmin,
        )
        # pre-filter based on energy only
        selectedcids_initial = pre_E_filter(
            sortedcids,
            cenergy,
            dup_data,
            dup_data_idx,
            self.log,
            self.program,
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
            self.log,
            self.program,
        )

        if self.program == "xtb":
            dup_data.at[dup_data_idx, "xTB-Initial-samples"] = len(self.mols)
        elif self.program == "ani":
            dup_data.at[dup_data_idx, "ANI-Initial-samples"] = len(self.mols)

        # write the filtered, ordered conformers to external file
        self.write_confs(
            outmols,
            selectedcids,
            name_mol,
            self.args,
            self.program,
            self.log
        )
        dup_data.at[dup_data_idx, "CMIN time (seconds)"] = round(
            time.time() - start_time, 2
        )
        return dup_data

    def xtb_calc(self, elements, coordinates, args, log, ase_metal, ase_metal_idx):
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
        ase_metal : [type]
                [description]
        ase_metal_idx : [type]
                [description]

        Returns
        -------
        tuple
                sqm_energy, coordinates
        """

        from xtb.ase.calculator import XTB

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
        if args.metal_complex:
            extension = os.path.splitext(args.input)[1]
            if extension in [".csv", ".cdx", ".smi"]:
                for i, atom in enumerate(ase_molecule):
                    if i in ase_metal:
                        ase_charge = args.charge[
                            args.metal_idx.index(ase_metal_idx[ase_metal.index(i)])
                        ]
                        # will update only for cdx, smi, and csv formats.
                        atom.charge = ase_charge
            else:
                atom.charge = args.charge
                if args.verbose:
                    log.write("o  The Overall charge is read from the input file ")

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
    def optimize(self, mol, args, program, log, dup_data, dup_data_idx):

        # Attempt an XTB import, if it fails log it and mock xtb_calc to delay the
        # system exit until it is used.
        try:
            from xtb.ase.calculator import XTB
        except (ModuleNotFoundError, AttributeError):
            log.write("\nx  xTB is not installed correctly - xTB is not available")
            print("\nx  xTB is not installed correctly - xTB is not available")
            sys.exit()

        # if large system increase stack size
        if args.stacksize != "1G":
            os.environ["OMP_STACKSIZE"] = args.stacksize

        if args.metal_complex and not args.program == "summ":
            set_metal_atomic_number(mol, args.metal_idx, args.metal_sym)

        elements = ""
        ase_metal = []
        ase_metal_idx = []
        for i, atom in enumerate(mol.GetAtoms()):
            if atom.GetIdx() in args.metal_idx:
                ase_metal.append(i)
                ase_metal_idx.append(atom.GetIdx())
            elements += atom.GetSymbol()

        extension = os.path.splitext(args.input)[1]
        if extension in [".cdx", ".smi", ".csv"]:
            args.charge = rules_get_charge(mol, args)
            # replace None values if there are metals that are not used
            for i, charge in enumerate(args.charge):
                if charge is None:
                    args.charge[i] = 0
            dup_data.at[dup_data_idx, "Overall charge"] = np.sum(args.charge)
        else:
            dup_data.at[dup_data_idx, "Overall charge"] = np.sum(args.charge)

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
                elements, coordinates, args, log, ase_metal, ase_metal_idx
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
    def write_confs(
        self, conformers, selectedcids, name, args, program, log
    ):
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


def rdkit_sdf_read(name, args, log):
    """
    Reads sdf files and stops the execution if the file was not accesible.

    Parameters
    ----------
    name : [type]
            [description]
    args : argparse.args
            [description]
    log : Logger
            [description]

    Returns
    -------
    list
            rdkit.Chem.Mol objects
    """
    inmols = Chem.SDMolSupplier(name + args.output, removeHs=False)
    if inmols is None:
        log.write(f"Could not open {name}{args.output}")
        sys.exit(-1)
    return inmols


def mult_min(name, args, program, charge, log, w_dir_main):
    """
    READ FILES FOR xTB AND ANI OPTIMIZATION, FILTER AND WRITING SDF FILES

    Parameters
    ----------
    name : [type]
            [description]
    args : [type]
            [description]
    program : [type]
            [description]
    log : [type]
            [description]
    dup_data : [type]
            [description]
    dup_data_idx : [type]
            [description]
    """

    if args.verbose:
        if args.cmin == "xtb":
            if args.xtb_solvent == "none":
                method = f"xTB ({args.xtb_method}"
            else:
                method = f"xTB ({args.xtb_method} in {args.xtb_solvent})"
        if args.cmin == "ani":
            method = f"ANI ({args.ani_method})"
        if args.cmin in ["xtb", "ani"]:
            filename = name + args.output
            log.write(f"\no  Multiple minimization of {filename} with {method}")

    # bar = IncrementalBar('o  Minimizing', max = len(inmols))

    # read SDF files from RDKit optimization
    inmols = rdkit_sdf_read(name, args, log)

    name_mol = os.path.basename(name).split("_" + args.program)[0]

    # bar.next()
    obj = cmin(mols=inmols, name=name_mol, w_dir_main=w_dir_main, varfile=args.varfile,
            charge_default=charge, program=program, options=args)
    total_data = obj.compute_cmin()

    return total_data
