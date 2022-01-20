#####################################################.
#      This file stores all the functions used       #
#    in conformer minimization with xTB and ANI      #
#####################################################.

from operator import itemgetter
import os
import sys
import numpy as np
from pathlib import Path
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.PropertyMol import PropertyMol
from rdkit.Geometry import Point3D
import time
from aqme.argument_parser import set_options
from aqme.filter import CompoundFilter, EnergyFilter, RMSDFilter
from aqme.utils import (
    set_metal_atomic_number,
    Logger,
    rules_get_charge,
    substituted_mol,
    creation_of_dup_csv_cmin,
    load_from_yaml,
)

hartree_to_kcal = 627.509


class cmin:
    """
    Representation of the neccesary information related with cmin.

    Parameters
    ----------
    mol : RDKit Mol object, needed
            SMILES string necessary for setting up csearch object
    name : Name of smiles, needed
            string representing the code name for the object.
    """

    def __init__(
        self,
        mols=None,
        name=None,
        w_dir_initial=os.getcwd(),
        varfile=None,
        charge_defualt=0,
        **kwargs,
    ):
        self.mols = mols
        self.name = name
        self.w_dir_initial = w_dir_initial
        self.charge_default = charge_defualt

        if "options" in kwargs:
            self.args = kwargs["options"]
        else:
            self.args = set_options(kwargs)
        self.args.varfile = varfile

        dat_dir = Path(self.w_dir_initial + "/CMIN/dat_files")
        dat_dir.mkdir(exist_ok=True, parents=True)
        self.log = Logger(dat_dir / self.name, self.args.output_name)

        if varfile is not None:
            self.args, self.log = load_from_yaml(self.args, self.log)

        self.args.charge_default = self.charge_default
        self.args.charge = rules_get_charge(self.mols[0], self.args)

        self.program = self.args.CMIN

        self.cmin_folder = Path(self.w_dir_initial).joinpath(f"CMIN/{self.args.CMIN}")
        self.cmin_folder.mkdir(exist_ok=True)
        self.cmin_all_file = self.cmin_folder.joinpath(
            f"{self.name}_{self.program}_all_confs{self.args.output}"
        )
        self.sdwriterall = Chem.SDWriter(str(self.cmin_all_file))

        self.cmin_file = self.cmin_folder.joinpath(self.name + "_" + self.program + self.args.output)
        self.sdwriter = Chem.SDWriter(str(self.cmin_file))

    def compute_cmin(self):

        dup_data = creation_of_dup_csv_cmin(self.args.CMIN)
        dup_data_idx = 0
        dup_data.at[dup_data_idx, "Molecule"] = self.name
        cenergy, outmols = [], []
        start_time = time.time()
        for i, mol in enumerate(self.mols):

            if mol is not None:
                # optimize this structure and record the energy
                if self.args.metal_complex:
                    # fill the lists with None for every metal in the option
                    for _ in args.metal:
                        self.args.metal_idx.append(None)
                        self.args.complex_coord.append(None)
                        self.args.metal_sym.append(None)

                    (
                        self.args.mol,
                        self.args.metal_idx,
                        self.args.complex_coord,
                        self.args.metal_sym,
                    ) = substituted_mol(mol, self.args)

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

        # cmin_folder = Path(self.w_dir_initial).joinpath(f"CMIN/{self.args.CMIN}")
        # cmin_folder.mkdir(exist_ok=True)
        # cmin_file = cmin_folder.joinpath(
        #     f"{name_mol}_{self.program}_all_confs{self.args.output}"
        # )
        # sdwriter = Chem.SDWriter(str(cmin_file))
        # writing all conformers to files after minimization
        # sdwriter = Chem.SDWriter(f'{name_mol}_{program}_all_confs{args.output}')

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
            f"\n\no  Applying filters to intial conformers after {self.program} minimization"
        )

        # THIS IS THE VERSION FROM RAUL, NEED TO ADAPT THIS!
        # # filter based on energy window ewin_csearch
        # Ewindow_filter = EnergyFilter(args.ewin_csearch,'window')
        # Ediff_filter1 = EnergyFilter(args.initial_energy_threshold,'difference')
        # Ediff_filter2 = EnergyFilter(args.energy_threshold,'difference')
        # RMSD_filter = RMSDFilter(threshold=args.rms_threshold,
        #                          maxmatches=args.max_matches_RMSD,
        #                          heavyonly=args.heavyonly,
        #                          is_rdkit=program=='rdkit')
        # Ediff_RMSD = CompoundFilter(Ediff_filter2,RMSD_filter)

        # items = [(cid,energy,mol) for cid,energy,mol in zip(sorted_all_cids,cenergy,outmols)]
        # keys = [itemgetter(1),
        #         itemgetter(1),
        #         [itemgetter(1),itemgetter(0,2)]]
        # total_filter = CompoundFilter(Ewindow_filter,Ediff_filter1,Ediff_RMSD)

        # total_filter.apply(items,keys=keys)

        # filter_to_pandas(total_filter,dup_data,dup_data_idx,program)

        # selectedcids, cenergy, outmols = zip(*total_filter.accepted)

        from aqme.filter import ewin_filter, pre_E_filter, RMSD_and_E_filter

        sortedcids = ewin_filter(
            sorted_all_cids,
            cenergy,
            self.args,
            dup_data,
            dup_data_idx,
            self.log,
            self.program,
            self.args.ewin_csearch,
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
            cenergy,
            selectedcids,
            name_mol,
            self.args,
            self.program,
            self.log,
            self.cmin_folder,
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

        import ase
        import ase.optimize
        from ase.units import Hartree
        import torch

        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        DEVICE = torch.device("cpu")

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
                    atom.charge = args.charge_default
                    if args.verbose:
                        log.write("o  The Overall charge is read from the .com file ")

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

        # REPEATED IMPORTS JUST AS PLACEHOLDERS!
        import ase
        import ase.optimize
        from ase.units import Hartree
        import torchani
        import torch

        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        DEVICE = torch.device("cpu")

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

        # imports for xTB and ANI
        try:
            import ase
            import ase.optimize
            from ase.units import Hartree

        except (ModuleNotFoundError, AttributeError):
            err_msg = "ASE is not installed correctly - xTB and ANI are not available"
            log.write(
                "\nx  ASE is not installed correctly - xTB and ANI are not available"
            )
            raise ModuleNotFoundError(err_msg)

        try:
            import torch

            os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
            DEVICE = torch.device("cpu")

        except (ModuleNotFoundError, AttributeError):
            err_msg = "TORCH is not installed correctly - xTB and ANI are not available"
            log.write(
                "\nx  TORCH is not installed correctly - xTB and ANI are not available"
            )
            raise ModuleNotFoundError(err_msg)

        # Attempt an XTB import, if it fails log it and mock xtb_calc to delay the
        # system exit until it is used.
        try:
            from xtb.ase.calculator import XTB
        except (ModuleNotFoundError, AttributeError):
            log.write("\nx  xTB is not installed correctly - xTB is not available")
            xtb_calc = lambda *x, **y: sys.exit()
        # Attempt a torchani import, if it fails log it and mock ani_calc function
        # to raise a sys.exit() if called
        try:
            import torchani

        except (ModuleNotFoundError, AttributeError):
            log.write("\nx  Torchani is not installed correctly - ANI is not available")
            ani_calc = lambda *x, **y: sys.exit()

        # if large system increase stack size
        if args.STACKSIZE != "1G":
            os.environ["OMP_STACKSIZE"] = args.STACKSIZE

        # removing the Ba atom if NCI complexes
        if args.nci_complex:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "I":
                    atom.SetAtomicNum(1)

        if args.metal_complex and not args.CSEARCH == "summ":
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
                # VERY DIRTY HACK! WE NEED TO FIX THE IMPORTS THROUGH CMIN!
                # from aqme.cmin_bug import ani_calc
                energy, coordinates = self.ani_calc(elements, coordinates, args)

            except KeyError:
                log.write(
                    f"\nx  {args.ani_method} could not optimize this molecule (i.e. check of atoms that are not compatible)"
                )
                ani_incompatible = True
                energy = 0
                coordinates = np.zeros((1, 3))

        elif program == "xtb":
            # VERY DIRTY HACK! WE NEED TO FIX THE IMPORTS THROUGH CMIN!
            # from aqme.cmin_bug import xtb_calc
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

    # FUNCTIONS TO FIX RAUL'S COMMIT, REORGANIZE THEM!
    # WRITE SDF FILES FOR xTB AND ANI1
    def write_confs(
        self, conformers, energies, selectedcids, name, args, program, log, cmin_folder
    ):
        if len(conformers) > 0:
            # name = name.split('_'+args.CSEARCH)[0]# a bit hacky
            # cmin_file2 = cmin_folder.joinpath(name + "_" + program + args.output)
            # sdwriter = Chem.SDWriter(str(cmin_file2))
            # sdwriter = Chem.SDWriter(name+'_'+program+args.output)

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

    def filter_to_pandas(compfilter, dataframe, row, program):
        """
        Writes the results of a filter in a dataframe at the specified row.

        Parameters
        ----------
        compfilter : CompoundFilter
                A filter with a dataset != None.
        dataframe : pd.Dataframe
                The dataframe where the results are to be added.
        row : int
                row of the dataframe where the results are to be stored.
        program : str
                program for the conformer search. ['rdkit','summ','ani','xtb']
        """
        columns = [
            "energy-window",
            "initial_energy_threshold",
            "RMSD-and-energy-duplicates",
        ]

        program2name = {"rdkit": "RDKit", "summ": "summ", "ani": "ANI", "xtb": "xTB"}

        prog = program2name[program]
        for i, col in enumerate(columns):
            duplicates = compfilter.discarded_from(i)
            dataframe.at[row, f"{prog}-{col}"] = len(duplicates)
        dataframe.at[row, f"{prog}-Unique-conformers"] = len(compfilter.accepted)


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


def mult_min(name, args, program, charge, log, w_dir_initial):
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
        if args.CMIN == "xtb":
            if args.xtb_solvent == "none":
                method = f"xTB ({args.xtb_method}"
            else:
                method = f"xTB ({args.xtb_method} in {args.xtb_solvent})"
        if args.CMIN == "ani":
            method = f"ANI ({args.ani_method})"
        if args.CMIN in ["xtb", "ani"]:
            filename = name + args.output
            log.write(f"\n\no  Multiple minimization of {filename} with {method}")

    # bar = IncrementalBar('o  Minimizing', max = len(inmols))

    # read SDF files from RDKit optimization
    inmols = rdkit_sdf_read(name, args, log)

    name_mol = os.path.basename(name).split("_" + args.CSEARCH)[0]

    # bar.next()
    obj = cmin(inmols, name_mol, w_dir_initial, args.varfile, charge, program)
    total_data = obj.compute_cmin()

    return total_data
