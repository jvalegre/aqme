"""
Parameters
----------

General
+++++++

   files : str or list of str, default=None
     Input files. Formats accepted: XYZ, SDF, GJF, COM and PDB. Also, lists can
     be used (i.e. [FILE1.sdf, FILE2.sdf] or \*.FORMAT such as \*.sdf).  
   program : str, default=None
     Program required in the conformational refining. 
     Current options: 'xtb', 'ani'
   w_dir_main : str, default=os.getcwd()
     Working directory  
   destination : str, default=None,
     Directory to create the output file(s)  
   varfile : str, default=None
     Option to parse the variables using a yaml file (specify the filename)  
   nprocs : int, default=None
     Number of processors used in the xTB optimizations  
   charge : int, default=None
     Charge of the calculations used in the xTB calculations. If charge isn't 
     defined, it automatically reads the charge from the input SDF files 
     (if the files come from CSEARCH, which adds the property "Real charge") 
     or calculates it from the generated mol object  
   mult : int, default=None
     Multiplicity of the calculations used in the xTB calculations. If charge 
     isn't defined, it automatically reads the charge from the input SDF files 
     (if the files come from CSEARCH, which adds the property "Mult") or 
     calculates it from the generated mol object. Be careful with the automated 
     calculation of mult from mol objects when using metals!  
   ewin_cmin : float, default=5.0
     Energy window in kcal/mol to discard conformers (i.e. if a conformer is 
     more than the E window compared to the most stable conformer)  
   initial_energy_threshold : float, default=0.0001
     Energy difference in kcal/mol between unique conformers for the first 
     filter of only E  
   energy_threshold : float, default=0.25
     Energy difference in kcal/mol between unique conformers for the second 
     filter of E + RMS  
   rms_threshold : float, default=0.25
     RMS difference between unique conformers for the second filter of E + RMS  
   stacksize : str, default='1G'
     Controls the stack size used (especially relevant for xTB/CREST 
     calculations of large systems, where high stack sizes are needed)
   prefix : str, default=''  
      Prefix added to all the names  
   suffix : str, default=''  
      Suffix added to all the names  

xTB only
++++++++

   xtb_keywords : str, default=None
     Define additional keywords to use in xTB that are not included in -c, 
     --uhf, -P and --input. For example: '--alpb ch2cl2 --gfn 1'
   constraints_atoms : list, default=[]
     Specify constrained atoms as [AT1,AT2,AT3]. An example of multiple constraints with
     atoms 1, 2 and 5 frozen: [1,2,5]
   constraints_dist : list of lists, default=[]
     Specify distance constraints as [AT1,AT2,DIST]. An example of multiple constraints with
     atoms 1 and 2 frozen at a distance of 1.8 Å, and atoms 4 and 5 with distance of 2.0 Å:
     [[1,2,1.8],[4,5,2.0]]
   constraints_angle : list of lists, default=[]
     Specify angle constraints as [AT1,AT2,AT3,ANGLE]. An example of multiple constraints with
     atoms 1, 2 and 3 frozen at an angle of 180 degrees, and atoms 4, 5 and 6 with an angle of 120:
     [[1,2,3,180],[4,5,6,120]]
   constraints_dihedral : list of lists, default=[]
     Specify dihedral constraints as [AT1,AT2,AT3,AT4,DIHEDRAL]. An example of multiple constraints
     with atoms 1, 2, 3 and 4 frozen at a dihedral angle of 180 degrees, and atoms 4, 5, 6 and 7
     with a dihedral angle of 120: [[1,2,3,4,180],[4,5,6,7,120]]

ANI only
++++++++

   opt_steps : int, default=1000
     Maximum number of steps used in the ase.optimize.BFGS optimizer.  
   opt_fmax : float, default=0.05
     Maximum force value to determine convergence in the ase.optimize.BFGS optimizer.  
   ani_method : str, default='ANI2x'
     ANI model used in the ase.optimize.BFGS optimizer.  
"""
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
import time
from aqme.utils import (
    load_variables,
    mol_from_sdf_or_mol_or_mol2,
    add_prefix_suffix,
    check_xtb,
    check_dependencies,
    set_destination
)
from aqme.filter import conformer_filters
from aqme.csearch.crest import xtb_opt_main
from aqme.csearch.utils import prepare_com_files

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

        # check whether dependencies are installed
        _ = check_dependencies(self)

        # most users employ xTB to refine RDKit structures
        if self.args.program is None:
            self.args.program = "xtb"
        elif self.args.program.lower() not in ["xtb", "ani"]:
            self.args.log.write('\nx  Program not supported for CMIN refinement! Specify: program="xtb" (or "ani")')
            self.args.log.finalize()
            sys.exit()

        # set number of processors
        if self.args.nprocs is None:
            self.args.nprocs = 1

        # check if xTB is installed
        if self.args.program.lower() == "xtb":
            _ = check_xtb(self)

        # retrieves the different files to run in CMIN
        if len(self.args.files) == 0:
            self.args.log.write('\nx  No files were found! Make sure you use quotation marks if you are using * (i.e. --files "*.sdf")')
            self.args.log.finalize()
            sys.exit()

        bar = IncrementalBar(
            "\no  Number of finished jobs from CMIN", max=len(self.args.files)
        )

        file_format = os.path.basename(Path(self.args.files[0])).split('.')[-1]

        if file_format.lower() in ['xyz', 'gjf', 'com']:
            for file in self.args.files:
                prepare_com_files(self.args, file)
            if file_format.lower() in ['gjf', 'com']:
                files_temp_extra = glob.glob('*.xyz')
            files_cmin = glob.glob('*.sdf')
        elif file_format.lower() == 'pdb':
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
        elif file_format.lower() == 'sdf':
            files_cmin = self.args.files
        else:
            self.args.log.write(f"\nx  The input format {file_format} is not supported for CMIN refinement! Formats allowed: SDF, XYZ, COM, GJF and PDB")
            self.args.log.finalize()
            sys.exit()

        for file in files_cmin:
            # load jobs for cmin minimization
            self.mols, self.name = self.load_jobs(file)
            self.name = add_prefix_suffix(self.name, self.args)

            self.args.log.write(f"\n\n   ----- {self.name} -----")

            self.cmin_folder = set_destination(self,'CMIN')

            # two destinations, one to store all the confs optimized and another to store
            # only the confs that passed the duplicate and energy filters
            self.cmin_folder.mkdir(exist_ok=True, parents=True)
            self.cmin_folder.joinpath('All_confs').mkdir(exist_ok=True, parents=True)

            self.cmin_all_file = self.cmin_folder.joinpath(
                f"All_confs/{self.name}_{self.args.program.lower()}_all_confs{self.args.output}"
            )
            self.sdwriterall = Chem.SDWriter(str(self.cmin_all_file))

            self.cmin_file = self.cmin_folder.joinpath(
                self.name + "_" + self.args.program.lower() + self.args.output
            )
            self.sdwriter = Chem.SDWriter(str(self.cmin_file))

            # run the optimizations
            _ = self.compute_cmin(file)

            bar.next()
        bar.finish()

        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime CMIN: {elapsed_time} seconds\n")
        self.args.log.finalize()

        # delete extra temporary files created when using XYZ, GJF, COM and PDB files
        if file_format.lower() in ['xyz', 'gjf', 'com', 'pdb']:
            if file_format.lower() in ['gjf', 'com']:
                files_cmin = files_cmin + files_temp_extra
            for temp_file in files_cmin:
                os.remove(temp_file)
        
        # removes systems that did not generate any conformers
        sdf_files_created = glob.glob(f'{self.cmin_folder}/*.sdf') + glob.glob(f'{self.cmin_folder.joinpath("All_confs")}/*.sdf')
        for sdf_file in sdf_files_created:
            if os.path.getsize(sdf_file) == 0:
                os.remove(sdf_file)

        # this is added to avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

    def load_jobs(self, file):

        # read SDF files
        try:
            inmols = mol_from_sdf_or_mol_or_mol2(file, 'cmin', self.args)
        except OSError:
            file_path = Path(self.args.initial_dir).joinpath(file)
            file_path = file_path.as_posix()
            inmols = mol_from_sdf_or_mol_or_mol2(file_path, 'cmin', self.args)

        name_mol = os.path.basename(file).split(".sdf")[0]

        return inmols, name_mol

    def compute_cmin(self, file):

        cenergy, outmols = [], []

        if self.args.program.lower() == "ani":
            if self.args.charge is not None:
                self.args.log.write("\nx  Charge is automatically calculated for ANI methods, do not use the charge option!")
                self.args.log.finalize()
                sys.exit()
            elif self.args.mult is not None:
                self.args.log.write("\nx  Multiplicity is automatically calculated for ANI methods, do not use the mult option!")
                self.args.log.finalize()
                sys.exit()
            charge,mult,final_mult = self.charge_mult_cmin()

        elif self.args.program.lower() == "xtb":
            # sets charge and mult
            file_format = os.path.basename(Path(file)).split('.')[-1]
            charge_input, mult_input, final_mult = None, None, None
            if file_format.lower() == 'sdf':
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
                charge = 0
                self.args.log.write(f'nx  No charge was assigned! Setting a value of 0, it can be changed with the charge option (or column in CSV inputs).')
            elif self.args.charge is None:
                charge = charge_input
            else:
                charge = self.args.charge
            if self.args.mult is None and mult_input is None:
                mult = 1
                self.args.log.write(f'nx  No multiplicity was assigned! Setting a value of 1, it can be changed with the mult option (or column in CSV inputs).')
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
                    name_init = mol.GetProp('_Name')
                    mol, energy, cmin_valid = xtb_opt_main(
                        f'{self.name}_conf_{i}',
                        self,
                        charge,
                        mult,
                        None,
                        self.args.constraints_atoms,
                        self.args.constraints_dist,
                        self.args.constraints_angle,
                        self.args.constraints_dihedral,
                        'xtb',
                        self.args.geom,
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

            for cid in sorted_all_cids:
                self.sdwriterall.write(outmols[cid])
            self.sdwriterall.close()

            self.args.log.write(f"\no  Applying filters to initial conformers after {self.args.program.lower()} minimization")
            selectedcids = conformer_filters(self,sorted_all_cids,cenergy,outmols)

            # write the filtered, ordered conformers to external file
            self.write_confs(
                outmols, selectedcids, self.args.log
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

    # ANI MAIN OPTIMIZATION PROCESS
    def ani_optimize(self, mol, charge, mult):
        
        import torch
        import ase
        self.args.log.write(f"\no  Starting ANI optimization")

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

        import torchani
        model = getattr(torchani.models,self.args.ani_method)()

        return model

    # write SDF files for xTB and ANI
    def write_confs(self, conformers, selectedcids, log):
        if len(conformers) > 0:
            for cid in selectedcids:
                self.sdwriter.write(conformers[cid])
            self.sdwriter.close()
        else:
            log.write("x  No conformers found!")


    def charge_mult_cmin(self):
        """
        Retrieves charge and multiplicity arrays for ANI optimizations.

        Parameters
        ----------
        mol_cmin : Mol object
            Mol used to analyze charges and multiplicity of each atom

        Returns
        -------
        charge : list
            List of atomic charges
        mult : list
            List of atomic unpaired electrons
        """

        mol_cmin = self.mols[0]


        charge = []
        mult = []
        for _, atom in enumerate(mol_cmin.GetAtoms()):
            charge.append(atom.GetFormalCharge())
            mult.append(atom.GetNumRadicalElectrons())
        TotalElectronicSpin = np.sum(mult) / 2
        final_mult = int((2 * TotalElectronicSpin) + 1)

        return charge,mult,final_mult