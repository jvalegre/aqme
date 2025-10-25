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

    def _validate_and_setup_program(self):
        """Validate program selection and set defaults.
        
        Sets program to 'xtb' if not specified, validates it's supported,
        and checks for xTB installation if needed.
        """
        if self.args.program is None:
            self.args.program = "xtb"
        elif self.args.program.lower() not in ["xtb", "ani"]:
            self.args.log.write('\nx  Program not supported for CMIN refinement! Specify: program="xtb" (or "ani")')
            self.args.log.finalize()
            sys.exit()
        
        # Set number of processors
        if self.args.nprocs is None:
            self.args.nprocs = 1
        
        # Check if xTB is installed
        if self.args.program.lower() == "xtb":
            _ = check_xtb(self)

    def _convert_input_files_to_sdf(self, file_format):
        """Convert input files to SDF format based on file type.
        
        Handles conversion for XYZ, GJF, COM, and PDB formats.
        Returns list of SDF files and any extra temporary files created.
        
        Args:
            file_format (str): Input file format extension
        
        Returns:
            tuple: (files_cmin, files_temp_extra) - SDF files and temporary files
        """
        files_temp_extra = []
        
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
            self.args.log.write(
                f"\nx  The input format {file_format} is not supported for CMIN refinement! "
                f"Formats allowed: SDF, XYZ, COM, GJF and PDB"
            )
            self.args.log.finalize()
            sys.exit()
        
        return files_cmin, files_temp_extra

    def _setup_output_directories(self):
        """Create output directories and initialize SDF writers.
        
        Creates CMIN folder and All_confs subfolder, then initializes
        SDWriter objects for both the full conformer set and filtered set.
        """
        self.cmin_folder = set_destination(self, 'CMIN')
        
        # Two destinations: all optimized confs and filtered confs
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

    def _cleanup_temporary_files(self, file_format, files_cmin, files_temp_extra):
        """Remove temporary files created during processing.
        
        Deletes conversion intermediates and empty SDF files.
        
        Args:
            file_format (str): Input file format
            files_cmin (list): List of SDF files created
            files_temp_extra (list): Additional temporary files
        """
        # Delete extra temporary files created when using XYZ, GJF, COM and PDB files
        if file_format.lower() in ['xyz', 'gjf', 'com', 'pdb']:
            if file_format.lower() in ['gjf', 'com']:
                files_cmin = files_cmin + files_temp_extra
            for temp_file in files_cmin:
                os.remove(temp_file)
        
        # Remove systems that did not generate any conformers
        sdf_files_created = (
            glob.glob(f'{self.cmin_folder}/*.sdf') + 
            glob.glob(f'{self.cmin_folder.joinpath("All_confs")}/*.sdf')
        )
        for sdf_file in sdf_files_created:
            if os.path.getsize(sdf_file) == 0:
                os.remove(sdf_file)

    def __init__(self, **kwargs):
        """Initialize CMIN conformer refinement workflow.
        
        Loads configuration, validates programs, converts input files,
        runs optimizations, and cleans up temporary files.
        
        Args:
            **kwargs: CMIN module arguments (see module documentation)
        """
        start_time_overall = time.time()
        
        # Load default and user-specified variables
        self.args = load_variables(kwargs, "cmin")
        
        # Check whether dependencies are installed
        _ = check_dependencies(self)
        
        # Validate and setup program
        self._validate_and_setup_program()
        
        # Validate input files exist
        if len(self.args.files) == 0:
            self.args.log.write(
                '\nx  No files were found! Make sure you use quotation marks '
                'if you are using * (i.e. --files "*.sdf")'
            )
            self.args.log.finalize()
            sys.exit()
        
        # Setup progress bar
        bar = IncrementalBar(
            "\no  Number of finished jobs from CMIN", max=len(self.args.files)
        )
        
        # Convert input files to SDF format
        file_format = os.path.basename(Path(self.args.files[0])).split('.')[-1]
        files_cmin, files_temp_extra = self._convert_input_files_to_sdf(file_format)
        
        # Process each file
        for file in files_cmin:
            # Load jobs for cmin minimization
            self.mols, self.name = self.load_jobs(file)
            self.name = add_prefix_suffix(self.name, self.args)
            
            self.args.log.write(f"\n\n   ----- {self.name} -----")
            
            # Setup output directories and writers
            self._setup_output_directories()
            
            # Run the optimizations
            _ = self.compute_cmin(file)
            
            bar.next()
        
        bar.finish()
        
        # Report timing
        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\nTime CMIN: {elapsed_time} seconds\n")
        self.args.log.finalize()
        
        # Cleanup temporary files
        self._cleanup_temporary_files(file_format, files_cmin, files_temp_extra)
        
        # Avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

    def load_jobs(self, file):
        """Load molecules from input file.
        
        Reads SDF file and extracts molecule objects and base name.
        Tries loading from current directory first, then from initial directory.
        
        Args:
            file (str): Path to SDF file
        
        Returns:
            tuple: (inmols, name_mol) - list of molecules and file basename
        """
        try:
            inmols = mol_from_sdf_or_mol_or_mol2(file, 'cmin', self.args)
        except OSError:
            file_path = Path(self.args.initial_dir).joinpath(file)
            file_path = file_path.as_posix()
            inmols = mol_from_sdf_or_mol_or_mol2(file_path, 'cmin', self.args)
        
        name_mol = os.path.basename(file).split(".sdf")[0]
        
        return inmols, name_mol

    def _read_charge_mult_from_sdf(self, file):
        """Extract charge and multiplicity from SDF file properties.
        
        Reads 'Real charge' and 'Mult' properties from SDF files
        (typically created by CSEARCH module).
        
        Args:
            file (str): Path to SDF file
        
        Returns:
            tuple: (charge_input, mult_input) or (None, None) if not found
        """
        charge_input, mult_input = None, None
        
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
        
        return charge_input, mult_input

    def _determine_charge_mult_for_xtb(self, file):
        """Determine charge and multiplicity for xTB calculations.
        
        Uses the following precedence:
        1. User-specified values (self.args.charge/mult)
        2. SDF file properties (from CSEARCH)
        3. Default values (0 for charge, 1 for mult)
        
        Args:
            file (str): Path to input SDF file
        
        Returns:
            tuple: (charge, mult, final_mult)
        """
        file_format = os.path.basename(Path(file)).split('.')[-1]
        charge_input, mult_input = None, None
        
        # Try to read from SDF file properties
        if file_format.lower() == 'sdf':
            if self.args.charge is None or self.args.mult is None:
                charge_input, mult_input = self._read_charge_mult_from_sdf(file)
        
        # Determine final charge
        if self.args.charge is None and charge_input is None:
            charge = 0
            self.args.log.write(
                'nx  No charge was assigned! Setting a value of 0, '
                'it can be changed with the charge option (or column in CSV inputs).'
            )
        elif self.args.charge is None:
            charge = charge_input
        else:
            charge = self.args.charge
        
        # Determine final multiplicity
        if self.args.mult is None and mult_input is None:
            mult = 1
            self.args.log.write(
                'nx  No multiplicity was assigned! Setting a value of 1, '
                'it can be changed with the mult option (or column in CSV inputs).'
            )
        elif self.args.mult is None:
            mult = mult_input
        else:
            mult = self.args.mult
        
        return charge, mult, None

    def _check_constraints_exist(self):
        """Check if any constraints are specified for xTB optimization.
        
        Returns:
            bool: True if any constraint type has at least one constraint
        """
        return (
            len(self.args.constraints_atoms) >= 1 or
            len(self.args.constraints_dist) >= 1 or
            len(self.args.constraints_angle) >= 1 or
            len(self.args.constraints_dihedral) >= 1
        )

    def _optimize_single_conformer(self, mol, i, charge, mult):
        """Optimize a single conformer using ANI or xTB.
        
        Args:
            mol: RDKit molecule object
            i (int): Conformer index
            charge: Charge value (int for xTB, list for ANI)
            mult: Multiplicity value (int for xTB, list for ANI)
        
        Returns:
            tuple: (mol, energy, cmin_valid)
        """
        if self.args.program.lower() == "ani":
            return self.ani_optimize(mol, charge, mult)
        
        elif self.args.program.lower() == "xtb":
            complex_ts = self._check_constraints_exist()
            name_init = mol.GetProp('_Name')
            
            return xtb_opt_main(
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
                self.args.sample,
                complex_ts=complex_ts,
                mol=mol,
                name_init=name_init
            )

    def _add_molecule_properties(self, outmols, sorted_all_cids, cenergy, charge, mult, final_mult):
        """Add energy and charge/mult properties to optimized molecules.
        
        Updates molecule properties with energy, charge, and multiplicity
        based on the optimization program used.
        
        Args:
            outmols (list): List of optimized molecule objects
            sorted_all_cids (list): Conformer IDs sorted by energy
            cenergy (list): Conformer energies
            charge: Charge value (list for ANI, int for xTB)
            mult: Multiplicity value (list for ANI, int for xTB)
            final_mult: Final multiplicity (for ANI)
        """
        for cid in sorted_all_cids:
            outmols[cid].SetProp(
                "_Name", 
                outmols[cid].GetProp("_Name") + " " + self.args.program.lower()
            )
            outmols[cid].SetProp("Energy", cenergy[cid])
            
            if self.args.program.lower() == "ani":
                outmols[cid].SetProp("Real charge", str(np.sum(charge)))
                outmols[cid].SetProp("Mult", str(final_mult))
            elif self.args.program.lower() == "xtb":
                outmols[cid].SetProp("Real charge", str(charge))
                outmols[cid].SetProp("Mult", str(mult))

    def _write_all_conformers(self, outmols, sorted_all_cids):
        """Write all optimized conformers to SDF file.
        
        Args:
            outmols (list): List of optimized molecule objects
            sorted_all_cids (list): Conformer IDs sorted by energy
        """
        for cid in sorted_all_cids:
            self.sdwriterall.write(outmols[cid])
        self.sdwriterall.close()

    def _cleanup_optimization_files(self):
        """Remove temporary files created during optimization."""
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

    def compute_cmin(self, file):
        """Run conformer minimization for all molecules in file.
        
        Optimizes each conformer using ANI or xTB, applies energy and
        RMSD filters, and writes results to SDF files.
        
        Args:
            file (str): Path to input SDF file
        """
        cenergy, outmols = [], []
        
        # Determine charge and multiplicity based on program
        if self.args.program.lower() == "ani":
            if self.args.charge is not None:
                self.args.log.write(
                    "\nx  Charge is automatically calculated for ANI methods, "
                    "do not use the charge option!"
                )
                self.args.log.finalize()
                sys.exit()
            elif self.args.mult is not None:
                self.args.log.write(
                    "\nx  Multiplicity is automatically calculated for ANI methods, "
                    "do not use the mult option!"
                )
                self.args.log.finalize()
                sys.exit()
            
            charge, mult, final_mult = self.charge_mult_cmin()
        
        elif self.args.program.lower() == "xtb":
            charge, mult, final_mult = self._determine_charge_mult_for_xtb(file)
        
        # Optimize each conformer
        for i, mol in enumerate(self.mols):
            if mol is not None:
                mol, energy, cmin_valid = self._optimize_single_conformer(
                    mol, i, charge, mult
                )
                
                if cmin_valid:
                    pmol = PropertyMol(mol)
                    outmols.append(pmol)
                    cenergy.append(energy)
        
        # Process results if optimizations succeeded
        if len(cenergy) >= 1:
            # Sort conformers by energy
            cids = list(range(len(outmols)))
            sorted_all_cids = sorted(cids, key=lambda cid: cenergy[cid])
            
            # Add properties to molecules
            self._add_molecule_properties(
                outmols, sorted_all_cids, cenergy, charge, mult, final_mult
            )
            
            # Write all conformers to file
            self._write_all_conformers(outmols, sorted_all_cids)
            
            # Apply filters
            self.args.log.write(
                f"\no  Applying filters to initial conformers after "
                f"{self.args.program.lower()} minimization"
            )
            selectedcids = conformer_filters(self, sorted_all_cids, cenergy, outmols)
            
            # Write filtered conformers
            self.write_confs(outmols, selectedcids, self.args.log)
        
        # Cleanup temporary files
        self._cleanup_optimization_files()

    def _setup_ani_environment(self):
        """Setup environment variables and device for ANI optimization.
        
        Returns:
            torch.device: CPU device for tensor operations
        """
        import torch
        
        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        os.environ["OMP_STACKSIZE"] = self.args.stacksize
        
        return torch.device("cpu")

    def _extract_molecular_elements(self, mol):
        """Extract element symbols from RDKit molecule.
        
        Args:
            mol: RDKit molecule object
        
        Returns:
            str: Concatenated element symbols
        """
        elements = ""
        for _, atom in enumerate(mol.GetAtoms()):
            elements += atom.GetSymbol()
        return elements

    def _prepare_ani_molecule(self, mol, charge, mult, model, DEVICE):
        """Prepare ASE molecule with ANI calculator for optimization.
        
        Args:
            mol: RDKit molecule object
            charge (list): Atomic charges
            mult (list): Atomic unpaired electrons
            model: ANI model
            DEVICE: Torch device
        
        Returns:
            tuple: (ase_molecule, coordinates, elements)
        """
        import torch
        import ase
        
        elements = self._extract_molecular_elements(mol)
        cartesians = mol.GetConformers()[0].GetPositions()
        coordinates = torch.tensor(
            [cartesians.tolist()], requires_grad=True, device=DEVICE
        )
        
        # Define ASE molecule using ANI calculator
        ase_molecule = ase.Atoms(
            elements, positions=coordinates.tolist()[0], calculator=model.ase()
        )
        
        # Adjust charge and multiplicity
        for i, atom in enumerate(ase_molecule):
            atom.charge = charge[i]
            atom.magmom = mult[i]
        
        return ase_molecule, coordinates, elements

    def _run_ani_optimization(self, ase_molecule):
        """Run BFGS optimization on ASE molecule.
        
        Args:
            ase_molecule: ASE Atoms object with ANI calculator
        
        Returns:
            bool: True if optimization succeeded, False otherwise
        """
        import ase
        
        optimizer = ase.optimize.BFGS(
            ase_molecule, trajectory="cmin_opt.traj", logfile="cmin.opt"
        )
        
        try:
            optimizer.run(fmax=self.args.opt_fmax, steps=self.args.opt_steps)
            return True
        except KeyError:
            self.args.log.write(
                f"\nx  {self.args.ani_method} could not optimize this molecule "
                f"(i.e. check if all the atoms used are compatible with ANI)"
            )
            return False

    def _compute_ani_energy(self, model, elements, coordinates, DEVICE):
        """Compute energy using ANI model.
        
        Args:
            model: ANI model
            elements (str): Concatenated element symbols
            coordinates: Torch tensor of atomic positions
            DEVICE: Torch device
        
        Returns:
            float: Energy in kcal/mol
        """
        species = model.species_to_tensor(elements).to(DEVICE).unsqueeze(0)
        _, ani_energy = model((species, coordinates))
        return ani_energy.item() * hartree_to_kcal

    def _update_molecule_coordinates(self, mol, coordinates):
        """Update RDKit molecule coordinates from optimized positions.
        
        Args:
            mol: RDKit molecule object to update
            coordinates: Torch tensor of optimized positions
        """
        cartesians = np.array(coordinates.tolist()[0])
        for j in range(mol.GetNumAtoms()):
            [x, y, z] = cartesians[j]
            mol.GetConformer().SetAtomPosition(j, Point3D(x, y, z))

    def ani_optimize(self, mol, charge, mult):
        """Perform ANI optimization on a single molecule.
        
        Uses ASE with ANI calculator to optimize molecular geometry
        and compute energy using semi-empirical neural network potentials.
        
        Args:
            mol: RDKit molecule object
            charge (list): List of atomic charges
            mult (list): List of atomic unpaired electrons
        
        Returns:
            tuple: (mol, energy, cmin_valid)
                - mol: Updated molecule with optimized coordinates
                - energy: Final energy in kcal/mol
                - cmin_valid: True if optimization succeeded
        """
        import torch
        import ase
        
        self.args.log.write(f"\no  Starting ANI optimization")
        
        # Setup environment and device
        DEVICE = self._setup_ani_environment()
        
        # Get ANI model
        model = self.get_cmin_model()
        
        # Prepare molecule for optimization
        ase_molecule, coordinates, elements = self._prepare_ani_molecule(
            mol, charge, mult, model, DEVICE
        )
        
        # Run optimization
        cmin_valid = self._run_ani_optimization(ase_molecule)
        
        if not cmin_valid:
            return mol, 0, False
        
        # Update coordinates if optimization didn't complete all steps
        if len(ase.io.Trajectory("cmin_opt.traj", mode="r")) != (self.args.opt_steps + 1):
            species_coords = ase_molecule.get_positions().tolist()
            coordinates = torch.tensor(
                [species_coords], requires_grad=True, device=DEVICE
            )
        
        # Compute final energy
        energy = self._compute_ani_energy(model, elements, coordinates, DEVICE)
        
        # Update molecule coordinates
        self._update_molecule_coordinates(mol, coordinates)
        
        return mol, energy, True

    def get_cmin_model(self):
        """Generate the ANI optimization model for CMIN.
        
        Loads the specified ANI model (e.g., ANI2x) from torchani.
        
        Returns:
            torchani model: ANI neural network potential model
        """
        import torchani
        model = getattr(torchani.models, self.args.ani_method)()
        return model

    def write_confs(self, conformers, selectedcids, log):
        """Write filtered conformers to SDF file.
        
        Writes only the conformers that passed filters to the output SDF.
        
        Args:
            conformers (list): List of molecule objects
            selectedcids (list): IDs of conformers to write
            log: Logger object for error messages
        """
        if len(conformers) > 0:
            for cid in selectedcids:
                self.sdwriter.write(conformers[cid])
            self.sdwriter.close()
        else:
            log.write("x  No conformers found!")

    def charge_mult_cmin(self):
        """Retrieve charge and multiplicity arrays for ANI optimizations.
        
        Extracts formal charges and radical electrons from each atom,
        then calculates total multiplicity.
        
        Returns:
            tuple: (charge, mult, final_mult)
                - charge: List of atomic formal charges
                - mult: List of unpaired electrons per atom
                - final_mult: Overall spin multiplicity (2S+1)
        """
        mol_cmin = self.mols[0]
        
        charge = []
        mult = []
        for _, atom in enumerate(mol_cmin.GetAtoms()):
            charge.append(atom.GetFormalCharge())
            mult.append(atom.GetNumRadicalElectrons())
        
        TotalElectronicSpin = np.sum(mult) / 2
        final_mult = int((2 * TotalElectronicSpin) + 1)
        
        return charge, mult, final_mult