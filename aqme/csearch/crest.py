#######################################################################
#             CREST functions from the CSEARCH module                 #
#######################################################################

from __future__ import print_function, absolute_import
import os
import glob
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdmolfiles
from rdkit import Geometry
import subprocess
import rdkit
from pathlib import Path
import shutil
from aqme.utils import read_file, run_command, mol_from_sdf_or_mol_or_mol2,set_destination,load_sdf
from aqme.filter import geom_filter,cluster_conformers
from rdkit.Chem import rdMolTransforms


def atompairs(mol, atom1, atom2, constraints):
    """Find non-constrained atom pairs with specific elements and their distances.
    
    Identifies pairs of bonded atoms that are not in the constraints list and 
    match the specified element types. Used to generate additional distance
    constraints for geometry optimization.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to analyze
        atom1 (str): First atom symbol ("X" for any, "H" for hydrogen)
        atom2 (str): Second atom symbol ("X" for any, "H" for hydrogen) 
        constraints (list): List of existing constraints, each containing
            at least [atom1_idx, atom2_idx] as first elements

    Returns:
        list: List of [atom1_idx, atom2_idx, distance] for matching pairs.
            Indices are 1-based, distances in Angstroms.

    Note:
        Special values for atom types:
        - "X": Matches any atom type
        - "H": Matches only hydrogen
        Otherwise matches exact element symbols
    """
    # Extract atom index pairs from constraints (convert to 1-based)
    active_pairs = []
    for constraint in constraints:
        idx1, idx2 = map(int, constraint[:2])
        active_pairs.append([idx1, idx2])

    # Find non-constrained bonds matching element types
    pairs = []
    for bond in mol.GetBonds():
        # Get 0-based indices and convert to 1-based
        idx1 = bond.GetBeginAtomIdx() + 1
        idx2 = bond.GetEndAtomIdx() + 1
        
        # Skip if bond is already constrained
        if [idx1, idx2] in active_pairs or [idx2, idx1] in active_pairs:
            continue
        
        # Get atoms and their symbols
        atom1_obj = mol.GetAtomWithIdx(idx1 - 1)
        atom2_obj = mol.GetAtomWithIdx(idx2 - 1)
        sym1 = atom1_obj.GetSymbol()
        sym2 = atom2_obj.GetSymbol()
        
        # Calculate bond distance (using 0-based indices)
        dist = round(rdMolTransforms.GetBondLength(
            mol.GetConformer(), idx1 - 1, idx2 - 1), 3)

        # Check atom type matches
        if atom1 == "X" and atom2 == "X":  # Any pair of atoms
            pairs.append([float(idx1), float(idx2), dist])
        elif atom1 == "X" and atom2 == "H":  # Any atom with hydrogen
            if "H" in (sym1, sym2):
                pairs.append([float(idx1), float(idx2), dist])
        else:  # Specific atom types
            if (sym1 == atom1 and sym2 == atom2) or \
               (sym1 == atom2 and sym2 == atom1):
                pairs.append([float(idx1), float(idx2), dist])
    
    return pairs


def get_constraint(mol, constraints):
    """Get all constraints for XTB optimization, including bond distances.
    
    Combines user-defined constraints with all existing bond distances in the
    molecule that are not already constrained. This ensures all bonds are
    fixed during initial optimization steps.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to analyze
        constraints (list): List of user-defined constraints, each containing
            at least [atom1_idx, atom2_idx, value] (indices are 1-based)
    
    Returns:
        list: Combined list of all constraints (user-defined + bond distances),
            each containing [atom1_idx, atom2_idx, value]
    """
    # Get all bond pairs not in constraints
    xx_pairs = atompairs(mol, "X", "X", constraints)
    
    # Combine user constraints and bond constraints
    all_constraints = [list(constraint) for constraint in constraints]
    all_constraints.extend(xx_pairs)
    
    return all_constraints


def xyzall_2_xyz(xyzin, name):
    """Split a multi-molecule XYZ file into individual XYZ files.
    
    Uses OpenBabel to convert a file containing multiple XYZ structures into 
    separate files, one for each conformer. Output files are named with 
    the pattern {name}_conf_.xyz.
    
    Args:
        xyzin (str): Path to input XYZ file containing multiple structures
        name (str): Base name for output files (will be appended with _conf_.xyz)
    
    Note:
        Requires OpenBabel (obabel) to be installed and accessible in PATH.
        Silently ignores any errors from OpenBabel execution.
    """
    # Convert multi-XYZ to separate files using OpenBabel
    command = ["obabel", xyzin, "-oxyz", f"-O{name}_conf_.xyz", "-m"]
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def xtb_opt_main(
    name,
    self,
    charge,
    mult,
    smi,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    method_opt,
    geom,
    sample,
    complex_ts=False,
    mol=None,
    name_init=None,
):
    """Run XTB optimization and/or CREST conformer sampling with optional constraints.
    
    This function handles both XTB geometry optimization and CREST conformer sampling,
    with support for various types of constraints (atoms, distances, angles, dihedrals).
    For complex systems or transition states, it performs staged optimizations with
    different constraint sets.
    
    Args:
        name (str): Base name for input/output files
        self (object): AQME instance containing program settings
        charge (int): Molecular charge
        mult (int): Spin multiplicity
        smi (str): SMILES string of the molecule, or None
        constraints_atoms (list): Atom indices to constrain
        constraints_dist (list): Distance constraints [[i,j,dist], ...]
        constraints_angle (list): Angle constraints [[i,j,k,angle], ...]
        constraints_dihedral (list): Dihedral constraints [[i,j,k,l,angle], ...]
        method_opt (str): Optimization method ('xtb' or 'crest')
        geom (bool): Whether to apply geometry filters
        sample (int): Number of conformers to generate
        complex_ts (bool, optional): Special handling for complexes/TS. Defaults to False.
        mol (rdkit.Mol, optional): RDKit molecule object. Defaults to None.
        name_init (str, optional): Initial molecule name. Defaults to None.
    
    Returns:
        For method_opt='crest': int (1 for success)
        For method_opt='xtb': tuple(rdkit.Mol, float, bool) 
            - Optimized molecule
            - Energy in kcal/mol  
            - Success flag
    
    Note:
        The function handles several complex scenarios:
        1. Pre-optimization with XTB before CREST
        2. Constrained optimizations for transition states
        3. CREGEN conformer sorting and clustering
        4. Error recovery with different settings
    """
    
    def _setup_environment_and_directories(self, name, method_opt, mol):
        """Set up working directories and environment variables."""
        name_no_path = os.path.basename(Path(name)).split(".xyz")[0]
        
        # Create destination directories
        if method_opt == 'crest':
            csearch_dir = set_destination(self, 'CSEARCH')
            dat_dir = csearch_dir / "crest_xyz"
            self.args.log.write("\no  Starting xTB pre-optimization before CREST sampling")
        else:  # xtb
            csearch_dir = set_destination(self, 'CMIN')
            dat_dir = csearch_dir / "xtb_xyz"
            rdmolfiles.MolToXYZFile(mol, f"{name}.xyz")
            self.args.log.write("\no  Starting xTB optimization")
        
        dat_dir.mkdir(exist_ok=True, parents=True)
        xyzin = f"{dat_dir}/{name_no_path}{'_xtb' if method_opt == 'xtb' else ''}.xyz"
        shutil.move(f"{name}.xyz", xyzin)
        
        # Set environment variables for parallel execution
        os.environ["OMP_STACKSIZE"] = self.args.stacksize
        os.environ["OMP_NUM_THREADS"] = str(self.args.nprocs)
        os.environ["MKL_NUM_THREADS"] = str(self.args.nprocs)
        
        os.chdir(dat_dir)
        
        return name_no_path, csearch_dir, dat_dir, xyzin
    
    def _has_constraints(constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral):
        """Check if any constraints are defined."""
        return (len(constraints_atoms) > 0 or len(constraints_dist) > 0 or 
                len(constraints_angle) > 0 or len(constraints_dihedral) > 0)
    
    def _build_xtb_command(xyzin, charge, mult, nprocs, xtb_keywords=None, input_file=None):
        """Build XTB command with options."""
        command = ["xtb", xyzin, "--opt"]
        
        if input_file:
            command.extend(["--input", input_file])
        
        command.extend([
            "-c", str(charge),
            "--uhf", str(int(mult) - 1),
            "-P", str(nprocs)
        ])
        
        if xtb_keywords:
            for keyword in xtb_keywords.split():
                if keyword == "--ohess":
                    command.remove("--opt")
                command.append(keyword)
        
        return command
    
    def _rename_xtb_output(dat_dir, output_file):
        """Rename XTB output file, trying both xtbopt.xyz and xtblast.xyz."""
        try:
            os.rename(f"{dat_dir}/xtbopt.xyz", output_file)
        except FileNotFoundError:
            os.rename(f"{dat_dir}/xtblast.xyz", output_file)
    
    def _cleanup_intermediate_files(keep_patterns=None):
        """Remove intermediate files that might interfere with subsequent calculations."""
        if keep_patterns is None:
            keep_patterns = []
        
        for file in glob.glob('*') + glob.glob('*.*') + glob.glob('.*'):
            if not os.path.exists(file):
                continue
            
            # Check if file matches any keep pattern
            should_keep = any(pattern in file for pattern in keep_patterns)
            if should_keep:
                continue
            
            try:
                os.remove(file)
            except IsADirectoryError:
                try:
                    shutil.rmtree(file)
                except OSError:  # Avoid problems on HPCs
                    pass
            except OSError:  # Avoid problems on HPCs
                pass
    
    # Initialize and setup
    name_no_path, csearch_dir, dat_dir, xyzin = _setup_environment_and_directories(
        self, name, method_opt, mol
    )
    opt_valid = True
    
    # Check if constraints are defined
    constrained_opt = _has_constraints(constraints_atoms, constraints_dist, 
                                      constraints_angle, constraints_dihedral)
    if constrained_opt:
        complex_ts = True
    
    # XTB pre-optimization paths
    xyzoutxtb1 = f"{dat_dir}/{name_no_path}_xtb1.xyz"
    
    if complex_ts:
        # First XTB optimization with all bonds frozen
        all_fix = get_constraint(mol, constraints_dist)
        
        create_xcontrol(
            self.args,
            constraints_atoms,
            all_fix,
            [],
            [],
            xyzin,
            "constrain1.inp",
        )
        
        command1 = _build_xtb_command(
            xyzin, charge, mult, self.args.nprocs,
            self.args.xtb_keywords, "constrain1.inp"
        )
        
        xtb_out1 = f'{os.path.dirname(Path(xyzoutxtb1))}/{os.path.basename(Path(xyzoutxtb1)).split(".xyz")[0]}'
        run_command(command1, f"{xtb_out1}.out")
        _rename_xtb_output(dat_dir, xyzoutxtb1)
        
        # Clean up intermediate files
        _cleanup_intermediate_files(keep_patterns=['_xtb2', '_xtb1', '.out'])
        
        if constrained_opt:
            # Second XTB optimization with user-defined constraints
            xyzoutxtb2 = f"{dat_dir}/{name_no_path}_xtb2.xyz"
            
            create_xcontrol(
                self.args,
                list(constraints_atoms),
                list(constraints_dist),
                list(constraints_angle),
                list(constraints_dihedral),
                xyzoutxtb1,
                "constrain2.inp",
            )
            
            command2 = _build_xtb_command(
                xyzoutxtb1, charge, mult, self.args.nprocs,
                self.args.xtb_keywords, "constrain2.inp"
            )
            
            xtb_out2 = f'{os.path.dirname(Path(xyzoutxtb2))}/{os.path.basename(Path(xyzoutxtb2)).split(".xyz")[0]}'
            run_command(command2, f"{xtb_out2}.out")
            _rename_xtb_output(dat_dir, xyzoutxtb2)
        else:
            xyzoutxtb2 = xyzoutxtb1
    
    else:
        # Unconstrained XTB pre-optimization
        try:
            command = _build_xtb_command(
                xyzin, charge, mult, self.args.nprocs,
                self.args.xtb_keywords
            )
            
            xtb_out1 = f'{os.path.dirname(Path(xyzin))}/{os.path.basename(Path(xyzin)).split(".xyz")[0]}'
            run_command(command, f"{xtb_out1}_xtb1.out")
            _rename_xtb_output(dat_dir, xyzoutxtb1)
            
        except FileNotFoundError:
            # Handle parallelization issues with metal complexes
            self.args.log.write(
                "\nx  There was an error during the xTB pre-optimization. This error might be "
                "related to parallelization of xTB jobs and is normally observed when using metal "
                "complexes in some operative systems/OpenMP versions. AQME is switching to using "
                "one processor (nprocs=1).\n"
            )
            self.args.nprocs = 1
            
            try:
                xtb_out1 = f'{os.path.dirname(Path(xyzin))}/{os.path.basename(Path(xyzin)).split(".xyz")[0]}'
                
                # Build shell command for single processor
                if self.args.xtb_keywords is None:
                    comm_xtb = (
                        f"export OMP_STACKSIZE={self.args.stacksize} && "
                        f"export OMP_NUM_THREADS={self.args.nprocs},1 && "
                        f"xtb {xyzin} --opt -c {charge} --uhf {int(mult) - 1} -P 1 >> {xtb_out1}_xtb1.out"
                    )
                else:
                    comm_xtb = (
                        f"export OMP_STACKSIZE={self.args.stacksize} && "
                        f"export OMP_NUM_THREADS={self.args.nprocs},1 && "
                        f"xtb {xyzin} --opt -c {charge} --uhf {int(mult) - 1} {self.args.xtb_keywords} "
                        f"-P 1 >> {xtb_out1}_xtb1.out"
                    )
                    if " --ohess " in comm_xtb:
                        comm_xtb = comm_xtb.replace("--opt", "")
                
                subprocess.call(comm_xtb, shell=False)
                _rename_xtb_output(dat_dir, xyzoutxtb1)
                
            except FileNotFoundError:
                if self.args.program.lower() == "crest":
                    self.args.log.write(
                        "\nx  There was another error during the xTB pre-optimization that could not be "
                        "fixed even with nprocs=1. Trying CREST directly with no xTB preoptimization.\n"
                    )
                else:
                    self.args.log.write(
                        "\nx  There was another error during the xTB pre-optimization that could not be "
                        "fixed even with nprocs=1 (this molecule will be skipped).\n"
                    )
                opt_valid = False
                mol_rd = None
        
        xyzoutxtb2 = xyzoutxtb1
    
    xyzoutall = f"{dat_dir}/{name_no_path}_conformers.xyz"
    
    def _build_crest_command(xyzfile, charge, mult, nprocs, ewin, keywords=None, cinp_file=None):
        """Build CREST command with options."""
        command = [
            "crest", xyzfile,
            "--chrg", str(charge),
            "--uhf", str(int(mult) - 1),
            "-T", str(nprocs),
            "--ewin", str(ewin)
        ]
        
        if cinp_file:
            command.extend(["--cinp", cinp_file])
        
        if keywords:
            command.extend(keywords.split())
        
        return command
    
    def _run_crest_with_fallbacks(self, command, dat_dir, name_no_path, constrained_opt, const_command=None):
        """Run CREST with multiple fallback strategies."""
        try:
            run_command(command, f"/{dat_dir}/{name_no_path}.out")
            natoms = open("crest_best.xyz").readlines()[0].strip()
            return natoms, True
        except FileNotFoundError:
            self.args.log.write(
                "\nx  CREST optimization failed! This might be caused by different reasons:\n"
                "   1) In metal complexes: using metal complexes without specifying any kind of template "
                "in the complex_type option (i.e. squareplanar).\n"
                "   2) In TSs: include the \"--noreftopo\" option in CREST with the crest_keywords option "
                "(i.e. crest_keywords=\"--noreftopo\").\n"
                "   3) In big systems: increase stacksize with the stacksize option (i.e. stacksize=\"4GB\")."
            )
            
            # Try with --noreftopo if constraints are present
            if constrained_opt and "--noreftopo" not in command and const_command:
                try:
                    self.args.log.write(
                        "\no  Constraints were detected, trying a new CREST run with --noreftopo. "
                        "WARNING! Check that your geometry doesn't isomerize!\n"
                    )
                    const_command.append('--noreftopo')
                    run_command(const_command, f"/{dat_dir}/{name_no_path}.out")
                    natoms = open("crest_best.xyz").readlines()[0].strip()
                    return natoms, True
                except FileNotFoundError:
                    self.args.log.write(
                        "\nx  CREST optimization failed again even with --noreftopo! "
                        "Contact the administrators to check this issue in more detail.\n"
                    )
            
            # Try with increased stacksize
            try:
                self.args.log.write("\no  Trying the CREST calculations with stacksize=\"4GB\".")
                os.environ["OMP_STACKSIZE"] = '4GB'
                run_command(command, f"/{dat_dir}/{name_no_path}.out")
                natoms = open("crest_best.xyz").readlines()[0].strip()
                return natoms, True
            except FileNotFoundError:
                self.args.log.write(
                    "\nx  CREST optimization failed again even with stacksize=\"4GB\"! "
                    "Contact the administrators to check this issue in more detail.\n"
                )
                return None, False
    
    def _run_cregen(self, dat_dir, name_no_path, natoms):
        """Run CREGEN conformer sorting."""
        try:
            if self.args.cregen and int(natoms) != 1:
                self.args.log.write("\no  Starting CREGEN sorting")
                command = ["crest", "crest_best.xyz", "--cregen", "crest_conformers.xyz", '--esort']
                
                if self.args.cregen_keywords:
                    command.extend(self.args.cregen_keywords.split())
                
                run_command(command, f"{dat_dir}/{name_no_path}_cregen.out")
        except UnboundLocalError:
            pass
    
    def _copy_final_crest_output(dat_dir, xyzoutall):
        """Copy final CREST output file."""
        output_files = [
            'crest_clustered.xyz',
            'crest_conformers.xyz.sorted',
            'crest_ensemble.xyz',
            'crest_conformers.xyz'
        ]
        
        for file_name in output_files:
            file_path = f'{dat_dir}/{file_name}'
            if os.path.exists(file_path):
                shutil.copy(file_path, xyzoutall)
                return True
        return False

    # CREST sampling
    if self.args.program.lower() == "crest":
        self.args.log.write("\no  Starting CREST sampling")
        
        # Create constraint file if needed
        if constrained_opt:
            create_xcontrol(
                self.args,
                list(constraints_atoms),
                list(constraints_dist),
                list(constraints_angle),
                list(constraints_dihedral),
                xyzoutxtb2,
                ".xcontrol.sample",
            )
        
        # Build CREST command
        command = _build_crest_command(
            xyzoutxtb2, charge, mult, self.args.nprocs,
            self.args.ewin_csearch, self.args.crest_keywords,
            ".xcontrol.sample" if constrained_opt else None
        )
        
        # Keep copy of command for fallback with constraints
        const_command = command.copy() if constrained_opt else None
        
        # Run CREST with fallback strategies
        natoms, opt_valid = _run_crest_with_fallbacks(
            self, command, dat_dir, name_no_path, constrained_opt, const_command
        )
        
        # Run CREGEN sorting if successful
        if opt_valid and natoms:
            _run_cregen(self, dat_dir, name_no_path, natoms)
        
        # Copy final output file
        try:
            if opt_valid:
                if not _copy_final_crest_output(dat_dir, xyzoutall):
                    raise FileNotFoundError
        except FileNotFoundError:
            self.args.log.write(
                "\nx  CREST conformer sampling failed! Please, try other options "
                "(i.e. include constrains, change the crest_keywords option, etc.)"
            )
            opt_valid = False
    
    def _convert_xyz_to_sdf(xyz_files):
        """Convert XYZ files to SDF format using OpenBabel."""
        for xyz_file in xyz_files:
            name_conf = os.path.basename(Path(xyz_file)).split(".xyz")[0]
            command = ["obabel", "-ixyz", xyz_file, "-osdf", f"-O{name_conf}.sdf"]
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    def _sort_sdf_files_by_number(sdf_files):
        """Sort SDF files numerically (handles 1, 2, ..., 10 correctly)."""
        try:
            return sorted(sdf_files, key=lambda x: int(x.split('_')[-1].split('.')[0]))
        except ValueError:  # For CMIN and QDESCP xTB optimizations
            return sdf_files
    
    def _process_sdf_file(file, self, charge, mult, smi, geom, name_init, sdwriter=None):
        """Process a single SDF file and extract molecular data."""

        supplier = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
        mol_supplier = [mol for mol in supplier]
        # IMPORTANT: release file handle
        del supplier

        # Handle invalid SDF files from GaussView
        if any(m is None for m in mol_supplier):
            mol_supplier = load_sdf(file)
        
        mol_rd = rdkit.Chem.RWMol(mol_supplier[0])
        file_nopath = os.path.basename(Path(file)).split(".sdf")[0]
        
        if self.args.program.lower() == "xtb":
            # Convert energy from Hartree to kcal/mol
            with open(f'{file_nopath}.xyz', "r") as F:
                energy_Eh = float(F.readlines()[1].split()[1])
            energy_kcal = energy_Eh * 627.5
            mol_rd.SetProp("_Name", name_init)
            os.remove(file)
            return mol_rd, energy_kcal
            
        elif self.args.program.lower() == "crest":
            # Convert energy from Hartree to kcal/mol
            try:
                energy_Eh = float(open(file, "r").readlines()[0])
            except ValueError:  # For single atom calculations
                energy_Eh = float(open(file, "r").readlines()[0].split()[1])
            
            energy_kcal = str(energy_Eh * 627.5)
            mol_rd.SetProp("_Name", '.'.join(file.split('.')[:-1]))
            mol_rd.SetProp("Energy", energy_kcal)
            mol_rd.SetProp("Real charge", str(charge))
            mol_rd.SetProp("Mult", str(int(mult)))
            
            if smi is not None:
                mol_rd.SetProp("SMILES", str(smi))
            
            # Apply geometry filter
            mol_ensemb = Chem.Mol(mol_rd)
            passing_geom = geom_filter(self, mol_ensemb, mol_rd, geom)
            
            if passing_geom and sdwriter:
                sdwriter.write(mol_rd)
            
            # Clean up files
            os.remove(file)
            os.remove(f'{file_nopath}.xyz')
            
            return mol_rd, energy_kcal
    
    def _sort_and_cluster_conformers(self, csearch_file, sample, name):
        """Sort conformers by energy and optionally cluster them."""
        suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(csearch_file, "csearch", self.args, keep_xyz=True)
        os.remove(csearch_file)
        
        # Sort by energy (CREGEN sometimes fails to do this)
        energies = [float(mol.GetProp('Energy')) for mol in suppl]
        suppl = [mol for _, mol in sorted(zip(energies, suppl), key=lambda pair: pair[0])]
        
        # Write sorted molecules
        sdwriter = Chem.SDWriter(csearch_file)
        for mol in suppl:
            sdwriter.write(mol)
        sdwriter.close()
        
        # Apply clustering if needed
        if len(suppl) > sample and self.args.auto_cluster:
            cluster_conformers(self, suppl, "rdkit", csearch_file, name, sample)

    # Post-processing: Convert XYZ to SDF and process
    if opt_valid:
        # Get XYZ files
        if self.args.program.lower() == "crest":
            xyzall_2_xyz(xyzoutall, name_no_path)
            xyz_files = glob.glob(f"{name_no_path}_conf_*.xyz")
        else:  # xtb
            xyz_files = [xyzoutxtb1]
        
        # Convert all XYZ files to SDF
        _convert_xyz_to_sdf(xyz_files)
        
        # Setup SDF writer for CREST
        if self.args.program.lower() == "crest":
            csearch_file = f"{csearch_dir}/{name_no_path}.sdf"
            sdwriter = Chem.SDWriter(csearch_file)
        else:
            sdwriter = None
        
        # Process all SDF files
        sdf_files = glob.glob(f"{name_no_path}*.sdf")
        sdf_files = _sort_sdf_files_by_number(sdf_files)
        
        for file in sdf_files:
            mol_rd, energy_kcal = _process_sdf_file(
                file, self, charge, mult, smi, geom, name_init, sdwriter
            )
        
        # Sort and cluster for CREST
        if self.args.program.lower() == "crest":
            sdwriter.close()
            _sort_and_cluster_conformers(self, csearch_file, sample, name)
    
    else:
        xyz_files = []
        energy_kcal = None
    
    # Final cleanup: Remove xTB/CREST intermediate files
    for file in glob.glob('*') + glob.glob('*.*') + glob.glob('.*'):
        if not os.path.exists(file):
            continue
        
        # Skip output files
        if '.out' in file:
            continue
        
        try:
            if self.args.program.lower() == "xtb":
                os.remove(file)
            elif self.args.program.lower() == "crest":
                # Keep XTB intermediate and output files, and final conformers
                if ('_xtb2' not in file and '_xtb1' not in file and 
                    os.path.basename(file) != os.path.basename(xyzoutall)):
                    os.remove(file)
        except IsADirectoryError:
            try:
                shutil.rmtree(file)
            except OSError:  # Avoid problems on HPCs
                pass
        except OSError:  # Avoid problems on HPCs
            pass
    
    # Return to main working directory
    os.chdir(self.args.w_dir_main)
    
    # Return results based on method
    if method_opt == 'crest':
        return 1
    else:  # xtb
        return mol_rd, energy_kcal, opt_valid


def create_xcontrol(
    args,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    xyzin,
    name_constraint,
):
    """Create an XTB control file for constrained optimizations.
    
    Generates a control file (.xcontrol) for XTB/CREST with various types of 
    constraints including fixed atoms, distances, angles, and dihedrals. For CREST
    sampling, also handles metadynamics settings for non-constrained atoms.
    
    Args:
        args: AQME arguments object containing settings (needs crest_force)
        constraints_atoms (list): List of atom indices to fix
        constraints_dist (list): List of [atom1, atom2, distance] constraints
        constraints_angle (list): List of [atom1, atom2, atom3, angle] constraints
        constraints_dihedral (list): List of [atom1, atom2, atom3, atom4, angle] constraints
        xyzin (str): Path to input XYZ file
        name_constraint (str): Output filename for constraint file
    
    Returns:
        bool: True if any constraints were written, False otherwise
    
    Note:
        For CREST sampling (name_constraint='.xcontrol.sample'), adds metadynamics
        settings for non-constrained atoms to enable efficient conformer sampling.
        Uses ranges in atom lists to handle large systems efficiently.
    """
    
    def _collect_unique_atoms(constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral):
        """Collect all unique atom indices from all constraint types."""
        unique_atoms = []
        
        # Add atoms from atom constraints
        for atom in constraints_atoms:
            if atom not in unique_atoms:
                unique_atoms.append(int(atom))
        
        # Add atoms from distance constraints (first 2 elements)
        for constraint in constraints_dist:
            for atom_idx in constraint[:2]:
                if atom_idx not in unique_atoms:
                    unique_atoms.append(int(atom_idx))
        
        # Add atoms from angle constraints (first 3 elements)
        for constraint in constraints_angle:
            for atom_idx in constraint[:3]:
                if atom_idx not in unique_atoms:
                    unique_atoms.append(int(atom_idx))
        
        # Add atoms from dihedral constraints (first 4 elements)
        for constraint in constraints_dihedral:
            for atom_idx in constraint[:4]:
                if atom_idx not in unique_atoms:
                    unique_atoms.append(int(atom_idx))
        
        return unique_atoms
    
    def _build_constraint_section(constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral):
        """Build the $constrain section of the xcontrol file."""
        content = "$constrain\n"
        
        # Add fixed atoms
        if constraints_atoms:
            content += "atoms: "
            content += ",".join(str(atom_idx) for atom_idx in constraints_atoms)
            content += "\n"
        
        # Add geometric constraints (distances, angles, dihedrals)
        constraint_types = [
            (constraints_dist, "distance", 2),
            (constraints_angle, "angle", 3),
            (constraints_dihedral, "dihedral", 4)
        ]
        
        for constraint_list, constraint_name, n_indices in constraint_types:
            if constraint_list:
                for constraint in constraint_list:
                    content += f"{constraint_name}: "
                    # Join atom indices
                    content += ",".join(str(int(val)) for val in constraint[:n_indices])
                    # Add constraint value
                    content += f",{constraint[n_indices]}\n"
        
        return content
    
    def _build_metadyn_section(n_atoms, unique_atoms):
        """Build the $metadyn section for CREST sampling."""
        content = "$metadyn\natoms: "
        
        if n_atoms == 1:
            content += "1"
            return content
        
        # Build atom ranges to avoid overly long lists
        ranges = []
        in_range = False
        range_start = None
        
        for atom_idx in range(1, n_atoms + 1):
            if atom_idx not in unique_atoms:
                if not in_range:
                    # Start new range
                    in_range = True
                    range_start = atom_idx
                
                if atom_idx == n_atoms:
                    # End of molecule - close range
                    if range_start == atom_idx:
                        ranges.append(str(range_start))
                    else:
                        ranges.append(f"{range_start}-{atom_idx}")
            else:
                if in_range:
                    # Close current range
                    if range_start == atom_idx - 1:
                        ranges.append(str(range_start))
                    else:
                        ranges.append(f"{range_start}-{atom_idx - 1}")
                    in_range = False
        
        content += ",".join(ranges)
        return content

    # Collect all unique atom indices from constraints
    unique_atoms = _collect_unique_atoms(
        constraints_atoms, 
        constraints_dist, 
        constraints_angle, 
        constraints_dihedral
    )
    
    if len(unique_atoms) == 0:
        return False
    
    # Generate coord.ref file using CREST
    subprocess.run(
        ["crest", xyzin, "--constrain", "1"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    
    # Build constraint section
    xcontrol_content = _build_constraint_section(
        constraints_atoms,
        constraints_dist,
        constraints_angle,
        constraints_dihedral
    )
    
    # Add force constant and reference
    xcontrol_content += f"force constant={args.crest_force}\n"
    xcontrol_content += "reference=coord.ref\n"
    
    # Add metadynamics section for CREST sampling
    if name_constraint == ".xcontrol.sample":
        outlines = read_file(os.getcwd(), os.getcwd(), xyzin)
        n_atoms = int(outlines[0])
        xcontrol_content += _build_metadyn_section(n_atoms, unique_atoms)
    
    xcontrol_content += "\n$end\n"
    
    # Write the control file
    with open(name_constraint, "w") as xcontrol_file:
        xcontrol_file.write(xcontrol_content)
    
    return True


def nci_ts_mol(
    smi,
    log,
    seed,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
):
    """Generate a 3D structure for non-covalent complexes or transition states.
    
    Creates a 3D structure from multiple SMILES fragments, handling atom mapping
    and constraints properly. Suitable for building non-covalent complexes and
    transition state geometries.
    
    Args:
        smi (list): List of SMILES strings for each fragment
        log: Logging object to write status messages
        seed (int): Random seed for conformer generation
        constraints_atoms (list): List of atom indices to constrain
        constraints_dist (list): List of [atom1, atom2, distance] constraints
        constraints_angle (list): List of [atom1, atom2, atom3, angle] constraints
        constraints_dihedral (list): List of [atom1, atom2, atom3, atom4, angle] constraints
    
    Returns:
        tuple: (
            rdkit.Mol: Combined 3D molecule,
            list: Adapted atom constraints,
            list: Adapted distance constraints,
            list: Adapted angle constraints,
            list: Adapted dihedral constraints
        )
    
    Note:
        - Handles atom mapping between input SMILES and output structure
        - Positions fragments with 5Å spacing to avoid clashes
        - Falls back to different embedding strategies if initial attempt fails
        - Preserves constraints when atom mapping is present in SMILES
    """
    
    def _convert_constraints_to_array(constraints):
        """Convert constraint list to numpy array of floats."""
        if constraints is not None and constraints != []:
            return np.array([[float(y) for y in x] for x in constraints])
        return constraints
    
    def _create_3d_fragments(smiles_list, with_hydrogens=False):
        """Generate 3D structures from SMILES strings."""
        fragments = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if with_hydrogens:
                mol = Chem.AddHs(mol)
            Chem.EmbedMultipleConfs(mol, numConfs=1, randomSeed=seed)
            fragments.append(mol)
        return fragments
    
    def _combine_fragments(fragments, spacing=5.0):
        """Combine molecular fragments with spatial offsets."""
        if not fragments:
            return None
            
        combined = fragments[0]
        coord = [0.0, 0.0, spacing]
        
        for fragment in fragments[1:]:
            offset = Geometry.Point3D(coord[0], coord[1], coord[2])
            combined = Chem.CombineMols(combined, fragment, offset)
            coord[1] += spacing
            Chem.SanitizeMol(combined)
            
        return combined
    
    def _perform_constrained_embed(mol, template, seed, log):
        """Try constrained embedding with fallback options."""
        try:
            return Chem.ConstrainedEmbed(mol, template, randomseed=seed)
        except ValueError:
            # Try again with hydrogens removed from template
            template_no_h = Chem.RemoveHs(template)
            try:
                return Chem.ConstrainedEmbed(mol, template_no_h, randomseed=seed)
            except ValueError:
                log.write("\nx  Constrained optimization failed due to an embedding problem with RDKit that could not be fixed!")
                return mol

    def _update_hydrogen_atom_mapping(mol):
        """Assign atom mapping numbers to hydrogens."""
        atom_map_nums = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
        max_map = max(atom_map_nums)
        
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                max_map += 1
                atom.SetAtomMapNum(int(max_map))

    def _check_smiles_mapping(smiles_list):
        """Check if all SMILES have atom mapping."""
        return all(':' in s and '[' in s for s in smiles_list)
                
    def _adapt_atom_constraints(mol, atom_constraints):
        """Map atom constraints from SMILES atom mapping to molecule indices."""
        if atom_constraints is None or not atom_constraints:
            return []
            
        adapted = []
        for atom in mol.GetAtoms():
            atom_map_num = atom.GetAtomMapNum()
            for constraint_value in atom_constraints:
                if constraint_value == atom_map_num:
                    adapted.append(float(atom.GetIdx()) + 1)
        
        return np.array(adapted) if adapted else []

    def _adapt_geometric_constraints(mol, geometric_constraints, n_atoms):
        """Map geometric constraints from SMILES atom mapping to molecule indices."""
        if geometric_constraints is None or len(geometric_constraints) == 0:
            return []
            
        adapted = []
        for constraint in geometric_constraints:
            new_constraint = []
            # Map the first n_atoms indices
            for atom_map_num in constraint[:n_atoms]:
                for atom in mol.GetAtoms():
                    if atom_map_num == atom.GetAtomMapNum():
                        new_constraint.append(float(atom.GetIdx()) + 1)
                        break
            # Add the constraint value (distance, angle, or dihedral)
            new_constraint.append(constraint[n_atoms])
            adapted.append(new_constraint)
        
        return np.array(adapted) if adapted else []
    
    # Convert all constraints to numpy arrays
    using_const = None
    if constraints_atoms is not None and constraints_atoms != []:
        using_const = constraints_atoms
        constraints_atoms = _convert_constraints_to_array(constraints_atoms)
    if constraints_dist is not None and constraints_dist != []:
        using_const = constraints_dist
        constraints_dist = _convert_constraints_to_array(constraints_dist)
    if constraints_angle is not None and constraints_angle != []:
        using_const = constraints_angle
        constraints_angle = _convert_constraints_to_array(constraints_angle)
    if constraints_dihedral is not None and constraints_dihedral != []:
        using_const = constraints_dihedral
        constraints_dihedral = _convert_constraints_to_array(constraints_dihedral)
    
    # Generate 3D structures for fragments (both with and without H)
    molsH = _create_3d_fragments(smi, with_hydrogens=True)
    mols = _create_3d_fragments(smi, with_hydrogens=False)
    
    # Combine fragments into single molecules
    molH = _combine_fragments(molsH)
    mol = _combine_fragments(mols)
    mol = Chem.AddHs(mol)
    
    # Perform constrained embedding
    mol = _perform_constrained_embed(mol, molH, seed, log)
    
    # Update atom mapping for hydrogens
    _update_hydrogen_atom_mapping(mol)
    
    # Adapt constraints based on atom mapping
    adapted_atoms = []
    adapted_dist = []
    adapted_angle = []
    adapted_dihedral = []
    
    if using_const is not None and using_const != []:
        # Check if SMILES have proper atom mapping
        if not _check_smiles_mapping(smi):
            log.write(f"\nx  Constraints were specified {using_const} but atoms might not be mapped in the SMILES input!")
            # Return original constraints without adaptation
            adapted_atoms = constraints_atoms
            adapted_dist = constraints_dist
            adapted_angle = constraints_angle
            adapted_dihedral = constraints_dihedral
        else:
            # Adapt constraints using atom mapping
            if constraints_atoms is not None:
                adapted_atoms = _adapt_atom_constraints(mol, constraints_atoms)
            
            if constraints_dist is not None:
                adapted_dist = _adapt_geometric_constraints(mol, constraints_dist, 2)
            
            if constraints_angle is not None:
                adapted_angle = _adapt_geometric_constraints(mol, constraints_angle, 3)
            
            if constraints_dihedral is not None:
                adapted_dihedral = _adapt_geometric_constraints(mol, constraints_dihedral, 4)
    
    return (
        mol,
        adapted_atoms,
        adapted_dist,
        adapted_angle,
        adapted_dihedral,
    )
