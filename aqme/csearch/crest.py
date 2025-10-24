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

    name_no_path = os.path.basename(Path(name)).split(".xyz")[0]

    # folder to create the files
    if method_opt == 'crest':
        csearch_dir = set_destination(self,'CSEARCH')
    elif method_opt == 'xtb':
        csearch_dir = set_destination(self,'CMIN')

    # create the initial xyz input
    if method_opt == 'crest':
        self.args.log.write(f"\no  Starting xTB pre-optimization before CREST sampling")
        dat_dir = csearch_dir / "crest_xyz"
        xyzin = f"{dat_dir}/{name_no_path}.xyz"
    elif method_opt == 'xtb':
        rdmolfiles.MolToXYZFile(mol, f"{name}.xyz")
        self.args.log.write(f"\no  Starting xTB optimization")
        dat_dir = csearch_dir / "xtb_xyz"
        xyzin = f"{dat_dir}/{name_no_path}_xtb.xyz"
    dat_dir.mkdir(exist_ok=True, parents=True)
    shutil.move(f"{name}.xyz", xyzin)

    os.environ["OMP_STACKSIZE"] = self.args.stacksize
    # to run xTB/CREST with more than 1 processor
    os.environ["OMP_NUM_THREADS"] = str(self.args.nprocs)
    os.environ["MKL_NUM_THREADS"] = str(self.args.nprocs)
    opt_valid = True

    os.chdir(dat_dir)
    
    # for systems that were created from 1D and 2D inputs (i.e. SMILES), this part includes two xTB
    # constrained optimizations to avoid geometry problems in noncovalent complexes and transition states

    # xTB optimization with all bonds frozen
    constrained_opt = False
    if len(constraints_atoms) > 0 or len(constraints_dist) > 0 or len(constraints_angle) > 0 or len(constraints_dihedral) > 0:
        constrained_opt = True
        complex_ts = True

    xyzoutxtb1 = str(dat_dir) + "/" + name_no_path + "_xtb1.xyz"
    if complex_ts:
        all_fix = get_constraint(mol, constraints_dist)

        _ = create_xcontrol(
            self.args,
            constraints_atoms,
            all_fix,
            [],
            [],
            xyzin,
            "constrain1.inp",
        )

        command1 = [
            "xtb",
            xyzin,
            "--opt",
            "--input",
            "constrain1.inp",
            "-c",
            str(charge),
            "--uhf",
            str(int(mult) - 1),
            "-P",
            str(self.args.nprocs),
        ]

        if self.args.xtb_keywords is not None:
            for keyword in self.args.xtb_keywords.split():
                if keyword == "--ohess":
                    command1.remove("--opt ")
                command1.append(keyword)

        xtb_out1 = f'{os.path.dirname(Path(xyzoutxtb1))}/{os.path.basename(Path(xyzoutxtb1)).split(".xyz")[0]}'
        run_command(command1, f"{xtb_out1}.out")
        try:
            os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb1)
        except FileNotFoundError:
            os.rename(str(dat_dir) + "/xtblast.xyz", xyzoutxtb1)

        # remove files that might interfere in subsequent calculations (i.e. wrong electron readings)
        for file in glob.glob('*') + glob.glob('*.*') + glob.glob('.*'):
            if os.path.exists(file):
                if file.find('_xtb2') == -1 and file.find('_xtb1') == -1 and file.find('.out') == -1:
                    try:
                        os.remove(file)
                    except IsADirectoryError:
                        try:
                            shutil.rmtree(file)
                        except OSError: # this avoids problems when running AQME in HPCs
                            pass

        if constrained_opt:
            xyzoutxtb2 = str(dat_dir) + "/" + name_no_path + "_xtb2.xyz"
            # xTB optimization with the user-defined constraints
            _ = create_xcontrol(
                self.args,
                list(constraints_atoms),
                list(constraints_dist),
                list(constraints_angle),
                list(constraints_dihedral),
                xyzoutxtb1,
                "constrain2.inp",
            )

            command2 = [
                "xtb",
                xyzoutxtb1,
                "--opt",
                "--input",
                "constrain2.inp",
                "-c",
                str(charge),
                "--uhf",
                str(int(mult) - 1),
                "-P",
                str(self.args.nprocs),
            ]

            if self.args.xtb_keywords is not None:
                for keyword in self.args.xtb_keywords.split():
                    if keyword == "--ohess":
                        command2.remove("--opt ")
                    command2.append(keyword)

            xtb_out2 = f'{os.path.dirname(Path(xyzoutxtb2))}/{os.path.basename(Path(xyzoutxtb2)).split(".xyz")[0]}'
            run_command(command2, f"{xtb_out2}.out")

            try:
                os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb2)
            except FileNotFoundError:
                os.rename(str(dat_dir) + "/xtblast.xyz", xyzoutxtb2)
        else:
            xyzoutxtb2 = xyzoutxtb1

    else:
        # Preoptimization with xTB to avoid issues from innacurate starting structures in CREST.
        # If you're dealing with a large system, increase the stack size
        try:
            command = [
                "xtb",
                xyzin,
                "--opt",
                "-c",
                str(charge),
                "--uhf",
                str(int(mult) - 1),
                "-P",
                str(self.args.nprocs), 
            ]

            if self.args.xtb_keywords is not None:
                for keyword in self.args.xtb_keywords.split():
                    if keyword == "--ohess":
                        command.remove("--opt ")
                    command.append(keyword)
            xtb_out1 = f'{os.path.dirname(Path(xyzin))}/{os.path.basename(Path(xyzin)).split(".xyz")[0]}'
            run_command(command, f"{xtb_out1}_xtb1.out")

            os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb1)
        except FileNotFoundError:
            self.args.log.write(f"\nx  There was an error during the xTB pre-optimization. This error might be related to parallelization of xTB jobs and is normally observed when using metal complexes in some operative systems/OpenMP versions. AQME is switching to using one processor (nprocs=1).\n")
            self.args.nprocs = 1
            try:
                xtb_out1 = f'{os.path.dirname(Path(xyzin))}/{os.path.basename(Path(xyzin)).split(".xyz")[0]}'
                if self.args.xtb_keywords is None:
                    comm_xtb = f"export OMP_STACKSIZE={self.args.stacksize} && export OMP_NUM_THREADS={self.args.nprocs},1 \
                    && xtb {xyzin} --opt -c {charge} --uhf {int(mult) - 1} -P 1 >> {xtb_out1}_xtb1.out"
                else:
                    comm_xtb = f"export OMP_STACKSIZE={self.args.stacksize} && export OMP_NUM_THREADS={self.args.nprocs},1 \
                    && xtb {xyzin} --opt -c {charge} --uhf {int(mult) - 1} {self.args.xtb_keywords} -P 1 >> {xtb_out1}_xtb1.out"
                    if " --ohess " in comm_xtb:
                        comm_xtb.remove("--opt ")
                subprocess.call(comm_xtb, shell=False)
                os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb1)
            except FileNotFoundError:
                if self.args.program.lower() == "crest":
                    self.args.log.write(f"\nx  There was another error during the xTB pre-optimization that could not be fixed even with nprocs=1. Trying CREST directly with no xTB preoptimization.\n")
                else:
                    self.args.log.write(f"\nx  There was another error during the xTB pre-optimization that could not be fixed even with nprocs=1 (this molecule will be skipped).\n")
                opt_valid = False
                mol_rd = None

        xyzoutxtb2 = xyzoutxtb1

    xyzoutall = str(dat_dir) + "/" + name_no_path + "_conformers.xyz"

    # CREST sampling
    if self.args.program.lower() == "crest":
        self.args.log.write(f"\no  Starting CREST sampling")
        if constrained_opt:
            _ = create_xcontrol(
                self.args,
                list(constraints_atoms),
                list(constraints_dist),
                list(constraints_angle),
                list(constraints_dihedral),
                xyzoutxtb2,
                ".xcontrol.sample",
            )

        command = [
            "crest",
            xyzoutxtb2,
            "--chrg",
            str(charge),
            "--uhf",
            str(int(mult) - 1),
            "-T",
            str(self.args.nprocs),
            "--ewin",
            str(self.args.ewin_csearch),
        ]

        if constrained_opt:
            command.append("--cinp")
            command.append(".xcontrol.sample")
            const_command = command.copy()
        if self.args.crest_keywords is not None:
            for keyword in self.args.crest_keywords.split():
                command.append(keyword)
        try:
            run_command(command, f"/{dat_dir}/{name_no_path}.out")
            natoms = open("crest_best.xyz").readlines()[0].strip()
        except FileNotFoundError:
            self.args.log.write(f"\nx  CREST optimization failed! This might be caused by different reasons:\n   1) In metal complexes: using metal complexes without specifying any kind of template in the complex_type option (i.e. squareplanar).\n   2) In TSs: include the \"--noreftopo\" option in CREST with the crest_keywords option (i.e. crest_keywords=\"--noreftopo\").\n   3) In big systems: increase stacksize with the stacksize option (i.e. stacksize=\"4GB\").")
            if constrained_opt and "--noreftopo" not in command:
                try:
                    self.args.log.write(f"\no  Constraints were detected, trying a new CREST run with --noreftopo. WARNING! Check that your geometry doesn't isomerize!\n")
                    if self.args.crest_keywords is not None:
                        for keyword in self.args.crest_keywords.split():
                            const_command.append(keyword)
                    const_command.append('--noreftopo')
                    run_command(const_command, f"/{dat_dir}/{name_no_path}.out")
                    natoms = open("crest_best.xyz").readlines()[0].strip()
                except FileNotFoundError:
                    self.args.log.write(f"\nx  CREST optimization failed again even with --noreftopo! Contact the administrators to check this issue in more detail.\n")
                    opt_valid = False
            else:
                opt_valid = False
            if not opt_valid:
                try:
                    self.args.log.write(f"\no  Trying the CREST calculations with stacksize=\"4GB\".")
                    os.environ["OMP_STACKSIZE"] = '4GB'
                    run_command(command, f"/{dat_dir}/{name_no_path}.out")
                    natoms = open("crest_best.xyz").readlines()[0].strip()
                except FileNotFoundError:
                    self.args.log.write(f"\nx  CREST optimization failed again even with stacksize=\"4GB\"! Contact the administrators to check this issue in more detail.\n")

        # CREGEN sorting
        try:
            if self.args.cregen and int(natoms) != 1 and opt_valid:
                cregen_text = f"\no  Starting CREGEN sorting"
                command = ["crest", "crest_best.xyz", "--cregen", "crest_conformers.xyz", '--esort']
                # we don't use this part, since the code runs a Butina clustering afterwards
                # if self.args.auto_cluster:
                #     cregen_text += ' and conformer selection through clustering (users can disable conformer selection with --auto_cluster False)'
                #     command = command + ['--cluster', f'{sample}']

                if self.args.cregen_keywords is not None:
                    for keyword in self.args.cregen_keywords.split():
                        command.append(keyword)
                self.args.log.write(cregen_text)
                run_command(command, f"{dat_dir}/{name_no_path}_cregen.out")

        except UnboundLocalError:
            pass

        # rename final XYZ file
        try:
            if opt_valid:
                for file_name in ['crest_clustered.xyz','crest_conformers.xyz.sorted','crest_ensemble.xyz','crest_conformers.xyz']:
                    if os.path.exists(f'{dat_dir}/{file_name}'):
                        shutil.copy(f'{dat_dir}/{file_name}', xyzoutall)
                        break
        except FileNotFoundError:
            self.args.log.write("\nx  CREST conformer sampling failed! Please, try other options (i.e. include constrains, change the crest_keywords option, etc.)")
            opt_valid = False

    if opt_valid:
        if self.args.program.lower() == "crest":
            xyzall_2_xyz(xyzoutall, name_no_path)
            xyz_files = glob.glob(name_no_path + "_conf_*.xyz")
        if self.args.program.lower() == "xtb":
            xyz_files = [xyzoutxtb1]
        for _, file in enumerate(xyz_files):
            name_conf = f'{os.path.basename(Path(file)).split(".xyz")[0]}'
            command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name_conf + ".sdf"]
            subprocess.run(command_xyz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        if self.args.program.lower() == "crest":
            csearch_file = str(f"{csearch_dir}/{name_no_path}.sdf")
            sdwriter = Chem.SDWriter(csearch_file)

        sdf_files = glob.glob(name_no_path + "*.sdf")
        # the next function is needed to keep the order (glob.glob sorts first 1, then 10 instead of 2)
        try:
            def func(x):
                return int(x.split('_')[-1].split('.')[0])
            sdf_files = sorted(sdf_files, key=func)
        except ValueError: # for CMIN and QDESCP xTB optimizations
            pass
        for file in sdf_files:
            mol = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
            # transform invalid SDF files created with GaussView into valid SDF from Open Babel
            if None in mol:
                mol = load_sdf(file)
            mol_rd = rdkit.Chem.RWMol(mol[0])
            file_nopath = f'{os.path.basename(Path(file)).split(".sdf")[0]}'
            if self.args.program.lower() == "xtb":
                # convert from hartree (default in xtb) to kcal
                energy_Eh = float(open(f'{file_nopath}.xyz', "r").readlines()[1].split()[1])
                energy_kcal = energy_Eh*627.5
                mol_rd.SetProp("_Name", name_init)
                os.remove(file)
            elif self.args.program.lower() == "crest":
                # convert from hartree (default in xtb) to kcal
                try:
                    energy_Eh = float(open(file, "r").readlines()[0])
                except ValueError: # for calcs with a single atom
                    energy_Eh = float(open(f'{file}', "r").readlines()[0].split()[1])
                energy_kcal = str(energy_Eh*627.5)
                mol_rd.SetProp("_Name", '.'.join(file.split('.')[:-1]))
                mol_rd.SetProp("Energy", energy_kcal)
                mol_rd.SetProp("Real charge", str(charge))
                mol_rd.SetProp("Mult", str(int(mult)))
                if smi is not None:
                    mol_rd.SetProp("SMILES", str(smi))
                mol_ensemb = Chem.Mol(mol_rd)
                passing_geom = geom_filter(self,mol_ensemb,mol_rd,geom)
                if passing_geom:
                    sdwriter.write(mol_rd)
                os.remove(file)
                os.remove(f'{file_nopath}.xyz')

        # sorting and clusterization
        if self.args.program.lower() == "crest":
            sdwriter.close()
            suppl, _, _, _ = mol_from_sdf_or_mol_or_mol2(f'{csearch_file}', "csearch", self.args)
            os.remove(f'{csearch_file}')
            # sort by energy (even though CREGEN should do that automatically, it fails to do so sometimes)
            allenergy = []
            for mol in suppl:
                allenergy.append(float(mol.GetProp('Energy')))
            suppl = [mol for _, mol in sorted(zip(allenergy, suppl), key=lambda pair: pair[0])]

            sdwriter = Chem.SDWriter(f'{csearch_file}')
            for mol in suppl:
                sdwriter.write(mol)
            sdwriter.close()
            if len(suppl) > sample and self.args.auto_cluster:
                _ = cluster_conformers(self,suppl,"rdkit",csearch_file,name,sample)
            
    else:
        xyz_files = []
        energy_kcal = None

    # remove xTB/CREST files to avoid wrong readings of molecular information
    for file in glob.glob('*') + glob.glob('*.*') + glob.glob('.*'):
        if os.path.exists(file):
            if file.find('.out') == -1:
                if self.args.program.lower() == "xtb":
                    try:
                        os.remove(file)
                    except OSError: # this avoids problems when running AQME in HPCs
                        pass
                elif self.args.program.lower() == "crest":
                    if file.find('_xtb2') == -1 and file.find('_xtb1') == -1 and file.find('.out') == -1:
                        try:
                            if os.path.basename(file) != os.path.basename(xyzoutall):
                                try:
                                    os.remove(file)
                                except OSError: # this avoids problems when running AQME in HPCs
                                    pass
                        except IsADirectoryError:
                            try:
                                shutil.rmtree(file)
                            except OSError: # this avoids problems when running AQME in HPCs
                                pass

    os.chdir(self.args.w_dir_main)

    if method_opt == 'crest':
        return 1

    if method_opt == 'xtb':
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

    constrained_sampling = False

    unique_atoms = []
    for atom in constraints_atoms:
        if atom not in unique_atoms:
            unique_atoms.append(int(atom))
    for x in constraints_dist:
        for i in x[:2]:
            if i not in unique_atoms:
                unique_atoms.append(int(i))
    for x in constraints_angle:
        for i in x[:3]:
            if i not in unique_atoms:
                unique_atoms.append(int(i))
    for x in constraints_dihedral:
        for i in x[:4]:
            if i not in unique_atoms:
                unique_atoms.append(int(i))

    if len(unique_atoms) > 0:

        constrained_sampling = True

        # call --constrain just fo create a coord.ref file
        subprocess.run(
            ["crest", xyzin, "--constrain", "1"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # add the constraints part
        edited_xcontrol = "$constrain\n"

        if constraints_atoms != []:
            edited_xcontrol += "atoms: "
            edited_xcontrol += ",".join(str(atom_idx) for atom_idx in constraints_atoms)
            edited_xcontrol += "\n"

        for constraint_type in [
            constraints_dist,
            constraints_angle,
            constraints_dihedral,
        ]:
            if constraint_type != []:
                for const in constraint_type:
                    if constraint_type == constraints_dist:
                        edited_xcontrol += "distance: "
                        edited_xcontrol += (
                            ",".join(str(int(val)) for val in const[:2])
                            + ","
                            + str(const[2])
                        )
                    elif constraint_type == constraints_angle:
                        edited_xcontrol += "angle: "
                        edited_xcontrol += (
                            ",".join(str(int(val)) for val in const[:3])
                            + ","
                            + str(const[3])
                        )
                    elif constraint_type == constraints_dihedral:
                        edited_xcontrol += "dihedral: "
                        edited_xcontrol += (
                            ",".join(str(int(val)) for val in const[:4])
                            + ","
                            + str(const[4])
                        )
                    edited_xcontrol += "\n"

        edited_xcontrol += f"force constant={args.crest_force}\n"
        edited_xcontrol += "reference=coord.ref\n"

        # metadyn part
        if name_constraint == ".xcontrol.sample":
            outlines = read_file(os.getcwd(), os.getcwd(), xyzin)
            n_atoms = int(outlines[0])
            edited_xcontrol += "$metadyn\n"
            edited_xcontrol += "atoms: "
            if n_atoms == 1: # just to avoid bugs when parsing single atoms
                edited_xcontrol += "atoms: 1"
            else:
                # if the list is too long, CREST doesn't read it when called from subprocess() in Python
                # I need to include ranges to shorten the lists of atoms for the $metadyn section
                new_cycle = True
                start = True
                for atom_idx in range(1, n_atoms + 1):
                    if new_cycle:
                        new_cycle = False
                        start_idx = atom_idx
                    if atom_idx not in unique_atoms:
                        if start: # just in case the first atom isn't part of the list
                            start = False
                            start_idx = atom_idx
                        elif atom_idx == n_atoms:
                            if start_idx == (atom_idx):
                                edited_xcontrol += f"{start_idx}"
                            else:
                                edited_xcontrol += f"{start_idx}-{atom_idx}"
                    elif atom_idx in unique_atoms and not start:
                        new_cycle = True
                        if start_idx == (atom_idx-1):
                            edited_xcontrol += f"{start_idx}"
                        else:
                            edited_xcontrol += f"{start_idx}-{atom_idx-1}"
                        if atom_idx != n_atoms:
                            edited_xcontrol += ','

        edited_xcontrol += "\n$end\n"

        # write the file
        xcontrol_file = open(name_constraint, "w")
        xcontrol_file.write(edited_xcontrol)
        xcontrol_file.close()

    return constrained_sampling


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
    using_const = None
    if constraints_atoms is not None and constraints_atoms != []:
        using_const = constraints_atoms
        constraints_atoms = [[float(y) for y in x] for x in constraints_atoms]
        constraints_atoms = np.array(constraints_atoms)
    if constraints_dist is not None and constraints_dist != []:
        using_const = constraints_dist
        constraints_dist = [[float(y) for y in x] for x in constraints_dist]
        constraints_dist = np.array(constraints_dist)
    if constraints_angle is not None and constraints_angle != []:
        using_const = constraints_angle
        constraints_angle = [[float(y) for y in x] for x in constraints_angle]
        constraints_angle = np.array(constraints_angle)
    if constraints_dihedral is not None and constraints_dihedral != []:
        using_const = constraints_dihedral
        constraints_dihedral = [[float(y) for y in x] for x in constraints_dihedral]
        constraints_dihedral = np.array(constraints_dihedral)

    molsH = []
    mols = []
    for m in smi:
        mols.append(Chem.MolFromSmiles(m))
        molsH.append(Chem.AddHs(Chem.MolFromSmiles(m)))

    for m in molsH:
        Chem.EmbedMultipleConfs(m, numConfs=1, randomSeed=seed)
    for m in mols:
        Chem.EmbedMultipleConfs(m, numConfs=1, randomSeed=seed)

    coord = [0.0, 0.0, 5.0]
    molH = molsH[0]
    for _, fragment in enumerate(molsH[1:]):
        offset_3d = Geometry.Point3D(coord[0], coord[1], coord[2])
        molH = Chem.CombineMols(molH, fragment, offset_3d)
        coord[1] += 5
        Chem.SanitizeMol(molH)

    coord = [0.0, 0.0, 5.0]
    mol = mols[0]
    for _, fragment in enumerate(mols[1:]):
        offset_3d = Geometry.Point3D(coord[0], coord[1], coord[2])
        mol = Chem.CombineMols(mol, fragment, offset_3d)
        coord[1] += 5
        Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)

    try:
        mol = Chem.ConstrainedEmbed(mol, molH, randomseed=seed)
    except ValueError: # removingH in the core (molH) avoids problems related to RDKit embedding
        molH = Chem.RemoveHs(molH)
        try:
            mol = Chem.ConstrainedEmbed(mol, molH, randomseed=seed)
        except ValueError:
            log.write(f"\nx  Constrained optimization failed due to an embedding problem with RDKit that could not be fixed!")
            pass

    atom_map = []
    for atom in mol.GetAtoms():
        atom_map.append(atom.GetAtomMapNum())

    max_map = max(atom_map)
    for a in mol.GetAtoms():
        if a.GetSymbol() == "H":
            max_map += 1
            a.SetAtomMapNum(int(max_map))

    adapted_atoms = []
    adapted_dist = []
    adapted_angle = []
    adapted_dihedral = []
    # assign constraints
    if using_const != [] and using_const is not None:
        for smi_part in smi:
            if ':' not in smi_part or '[' not in smi_part: # for SMILES that are not mapped
                log.write(f"\nx  Constraints were specified {using_const} but atoms might not be mapped in the SMILES input!")
                adapted_atoms = constraints_atoms
                adapted_dist = constraints_dist
                adapted_angle = constraints_angle
                adapted_dihedral = constraints_dihedral
                break

        if constraints_atoms is not None:
            for _, ele in enumerate(constraints_atoms):
                if adapted_atoms != []: # for mapped SMILES
                    for atom in mol.GetAtoms():
                        if ele == atom.GetAtomMapNum():
                            adapted_atoms.append(float(atom.GetIdx()) + 1)
                else:
                    adapted_atoms = constraints_atoms
            adapted_atoms = np.array(adapted_atoms)

        constraints = [constraints_dist, constraints_angle, constraints_dihedral]
        adapted_consts = [adapted_dist, adapted_angle, adapted_dihedral]
        n_consts= [2, 3, 4]
        for const,adapted_const,n_const in zip(constraints,adapted_consts,n_consts):
            if adapted_const == []: # for mapped SMILES
                for _, r in enumerate(const):
                    nr = []
                    for _, ele in enumerate(r[:n_const]):
                        if const is not None:
                            for atom in mol.GetAtoms():
                                if ele == atom.GetAtomMapNum():
                                    nr.append(float(atom.GetIdx()) + 1)
                    nr.append(r[-1])
                    adapted_const.append(nr)
                adapted_const = np.array(adapted_const)

    return (
        mol,
        adapted_atoms,
        adapted_dist,
        adapted_angle,
        adapted_dihedral,
    )
