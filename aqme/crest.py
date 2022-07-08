#!/usr/bin/python

#######################################################################
# Runs crest on an xyz file, can add in options to the code if needed #
#######################################################################

from __future__ import print_function, absolute_import
import os
import glob
from rdkit.Chem import AllChem as Chem
import subprocess
import rdkit
from pathlib import Path
import shutil
from aqme.utils import read_file, run_command
from rdkit.Chem import rdMolTransforms


def atompairs(mol, atom1, atom2, constraints):
    active = []
    for x in constraints:
        active.append(x[:2])

    for i, x in enumerate(active):
        active[i] = [int(j) for j in x]

    pairs = []
    bonds = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in mol.GetBonds()]
    for [a, b] in bonds:
        if [a + 1, b + 1] not in active and [b + 1, a + 1] not in active:
            dist = round(rdMolTransforms.GetBondLength(mol.GetConformer(), a, b), 3)
            at_a, at_b = mol.GetAtoms()[a].GetSymbol(), mol.GetAtoms()[b].GetSymbol()
            if atom1 == "X" and atom2 == "X":
                pairs.append([float(a + 1), float(b + 1), dist])
            elif atom1 == "X" and atom2 == "H":
                if at_a == "H" or at_b == "H":
                    pairs.append([float(a + 1), float(b + 1), dist])
            else:
                if (at_a == atom1 and at_b == atom2) or (
                    at_a == atom2 and at_b == atom1
                ):
                    pairs.append([float(a + 1), float(b + 1), dist])
    return pairs


def get_constraint(mol, constraints):
    # constrained optimization with xtb
    xx_pairs = atompairs(mol, "X", "X", constraints)
    all_fix = []
    for x in constraints:
        all_fix.append(list(x))
    for x in xx_pairs:
        all_fix.append(x)
    return all_fix


def xyzall_2_xyz(xyzin, name):
    # converting multiple xyz to single
    command_run_1 = ["obabel", xyzin, "-oxyz", "-O" + name + "_conf_.xyz", "-m"]
    subprocess.run(command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def crest_opt(
    name,
    dup_data,
    dup_data_idx,
    args,
    charge,
    mult,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    mol=None,
):

    """
    Run xTB using subprocess to perform CREST/CREGEN conformer sampling
    """

    nci_ts_complex = False
    if (
        (len(constraints_atoms) != 0)
        or (len(constraints_dist) != 0)
        or (len(constraints_angle) != 0)
        or (len(constraints_dihedral) != 0)
    ):
        nci_ts_complex = True

    name_no_path = name.replace("/", "\\").split("\\")[-1].split(".")[0]
    csearch_dir = Path(args.w_dir_main)
    dat_dir = csearch_dir / "CSEARCH" / "crest_xyz" / name_no_path
    dat_dir.mkdir(exist_ok=True, parents=True)

    xyzin = f"{dat_dir}/{name_no_path}.xyz"
    sdwriter = Chem.SDWriter(str(f"{csearch_dir}/CSEARCH/crest/{name_no_path}.sdf"))

    shutil.move(f"{name}.xyz", xyzin)

    os.chdir(dat_dir)
    # for systems that were created from 1D and 2D inputs (i.e. SMILES), this part includes two xTB
    # constrained optimizations to avoid geometry problems in noncovalent complexes and transition states
    if nci_ts_complex:
        if len(constraints_dist) != 0:
            # xTB optimization with all bonds frozen
            xyzoutxtb1 = str(dat_dir) + "/" + name_no_path + "_xtb1.xyz"

            all_fix = get_constraint(mol, constraints_dist)

            _ = create_xcontrol(
                args,
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
                str(args.nprocs),
            ]

            run_command(command1, "{}.out".format(xyzoutxtb1.split(".xyz")[0]))
            os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb1)

        else:
            xyzoutxtb1 = xyzin

        xyzoutxtb2 = str(dat_dir) + "/" + name_no_path + "_xtb2.xyz"
        # xTB optimization with the user-defined constraints
        _ = create_xcontrol(
            args,
            list(constraints_atoms),
            list(constraints_dist),
            list(constraints_angle),
            list(constraints_dihedral),
            xyzin,
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
            str(args.nprocs),
        ]
        run_command(command2, "{}.out".format(xyzoutxtb2.split(".xyz")[0]))
        os.rename(str(dat_dir) + "/xtbopt.xyz", xyzoutxtb2)

    else:
        xyzoutxtb2 = xyzin

    xyzoutall = str(dat_dir) + "/" + name_no_path + "_conformers.xyz"

    constrained_sampling = False
    if nci_ts_complex:
        constrained_sampling = create_xcontrol(
            args,
            list(constraints_atoms),
            list(constraints_dist),
            list(constraints_angle),
            list(constraints_dihedral),
            xyzin,
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
        str(args.nprocs),
        "--ewin",
        str(args.ewin_csearch),
    ]

    if constrained_sampling:
        command.append("-cinp")
        command.append(".xcontrol.sample")

    if args.crest_keywords is not None:
        for keyword in args.crest_keywords.split():
            command.append(keyword)

    run_command(command, str(dat_dir) + "/crest.out")

    # get number of n_atoms
    natoms = open("crest_best.xyz").readlines()[0].strip()

    if args.cregen and int(natoms) != 1:
        command = ["crest", "crest_best.xyz", "--cregen", "crest_conformers.xyz"]

        if args.cregen_keywords is not None:
            for keyword in args.cregen_keywords.split():
                command.append(keyword)

        run_command(command, str(dat_dir) + "/cregen.out")

    try:
        if os.path.exists(str(dat_dir) + "/crest_clustered.xyz"):
            shutil.copy(str(dat_dir) + "/crest_clustered.xyz", xyzoutall)

        elif os.path.exists(str(dat_dir) + "/crest_ensemble.xyz"):
            shutil.copy(str(dat_dir) + "/crest_ensemble.xyz", xyzoutall)
        else:
            shutil.copy(str(dat_dir) + "/crest_conformers.xyz", xyzoutall)
    except FileNotFoundError:
        args.log.write(
            "x   CREST conformer sampling failed! Please, try other options (i.e. include constrains, change the crest_keywords option, etc.)"
        )

    xyzall_2_xyz(xyzoutall, name_no_path)

    xyz_files = glob.glob(name_no_path + "_conf_*.xyz")
    for _, file in enumerate(xyz_files):
        name_conf = file.split(".xyz")[0]
        command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name_conf + ".sdf"]
        subprocess.run(
            command_xyz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )

    sdf_files = glob.glob(name_no_path + "*.sdf")
    for file in sdf_files:
        mol = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
        mol_rd = rdkit.Chem.RWMol(mol[0])
        energy = str(open(file, "r").readlines()[0])
        mol_rd.SetProp("Energy", energy)
        mol_rd.SetProp("Real charge", str(charge))
        mol_rd.SetProp("Mult", str(int(mult)))
        sdwriter.write(mol_rd)

    dup_data.at[dup_data_idx, "crest-conformers"] = len(xyz_files)

    os.chdir(args.w_dir_main)

    return 1


def create_xcontrol(
    args,
    constraints_atoms,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    xyzin,
    name_constraint,
):
    """
    Function to create the .xcontrol.sample if constraints are defined
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
        os.remove(".xcontrol.sample")

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
            for atom_idx in range(1, n_atoms + 1):
                if atom_idx not in unique_atoms:
                    if atom_idx == n_atoms:
                        edited_xcontrol += f"{atom_idx}\n"
                    else:
                        edited_xcontrol += f"{atom_idx},"
        edited_xcontrol += "$end\n"

        # write the file
        xcontrol_file = open(name_constraint, "w")
        xcontrol_file.write(edited_xcontrol)
        xcontrol_file.close()

    return constrained_sampling
