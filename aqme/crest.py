#!/usr/bin/python
from __future__ import print_function, absolute_import

#######################################################################
# Runs crest on an xyz file, can add in options to the code if needed #
#######################################################################

# Python Libraries
import os
from rdkit.Chem import rdMolTransforms
import glob
import subprocess
import rdkit
from rdkit.Chem import AllChem
from pathlib import Path


def atompairs(mol, atom1, atom2, constraints):
    """
    Returns addtional constraints

    Parameters
    ----------
    mol : RDKit mol object
        Mol object
    atom1 : str
        type of atom 1
    atom2 : str
        type of atom 2

    Returns
    -------
    [int,int] [int,int]
        atom pairs, constraint value
    """
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
    # xh_pairs, xh_vals = atompairs(rdmol, "X", "H", active)
    # print(active, template)
    all_fix = []
    for x in constraints:
        all_fix.append(x)
    for x in xx_pairs:
        all_fix.append(x)
    return all_fix


def xyzall_2_xyz(xyzin, name):
    # converting multiple xyz to single
    command_run_1 = ["obabel", xyzin, "-oxyz", "-O" + name + "_conf_.xyz", "-m"]
    subprocess.run(command_run_1)


def run_xtb(xyzin, xyzoutxtb, constraints_all, charge, mult):
    # check whether job has already been run
    if os.path.exists(xyzoutxtb):
        print("   {0} already exists: skipping xTB optimization".format(xyzoutxtb))
        # pass
    else:
        command = [
            os.path.abspath(os.path.dirname(__file__)) + "/run_xtb.sh",
            xyzin,
            "--xyzout",
            str(xyzoutxtb),
            "--charge",
            str(charge),
            "--uhf",
            str(mult - 1),
        ]

        for i, _ in enumerate(constraints_all["Constraint_values"]):
            if constraints_all["Constraint_values"][i] is not None:
                for _, constraint in enumerate(constraints_all["Constraint_values"][i]):
                    command.append(constraints_all["Constraint_types"][i])
                    command.append(",".join([str(int(elem)) for elem in constraint]))

        subprocess.run(command)


def crest_opt(mol, name, dup_data, dup_data_idx, sdwriter, args, charge, mult):

    """
    Run xtb using shell script and args to perform crest/cregen conformer search

    Parameters
    ----------
    mol : RDKit Mol
        RDkit mol object
    name : str
        Name of molecule
    dup_data : pandas df
        csv
    dup_data_idx : int
        index of pandas df
    sdwriter : Rdkit writer
        sdwriter to write sdf file
    args : pyconfot args
        arguments for crest
    charge : charge
        charge
    mult : mult
        mult
    """

    csearch_dir = Path(args.w_dir_main)
    dat_dir = csearch_dir / "CSEARCH" / "crest_xyz" / name
    dat_dir.mkdir(exist_ok=True, parents=True)
    os.chdir(dat_dir)

    xyzin = str(dat_dir) + "/" + name + ".xyz"

    if not args.complex:
        AllChem.EmbedMolecule(mol)
        rdkit.Chem.rdmolfiles.MolToXYZFile(mol, xyzin)

    cbonds = args.cbonds
    constraints_dist = args.constraints_dist
    constraints_angle = args.constraints_angle
    constraints_dihedral = args.constraints_dihedral
    cregen = args.cregen
    cregen_ethr = args.cregen_ethr
    cregen_rthr = args.cregen_rthr
    cregen_bthr = args.cregen_bthr
    cregen_ewin = args.cregen_ewin

    xyzoutall = str(dat_dir) + "/" + name + "_conformers.xyz"
    xyzoutbest = str(dat_dir) + "/" + name + "_best.xyz"

    if args.complex:
        xyzoutxtb1 = str(dat_dir) + "/" + name + "_xtb1.xyz"
        xyzoutxtb2 = str(dat_dir) + "/" + name + "_xtb2.xyz"
        constraints_dist = get_constraint(mol, constraints_dist)
        constraints_all = {
            "Constraint_values": [
                constraints_dist,
                constraints_angle,
                constraints_dihedral,
            ],
            "Constraint_types": ["--dist", "--angle", "--dihedral"],
        }
        run_xtb(xyzin, xyzoutxtb1, constraints_all, charge, mult)
        run_xtb(xyzoutxtb1, xyzoutxtb2, constraints_all, charge, mult)
    else:
        xyzoutxtb2 = xyzin

    if os.path.exists(xyzoutall):
        print("   {0} already exists: skipping crest search".format(xyzoutall))
        # pass
    else:
        unique_atoms = []
        # list of list with types of contraints and number of elements (as indexes) that define the constrain
        # i.e. for a distance constrain, 3 elements are needed atom1, atom2 and distance
        contraint_types = [
            [constraints_dist, 2],
            [constraints_angle, 3],
            [constraints_dihedral, 4],
        ]
        for contraint_type in contraint_types:
            if contraint_type[0] is not None:
                for x in contraint_type[0]:
                    for i in x[: contraint_type[1]]:
                        if i not in unique_atoms:
                            unique_atoms.append(int(i))
        command = [
            os.path.abspath(os.path.dirname(__file__)) + "/run_crest.sh",
            xyzoutxtb2,
            "--xyzoutall",
            str(xyzoutall),
            "--xyzoutbest",
            str(xyzoutbest),
            "--charge",
            str(charge),
            "--uhf",
            str(mult - 1),
        ]
        if len(unique_atoms) != 0:
            command.append("--constraint")
            toadd = ",".join([str(elem) for elem in unique_atoms])
            command.append(toadd)
        if cbonds is not None:
            command.append("--cbonds")
            command.append(str(cbonds))

        subprocess.run(command)

    if cregen:
        xyzcregenensemble = xyzoutall.split(".xyz")[0] + "_cregen.xyz"
        xyzcregenbest = xyzoutbest.split(".xyz")[0] + "_cregen.xyz"

        # check whether job has already been run
        if os.path.exists(xyzcregenensemble):
            print(
                "   {0} already exists: skipping cregen search".format(
                    xyzcregenensemble
                )
            )
            # pass
        else:
            command = [
                os.path.abspath(os.path.dirname(__file__)) + "/run_cregen.sh",
                xyzoutall,
                xyzoutbest,
                "--xyzout",
                str(xyzcregenensemble),
                "--ethr",
                str(cregen_ethr),
                "--rthr",
                str(cregen_rthr),
                "--bthr",
                str(cregen_bthr),
                "--ewin",
                str(cregen_ewin),
            ]
            subprocess.run(command)

            # converting multiple xyz to single
            command_run_1 = [
                "obabel",
                xyzcregenensemble,
                "-oxyz",
                "-O" + xyzcregenbest,
                "-l",
                "1",
            ]
            subprocess.run(command_run_1)

    if cregen:
        xyzall_2_xyz(xyzcregenensemble, name)
    else:
        xyzall_2_xyz(xyzoutall, name)

    xyz_files = glob.glob(name + "_conf_*.xyz")
    for i, file in enumerate(xyz_files):
        name1 = file.split(".xyz")[0]
        command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name1 + ".sdf"]
        subprocess.call(command_xyz)
        # os.remove(file)

    sdf_files = glob.glob(name + "*.sdf")
    for file in sdf_files:
        mol = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        mol_rd = rdkit.Chem.RWMol(mol[0])
        energy = str(open(file, "r").readlines()[0])
        mol_rd.SetProp("Energy", energy)
        mol_rd.SetProp("Real charge", str(charge))
        mol_rd.SetProp("Mult", str(mult))
        sdwriter.write(mol_rd)
        # os.remove(file)

    dup_data.at[dup_data_idx, "crest-conformers"] = len(xyz_files)

    os.chdir(args.w_dir_main)

    return 1
