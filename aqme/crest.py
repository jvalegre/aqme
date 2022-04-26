#!/usr/bin/python
from __future__ import print_function, absolute_import
from rdkit.Chem import AllChem as Chem

#######################################################################
# Runs crest on an xyz file, can add in options to the code if needed #
#######################################################################

# Python Libraries
import os
import glob
import subprocess
import rdkit
from pathlib import Path
from aqme.utils import read_file


def xyzall_2_xyz(xyzin, name):
    # converting multiple xyz to single
    command_run_1 = ["obabel", xyzin, "-oxyz", "-O" + name + "_conf_.xyz", "-m"]
    subprocess.run(command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def crest_opt(name, dup_data, dup_data_idx, args, charge, mult, constraints_atoms, constraints_dist, constraints_angle, constraints_dihedral):

    """
    Run xTB using subprocess to perform CREST/CREGEN conformer sampling
    """

    name_no_path = name.replace('/','\\').split("\\")[-1].split('.')[0]
    csearch_dir = Path(args.w_dir_main)
    dat_dir = csearch_dir / "CSEARCH" / "crest_xyz" / name_no_path
    dat_dir.mkdir(exist_ok=True, parents=True)

    xyzin = f'{dat_dir}/{name_no_path}.xyz'
    sdwriter = Chem.SDWriter(str(f'{dat_dir}/{name_no_path}'))

    os.rename(f'{name}.xyz', xyzin)

    os.chdir(dat_dir)
    xyzoutall = str(dat_dir) + "/" + name_no_path + "_conformers.xyz"

    constrained_sampling = create_xcontrol(args,constraints_atoms,constraints_dist,constraints_angle,constraints_dihedral,xyzin)

    command = [
        "crest",
        xyzin,
        "--chrg",
        str(charge),
        "--uhf",
        str(mult - 1),
        "-T",
        str(args.nprocs)
    ]

    if constrained_sampling:
        command.append('-cinp')
        command.append('.xcontrol.sample')

    if args.crest_keywords is not None:
        for keyword in args.crest_keywords.split():
            if keyword not in command:
                command.append(keyword)

    subprocess.run(command)

    if args.cregen:
        command = [
        "crest",
        xyzin,
        "--cregen",
        'crest_conformers.xyz']

        if args.cregen_keywords is not None:
            for keyword in args.cregen_keywords.split():
                if keyword not in command:
                    command.append(keyword)

        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    try:
        if args.cregen:
            os.rename(str(dat_dir) + "/crest_ensemble.xyz", xyzoutall)
        else:
            os.rename(str(dat_dir) + "/crest_conformers.xyz", xyzoutall)
    except FileNotFoundError:
        args.log.write('x   CREST conformer sampling failed! Please, try other options (i.e. include constrains, change the crest_keywords option, etc.)')

    xyzall_2_xyz(xyzoutall, name_no_path)

    xyz_files = glob.glob(name_no_path + "_conf_*.xyz")
    for _, file in enumerate(xyz_files):
        name_conf = file.split(".xyz")[0]
        command_xyz = ["obabel", "-ixyz", file, "-osdf", "-O" + name_conf + ".sdf"]
        subprocess.run(command_xyz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    sdf_files = glob.glob(name_no_path + "*.sdf")
    for file in sdf_files:
        mol = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
        mol_rd = rdkit.Chem.RWMol(mol[0])
        energy = str(open(file, "r").readlines()[0])
        mol_rd.SetProp("Energy", energy)
        mol_rd.SetProp("Real charge", str(charge))
        mol_rd.SetProp("Mult", str(mult))
        sdwriter.write(mol_rd)

    dup_data.at[dup_data_idx, "crest-conformers"] = len(xyz_files)

    os.chdir(args.w_dir_main)

    return 1


def create_xcontrol(args,constraints_atoms,constraints_dist,constraints_angle,constraints_dihedral,xyzin):
    '''
    Function to create the .xcontrol.sample if constraints are defined
    '''

    constrained_sampling = False

    # this avoids problems when running AQME through command lines
    if not isinstance(constraints_atoms, list):
        constraints_atoms = constraints_atoms.strip('][').split(',')

    if not isinstance(constraints_dist, list):
        constraints_dist = [constraints_dist.strip('][').split(',')]

    if not isinstance(constraints_angle, list):
        constraints_angle = [constraints_angle.strip('][').split(',')]

    if not isinstance(constraints_dihedral, list):
        constraints_dihedral = [constraints_dihedral.strip('][').split(',')]

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
        subprocess.run(['crest', xyzin, '--constrain', '1'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.remove('.xcontrol.sample')

        # add the constraints part
        edited_xcontrol = '$constrain\n'

        if constraints_atoms != []:
            edited_xcontrol += 'atoms: '
            edited_xcontrol += ','.join(atom_idx for atom_idx in constraints_atoms)

        for constraint_type in [constraints_dist,constraints_angle,constraints_dihedral]:
            if constraint_type != []:
                for const in constraint_type: 
                    if constraint_type == constraints_dist:
                        edited_xcontrol += 'distance: '
                    elif constraint_type == constraints_angle:
                        edited_xcontrol += 'angle: '
                    elif constraint_type == constraints_dihedral:  
                        edited_xcontrol += 'dihedral: '
                    edited_xcontrol += ','.join(val for val in const)
                    edited_xcontrol += '\n'

        edited_xcontrol += f'force constant={args.crest_force}\n'
        edited_xcontrol += 'reference=coord.ref\n'

        # metadyn part
        outlines = read_file(os.getcwd(),os.getcwd(),xyzin)
        n_atoms = int(outlines[0])
        edited_xcontrol += '$metadyn\n'
        edited_xcontrol += 'atoms: '
        for atom_idx in range(1,n_atoms+1):
            if atom_idx not in unique_atoms:
                if atom_idx == n_atoms:
                    edited_xcontrol += f'{atom_idx}\n'
                else:
                    edited_xcontrol += f'{atom_idx},'
        edited_xcontrol += '$end\n'

        # write the file
        xcontrol_file = open('.xcontrol.sample', 'w')
        xcontrol_file.write(edited_xcontrol)
        xcontrol_file.close()
    
    return constrained_sampling