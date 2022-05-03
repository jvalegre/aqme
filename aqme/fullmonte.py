#####################################################.
#        This file stores all the functions         #
#     used for genrating details for fullmonte      #
#####################################################.

import sys
import numpy as np
import math
import random
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, rdMolAlign

from aqme.utils import set_metal_atomic_number, get_conf_RMS
from aqme.csearch_utils import minimize_rdkit_energy


def realign_mol(
    mol, conf, coord_Map, alg_Map, mol_template, maxsteps
):  # RAUL: This function requires a clear separation between minimization and alignment.
    """
    Minimizes and aligns the molecule provided freezing the atoms that match the mol_template

    Parameters
    ----------
    mol : RDKit mol object
        Molecule to be minimized and aligned
    conf : int
        Number that indicates which conformation of the molecule will be minimized and aligned
    coord_Map : [type]
        [description]
    alg_Map : [type]
        [description]
    mol_template : [type]
        [description]
    maxsteps : int
        Maximum number of iterations in FF minimization

    Returns
    -------
    mol,energy
        The updated mol object and the final forcefield energy.
    """

    num_atom_match = mol.GetSubstructMatch(mol_template)
    forcefield = Chem.UFFGetMoleculeForceField(mol, confId=conf)
    for i, idxI in enumerate(num_atom_match):
        for idxJ in num_atom_match[i + 1 :]:
            d = coord_Map[idxI].Distance(coord_Map[idxJ])
            forcefield.AddDistanceConstraint(idxI, idxJ, d, d, 10000)
    forcefield.Initialize()
    forcefield.Minimize(maxIts=maxsteps)
    # rotate the embedded conformation onto the core_mol:
    rdMolAlign.AlignMol(
        mol,
        mol_template,
        prbCid=conf,
        refCid=-1,
        atomMap=alg_Map,
        reflect=True,
        maxIters=100,
    )
    energy = float(forcefield.CalcEnergy())
    return mol, energy


def rotate_dihedrals(conformer, dihedrals, seed, stepsize):
    """
    Applies a random rotation to all the dihedrals

    Parameters
    ----------
    conformer : rdkit.Chem.rdchem.Conformer
        The conformer whose angles are going to be rotated (conformer = mol.GetConformer(cid))
    dihedrals : list
        A list of tuples of all the dihedrals that are going to be rotated.
    seed : int
        seed for the random module
    stepsize : float
        Angle in Degrees to do the steps between 0.0 and 360.0
    """

    rad_range = np.arange(stepsize, 360.0, stepsize)
    for dihedral in dihedrals:
        random.seed(seed)  # RAUL: Any good reason to keep reseting the seed?
        rad_ang = random.choice(rad_range)
        rad = math.pi * rad_ang / 180.0
        rdMolTransforms.SetDihedralRad(conformer, *dihedral, value=rad)


def generating_conformations_fullmonte(
    name,
    args,
    rotmatches,
    selectedcids_rdkit,
    outmols,
    sdwriter,
    dup_data,
    dup_data_idx,
    coord_Map,
    alg_Map,
    mol_template,
    ff,
):

    ##working with fullmonte
    n_unique_conformers = len(selectedcids_rdkit)
    args.log.write(
        f"\no  Generation of confomers using FULLMONTE using "
        f"{n_unique_conformers} unique conformer(s) as starting point(s)"
    )

    # Writing the conformers as mol objects to sdf
    sdtemp = Chem.SDWriter(name + "_" + "rdkit" + args.output)
    for conf in selectedcids_rdkit:
        sdtemp.write(outmols[conf], conf)
    sdtemp.close()

    fmmols = Chem.SDMolSupplier(name + "_" + "rdkit" + args.output, removeHs=False)
    if fmmols is None:
        args.log.write("Could not open " + name + args.output)
        args.log.finalize()
        sys.exit()

    # array for each each unique from rdkit
    unique_mol, c_energy, unique_mol_sample = [], [], []

    # STEP 1: Use start conformation for and append to unique list
    nsteps = 1
    for mol_fm in fmmols:
        unique_mol.append(mol_fm)
        c_energy.append(float(mol_fm.GetProp("Energy")))

    # defining unique mol sample for choosing
    globmin = min(c_energy)
    for ene in reversed(c_energy):
        if abs(globmin - ene) < args.ewin_sample_fullmonte:
            unique_mol_sample.append(unique_mol[c_energy.index(ene)])

    while nsteps < args.nsteps_fullmonte + 1:
        seed = nsteps

        # STEP 2: Choose mol object form unique_mol:
        random.seed(seed)
        mol_rot = random.choices(unique_mol_sample, k=1)[0]

        # updating the location of mol object i.e., the hexadecimal locaiton to a new one so the older one isnt affected
        mol = Chem.RWMol(mol_rot)
        rot_mol = mol.GetMol()

        # STEP 3: Choose random subset of dihedral from rotmatches
        random.seed(seed)  # RAUL: Any good reason to keep reseting the seed ?
        k = min(len(rotmatches), args.nrot_fullmonte)
        mutable_dihedrals = random.choices(rotmatches, k=k)

        # STEP 4: for the given conformation, then apply a random rotation to each torsion in the subset
        conformer = rot_mol.GetConformer()
        rotate_dihedrals(conformer, mutable_dihedrals, seed, args.ang_fullmonte)

        # STEP 5: Optimize geometry rot_mol
        if (coord_Map, alg_Map, mol_template) == (None, None, None):
            energy = minimize_rdkit_energy(
                rot_mol, -1, args.log, ff, args.opt_steps_rdkit
            )
        else:
            mol, energy = realign_mol(
                rot_mol, -1, coord_Map, alg_Map, mol_template, args.opt_steps_rdkit
            )

        # STEP 6 : Check for DUPLICATES - energy and rms filter (reuse)
        #  if the conformer is unique then save it the list
        exclude_conf = False
        # compare against allprevious conformers located
        for j, seenmol in enumerate(unique_mol):
            if abs(energy - c_energy[j]) < args.initial_energy_threshold:
                exclude_conf = True
                break
            if abs(energy - c_energy[j]) < args.energy_threshold:
                rms = get_conf_RMS(
                    rot_mol, seenmol, -1, -1, args.heavyonly, args.max_matches_rmsd
                )
                if rms < args.rms_threshold:
                    exclude_conf = True
                    break
        if not exclude_conf:
            unique_mol.append(rot_mol)
            c_energy.append(energy)
            unique_mol[c_energy.index(energy)].SetProp("Energy", str(energy))

        unique_mol_sample = []
        # STEP 7: ANALYSE THE UNIQUE list for lowest energy, reorder the uniques if greater the given thershold remove
        globmin = min(c_energy)
        for ene in reversed(c_energy):
            indx = c_energy.index(ene)
            if abs(globmin - ene) > args.ewin_fullmonte:
                unique_mol.pop(indx)
                c_energy.pop(indx)
            if abs(globmin - ene) < args.ewin_sample_fullmonte:
                unique_mol_sample.append(unique_mol[indx])

        nsteps += 1

    dup_data.at[dup_data_idx, "FullMonte-Unique-conformers"] = len(unique_mol)

    if args.verbose:
        args.log.write("o  " + str(len(unique_mol)) + " unique conformers remain")

    cids = list(range(len(unique_mol)))
    sorted_all_cids = sorted(cids, key=lambda cid: c_energy[cid])

    # STEP 9: WRITE FINAL uniques to sdf for xtb or ani
    for i, cid in enumerate(sorted_all_cids):
        unique_mol[cid].SetProp("_Name", name + " " + str(i))
        if coord_Map is None and alg_Map is None and mol_template is None:
            if args.metal_complex:
                set_metal_atomic_number(unique_mol[cid], args.metal_idx, args.metal_sym)
            sdwriter.write(unique_mol[cid])
        else:
            mol_realigned, _ = realign_mol(
                unique_mol[cid],
                -1,
                coord_Map,
                alg_Map,
                mol_template,
                args.opt_steps_rdkit,
            )
            if args.metal_complex:
                set_metal_atomic_number(mol_realigned, args.metal_idx, args.metal_sym)
            sdwriter.write(mol_realigned)

    status = 1

    return status
