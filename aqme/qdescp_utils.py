######################################################.
#        This file stores QDESCP functions           #
######################################################.

import json
import numpy as np
import pandas as pd
import ast
import json
import math
from aqme.xtb_to_json import read_json
import rdkit

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15


def get_boltz(energy):
    """
    Calculates the Boltzmann weights for a list of energies
    """

    energ = [number - min(energy) for number in energy]

    boltz_sum = 0.0
    for e in energ:
        boltz_sum += math.exp(-e * J_TO_AU / GAS_CONSTANT / T)
    weights = []

    for e in energ:
        weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / T) / boltz_sum
        weights.append(weight)

    return weights


def get_boltz_avg_properties_xtb(
    json_files,
    name,
    boltz_dir,
    type,
    nmr_atoms=None,
    nmr_slope=None,
    nmr_intercept=None,
    nmr_experim=None,
    mol=None,
):
    """
    Retrieves the properties from json files and gives Boltzmann averaged properties
    """

    if type == "xtb":
        mol_prop = [
            "total energy",
            "HOMO-LUMO gap/eV",
            "electronic energy",
            "Dipole module/D",
            "Total charge",
            "HOMO",
            "LUMO",
            "Fermi-level/eV",
            "Total dispersion C6",
            "Total dispersion C8",
            "Total polarizability alpha",
            "Total FOD",
        ]
        atom_prop = [
            "dipole",
            "partial charges",
            "mulliken charges",
            "cm5 charges",
            "FUKUI+",
            "FUKUI-",
            "FUKUIrad",
            "s proportion",
            "p proportion",
            "d proportion",
            "Coordination numbers",
            "Dispersion coefficient C6",
            "Polarizability alpha",
            "FOD",
            "FOD s proportion",
            "FOD p proportion",
            "FOD d proportion",
        ]
    elif type == "nmr":
        atom_prop = [
            "NMR Chemical Shifts",
        ]
        if nmr_experim is not None:
            exp_data = pd.read_csv(nmr_experim)

    energy = []

    for k, json_file in enumerate(json_files):
        json_data = read_json(json_file)
        if type == "xtb":
            energy.append(json_data["total energy"])
        elif type == "nmr":
            energy.append(json_data["optimization"]["scf"]["scf energies"][-1])

            json_data["properties"]["NMR"]["NMR Chemical Shifts"] = get_chemical_shifts(
                json_data, nmr_atoms, nmr_slope, nmr_intercept
            )
            if nmr_experim is not None:
                list_shift = json_data["properties"]["NMR"]["NMR Chemical Shifts"]
                df = pd.DataFrame(
                    list_shift.items(),
                    columns=["atom_idx", "conf_{}".format(k + 1)],
                )
                df["atom_idx"] = df["atom_idx"] + 1
                exp_data = exp_data.merge(df, on=["atom_idx"])
        with open(json_file, "w") as outfile:
            json.dump(json_data, outfile)

    boltz = get_boltz(energy)

    avg_json_data = {}
    for prop in atom_prop:
        prop_list = []
        for json_file in json_files:
            json_data = read_json(json_file)
            if type == "xtb":
                prop_list.append(json_data[prop])
            if type == "nmr":
                prop_list.append(json_data["properties"]["NMR"][prop].values())
        avg_prop = average_prop_atom(boltz, prop_list)

        if type == "nmr":
            dictavgprop = {}
            for i, key in enumerate(json_data["properties"]["NMR"][prop].keys()):
                dictavgprop[key] = avg_prop[i]
            avg_json_data[prop] = dictavgprop

            if nmr_experim is not None:
                list_shift = avg_json_data[prop]
                df = pd.DataFrame(
                    list_shift.items(),
                    columns=["atom_idx", "boltz_avg"],
                )
                df["atom_idx"] = df["atom_idx"].astype(int) + 1
                exp_data = exp_data.merge(df, on=["atom_idx"])
                exp_data["error_boltz"] = abs(
                    exp_data["experimental_ppm"] - exp_data["boltz_avg"]
                )
                exp_data.round(2).to_csv(
                    nmr_experim.split(".csv")[0] + "_predicted.csv", index=False
                )
        elif type == "xtb":
            avg_json_data[prop] = avg_prop.tolist()

    if type == "xtb":
        for prop in mol_prop:
            prop_list = []
            for json_file in json_files:
                json_data = read_json(json_file)
                prop_list.append(json_data[prop])
            avg_prop = average_prop_mol(boltz, prop_list)
            avg_json_data[prop] = avg_prop

    final_boltz_file = str(boltz_dir) + "/" + name + "_boltz.json"
    if mol is not None:
        avg_json_data = get_rdkit_properties(avg_json_data, mol)
    with open(final_boltz_file, "w") as outfile:
        json.dump(avg_json_data, outfile)


def get_chemical_shifts(json_data, nmr_atoms, nmr_slope, nmr_intercept):
    """
    Retrieves and scales NMR shifts from json files
    """

    if not isinstance(nmr_atoms, list):
        nmr_atoms = ast.literal_eval(nmr_atoms)
    if not isinstance(nmr_slope, list):
        nmr_slope = ast.literal_eval(nmr_slope)
    if not isinstance(nmr_intercept, list):
        nmr_intercept = ast.literal_eval(nmr_intercept)

    atoms = json_data["atoms"]["elements"]["number"]
    tensor = json_data["properties"]["NMR"]["NMR isotopic tensors"]
    shifts = {}
    i = 0
    for atom, ten in zip(atoms, tensor):
        if atom in nmr_atoms:
            # assigning values from arrays
            index = nmr_atoms.index(atom)
            slope_nuc = nmr_slope[index]
            intercept_nuc = nmr_intercept[index]

            scaled_nmr = (intercept_nuc - ten) / (-slope_nuc)
            shifts[i] = scaled_nmr
        else:
            pass
        i += 1

    return shifts


def average_prop_atom(weights, prop):
    """
    Returns Boltzmann averaged atomic properties
    """

    boltz_avg = []
    for i, p in enumerate(prop):
        boltz_avg.append([number * weights[i] for number in p])
    boltz_res = np.sum(boltz_avg, 0)
    return boltz_res


def average_prop_mol(weights, prop):
    """
    Returns Boltzmann averaged molecular properties
    """

    boltz_avg = 0.0
    for i, p in enumerate(prop):
        boltz_avg += p * weights[i]
    return boltz_avg


def get_rdkit_properties(avg_json_data, mol):
    """
    Calculates RDKit molecular descriptors
    """

    avg_json_data["NHOHCount"] = rdkit.Chem.Lipinski.NHOHCount(mol)
    avg_json_data["FractionCSP3"] = rdkit.Chem.Lipinski.FractionCSP3(mol)
    avg_json_data["NOCount"] = rdkit.Chem.Lipinski.NOCount(mol)
    avg_json_data["NumAliphaticRings"] = rdkit.Chem.Lipinski.NumAliphaticRings(mol)
    avg_json_data["NumAromaticRings"] = rdkit.Chem.Lipinski.NumAromaticRings(mol)
    avg_json_data["NumHAcceptors"] = rdkit.Chem.Lipinski.NumHAcceptors(mol)
    avg_json_data["NumHDonors"] = rdkit.Chem.Lipinski.NumHDonors(mol)
    avg_json_data["NumHeteroatoms"] = rdkit.Chem.Lipinski.NumHeteroatoms(mol)
    avg_json_data["NumRotatableBonds"] = rdkit.Chem.Lipinski.NumRotatableBonds(mol)

    avg_json_data["TPSA"] = rdkit.Chem.Descriptors.TPSA(mol)
    avg_json_data["MolLogP"] = rdkit.Chem.Descriptors.MolLogP(mol)
    # avg_json_data["NumAmideBonds"] = rdkit.Chem.Descriptors.NumAmideBonds(mol)

    return avg_json_data