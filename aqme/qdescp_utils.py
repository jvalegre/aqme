######################################################.
#        This file stores QDESCP functions           #
######################################################.

import json
import sys
import numpy as np
import pandas as pd
import ast
import math
import rdkit
import warnings
warnings.filterwarnings('ignore')

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


def get_boltz_props(
    json_files,
    name,
    boltz_dir,
    type,
    self,
    mol_props,
    atom_props,
    nmr_atoms=None,
    nmr_slope=None,
    nmr_intercept=None,
    nmr_experim=None,
    mol=None,
):
    """
    Retrieves the properties from json files and gives Boltzmann averaged properties
    """

    if type == "nmr":
        if nmr_experim is not None:
            try:
                exp_data = pd.read_csv(nmr_experim)
            except:
                self.args.log.write(f'\nx  The CSV file with experimental NMR shifts specified ({nmr_experim}) was not found!')
                self.args.log.finalize()
                sys.exit()

    energy = []
    for k, json_file in enumerate(json_files):
        json_data = read_json(json_file)
        if type == "xtb":
            # filter off molecules with no atomic properties found when using the qdescp_atoms option
            for prop in atom_props:
                if prop not in json_data:
                    return None
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
    for prop in atom_props:
        prop_list = []
        for json_file in json_files:
            json_data = read_json(json_file)
            if len(self.args.qdescp_atoms) == 0:
                json_data['DBSTEP_Vbur'] = 'NaN'
            if type == "xtb":
                prop_list.append(json_data[prop])
            if type == "nmr":
                prop_list.append(json_data["properties"]["NMR"][prop].values())
        if len(self.args.qdescp_atoms) == 0:
            avg_prop = average_prop_atom(boltz, prop_list)
        else:
            avg_prop = average_prop_mol(boltz, prop_list)

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
                qdescp_nmr = nmr_experim.split(".csv")[0] + "_predicted.csv"
                exp_data.round(2).to_csv(qdescp_nmr, index=False)
                self.args.log.write(f"o  The {qdescp_nmr} file containing Boltzmann weighted NMR shifts was successfully created in {self.args.initial_dir}")

        elif type == "xtb":
            if len(self.args.qdescp_atoms) > 0 or avg_prop == 'NaN':
                avg_json_data[prop] = avg_prop
            else:
                avg_json_data[prop] = avg_prop.tolist()                

    if type == "xtb":
        for prop in mol_props:
            prop_list = []
            for json_file in json_files:
                json_data = read_json(json_file)
                prop_list.append(json_data[prop])
            avg_prop = average_prop_mol(boltz, prop_list)
            avg_json_data[prop] = avg_prop

    final_boltz_file = str(boltz_dir) + "/" + name + "_boltz.json"
    
    # calculate RDKit descriptors
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
        if p == 'NaN':
            boltz_avg = 'NaN'
            break
        boltz_avg.append([number * weights[i] for number in p])
    if boltz_avg == 'NaN':
        boltz_res = 'NaN'
    else:
        boltz_res = np.sum(boltz_avg, 0)
    return boltz_res


def average_prop_mol(weights, prop):
    """
    Returns Boltzmann averaged molecular properties
    """

    boltz_avg = 0.0
    for i, p in enumerate(prop):
        if p == 'NaN':
            boltz_avg = 'NaN'
            break
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


def read_fukui(file):
    """
    Read fukui output file created from XTB option. Return data.
    """
    f = open(file, "r")
    data = f.readlines()
    f.close()

    f_pos, f_negs, f_neutrals = [], [], []
    for i in range(0, len(data)):
        if data[i].find("f(+)") > -1:
            start = i + 1
            break
    for j in range(start, len(data)):
        if data[j].find("      -------------") > -1:
            end = j
            break

    fukui_data = data[start:end]

    for line in fukui_data:
        item = line.split()
        f_po = float(item[-3])
        f_neg = float(item[-2])
        f_neutral = float(item[-1])
        f_pos.append(f_po)
        f_negs.append(f_neg)
        f_neutrals.append(f_neutral)

    return f_pos, f_negs, f_neutrals


def read_gfn1(file):
    """
    Read fukui output file created from xTB. Return data.
    """

    if file.find(".gfn1") > -1:
        f = open(file, "r")
        data = f.readlines()
        f.close()

        for i in range(0, len(data)):
            if data[i].find("Mulliken/CM5 charges") > -1:
                start = i + 1
                break
        for j in range(start, len(data)):
            if (
                data[j].find("Wiberg/Mayer (AO) data") > -1
                or data[j].find("generalized Born model") > -1
            ):
                end = j - 1
                break

        pop_data = data[start:end]
        mulliken, cm5, s_prop, p_prop, d_prop = [], [], [], [], []
        for line in pop_data:
            item = line.split()
            q_mull = float(item[-5])
            q_cm5 = float(item[-4])
            s_prop_ind = float(item[-3])
            p_prop_ind = float(item[-2])
            d_prop_ind = float(item[-1])
            mulliken.append(q_mull)
            cm5.append(q_cm5)
            s_prop.append(s_prop_ind)
            p_prop.append(p_prop_ind)
            d_prop.append(d_prop_ind)

        return mulliken, cm5, s_prop, p_prop, d_prop


def read_wbo(file):
    """
    Read wbo output file created from xTB. Return data.
    """

    if file.find(".wbo") > -1:
        f = open(file, "r")
        data = f.readlines()
        f.close()

        bonds, wbos = [], []
        for line in data:
            item = line.split()
            bond = [int(item[0]), int(item[1])]
            wbo = float(item[2])
            bonds.append(bond)
            wbos.append(wbo)
        return bonds, wbos


def read_xtb(file):
    """
    Read xtb.out file. Return data.
    """

    f = open(file, "r")
    data = f.readlines()
    f.close()

    energy, homo_lumo, homo, lumo, atoms, numbers, chrgs, wbos = (
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    )
    dipole_module, Fermi_level, transition_dipole_moment = np.nan, np.nan, np.nan
    total_charge, total_SASA = np.nan, np.nan
    total_C6AA, total_C8AA, total_alpha = np.nan, np.nan, np.nan

    for i in range(0, len(data)):
        if data[i].find("SUMMARY") > -1:
            energy = float(data[i + 2].split()[3])
        if data[i].find("total charge") > -1:
            total_charge = int(float(data[i].split()[3]))
        if data[i].find("(HOMO)") > -1:
            if data[i].split()[3] != "(HOMO)":
                homo = float(data[i].split()[3])
                homo_occ = float(data[i].split()[1])
            else:
                homo = float(data[i].split()[2])
                homo_occ = 0
        if data[i].find("(LUMO)") > -1:
            if data[i].split()[3] != "(LUMO)":
                lumo = float(data[i].split()[3])
                lumo_occ = float(data[i].split()[1])
            else:
                lumo = float(data[i].split()[2])
                lumo_occ = 0
        homo_lumo = float(lumo - homo)
        if data[i].find("molecular dipole:") > -1:
            dipole_module = float(data[i + 3].split()[-1])
        if data[i].find("transition dipole moment") > -1:
            transition_dipole_moment = float(data[i + 2].split()[-1])
        if data[i].find("Fermi-level") > -1:
            Fermi_level = float(data[i].split()[-2])

    # get atomic properties related to charges, dispersion, etc
    start, end = 0, 0
    for j in range(0, len(data)):
        if data[j].find("#   Z          covCN") > -1:
            start = j + 1
            break
    for k in range(start, len(data)):
        if data[k].find("Mol. ") > -1:
            end = k - 1
            total_C6AA = float(data[k].split()[-1])
            total_C8AA = float(data[k + 1].split()[-1])
            total_alpha = float(data[k + 2].split()[-1])
            break

    chrg_data = data[start:end]
    atoms, numbers, chrgs = [], [], []
    covCN, C6AA, alpha = [], [], []
    for line in chrg_data:
        item = line.split()
        numbers.append(int(item[0]))
        atoms.append(item[2])
        covCN.append(float(item[3]))
        chrgs.append(float(item[4]))
        C6AA.append(float(item[5]))
        alpha.append(float(item[6]))

    # get atomic properties related to solvent
    start_solv, end_solv = 0, 0
    for j in range(0, len(data)):
        if data[j].find("#   Z     Born rad") > -1:
            start_solv = j + 1
            break
    for k in range(start_solv, len(data)):
        if data[k].find("total SASA ") > -1:
            end_solv = k - 1
            total_SASA = float(data[k].split()[-1])
            break

    solv_data = data[start_solv:end_solv]
    born_rad, SASA, h_bond = [], [], []
    for line in solv_data:
        item = line.split()
        born_rad.append(float(item[3]))
        SASA.append(float(item[4]))
        # in apolar solvents such as CH2Cl2, xTB doesn't return any H bond parameters
        try:
            h_bond.append(float(item[5]))
        except IndexError:
            h_bond.append(float(0))

    return (
        energy,
        total_charge,
        homo_lumo,
        homo,
        lumo,
        atoms,
        numbers,
        chrgs,
        dipole_module,
        Fermi_level,
        transition_dipole_moment,
        covCN,
        C6AA,
        alpha,
        homo_occ,
        lumo_occ,
        born_rad,
        SASA,
        h_bond,
        total_SASA,
        total_C6AA,
        total_C8AA,
        total_alpha,
    )


def read_json(file):
    """
    Takes json files and parses data into pandas table. Returns data.
    """

    if file.find(".json") > -1:
        f = open(file, "r")  # Opening JSON file
        data = json.loads(f.read())  # read file
        f.close()
        return data
    else:
        pass


def read_fod(file):
    """
    Read xtb.fod files. Return FOD-related properties.
    """

    f = open(file, "r")
    data = f.readlines()
    f.close()

    # get fractional occupation density (FOD)
    for j in range(0, len(data)):
        if data[j].find("Loewdin FODpop") > -1:
            start_fod = j + 1
            total_fod = float(data[j - 2].split()[-1])
            break
    for k in range(start_fod, len(data)):
        if data[k].find("Wiberg/Mayer") > -1:
            end_fod = k - 1
            break

    fod_data = data[start_fod:end_fod]
    fod, s_prop_fod, p_prop_fod, d_prop_fod = [], [], [], []
    for line in fod_data:
        item = line.split()
        fod.append(float(item[1]))
        s_prop_fod.append(float(item[2]))
        p_prop_fod.append(float(item[3]))
        d_prop_fod.append(float(item[4]))

    return total_fod, fod, s_prop_fod, p_prop_fod, d_prop_fod