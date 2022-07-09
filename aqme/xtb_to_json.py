#!/opt/anaconda/envs/DL_CPU/bin/python

##########################################################################################
####            Authors: Liliana C. Gallegos and Juan V. Alegre-Requena               ####
####                        For any questions, contact:                               ####
####             LilianaC.Gallegos@colostate.edu or jvalegre@unizar.es                ####
##########################################################################################

from __future__ import print_function
import numpy as np
import json

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
    Read fukui output file created from XTB option. Return data.
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
    Read wbo output file created from XTB option. Return data.
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
        if data[i].find("charge                     :") > -1:
            total_charge = int(data[i].split()[-1])
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
    Takes json file and parses data into pandas table. Returns data.
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
