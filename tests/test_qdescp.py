#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#                   QDESCP module                    #
######################################################.

import os
import pytest
import pandas as pd
import numpy as np
from aqme.qdescp import qdescp
from aqme.csearch import csearch
from aqme.qdescp_utils import read_json
import glob
import math

# saves the working directory
w_dir_main = os.getcwd()
qdescp_input_dir = w_dir_main + "/tests/qdescp_inputs"

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15

# tests for QDESCP-xTB
@pytest.mark.parametrize(
    "file",
    [
        ("test.csv")
    ]
)

def test_qdescp_xtb(file):

    os.chdir(qdescp_input_dir)
    # conf gen
    csearch(program='rdkit',input=file)
    sdf_rdkit_files = glob.glob(f'CSEARCH/*.sdf')
    # QDESCP-xTB workflow
    qdescp(program='xtb',files=sdf_rdkit_files)

    # check if various parameters are stored correctly from xTB calculations to the generated json files
    folder_qdescp = f'{qdescp_input_dir}/QDESCP'
    os.chdir(folder_qdescp)

    file_1 = "mol_1_rdkit_conf_1"
    file_2 = "mol_1_rdkit_conf_2"
    files_xtb = [file_1, file_2]
    energies_xtb, Fermi_lvls_xtb = [], []
    energies_json, Fermi_lvls_json = [], []

    for file in files_xtb:
        # get data from the raw xTB calculations
        f = open(f'{qdescp_input_dir}/QDESCP/xtb_data/{file}/{file}.out', "r")
        data = f.readlines()
        f.close()
        energy_found, fermi_found = False, False
        for i in range(0, len(data)):
            if data[i].find("SUMMARY") > -1 and not energy_found:
                energies_xtb.append(float(data[i + 2].split()[3]))
                energy_found = True
            if data[i].find("Fermi-level") > -1:
                Fermi_lvls_xtb.append(float(data[i].split()[-2]))
                fermi_found = True
            if energy_found and fermi_found:
                break
        # get data from the generated json files
        json_data = read_json(f'{file}.json')
        energies_json.append(json_data["total energy"])
        Fermi_lvls_json.append(json_data["Fermi-level/eV"])
    
    for i,_ in enumerate(energies_xtb):
        assert round(energies_xtb[i],4) == round(energies_json[i],4)
        assert round(Fermi_lvls_xtb[i],2) == round(Fermi_lvls_json[i],2)

    # calculate Boltzmann averaged values
    energ = [number - min(energies_xtb) for number in energies_xtb]
    boltz_sum = 0.0
    for e in energ:
        boltz_sum += math.exp(-e * J_TO_AU / GAS_CONSTANT / T)
    weights = []
    for e in energ:
        weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / T) / boltz_sum
        weights.append(weight)

    energy_boltz_calc = 0.0
    for i, p in enumerate(energies_xtb):
        energy_boltz_calc += p * weights[i]

    Fermi_lvl_boltz_calc = 0.0
    for i, p in enumerate(Fermi_lvls_xtb):
        Fermi_lvl_boltz_calc += p * weights[i]

    # retrieve Boltzman avg values from json files
    folder_boltz = f'{folder_qdescp}/boltz'
    os.chdir(folder_boltz)
    json_data = read_json('mol_1_rdkit_boltz.json')
    energy_boltz_file = json_data["total energy"]
    Fermi_lvl_boltz_file = json_data["Fermi-level/eV"]

    assert round(energy_boltz_calc,4) == round(energy_boltz_file,4)
    assert round(Fermi_lvl_boltz_calc,2) == round(Fermi_lvl_boltz_file,2)

    # retrieve Boltzman avg values and RDKit descriptors from the generated csv file
    os.chdir(qdescp_input_dir)
    pd_boltz = pd.read_csv('QDESCP_boltz_descriptors.csv')
    energy_boltz_csv = pd_boltz["total energy"][0]
    Fermi_lvl_boltz_csv = pd_boltz["Fermi-level/eV"][0]

    assert round(energy_boltz_calc,4) == round(energy_boltz_csv,4)
    assert round(Fermi_lvl_boltz_calc,2) == round(Fermi_lvl_boltz_csv,2)
    assert pd_boltz["NumRotatableBonds"][0] == 3
    assert pd_boltz["NumRotatableBonds"][1] == 4

    os.chdir(w_dir_main)

# tests for QDESCP-NMR
@pytest.mark.parametrize(
    "json_files",
    [
        ("*.json")
    ]
)

def test_qdescp_nmr(json_files):

    os.chdir(qdescp_input_dir)
    json_files = glob.glob(json_files)
    nmr_atoms = [6, 1]  # [C,H]
    nmr_slope=[-1.0537, -1.0784]
    nmr_intercept=[181.7815,31.8723]

    # QDESCP-NMR workflow
    qdescp(program='nmr',files=json_files,nmr_slope=nmr_slope,
            nmr_intercept=nmr_intercept,nmr_experim="Experimental_NMR_shifts.csv")

    # read json files and calculate Boltzmann shifts
    energies_json,nmr_json = [],[]
    for json_file in json_files:
        json_data = read_json(json_file)
        energies_json.append(json_data["optimization"]["scf"]["scf energies"][-1])
        # retrieves and scales NMR shifts from json files
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
        nmr_json.append(shifts)

    # calculate Boltzmann averaged values
    energ = [number - min(energies_json) for number in energies_json]
    boltz_sum = 0.0
    for e in energ:
        boltz_sum += math.exp(-e * J_TO_AU / GAS_CONSTANT / T)
    weights = []
    for e in energ:
        weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / T) / boltz_sum
        weights.append(weight)

    boltz_avg = []
    for i, p in enumerate(nmr_json):
        boltz_avg.append([number * weights[i] for number in p.values()])
    nmr_json_calc = np.sum(boltz_avg, 0)

    # check that the averaged NMR shifts are the same as the values included in the json file
    folder_boltz = f'{qdescp_input_dir}/QDESCP/boltz'
    os.chdir(folder_boltz)
    json_data_boltz = read_json('test_boltz.json')
    nmr_json_file = list(json_data_boltz["NMR Chemical Shifts"].values())
    for i,_ in enumerate(nmr_json_calc):
        assert round(nmr_json_calc[i],2) == round(nmr_json_file[i],2)

    # check that the averaged NMR shifts are the same as the values included in the csv file
    folder_boltz = f'{qdescp_input_dir}/QDESCP/boltz'
    os.chdir(qdescp_input_dir)
    pd_boltz = pd.read_csv('Experimental_NMR_shifts_predicted.csv')
    nmr_boltz_csv = pd_boltz["boltz_avg"]
    nmr_json_calc_1H = nmr_json_calc[21:]
    for i,_ in enumerate(nmr_json_calc_1H):
        assert round(nmr_json_calc_1H[i],2) == round(nmr_boltz_csv[i],2)

    # check that the shift errors in the csv file are correct
    nmr_experim_csv = pd_boltz["experimental_ppm"]
    error_calc = abs(nmr_boltz_csv - nmr_experim_csv)
    error_csv = pd_boltz["error_boltz"]
    print('nmr_json_calc',error_calc)
    print('nmr_boltz_csv',error_csv)
    for i,_ in enumerate(error_calc):
        if str(error_calc[i]) not in ['nan']:
            assert round(error_calc[i],2) == round(error_csv[i],2)
