#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#                   QDESCP module                    #
######################################################.

import os
import pytest
import pandas as pd
import numpy as np
import subprocess
from aqme.qdescp_utils import read_json
import glob
import math
import shutil

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
        ("test.csv"), # standard test
        ("test_atom.csv"), # test with qdescp_atoms using an atom
        ("test_group.csv"), # test with qdescp_atoms using a functional group
        ("test_multigroup.csv"), # test with qdescp_atoms using a multiple atoms and functional groups
        ("test_robert.csv") # test for the AQME-ROBERT workflow
    ]
)

def test_qdescp_xtb(file):

    # reset folder and files
    folder_csearch = f'{qdescp_input_dir}/CSEARCH'
    folder_qdescp = f'{qdescp_input_dir}/QDESCP'
    folder_boltz = f'{folder_qdescp}/boltz'
    for folder in [folder_csearch,folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    if file in ['test_multigroup.csv','test_robert.csv']:
        file_csearch = 'test_atom.csv'
    else:
        file_csearch = file
    if file == 'test_robert.csv':
        file_descriptors = f'{w_dir_main}/AQME-ROBERT_{file_csearch}'
    else:
        file_descriptors = f'{w_dir_main}/QDESCP_boltz_descriptors.csv'
    if os.path.exists(file_descriptors):
        os.remove(file_descriptors)

    # CSEARCH conformer generation
    cmd_csearch = [
    "python",
    "-m",
    "aqme",
    "--csearch",
    "--program",
    "rdkit",
    "--input",
    f'{qdescp_input_dir}/{file_csearch}',
    "--destination",
    f'{folder_csearch}',
    ]
    subprocess.run(cmd_csearch)

    # QDESCP-xTB workflow
    cmd_qdescp = [
    "python",
    "-m",
    "aqme",
    "--qdescp",
    "--program",
    "xtb",
    "--files",
    f'{folder_csearch}/*.sdf',
    "--destination",
    f'{folder_qdescp}',
    ]

    if file == 'test_atom.csv':
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[P]", "--dbstep_calc"]
    elif file == 'test_group.csv':
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[C=O]"]
    elif file in ['test_multigroup.csv','test_robert.csv']:
        # Pd is included to check the try/except in the SMARTS pattern match,
        # and to check if the code works even if there are atoms that aren't used
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[P,CC,Pd]"]
        if file == 'test_robert.csv':
            cmd_qdescp = cmd_qdescp + ["--csv_name", f'{qdescp_input_dir}/{file_csearch}']

    subprocess.run(cmd_qdescp)

    # check if various parameters are stored correctly from xTB calculations to the generated json files
    if file == 'test.csv':
        file_1 = "mol_1_rdkit_conf_1"
        file_2 = "mol_1_rdkit_conf_2"
        files_xtb = [file_1, file_2]
        energies_xtb, Fermi_lvls_xtb = [], []
        energies_json, Fermi_lvls_json = [], []

        for file in files_xtb:
            # get data from the raw xTB calculations
            f = open(f'{folder_qdescp}/xtb_data/{file}/{file}.out', "r")
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
            json_data = read_json(f'{folder_qdescp}/{file}.json')
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
        json_data = read_json(f'{folder_boltz}/mol_1_rdkit_boltz.json')
        energy_boltz_file = json_data["total energy"]
        Fermi_lvl_boltz_file = json_data["Fermi-level/eV"]

        assert round(energy_boltz_calc,4) == round(energy_boltz_file,4)
        assert round(Fermi_lvl_boltz_calc,2) == round(Fermi_lvl_boltz_file,2)

        # retrieve Boltzman avg values and RDKit descriptors from the generated csv file
        pd_boltz = pd.read_csv(file_descriptors)
        energy_boltz_csv = pd_boltz["total energy"][0]
        Fermi_lvl_boltz_csv = pd_boltz["Fermi-level/eV"][0]

        assert round(energy_boltz_calc,4) == round(energy_boltz_csv,4)
        assert round(Fermi_lvl_boltz_calc,2) == round(Fermi_lvl_boltz_csv,2)
        assert pd_boltz["NumRotatableBonds"][0] == 3
        assert pd_boltz["NumRotatableBonds"][1] == 4
        assert 'DBSTEP_Vbur' not in pd_boltz

    elif file in ['test_atom.csv','test_group.csv','test_multigroup.csv','test_robert.csv']:
        pd_boltz = pd.read_csv(file_descriptors)

        # mol_1 is methane and it doesn't have any P/O atoms or CC/C=O groups so the boltz json file shouldn't appear
        assert os.path.exists(f'{folder_qdescp}/mol_1_rdkit_conf_1.json')
        assert not os.path.exists(f'{folder_boltz}/mol_1_rdkit_boltz.json')
        assert 'DBSTEP_Vbur' not in pd_boltz
        if file in ['test_atom.csv','test_group.csv']:
            assert len(pd_boltz["NumRotatableBonds"]) == 3
        if file in ['test_multigroup.csv','test_robert.csv']:
            assert len(pd_boltz["NumRotatableBonds"]) == 2

        # check variables and X_ prefixes in variable names
        if file in ['test_atom.csv','test_multigroup.csv','test_robert.csv']:
            if file == 'test_atom.csv':
                assert 'P_DBSTEP_Vbur' in pd_boltz
            else:
                assert 'P_DBSTEP_Vbur' not in pd_boltz
            assert 'P_FUKUI+' in pd_boltz
            
            if file == 'test_atom.csv':
                assert pd_boltz["NumRotatableBonds"][0] == 1

            elif file in ['test_multigroup.csv','test_robert.csv']:
                assert pd_boltz["NumRotatableBonds"][0] == 2
                assert 'CC_C1_DBSTEP_Vbur' not in pd_boltz
                assert 'CC_C2_DBSTEP_Vbur' not in pd_boltz
                assert 'CC_C1_FUKUI+' in pd_boltz
                assert 'CC_C2_FUKUI+' in pd_boltz
                # max and min values
                assert 'CC_max_DBSTEP_Vbur' not in pd_boltz
                assert 'CC_min_FUKUI+' in pd_boltz
                if file == 'test_robert.csv':
                    assert 'Name' not in pd_boltz

        elif file == 'test_group.csv':
            assert 'C=O_C_DBSTEP_Vbur' not in pd_boltz
            assert 'C=O_O_DBSTEP_Vbur' not in pd_boltz
            assert 'C=O_C_FUKUI+' in pd_boltz
            assert 'C=O_O_FUKUI+' in pd_boltz


# tests for QDESCP-NMR
@pytest.mark.parametrize(
    "json_files",
    [
        ("*.json")
    ]
)

def test_qdescp_nmr(json_files):

    # reset folder
    folder_csearch = f'{qdescp_input_dir}/CSEARCH'
    folder_qdescp = f'{qdescp_input_dir}/QDESCP'
    folder_boltz = f'{qdescp_input_dir}/boltz'
    for folder in [folder_csearch,folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    json_files = f'{qdescp_input_dir}/test_conf_*.json'
    nmr_atoms = [6, 1]  # [C,H]
    nmr_slope=[-1.0537, -1.0784]
    nmr_intercept=[181.7815,31.8723]

    # QDESCP-NMR workflow
    cmd_qdescp = [
    "python",
    "-m",
    "aqme",
    "--qdescp",
    "--program",
    "nmr",
    "--files",
    f'{json_files}',
    "--destination",
    f'{qdescp_input_dir}',
    "--nmr_slope",
    f'{nmr_slope}',
    "--nmr_intercept",
    f'{nmr_intercept}',
    "--nmr_experim",
    f'{qdescp_input_dir}/Experimental_NMR_shifts.csv'
    ]

    subprocess.run(cmd_qdescp)

    # read json files and calculate Boltzmann shifts
    energies_json,nmr_json = [],[]
    for json_file in glob.glob(json_files):
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
    folder_boltz = f'{qdescp_input_dir}/boltz'
    json_data_boltz = read_json(f'{folder_boltz}/test_boltz.json')
    nmr_json_file = list(json_data_boltz["NMR Chemical Shifts"].values())
    for i,_ in enumerate(nmr_json_calc):
        assert round(nmr_json_calc[i],2) == round(nmr_json_file[i],2)

    # check that the averaged NMR shifts are the same as the values included in the csv file
    pd_boltz = pd.read_csv(f'{qdescp_input_dir}/Experimental_NMR_shifts_predicted.csv')
    nmr_boltz_csv = pd_boltz["boltz_avg"]
    nmr_json_calc_1H = nmr_json_calc[21:]
    for i,_ in enumerate(nmr_json_calc_1H):
        assert round(nmr_json_calc_1H[i],2) == round(nmr_boltz_csv[i],2)

    # check that the shift errors in the csv file are correct
    nmr_experim_csv = pd_boltz["experimental_ppm"]
    error_calc = abs(nmr_boltz_csv - nmr_experim_csv)
    error_csv = pd_boltz["error_boltz"]
    for i,_ in enumerate(error_calc):
        if str(error_calc[i]) not in ['nan']:
            assert round(error_calc[i],2) == round(error_csv[i],2)
