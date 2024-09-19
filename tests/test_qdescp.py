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
from aqme.qdescp_utils import read_json,get_descriptors, get_descriptors
import glob
import math
import shutil

# saves the working directory
w_dir_main = os.getcwd()
qdescp_input_dir = w_dir_main + "/tests/qdescp_inputs"

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15

# import descriptors
denovo_descriptors = get_descriptors('denovo')
interpret_descriptors = get_descriptors('interpret')
full_descriptors = get_descriptors('full')


# tests for QDESCP-xTB
@pytest.mark.parametrize(
    "file",
    [
        ("test.csv"), # standard test
        ("test_atom.csv"), # test with qdescp_atoms using an atom
        ("test_idx.csv"), # test with qdescp_atoms using an atom index mapped
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
        file_descriptors_interpret = f'{w_dir_main}/AQME-ROBERT_interpret_{file_csearch}'
        file_descriptors_full = f'{w_dir_main}/AQME-ROBERT_full_{file_csearch}'
        file_descriptors_denovo = f'{w_dir_main}/AQME-ROBERT_denovo_{file_csearch}'
    else:
        file_descriptors_interpret = f'{w_dir_main}/QDESCP_interpret_boltz_descriptors.csv'
        file_descriptors_full = f'{w_dir_main}/QDESCP_full_boltz_descriptors.csv'
        file_descriptors_denovo = f'{w_dir_main}/QDESCP_denovo_boltz_descriptors.csv'
    if os.path.exists(file_descriptors_interpret):
        os.remove(file_descriptors_interpret)

    # CSEARCH conformer generation
    cmd_csearch = ["python","-m","aqme","--csearch","--program","rdkit","--input",f'{qdescp_input_dir}/{file_csearch}',"--destination",f'{folder_csearch}',]
    subprocess.run(cmd_csearch)

    # QDESCP-xTB workflow
    cmd_qdescp = ["python","-m","aqme","--qdescp","--program","xtb","--files",f'{folder_csearch}/*.sdf', "--destination",f'{folder_qdescp}',]

    if file == 'test_atom.csv':
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[P]"]
    
    if file == "test_idx.csv":
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[1]"]

    elif file == 'test_group.csv':
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[C=O]"]

    elif file in ['test_multigroup.csv','test_robert.csv']:
        # Pd is included to check the try/except in the SMARTS pattern match,
        # and to check if the code works even if there are atoms that aren't used
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[P,CC,Pd]"]
        if file == 'test_robert.csv':
            cmd_qdescp = cmd_qdescp + ["--csv_name", f'{qdescp_input_dir}/{file_csearch}']

    subprocess.run(cmd_qdescp)

    #TESTING test.csv
    # 1) check if various parameters are stored correctly from xTB calculations to the generated json files
    if file == 'test.csv':
        file_1 = "mol_1_rdkit_conf_1"
        file_2 = "mol_1_rdkit_conf_2"
        files_xtb = [file_1, file_2]
        energies_xtb, Fermi_lvls_xtb, SCC_IP_lvls_xtb = [], [], []
        energies_json, Fermi_lvls_json, SCC_IP_lvls_json = [], [], []

        for file in files_xtb:
            # Build the path for the .out file
            xtb_file_path = f'{folder_qdescp}/xtb_data/{file}/{file}.out' # mol_1_rdkit_conf_1.out
            #print(f"Reading xTB file: {xtb_file_path}")
            xtb_file_path_gn1 = f'{folder_qdescp}/xtb_data/{file}/{file}.gfn1' #mol_1_rdkit_conf_1.gfn1
            
            # Get data from the raw xTB calculations
            f = open(xtb_file_path, "r")
            data = f.readlines()
            f.close()
            f1 = open(xtb_file_path_gn1, "r")
            data1 = f1.readlines()
            f1.close()

            energy_xtb, fermi_xtb, delta_SCC_IP = None, None, None
            for i in range(0, len(data)):
                if data[i].find("TOTAL ENERGY") > -1:
                    energy_xtb = float(data[i].split()[3])
                    print(f"Energy found in {xtb_file_path}: {energy_xtb}")
                
                if data[i].find("Fermi-level") > -1:
                    fermi_xtb = float(data[i].split()[-2])  # Store the latest Fermi level found
                    print(f"Fermi-level found in {xtb_file_path}: {fermi_xtb}")

                if data1[i].find("delta SCC IP (eV):") > -1:
                    delta_SCC_IP = float(data1[i].split()[-1])
                    print(f"IP found in {xtb_file_path_gn1}: {delta_SCC_IP}")

            # Append the last found values to the lists
            if energy_xtb is not None:
                energies_xtb.append(energy_xtb)
            if fermi_xtb is not None:
                Fermi_lvls_xtb.append(fermi_xtb)
            if delta_SCC_IP is not None:
                SCC_IP_lvls_xtb.append(delta_SCC_IP)


            # Build the path for the JSON file (mol_1_rdkit_conf_1.json and mol_1_rdkit_conf_2.json)
            json_file_path = f'{folder_qdescp}/{file}.json'
            print(f"Reading JSON file: {json_file_path}")

            # Get data from the generated JSON files (mol_1_rdkit_conf_1.json and mol_1_rdkit_conf_2.json)
            json_data = read_json(json_file_path)
            energy_json = json_data["Total energy"]
            fermi_json = json_data["Fermi-level"]
            SCC_IP_json = json_data["IP"]

            energies_json.append(energy_json)
            Fermi_lvls_json.append(fermi_json)
            SCC_IP_lvls_json.append(SCC_IP_json)


        # Compare energies, Fermi level and IP values
        for i, _ in enumerate(energies_xtb):
            assert round(energies_xtb[i], 4) == round(energies_json[i], 4), \
                f"Energy mismatch at index {i}: xTB = {energies_xtb[i]}, JSON = {energies_json[i]}"
            assert round(Fermi_lvls_xtb[i], 2) == round(Fermi_lvls_json[i], 2), \
                f"Fermi level mismatch at index {i}: xTB = {Fermi_lvls_xtb[i]}, JSON = {Fermi_lvls_json[i]}"
        
        for i, _ in enumerate(SCC_IP_lvls_xtb):
            assert round(SCC_IP_lvls_xtb[i], 2) == round(SCC_IP_lvls_json[i], 2), \
                f"IP mismatch at index {i}: xTB = {SCC_IP_lvls_xtb[i]}, JSON = {SCC_IP_lvls_json[i]}"

        # 2) Calculate Boltzmann-averages
        # Calculate relative energies
        energ = [number - min(energies_xtb) for number in energies_xtb]

        # Calculate Boltzmann sum 
        boltz_sum = 0.0
        for e in energ:
            boltz_term = math.exp(-e * J_TO_AU / GAS_CONSTANT / T)
            boltz_sum += boltz_term
            print(f"Boltzmann term for energy {e}: {boltz_term}")

        # Calculate weights
        weights = []
        for e in energ:
            weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / T) / boltz_sum
            weights.append(weight)

        # Calculate Boltzmann-averaged energy
        energy_boltz_calc = 0.0
        for i, p in enumerate(energies_xtb):
            energy_boltz_calc += p * weights[i]
        print(f"Boltzmann-averaged energy XTB: {energy_boltz_calc}")

        # Calculate Boltzmann-averaged Fermi level
        Fermi_lvl_boltz_calc = 0.0
        for i, p in enumerate(Fermi_lvls_xtb):
            Fermi_lvl_boltz_calc += p * weights[i]
        print(f"Boltzmann-averaged Fermi level XTB: {Fermi_lvl_boltz_calc}")

        # Calculate Boltzmann-averaged IP
        IP_boltz_calc = 0.0
        for i, p in enumerate(SCC_IP_lvls_xtb):
            IP_boltz_calc += p * weights[i]

        # Retrieve Boltzmann-averaged values from the JSON file
        json_data = read_json(f'{folder_boltz}/mol_1_rdkit_full_boltz.json')
        json_data1 = read_json(f'{folder_boltz}/mol_1_rdkit_interpret_boltz.json')
        energy_boltz_file = json_data["Total energy"]
        Fermi_lvl_boltz_file = json_data["Fermi-level"]
        IP_lvl_boltz_file = json_data1["IP"]

        # Assertions with detailed error messages
        assert round(energy_boltz_calc, 4) == round(energy_boltz_file, 4), \
            f"Energy mismatch: calculated {energy_boltz_calc} vs file {energy_boltz_file}"

        assert round(Fermi_lvl_boltz_calc, 2) == round(Fermi_lvl_boltz_file, 2), \
            f"Fermi level mismatch: calculated {Fermi_lvl_boltz_calc} vs file {Fermi_lvl_boltz_file}"
        
        assert round(IP_boltz_calc, 2) == round(IP_lvl_boltz_file, 2), \
            f"IP mismatch: calculated {IP_boltz_calc} vs file {IP_lvl_boltz_file}"

        # 3) checking csv file
        # retrieve Boltzman avg values and RDKit descriptors from the generated csv file
        pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)
        IP_boltz_csv = pd_boltz_interpret["IP"][0]
        Fermi_lvl_boltz_csv = pd_boltz_interpret["Fermi-level"][0]

        assert round(IP_boltz_calc,4) == round(IP_boltz_csv,4), \
            f"IP mismatch: calculated {IP_boltz_calc} vs file {IP_boltz_csv} from CSV"
        assert round(Fermi_lvl_boltz_calc,2) == round(Fermi_lvl_boltz_csv,2), \
            f"Fermi level mismatch: calculated {Fermi_lvl_boltz_calc} vs file {Fermi_lvl_boltz_csv} from CSV"
        assert round(pd_boltz_interpret["MolLogP"][0],1) == 1.8
        assert round(pd_boltz_interpret["MolLogP"][1],1)== 2.2

    #Checking test_atom.csv,test_group.csv, test_multigroup.csv and test_robert.csv
    elif file in ['test_atom.csv','test_group.csv','test_multigroup.csv','test_robert.csv']:
        pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)

        # mol_1 is methane and it doesn't have any P/O atoms or CC/C=O groups so the boltz json file shouldn't appear
        assert os.path.exists(f'{folder_qdescp}/mol_1_rdkit_conf_1.json') # check if exist 
        assert not os.path.exists(f'{folder_boltz}/mol_1_rdkit_interpret_boltz.json') # check if exist
        assert len(pd_boltz_interpret["MolLogP"]) == 3

        # check variables and X_ prefixes in variable names
        if file in ['test_atom.csv','test_multigroup.csv','test_robert.csv']:
            assert 'P_Electrophil.' in pd_boltz_interpret
            
            if file == 'test_atom.csv':
                assert round(pd_boltz_interpret["MolLogP"][0],1) == 0.5

            if file == 'test_robert.csv':
                assert 'Name' not in pd_boltz_interpret

        elif file == 'test_group.csv':
            assert 'C=O_C_cm5 charges' in pd_boltz_interpret
            assert 'C=O_O_cm5 charges' in pd_boltz_interpret


    elif file == 'test_idx.csv':
        pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)
        assert 'C1_cm5 charges' in pd_boltz_interpret
        assert round(pd_boltz_interpret['C1_cm5 charges'][1],1) == -0.2

#Checking molecular and atomic descriptors
    def check_descriptors(pd_boltz, descriptors, excluded_descriptors, desc_type):
        """
        Function to check the presence and absence of descriptors in the DataFrame.
        pd_boltz: Pandas DataFrame with the calculated descriptors.
        descriptors: List of descriptors that should be present.
        excluded_descriptors: List of descriptors that should not be present.
        desc_type: Type of descriptors ('mol' or 'atoms') for printing in messages.
        """
        # Check for the presence of descriptors
        for descp in descriptors:
            assert descp in pd_boltz.columns, f"{desc_type.capitalize()} descriptor {descp} is missing from columns!"

        # Check for the absence of descriptors that should not be present
        for descp in excluded_descriptors:
            assert descp not in pd_boltz.columns, f"{desc_type.capitalize()} descriptor {descp} should not be present in columns!"

    # Checking molecular descriptors
    descp_denovo_mol = denovo_descriptors['mol'] 
    descp_denovo_atoms = denovo_descriptors['atoms']

    descp_interpret_mol = descp_denovo_mol + interpret_descriptors['mol']
    descp_interpret_atoms = descp_denovo_atoms + interpret_descriptors['atoms']

    descp_full_mol = descp_interpret_mol + full_descriptors['mol']
    descp_full_atoms = descp_interpret_atoms + full_descriptors['atoms']

    # Read the CSV files
    file_descriptors_denovo = f'{w_dir_main}/QDESCP_denovo_boltz_descriptors.csv'
    file_descriptors_interpret = f'{w_dir_main}/QDESCP_interpret_boltz_descriptors.csv'
    file_descriptors_full = f'{w_dir_main}/QDESCP_full_boltz_descriptors.csv'
    pd_boltz_denovo = pd.read_csv(file_descriptors_denovo)
    pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)
    pd_boltz_full = pd.read_csv(file_descriptors_full)

    # Check molecular and atomic descriptors
    if file == 'test.csv':
        check_descriptors(pd_boltz_denovo, descp_denovo_mol, [d for d in descp_full_mol if d not in descp_denovo_mol], 'mol')
        check_descriptors(pd_boltz_interpret, descp_interpret_mol, [d for d in descp_full_mol if d not in descp_interpret_mol], 'mol')
        check_descriptors(pd_boltz_full, descp_full_mol, [], 'mol')

        check_descriptors(pd_boltz_denovo, descp_denovo_atoms, [d for d in descp_full_atoms if d not in descp_denovo_atoms], 'atoms')
        check_descriptors(pd_boltz_interpret, descp_interpret_atoms, [d for d in descp_full_atoms if d not in descp_interpret_atoms], 'atoms')
        check_descriptors(pd_boltz_full, descp_full_atoms, [], 'atoms')
    



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
    pd_boltz_interpret = pd.read_csv(f'{qdescp_input_dir}/Experimental_NMR_shifts_predicted.csv')
    nmr_boltz_csv = pd_boltz_interpret["boltz_avg"]
    nmr_json_calc_1H = nmr_json_calc[21:]
    for i,_ in enumerate(nmr_json_calc_1H):
        assert round(nmr_json_calc_1H[i],2) == round(nmr_boltz_csv[i],2)

    # check that the shift errors in the csv file are correct
    nmr_experim_csv = pd_boltz_interpret["experimental_ppm"]
    error_calc = abs(nmr_boltz_csv - nmr_experim_csv)
    error_csv = pd_boltz_interpret["error_boltz"]
    for i,_ in enumerate(error_calc):
        if str(error_calc[i]) not in ['nan']:
            assert round(error_calc[i],2) == round(error_csv[i],2)
