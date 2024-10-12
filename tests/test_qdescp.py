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
from aqme.qdescp_utils import read_json,get_descriptors, get_descriptors, load_file_formats
import glob
import math
import shutil

# saves the working directory
w_dir_main = os.getcwd()
qdescp_input_dir = w_dir_main + "/tests/qdescp_inputs"
qdescp_empty_dir = w_dir_main + "/tests/qdescp_empty"
qdescp_sdf_dir = w_dir_main + "/tests/qdescp_sdf"

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15

# import descriptors and file formats
denovo_descriptors = get_descriptors('denovo')
interpret_descriptors = get_descriptors('interpret')
full_descriptors = get_descriptors('full')
file_formats = load_file_formats()


# tests for QDESCP-xTB
@pytest.mark.parametrize(
    "file",
    [
        ("test.csv"), # standard test
        ("test_atom.csv"), # test with qdescp_atoms using an atom
        ("test_idx.csv"), # test with qdescp_atoms using an atom index mapped
        ("test_group.csv"), # test with qdescp_atoms using a functional group
        ("test_multigroup.csv"), # test with qdescp_atoms using a multiple atoms and functional groups
        ("test_robert_atom.csv"), # test for the AQME-ROBERT workflow with atomic descriptors
        ("test_robert_mol.csv") # test for the AQME-ROBERT workflow with NO atomic descriptors
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

    if file in ['test_multigroup.csv','test_robert_atom.csv','test_robert_mol.csv']:
        file_csearch = 'test_atom.csv'
    else:
        file_csearch = file
    if file in ['test_robert_atom.csv','test_robert_mol.csv']:
        file_descriptors_interpret = f'{w_dir_main}/AQME-ROBERT_interpret_{file_csearch}'
        file_descriptors_full = f'{w_dir_main}/AQME-ROBERT_full_{file_csearch}'
        file_descriptors_denovo = f'{w_dir_main}/AQME-ROBERT_denovo_{file_csearch}'
        file2_descriptors_interpret = f'{w_dir_main}/QDESCP_interpret_descriptors.csv'
        file2_descriptors_full = f'{w_dir_main}/QDESCP_full_descriptors.csv'
        file2_descriptors_denovo = f'{w_dir_main}/QDESCP_denovo_descriptors.csv'

    else:
        file_descriptors_interpret = f'{w_dir_main}/QDESCP_interpret_descriptors.csv'
        file_descriptors_full = f'{w_dir_main}/QDESCP_full_descriptors.csv'
        file_descriptors_denovo = f'{w_dir_main}/QDESCP_denovo_descriptors.csv'
    if os.path.exists(file_descriptors_denovo): 
        os.remove(file_descriptors_denovo)
    if os.path.exists(file_descriptors_interpret): 
        os.remove(file_descriptors_interpret)
    if os.path.exists(file_descriptors_full): 
        os.remove(file_descriptors_full)

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

    elif file in ['test_multigroup.csv','test_robert_atom.csv']:
        # Pd is included to check the try/except in the SMARTS pattern match,
        # and to check if the code works even if there are atoms that aren't used
        cmd_qdescp = cmd_qdescp + ["--qdescp_atoms", "[P,CC,Pd]"]
        if file == 'test_robert_atom.csv':
            cmd_qdescp = cmd_qdescp + ["--csv_name", f'{qdescp_input_dir}/{file_csearch}']
    
    elif file == 'test_robert_mol.csv':
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

        for file_xtb in files_xtb:
            # only two files should remain inside the xtb_data folders
            xtb_data_path = f'{folder_qdescp}/xtb_data/{file_xtb}'
            assert len(glob.glob(f'{xtb_data_path}/*')) == 2
            assert os.path.exists(f'{xtb_data_path}/{file_xtb}.xyz')
            assert os.path.exists(f'{xtb_data_path}/{file_xtb}_All_Calcs.out')

            # Build the path for the out file
            xtb_file_path = f'{xtb_data_path}/{file_xtb}_All_Calcs.out' # mol_1_rdkit_conf_1.out
             
            # Get data from the raw xTB calculations
            f = open(xtb_file_path, "r")
            data = f.readlines()
            f.close()

            energy_xtb, fermi_xtb, delta_SCC_IP = None, None, None
            E_section = False
            IP_section = False
            for i in range(0, len(data)):
                if '----- Single-point -----' in data[i]:
                    E_section = True
                elif '----- FOD -----' in data[i]:
                    E_section = False

                if '----- GFN1 -----' in data[i]:
                    IP_section = True
                elif '----- Nminus1 -----' in data[i]:
                    IP_section = False

                if data[i].find("TOTAL ENERGY") > -1 and E_section:
                    energy_xtb = float(data[i].split()[3])
                
                if data[i].find("Fermi-level") > -1 and E_section:
                    fermi_xtb = float(data[i].split()[-2])  # Store the latest Fermi level found

                if data[i].find("delta SCC IP (eV):") > -1 and IP_section:
                    delta_SCC_IP = float(data[i].split()[-1])
            
            # check that all the calcs are inside the All_calcs.out file
            start_line,find_format = 0,0
            for file_format in file_formats:
                for i in range(start_line, len(data)):
                    if f'----- {file_formats[file_format]} -----' in data[i]:
                        find_format += 1
                        start_line = i
                        break
            assert find_format == len(file_formats)


            # Append the last found values to the lists
            if energy_xtb is not None:
                energies_xtb.append(energy_xtb)
            if fermi_xtb is not None:
                Fermi_lvls_xtb.append(fermi_xtb)
            if delta_SCC_IP is not None:
                SCC_IP_lvls_xtb.append(delta_SCC_IP)


            # Build the path for the JSON file (mol_1_rdkit_conf_1.json and mol_1_rdkit_conf_2.json)
            json_file_path = f'{folder_qdescp}/{file_xtb}.json'

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
            assert round(energies_xtb[i], 2) == round(energies_json[i], 2), \
                f"Energy mismatch at index {i}: xTB = {energies_xtb[i]}, JSON = {energies_json[i]}"
            assert round(Fermi_lvls_xtb[i], 1) == round(Fermi_lvls_json[i], 1), \
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

        # Calculate weights
        weights = []
        for e in energ:
            weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / T) / boltz_sum
            weights.append(weight)

        # Calculate Boltzmann-averaged energy
        energy_boltz_calc = 0.0
        for i, p in enumerate(energies_xtb):
            energy_boltz_calc += p * weights[i]

        # Calculate Boltzmann-averaged Fermi level
        Fermi_lvl_boltz_calc = 0.0
        for i, p in enumerate(Fermi_lvls_xtb):
            Fermi_lvl_boltz_calc += p * weights[i]

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

        assert round(Fermi_lvl_boltz_calc, 1) == round(Fermi_lvl_boltz_file, 1), \
            f"Fermi level mismatch: calculated {Fermi_lvl_boltz_calc} vs file {Fermi_lvl_boltz_file}"
        
        assert round(IP_boltz_calc, 1) == round(IP_lvl_boltz_file, 1), \
            f"IP mismatch: calculated {IP_boltz_calc} vs file {IP_lvl_boltz_file}"

        # 3) checking csv file
        # retrieve Boltzman avg values and RDKit descriptors from the generated csv file
        pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)
        IP_boltz_csv = pd_boltz_interpret["IP"][0]
        Fermi_lvl_boltz_csv = pd_boltz_interpret["Fermi-level"][0]

        assert round(IP_boltz_calc,4) == round(IP_boltz_csv,4), \
            f"IP mismatch: calculated {IP_boltz_calc} vs file {IP_boltz_csv} from CSV"
        assert round(Fermi_lvl_boltz_calc,1) == round(Fermi_lvl_boltz_csv,1), \
            f"Fermi level mismatch: calculated {Fermi_lvl_boltz_calc} vs file {Fermi_lvl_boltz_csv} from CSV"
        assert round(pd_boltz_interpret["MolLogP"][0],1) == 1.8
        assert round(pd_boltz_interpret["MolLogP"][1],1)== 2.2

        pd_boltz_full = pd.read_csv(file_descriptors_full)
        assert pd_boltz_full["NumRotatableBonds"][0] == 3
        assert pd_boltz_full["NumRotatableBonds"][1] == 4

    #Checking test_atom.csv,test_group.csv, test_multigroup.csv and test_robert_atom.csv
    elif file in ['test_atom.csv','test_group.csv','test_multigroup.csv','test_robert_atom.csv']:
        pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)

        # mol_1 is methane and it doesn't have any P/O atoms or CC/C=O groups so atomic descriptors shouldn't appear
        assert os.path.exists(f'{folder_qdescp}/mol_1_rdkit_conf_1.json') # check if exist 
        assert os.path.exists(f'{folder_boltz}/mol_1_rdkit_interpret_boltz.json') # check if exist
        assert len(pd_boltz_interpret["MolLogP"]) == 4
        
        count_nan = 0 # dirty hack that account for different sortings of the calcs within the CSV files
        for val in pd_boltz_interpret["MolLogP"]:
            if str(val).lower() == 'nan':
                count_nan += 1
        assert count_nan == 0

        if file == 'test_group.csv':
            assert len(pd_boltz_interpret["C=O_C_FOD"]) == 4
            count_nan = 0 # dirty hack that account for different sortings of the calcs within the CSV files
            for val in pd_boltz_interpret["C=O_C_FOD"]:
                if str(val).lower() == 'nan':
                    count_nan += 1
            assert count_nan == 1
        else:
            # from xTB
            assert len(pd_boltz_interpret["P_FOD"]) == 4
            count_nan = 0 # dirty hack that account for different sortings of the calcs within the CSV files
            for val in pd_boltz_interpret["P_FOD"]:
                if str(val).lower() == 'nan':
                    count_nan += 1
            assert count_nan == 1

            # from MORFEUS
            assert len(pd_boltz_interpret["P_Buried volume"]) == 4
            count_nan = 0 # dirty hack that account for different sortings of the calcs within the CSV files
            for val in pd_boltz_interpret["P_Buried volume"]:
                if str(val).lower() == 'nan':
                    count_nan += 1
            assert count_nan == 1

        # check variables and X_ prefixes in variable names
        if file in ['test_atom.csv','test_multigroup.csv','test_robert_atom.csv']:
            assert 'P_Electrophil.' in pd_boltz_interpret
            
            if file == 'test_atom.csv':
                assert round(pd_boltz_interpret["MolLogP"][0],1) == 0.6

            if file == 'test_robert_atom.csv':
                assert 'Name' not in pd_boltz_interpret

        elif file == 'test_group.csv':
            assert 'C=O_C_cm5 charges' in pd_boltz_interpret
            assert 'C=O_O_cm5 charges' in pd_boltz_interpret


    elif file == 'test_idx.csv':
        pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)
        assert 'C1_cm5 charges' in pd_boltz_interpret
        assert round(pd_boltz_interpret['C1_cm5 charges'][1],1) == -0.2

    # Checking molecular and atomic descriptors
    def check_descriptors(pd_boltz, descriptors, excluded_descriptors, desc_type, file_test):
        """
        Function to check the presence and absence of descriptors in the DataFrame.
        pd_boltz: Pandas DataFrame with the calculated descriptors.
        descriptors: List of descriptors that should be present.
        excluded_descriptors: List of descriptors that should not be present.
        desc_type: Type of descriptors ('mol' or 'atoms') for printing in messages.
        """
        # Check for the presence of descriptors
        for descp in descriptors:
            for val in pd_boltz[descp]:
                if file_test in ['test.csv','test_robert_mol.csv']:
                    assert str(val).lower() != 'nan'
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

    # add P_ prefix for SMARTS pattern
    descp_denovo_atoms = [f'P_{descp}' for descp in descp_denovo_atoms]
    descp_interpret_atoms = [f'P_{descp}' for descp in descp_interpret_atoms]
    descp_full_atoms = [f'P_{descp}' for descp in descp_full_atoms]

    # Read the CSV files   
    pd_boltz_denovo = pd.read_csv(file_descriptors_denovo)
    pd_boltz_interpret = pd.read_csv(file_descriptors_interpret)
    pd_boltz_full = pd.read_csv(file_descriptors_full)

    # Check molecular and atomic descriptors in the QDESCP_ files
    if file == 'test.csv':
        assert 'code_name' in pd_boltz_interpret.columns
        check_descriptors(pd_boltz_denovo, descp_denovo_mol, [d for d in descp_full_mol if d not in descp_denovo_mol], 'mol', file)
        check_descriptors(pd_boltz_interpret, descp_interpret_mol, [d for d in descp_full_mol if d not in descp_interpret_mol], 'mol', file)
        check_descriptors(pd_boltz_full, descp_full_mol, [], 'mol', file)
        assert len(pd_boltz_denovo.columns) == 12 == len(descp_denovo_mol)+2 # descps + 2 extra: code_name and MolLogP
        assert len(pd_boltz_interpret.columns) == 26 == len(descp_interpret_mol)+2
        assert len(pd_boltz_full.columns) == 244 # this might change in future RDKit versions

        # check whether the QDESCP original and raw files were moved to the raw_csv_databases folder
        for file_csv in [file_descriptors_denovo, file_descriptors_interpret,file_descriptors_full]:
            raw_csv = f'{folder_qdescp}/raw_csv_databases/Raw_{os.path.basename(file_csv)}'
            assert os.path.exists(raw_csv)
            raw_df = pd.read_csv(raw_csv)
            assert 'code_name' in raw_df.columns
            for raw_atom_val in ['cm5 charges','Electrophil.','Nucleophil.','Radical attack','SASA','Buried volume','Cone angle','H bond H2O']:
                assert raw_df[raw_atom_val][0][0] == '['
            if file_csv in [file_descriptors_interpret,file_descriptors_full]:
                for raw_atom_val in ["s proportion", "p proportion", "d proportion", "Coord. numbers","Polariz. alpha", "FOD", "FOD s proportion", "FOD p proportion", "FOD d proportion","Solid angle", "Pyramidaliz. P", "Pyramidaliz. Vol", "Dispersion"]:
                    assert raw_df[raw_atom_val][0][0] == '['

    # Check molecular and atomic descriptors in the AQME-ROBERT files
    elif file in ['test_robert_mol.csv','test_robert_atom.csv']:
        assert 'SMILES' in pd_boltz_interpret.columns
        # molecular descriptors must be in both cases
        check_descriptors(pd_boltz_denovo, descp_denovo_mol, [d for d in descp_full_mol if d not in descp_denovo_mol], 'mol', file)
        check_descriptors(pd_boltz_interpret, descp_interpret_mol, [d for d in descp_full_mol if d not in descp_interpret_mol], 'mol', file)
        check_descriptors(pd_boltz_full, descp_full_mol, [], 'mol', file)

        # atomic descriptors must be here
        if file == 'test_robert_atom.csv':
            check_descriptors(pd_boltz_denovo, descp_denovo_atoms, [d for d in descp_full_atoms if d not in descp_denovo_atoms], 'atoms', file)
            check_descriptors(pd_boltz_interpret, descp_interpret_atoms, [d for d in descp_full_atoms if d not in descp_interpret_atoms], 'atoms', file)
            check_descriptors(pd_boltz_full, descp_full_atoms, [], 'atoms', file)
            assert pd_boltz_denovo.columns[:2].tolist() == ['SMILES','code_name']
            assert len(pd_boltz_denovo.columns) == 21 == len(descp_denovo_mol)+len(descp_denovo_atoms)+3 # 3 extra: SMILES, code_name and MolLogP
            assert pd_boltz_interpret.columns[:2].tolist() == ['SMILES','code_name']
            assert len(pd_boltz_interpret.columns) == 48 == len(descp_interpret_mol)+len(descp_interpret_atoms)+3
            assert pd_boltz_full.columns[:2].tolist() == ['SMILES','code_name']
            assert len(pd_boltz_full.columns) == 282 # bunch of RDKit descps

        # atomic descriptors must not be here
        elif file == 'test_robert_mol.csv':
            assert pd_boltz_denovo.columns[:2].tolist() == ['SMILES','code_name']
            assert len(pd_boltz_denovo.columns) == 13 == len(descp_denovo_mol)+3 # 3 extra: SMILES, code_name and MolLogP
            assert pd_boltz_interpret.columns[:2].tolist() == ['SMILES','code_name']
            assert len(pd_boltz_interpret.columns) == 27 == len(descp_interpret_mol)+3
            assert pd_boltz_full.columns[:2].tolist() == ['SMILES','code_name']
            assert len(pd_boltz_full.columns) == 245 # this might change in future RDKit versions

        # check whether the QDESCP original and raw files were moved to the raw_csv_databases folder
        for file_csv in [file_descriptors_denovo, file_descriptors_interpret,file_descriptors_full]:
            assert os.path.exists(file_csv)
        for file_csv2 in [file2_descriptors_denovo, file2_descriptors_interpret, file2_descriptors_full]:
            assert os.path.exists(f'{folder_qdescp}/raw_csv_databases/{os.path.basename(file_csv2)}')
            raw_csv = f'{folder_qdescp}/raw_csv_databases/Raw_{os.path.basename(file_csv2)}'
            assert os.path.exists(raw_csv)
            raw_df = pd.read_csv(raw_csv)
            assert 'code_name' in raw_df.columns
            if file == 'test_robert_mol.csv':
                for raw_atom_val in ['cm5 charges','Electrophil.','Nucleophil.','Radical attack','SASA','Buried volume','Cone angle','H bond H2O']:
                    assert raw_df[raw_atom_val][0][0] == '['
                if file_csv2 in [file2_descriptors_interpret,file2_descriptors_full]:
                    for raw_atom_val in ["s proportion", "p proportion", "d proportion", "Coord. numbers","Polariz. alpha", "FOD", "FOD s proportion", "FOD p proportion", "FOD d proportion","Solid angle", "Pyramidaliz. P", "Pyramidaliz. Vol", "Dispersion"]:
                        assert raw_df[raw_atom_val][0][0] == '['

@pytest.mark.parametrize(
    "test",
    [
        # tests for ignoring not working optimizations (i.e. qdescp keeps going even if an xTB calc fails)
        ("empty_values"),
    ],
)
def test_qdescp_missing(
    test
):

    # reset folder and files
    folder_qdescp = f'{qdescp_empty_dir}/QDESCP'
    folder_boltz = f'{folder_qdescp}/boltz'
    for folder in [folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    file_descriptors_interpret = f'{w_dir_main}/QDESCP_interpret_descriptors.csv'
    file_descriptors_full = f'{w_dir_main}/QDESCP_full_descriptors.csv'
    file_descriptors_denovo = f'{w_dir_main}/QDESCP_denovo_descriptors.csv'
    if os.path.exists(file_descriptors_denovo): 
        os.remove(file_descriptors_denovo)
    if os.path.exists(file_descriptors_interpret): 
        os.remove(file_descriptors_interpret)
    if os.path.exists(file_descriptors_full): 
        os.remove(file_descriptors_full)

    # QDESCP-xTB workflow
    cmd_qdescp = ["python","-m","aqme","--qdescp","--program","xtb","--files",f'{qdescp_empty_dir}/*.sdf', "--destination",f'{folder_qdescp}',]

    subprocess.run(cmd_qdescp)

    # checking csv file
    df_interpret = pd.read_csv(file_descriptors_interpret)
    assert len(df_interpret['code_name']) == 2
    assert 'a' in df_interpret['code_name'][0]
    assert 'c' in df_interpret['code_name'][1]


@pytest.mark.parametrize(
    "test",
    [
        # tests for using SDF as inputs and automated detection of common atom patterns
        ("sdf_input_n_auto"),
    ],
)
def test_qdescp_sdf(
    test
):

    # reset folder and files
    folder_qdescp = f'{qdescp_sdf_dir}/QDESCP'
    folder_boltz = f'{folder_qdescp}/boltz'
    for folder in [folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    file_descriptors_interpret = f'{w_dir_main}/QDESCP_interpret_descriptors.csv'
    file_descriptors_full = f'{w_dir_main}/QDESCP_full_descriptors.csv'
    file_descriptors_denovo = f'{w_dir_main}/QDESCP_denovo_descriptors.csv'
    if os.path.exists(file_descriptors_denovo): 
        os.remove(file_descriptors_denovo)
    if os.path.exists(file_descriptors_interpret): 
        os.remove(file_descriptors_interpret)
    if os.path.exists(file_descriptors_full): 
        os.remove(file_descriptors_full)

    # QDESCP-xTB workflow
    cmd_qdescp = ["python","-m","aqme","--qdescp","--program","xtb","--files",f'{qdescp_sdf_dir}/*.sdf', "--destination",f'{folder_qdescp}',]

    subprocess.run(cmd_qdescp)

    # checking csv file
    df_interpret = pd.read_csv(file_descriptors_interpret)
    assert len(df_interpret['code_name']) == 2
    assert 'mol1' == df_interpret['code_name'][0]
    assert 'mol_2' == df_interpret['code_name'][1]
    assert len(df_interpret.columns) == 47

    # check if the automated detection of common pattern works
    count_nan = 0 # dirty hack that account for different sortings of the calcs within the CSV files
    for val in df_interpret["P_FOD"]:
        if str(val).lower() == 'nan':
            count_nan += 1
    assert count_nan == 0


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
