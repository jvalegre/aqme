#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#              Miscellaneous tests 2                 #
######################################################.

import os
import subprocess
import pytest
from definitions_testing import get_not_empty_files,calc_genecp,remove_data

# saves the working directory
path_misc3 = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, params_file, com_files",
[
    # tests with the options of Fullmonte
    ('QPREP_misc3', 'params_test_qprep_misc3.yaml', 2), # checks Fullmonte energies to original RDKit energies
])

def test_confgen_misc3(folder, params_file, com_files):
    # runs the program with the different tests
    cmd_misc3 = ['python', '-m', 'aqme', '--varfile', params_file]

    os.chdir(path_misc3+'/'+folder)
    subprocess.call(cmd_misc3)

    # check number of com files
    com_files_folder = path_misc3+'/'+folder+'/QMCALC/G16/wb97xd-6-31g(d)'
    test_com_files = 0
    test_com_files = get_not_empty_files(com_files_folder,test_com_files,'com')

    assert str(com_files) == str(test_com_files)

    # check that QPREP reads the charges from the SDF files
    _,_,_,_,charge_com1,multiplicity_com1 = calc_genecp('Ir_1_8.com', [])
    _,_,_,_,charge_com2,multiplicity_com2 = calc_genecp('Ir_2_5.com', [])

    assert str(charge_com1) == '0'
    assert str(multiplicity_com1) == '1'

    assert str(charge_com2) == '1'
    assert str(multiplicity_com2) == '1'

    remove_data(path_misc3, folder, smiles=False)
