#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#     Conformer generation of organic molecules      #
######################################################.

import os
import pytest
from definitions_testing import calc_energy, find_coordinates

# saves the working directory
path_xtb = os.getcwd()
# decimal digits for comparing E
precision_xtb = 5

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("params_file, coordinates, E_xtb",
[
    # tests for conformer generation with RDKit, xTB and ANI1
    ('params_xtb_test1.yaml', 'H           1.56490         0.21230         0.11960', [-5161.96238]), # test with standard options
    ('params_xtb_test2.yaml', 'H           1.57710         0.21230         0.11930', [-5623.22497]), # test with method = GFN1-xTB
    ('params_xtb_test3.yaml', 'H           1.59080         0.22530         0.12790', [-477.47884]), # test with method = GFN-FF
    ('params_xtb_test4.yaml', 'H           1.56550         0.21490         0.12170', [-5165.02187]), # test with solvent = water
    ('params_xtb_test5.yaml', 'H           1.56490         0.21230         0.11960', [-5161.96243]), # test with accuracy = 100000
    ('params_xtb_test6.yaml', 'H           1.56380         0.21220         0.11950', [-5161.96111]), # test with max_iterations = 1
])

def test_confgen_organic(params_file, coordinates, E_xtb):
    # runs the program with the different tests
    cmd_xtb = ['python', '-m', 'pyconfort', '--varfile', params_file]

    os.chdir(path+'/Organic_molecules/MeOH/')
    subprocess.call(cmd_xtb)

    # Gets energy from the SDF file
    os.chdir(path+'/Organic_molecules/MeOH/xtb_minimised_generated_sdf_files')
    test_E_xtb = calc_energy('MeOH_xtb.sdf')
    assert E_xtb == test_E_xtb

    # This returns a 1 if the coordinates are found inside the COM file
    os.chdir(path+'/Organic_molecules/MeOH/generated_gaussian_files/wb97xd-def2svp')
    coordinates_found = find_coordinates('MeOH_conformer_1.com',coordinates)
    assert 1 == coordinates_found
