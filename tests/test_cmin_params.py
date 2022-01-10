#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#     Conformer generation of organic molecules      #
#           designed to test CMIN options            #
######################################################.

import os
import glob
import subprocess
import pytest
from definitions_testing import calc_energy, find_coordinates, remove_data

# saves the working directory
path_cmin = os.getcwd()
# decimal digits for comparing E
precision_cmin = 2

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("params_file, coordinates, E_cmin, cmin_type",
[
    # reference test for conformer generation with RDKit
    ('params_rdkit_test1.yaml', 'H           1.54450         0.23880         0.13760', [0.66577], 'rdkit'), # test with standard options
    # tests for conformer generation with xTB
    ('params_xtb_test1.yaml', 'H           1.56490         0.21230         0.11960', [-5161.96238], 'xtb'), # test with standard options
    ('params_xtb_test2.yaml', 'H           1.57710         0.21230         0.11930', [-5623.22497], 'xtb'), # test with method = GFN1-xTB
    ('params_xtb_test3.yaml', 'H           1.59080         0.22530         0.12790', [-477.47884], 'xtb'), # test with method = GFN-FF
    ('params_xtb_test4.yaml', 'H           1.56550         0.21490         0.12170', [-5165.02187], 'xtb'), # test with solvent = water
    ('params_xtb_test5.yaml', 'H           1.56490         0.21230         0.11960', [-5161.96243], 'xtb'), # test with accuracy = 100000
    ('params_xtb_test6.yaml', 'H           1.56380         0.21220         0.11950', [-5161.96111], 'xtb'), # test with max_iterations = 1
    # tests for conformer generation with ANI
    ('params_ani_test1.yaml', 'H           1.57430         0.20860         0.11690', [-72590.69313], 'ani'), # test with ani_method = ANI1x
    ('params_ani_test2.yaml', 'H           1.57220         0.20940         0.11750', [-72554.92451], 'ani'), # test with ani_method = ANI1ccx
    ('params_ani_test3.yaml', 'H           1.56290         0.21390         0.12070', [-72590.57106], 'ani'), # test with ani_method = ANI2x
])

def test_confgen_organic(params_file, coordinates, E_cmin, cmin_type):
    # runs the program with the different tests
    cmd_cmin = ['python', '-m', 'aqme', '--varfile', params_file]

    os.chdir(path_cmin+'/Organic_molecules/MeOH/')
    subprocess.call(cmd_cmin)

    # Gets energy from the SDF file
    if cmin_type == 'rdkit':
        os.chdir(path_cmin+'/Organic_molecules/MeOH/CSEARCH/rdkit')
        test_E_cmin = calc_energy('MeOH_rdkit.sdf')
    if cmin_type == 'xtb':
        os.chdir(path_cmin+'/Organic_molecules/MeOH/CSEARCH/xtb')
        test_E_cmin = calc_energy('MeOH_xtb.sdf')
    elif cmin_type == 'ani':
        os.chdir(path_cmin+'/Organic_molecules/MeOH/CSEARCH/ani')
        test_E_cmin = calc_energy('MeOH_ani.sdf')

    test_E_cmin = [round(num, precision_cmin) for num in test_E_cmin]
    assert E_cmin == test_E_cmin

    # This returns a 1 if the coordinates are found inside the COM file
    os.chdir(path_cmin+'/Organic_molecules/MeOH/QMCALC/G16/wb97xd-6-31g(d)')
    file_cmin = glob.glob('MeOH_*.com')[0]
    coordinates_found = find_coordinates(file_cmin,coordinates)
    assert 1 == coordinates_found

    remove_data(path_cmin,"Organic_molecules","MeOH")
