#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#     Conformer generation of organic molecules      #
######################################################.

import os
import pytest
from definitions_testing import conf_gen

# saves the working directory
path_organic = os.getcwd()
# decimal digits for comparing E
precision_organic = 5

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1",
[
    # tests for conformer generation with RDKit, xTB and ANI1
    ('Organic_molecules', 'pentane.smi', 'params_test1.yaml', 240, 236, 0, [-5.27175,-4.44184,-3.84858,-1.57172], 0, False, False), # test sample = 'auto', auto_sample = 20
    ('Organic_molecules', 'pentane.smi', 'params_test2.yaml', 20, 16, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 0, False, False), # test sample = 20
    ('Organic_molecules', 'pentane.smi', 'params_test3.yaml', 20, 2, 13, [-5.27175, -4.44184, -4.44184, -3.84858, -1.57172], 0, False, False), # test initial_energy_threshold = 1E-10
    ('Organic_molecules', 'pentane.smi', 'params_test4.yaml', 20, 10, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858, -1.57172], 0, False, False), # test energy_threshold = 1E-15
    ('Organic_molecules', 'pentane.smi', 'params_test5.yaml', 20, 10, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858, -1.57172], 0, False, False), # test rms_threshold = 1E-15
    ('Organic_molecules', 'pentane.smi', 'params_test6.yaml', 20, 2, 11, [-5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858, -1.57172], 0, False, False),
    ('Organic_molecules', 'pentane.smi', 'params_test7.yaml', 60, 56, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 0, False, False), # test sample = 'auto', auto_sample = 5
    ('Organic_molecules', 'pentane.smi', 'params_test8.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test max_torsions = 1
    ('Organic_molecules', 'pentane.smi', 'params_test9.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test max_MolWt = 1
    ('Organic_molecules', 'pentane.smi', 'params_test10.yaml', 20, 16, 0, [2.52059, 3.68961, 4.94318, 6.51778], 0, False, False), # test ff = 'UFF'
    ('Organic_molecules', 'pentane.smi', 'params_test11.yaml', 20, 0, 8, [-5.26093, -4.41687, -4.39313, -4.10961, -3.93585, -2.95568, -2.43353, -2.03709, -1.51856, -1.45757, -0.22202, 0.46406], 0, False, False), # test opt_steps_RDKit = 40
    ('Organic_molecules', 'pentane.smi', 'params_test12.yaml', 20, 16, 0, [-5.27175,-4.44184,-3.84858,-1.57172], 0, False, True), # test xTB = True
    ('Organic_molecules', 'pentane.smi', 'params_test13.yaml', 20, 16, 0, [-5.27175,-4.44184,-3.84858,-1.57172], 0, False, True), # test ANI1ccx = True
    ('Organic_molecules', 'pentane.smi', 'params_test14.yaml', 20, 16, 0, [-5.27175, -4.44184], 0, False, False), # ewin = 1
    ('Organic_molecules', 'pentane.smi', 'params_test15.yaml', 36, 'nan', 4, [-5.27175,-4.44184,-3.84858,-1.57172], 0, True, False), # test dihedral scan
])

def test_confgen_organic(folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1):
    # runs the program with the different tests
    cmd_organic = ['python', '-m', 'pyconfort', '--varfile', params_file]

    conf_gen(path_organic, precision_organic, cmd_organic, folder, smiles, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1)
