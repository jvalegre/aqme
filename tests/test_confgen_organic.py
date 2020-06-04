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
@pytest.mark.parametrize("folder, smiles, params_file, n_confs_organic, prefilter_confs_rdkit_organic_organic, filter_confs_rdkit_organic, E_confs, charge_organic, multiplicity_organic, dihedral, xTB_ANI1",
[
    # tests for conformer generation with RDKit, xTB and ANI1
    ('Organic_molecules', 'pentane.smi', 'params_test1.yaml', 240, 236, 0, [-5.27175,-4.44184,-3.84858,-1.57172], 0, 1, False, False), # test sample = 'auto', auto_sample = 20
    ('Organic_molecules', 'pentane.smi', 'params_test2.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 1, False, False), # test sample = 20
    ('Organic_molecules', 'pentane.smi', 'params_test3.yaml', 20, 5, 11, [-5.27175, -4.44184, -4.44184, -3.84858], 0, 1, False, False), # test initial_energy_threshold = 1E-10
    ('Organic_molecules', 'pentane.smi', 'params_test4.yaml', 20, 11, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -3.84858], 0, 1, False, False), # test energy_threshold = 1E-15
    ('Organic_molecules', 'pentane.smi', 'params_test5.yaml', 20, 11, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -3.84858], 0, 1, False, False), # test rms_threshold = 1E-15
    ('Organic_molecules', 'pentane.smi', 'params_test6.yaml', 20, 5, 9, [-5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858], 0, 1, False, False),
    ('Organic_molecules', 'pentane.smi', 'params_test7.yaml', 60, 56, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 0, 1, False, False), # test sample = 'auto', auto_sample = 5
    ('Organic_molecules', 'pentane.smi', 'params_test8.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test max_torsions = 1
    ('Organic_molecules', 'pentane.smi', 'params_test9.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test max_MolWt = 1
    ('Organic_molecules', 'pentane.smi', 'params_test10.yaml', 20, 17, 0, [2.52059, 3.68961, 4.94318], 0, 1, False, False), # test ff = 'UFF'
    ('Organic_molecules', 'pentane.smi', 'params_test11.yaml', 20, 1, 12, [-5.27113, -4.44046, -4.43598, -4.06762, -3.90769, -3.81966, -2.53933], 0, 1, False, False), # test opt_steps_RDKit = 40
    ('Organic_molecules', 'pentane.smi', 'params_test12.yaml', 3, 0, 0, [-10560.62152, -10560.03876, -10559.41773], 0, 1, False, 'xTB'), # test xTB = True
    ('Organic_molecules', 'pentane.smi', 'params_test13.yaml', 3, 0, 0, [-5.27175,-4.44184,-3.84858,-1.57172], 0, 1, False, 'AN1ccx'), # test ANI1ccx = True
    ('Organic_molecules', 'pentane.smi', 'params_test14.yaml', 20, 17, 0, [-5.27175, -4.44184], 0, 1, False, False), # ewin = 1
    ('Organic_molecules', 'pentane.smi', 'params_test15.yaml', 27, 'nan', 4, [-5.27175,-4.44184,-3.84858,-1.57172], 0, 1, True, False), # test dihedral scan
    ('Organic_molecules', 'pentane.smi', 'params_test16.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 3, False, False), # test multiplicity = 3
])

def test_confgen_organic(folder, smiles, params_file, n_confs_organic, prefilter_confs_rdkit_organic_organic, filter_confs_rdkit_organic, E_confs, charge_organic, multiplicity_organic, dihedral, xTB_ANI1):
    # runs the program with the different tests
    cmd_organic = ['python', '-m', 'pyconfort', '--varfile', params_file]

    test_init_rdkit_confs_organic,test_prefilter_rdkit_confs_organic,test_filter_rdkit_confs_organic,round_confs_organic,test_round_confs_organic,test_charge_organic,test_unique_confs_organic,_,charge_organic_com,multiplicity_com_organic = conf_gen(path_organic, precision_organic, cmd_organic, folder, smiles, E_confs, dihedral, xTB_ANI1, metal=False, template=False)

    # the assert statements are placed here, otherwise pytest doesn't explain the AssertionError
    # first, dicard tests 8 and 9 since they are designed to fail
    if n_confs_organic != 'nan':
        # dihedral vs no dihedral scans
        if not dihedral:
            assert str(n_confs_organic) == str(test_init_rdkit_confs_organic[0])
            assert str(prefilter_confs_rdkit_organic_organic) == str(test_prefilter_rdkit_confs_organic[0])
            assert str(filter_confs_rdkit_organic) == str(test_filter_rdkit_confs_organic[0])
        else:
            assert str(n_confs_organic) == str(test_init_rdkit_confs_organic[0])
            # I use the filter_confs_rdkit_organic variable to assert for unique confs in dihedral scan
            assert str(filter_confs_rdkit_organic) == str(test_unique_confs_organic[0])

        assert str(round_confs_organic) == str(test_round_confs_organic)
        assert str(charge_organic) == str(test_charge_organic[0])

        # make sure the COM files have the right charge_organic and multiplicity
        assert str(charge_organic_com) == str(charge_organic)
        assert str(multiplicity_com_organic) == str(multiplicity_organic)

    elif params_file == 'params_test8.yaml' or params_file == 'params_test9.yaml':
        assert str(test_filter_rdkit_confs_organic) == 'nan'
        assert str(test_round_confs_organic) == 'nan'

    else:
        assert 3 ==2
