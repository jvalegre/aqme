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
precision_organic = 2

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, n_confs_organic, prefilter_confs_rdkit_organic_organic, filter_confs_rdkit_organic, E_confs, charge_organic, multiplicity_organic, dihedral, xTB_ANI",
[
    # tests for conformer generation with RDKit, xTB and ANI1
    ('Organic_molecules', 'pentane.smi', 'params_test1.yaml', 240, 236, 0, [-5.27175,-4.44184,-3.84858,-1.57172], 0, 1, False, False), # test sample = 'auto', auto_sample = 20
    ('Organic_molecules', 'pentane.smi', 'params_test2.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 1, False, False), # test sample = 20
    ('Organic_molecules', 'pentane.smi', 'params_test3.yaml', 20, 5, 11, [-5.27175, -4.44184, -4.44184, -3.84858], 0, 1, False, False), # test initial_energy_threshold = 1E-10
    ('Organic_molecules', 'pentane.smi', 'params_test4.yaml', 20, 11, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -3.84858], 0, 1, False, False), # test energy_threshold = 1E-15
    ('Organic_molecules', 'pentane.smi', 'params_test5.yaml', 20, 11, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -3.84858], 0, 1, False, False), # test rms_threshold = 1E-15
    ('Organic_molecules', 'pentane.smi', 'params_test6.yaml', 20, 5, 9, [-5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858], 0, 1, False, False), # max_matches_RMSD = 1
    ('Organic_molecules', 'pentane.smi', 'params_test7.yaml', 60, 56, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 0, 1, False, False), # test sample = 'auto', auto_sample = 5
    ('Organic_molecules', 'pentane.smi', 'params_test8.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test max_torsions = 1
    ('Organic_molecules', 'pentane.smi', 'params_test9.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test max_MolWt = 1
    ('Organic_molecules', 'pentane.smi', 'params_test10.yaml', 20, 17, 0, [2.52059, 3.68961, 4.94318], 0, 1, False, False), # test ff = 'UFF'
    ('Organic_molecules', 'pentane.smi', 'params_test11.yaml', 20, 1, 12, [-5.27113, -4.44046, -4.43598, -4.06762, -3.90769, -3.81966, -2.53933], 0, 1, False, False), # test opt_steps_RDKit = 40
    ('Organic_molecules', 'pentane.smi', 'params_test12.yaml', 3, 0, 0, [-10560.62152, -10560.03876, -10559.41775], 0, 1, False, 'xTB'), # test xTB = True
    ('Organic_molecules', 'pentane.smi', 'params_test13.yaml', 3, 0, 0, [-123942.60358,-123941.72746,-123940.80475], 0, 1, False, 'ANI'), # test ANI = True
    ('Organic_molecules', 'pentane.smi', 'params_test14.yaml', 20, 17, 0, [-5.27175, -4.44184], 0, 1, False, False), # ewin_rdkit = 1
    ('Organic_molecules', 'pentane.smi', 'params_test15.yaml', 27, 'nan', 4, [-5.27175,-4.44184,-3.84858,-1.57172], 0, 1, True, False), # test dihedral scan
    ('Organic_molecules', 'pentane.smi', 'params_test16.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 3, False, False), # test multiplicity = 3
    # these two tests go together to ensure that xTB and ANI work sequentially
# MISSING        # ('Organic_molecules', 'pentane.smi', 'params_test17.yaml', 3, 0, 0, [-10560.62152, -10560.03876, -10559.41775], 0, 1, False, 'xTB'), # test xTB = True with ANI = True
# MISSING        # ('Organic_molecules', 'pentane.smi', 'params_test17.yaml', 3, 0, 0, [-123942.60358,-123941.72746,-123940.80475], 0, 1, False, 'ANI'), # test xTB = True with ANI = True
    # xTB and ANI with dihedral
    ('Organic_molecules', 'pentane.smi', 'params_test18.yaml', 3, 0, 0, [-10560.62152, -10560.03876, -10559.41775], 0, 1, True, 'xTB'), # test xTB = True with dihedral
    ('Organic_molecules', 'pentane.smi', 'params_test19.yaml', 27, 'nan', 4, [-123942.60361, -123941.72741, -123940.80485, -123939.14739], 0, 1, True, 'ANI'), # test ANI = True with dihedral
    # E filters and windows for xTB and ANI
    ('Organic_molecules', 'pentane.smi', 'params_test20.yaml', 3, 0, 0, [-10560.62152, -10560.03876], 0, 1, False, 'xTB'), # test xTB = True with E window = 1
    ('Organic_molecules', 'pentane.smi', 'params_test21.yaml', 3, 0, 0, [-123942.60358,-123941.72746], 0, 1, False, 'ANI'), # test ANI = True with E window = 1
    ('Organic_molecules', 'pentane.smi', 'params_test22.yaml', 9, 0, 3, [-10560.62134,-10560.46555,-10560.03321,-10560.02739,-10559.88780,-10559.38807], 0, 1, False, 'xTB'), # test xTB = True, opt_steps_RDKit = 40 and E filter = 0.1
    ('Organic_molecules', 'pentane.smi', 'params_test23.yaml', 7, 0, 3, [-10560.62134,-10560.03321,-10560.02739,-10559.38807], 0, 1, False, 'xTB'), # test xTB = True, opt_steps_RDKit = 40 and RMSD filter = 0.1
    ('Organic_molecules', 'pentane.smi', 'params_test24.yaml', 6, 2, 1, [-10560.62134,-10560.03321,-10559.88780], 0, 1, False, 'xTB'), # test xTB = True, opt_steps_RDKit = 40 and E pre-filter = 0.1
    ('Organic_molecules', 'pentane.smi', 'params_test25.yaml', 9, 0, 2, [-123942.60309,-123942.50243,-123941.72166,-123941.71528,-123941.60528,-123941.49161,-123940.79906], 0, 1, False, 'ANI'), # test ANI = True, opt_steps_RDKit = 40 and E filter = 0.1
    ('Organic_molecules', 'pentane.smi', 'params_test26.yaml', 7, 0, 3, [-123942.60309,-123941.72166,-123941.71528,-123940.79906], 0, 1, False, 'ANI'), # test ANI = True, opt_steps_RDKit = 40 and RMSD filter = 0.1
    ('Organic_molecules', 'pentane.smi', 'params_test27.yaml', 6, 1, 2, [-123942.60309,-123941.72166,-123941.60528], 0, 1, False, 'ANI'), # test ANI = True, opt_steps_RDKit = 40 and E pre-filter = 0.1
])

def test_confgen_organic(folder, smiles, params_file, n_confs_organic, prefilter_confs_rdkit_organic_organic, filter_confs_rdkit_organic, E_confs, charge_organic, multiplicity_organic, dihedral, xTB_ANI):
    # runs the program with the different tests
    cmd_organic = ['python', '-m', 'aqme', '--varfile', params_file]

    test_init_rdkit_confs_organic,test_prefilter_rdkit_confs_organic,test_filter_rdkit_confs_organic,round_confs_organic,test_round_confs_organic,test_charge_organic,test_unique_confs_organic,_,charge_organic_com,multiplicity_com_organic = conf_gen(path_organic, precision_organic, cmd_organic, folder, smiles, E_confs, dihedral, xTB_ANI, metal=False, template=False)

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
