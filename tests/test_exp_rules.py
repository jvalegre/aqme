#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#   Experimental rules applied to metal complexes    #
######################################################.

import os
import pytest
from definitions_testing import conf_gen_exp_rules

# saves the working directory
path_exp_rules = os.getcwd()
# decimal digits for comparing E
precision_exp_rules = 5
params_file = 'params_Ir_exp_rules_test.yaml'

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, E_confs_no_rules, E_confs_rules, com_files, charge",
[
    # tests for conformer generation with experimental rules
    ('Ir_exp_rules', 'Ir_1', [68.86185, 69.01720], [68.86185], 1, 0), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_2', [69.54124, 69.60127, 69.63704], [69.60127], 1, 1), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_3', [48.28081, 48.3819, 48.64117, 48.92291], [48.38190], 1, 1), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_4', [57.24954, 57.4422, 57.48444, 57.84497], [57.24954], 1, 0), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_5', [122.76353, 123.88455, 124.00555, 125.27744, 125.43099], [125.43099], 1, 0), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_6', [129.91109, 130.12558, 130.49653], [129.91109], 1, 1), # test single metal with genecp
])

def test_confgen_exp_rules(folder, smiles, E_confs_no_rules, E_confs_rules, com_files, charge):
    # runs the program with the different tests
    cmd_exp_rules = ['python', '-m', 'pyconfort', '--varfile', params_file]

    round_E_confs_no_rules,round_E_confs_rules,test_round_E_confs_no_rules,test_round_E_confs_rules,test_charge,test_com_files = conf_gen_exp_rules(path_exp_rules, folder, precision_exp_rules, cmd_exp_rules, smiles, E_confs_no_rules, E_confs_rules)

    assert str(test_round_E_confs_no_rules) == str(round_E_confs_no_rules)
    assert str(test_round_E_confs_rules) == str(round_E_confs_rules)
    assert str(charge) == str(test_charge)
    assert str(com_files) == str(test_com_files)

# MISSING CHECKS:
# experimental rules for confgen SDF to COM, individually assigned
