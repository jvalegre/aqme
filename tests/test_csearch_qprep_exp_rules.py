#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#   Experimental rules applied to metal complexes    #
######################################################.

import os
import pytest
from definitions_testing import Ir_exp_rules,Pd_exp_rules

# saves the working directory
path_exp_rules = os.getcwd()
# decimal digits for comparing E
precision_exp_rules = 5

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, E_confs_no_rules, E_confs_rules, com_files, charge, format",
[
    # tests for conformer generation with experimental Ir_bidentate_x3 rules from SMILES
    ('Ir_exp_rules', 'Ir_1', 'params_Ir_exp_rules_test.yaml', [68.86185, 69.01720], [68.86185], 1, 0, None), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_2', 'params_Ir_exp_rules_test.yaml', [69.54124, 69.60127, 69.63704], [69.60127], 1, 1, None), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_3', 'params_Ir_exp_rules_test.yaml', [48.28081, 48.3819, 48.64117, 48.92291], [48.38190], 1, 1, None), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_4', 'params_Ir_exp_rules_test.yaml', [57.24954, 57.4422, 57.48444, 57.84497], [57.24954], 1, 0, None), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_5', 'params_Ir_exp_rules_test.yaml', [122.76353, 123.88455, 124.00555, 125.27744, 125.43099], [125.43099], 1, 0, None), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_6', 'params_Ir_exp_rules_test.yaml', [129.91109, 130.12558, 130.49653], [129.91109], 1, 1, None), # test single metal with genecp
    ('Ir_exp_rules', 'Ir_7', 'params_Ir_exp_rules_test.yaml', [124.00469, 124.09251, 124.10256], [124.09251], 1, 0, None), # test single metal with genecp
    # tests from other input formats with Ir_bidentate_x3 rules
# MISSING - charges are messed up    ('Ir_exp_rules2', 'Ir_4', 'params_Ir_exp_rules_test_sdf.yaml', [57.24954, 57.4422, 57.48444, 57.84497], [57.24954], 1, 0, 'sdf'), # test from sdf
# MISSING - charges are messed up    ('Ir_exp_rules2', 'Ir_4', 'params_Ir_exp_rules_test_csv.yaml', [57.24954, 57.4422, 57.48444, 57.84497], [57.24954], 1, 0, 'csv'), # test from csv
# MISSING - charges are messed up    ('Ir_exp_rules2', 'Ir_4', 'params_Ir_exp_rules_test_cdx.yaml', [57.24954, 57.4422, 57.48444, 57.84497], [57.24954], 1, 0, 'cdx'), # test from cdx
    # tests for manually inputted rules (i.e. 'I-Pd-Cl, 180'). The 5th parameter in the list is the number of com files
    # ('Pd_exp_rules', 'Pd_exp_rules.smi', 'params_Pd_test1.yaml', 3, 1, 1, None, None), # test 2 manual rules that apply to the complex
    # ('Pd_exp_rules','Pd_exp_rules.smi', 'params_Pd_test2.yaml', 3, 3, 4, None, None), # test 1 manual rule that doesn't apply
    # ('Pd_exp_rules2', 'Pd_exp_rules2.smi', 'params_Pd_test3.yaml', 3, 3, 17, None, None), # test 1 manual rule with too many optins (3 C-Pd bonds for C-Pd-C)
    # ('Pd_exp_rules3', 'Pd_exp_rules3.smi', 'params_Pd_test4.yaml', 3, 2, 9, None, None), # test 1 manual rule with 2 identical atoms in each side (C-Pd-C)
])

def test_confgen_exp_rules(folder, smiles, params_file, E_confs_no_rules, E_confs_rules, com_files, charge, format):
    # runs the program with the different tests
    cmd_exp_rules = ['python', '-m', 'pyconfort', '--varfile', params_file]

    Ir_rules = ['Ir_exp_rules','Ir_exp_rules2']
    Pd_rules = ['Pd_exp_rules','Pd_exp_rules2','Pd_exp_rules3']

    if folder in Ir_rules:
        round_E_confs_no_rules,round_E_confs_rules,test_round_E_confs_no_rules,test_round_E_confs_rules,test_charge,test_com_files = Ir_exp_rules(path_exp_rules, folder, precision_exp_rules, cmd_exp_rules, smiles, E_confs_no_rules, E_confs_rules)

        assert str(test_round_E_confs_no_rules) == str(round_E_confs_no_rules)
        assert str(test_round_E_confs_rules) == str(round_E_confs_rules)
        assert str(charge) == str(test_charge)
        assert str(com_files) == str(test_com_files)

    if folder in Pd_rules:
        # readjusting the names of the variables
        sdf_created = E_confs_no_rules
        sdf_final = E_confs_rules

        test_sdf_created,test_sdf_final,test_com_files = Pd_exp_rules(path_exp_rules, folder, precision_exp_rules, cmd_exp_rules, smiles)

        assert str(sdf_created) == str(test_sdf_created)
        assert str(sdf_final) == str(test_sdf_final)
        assert str(com_files) == str(test_com_files)
