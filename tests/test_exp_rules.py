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
params_file = 'params_Ir_exp_rules_test'

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("smiles, E_confs_no_rules, E_confs_rules, charge",
[
    # tests for conformer generation with experimental rules
    ('Ir_1', [68.86185, 69.01720], [68.86185], 0), # test single metal with genecp
    ('Ir_2', [], [69.60127], 1), # test single metal with genecp
    ('Ir_3', [], [47.91641, 48.05050], 1), # test single metal with genecp
    ('Ir_4', [], [56.25232, 56.75803], 0), # test single metal with genecp
    ('Ir_5', [], [124.67040], 0), # test single metal with genecp
    ('Ir_6', [129.91109, 130.12558], [129.91109], 1), # test single metal with genecp
])

def test_confgen_exp_rules(smiles, E_confs_no_rules, E_confs_rules, charge):
    # runs the program with the different tests
    cmd_exp_rules = ['python', '-m', 'pyconfort', '--varfile', params_file]

    conf_gen_exp_rules(path_exp_rules, precision_exp_rules, cmd_exp_rules, smiles, E_confs_no_rules, E_confs_rules, charge)

# MISSING CHECKS:
# experimental rules for confgen SDF to COM
