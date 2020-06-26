#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#       Conformer generation of other molecules      #
######################################################.

import os
import pytest
from definitions_testing import conf_gen

# saves the working directory
path_others = os.getcwd()
# decimal digits for comparing E
precision_others = 5

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, n_confs_others, prefilter_confs_rdkit_others, filter_confs_rdkit, E_confs, charge, multiplicity_others, dihedral, xTB_ANI1, metal, template",
[
    # tests of input files with different formats and charges. I included made up .smi names just to be coherent with the other tests
    ('Input_files', 'charged.csv', 'params_format_test1.yaml', 20, 19, 0, [-252.7254], 1, 1, False, False, False, False), # test smi with auto detecting of charges
    ('Input_files', 'charged.cdx', 'params_format_test2.yaml', 20, 19, 0, [-252.7254], 1, 1, False, False, False, False), # test cdx
    ('Input_files', 'pentane.com', 'params_format_test3.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 1, 2, False, False, False, False), # test com with charge 1 and multiplicity 2
    ('Input_files', 'pentane.gjf', 'params_format_test4.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], -1, 2, False, False, False, False), # test gjf with charge -1 and multiplicity 2
    ('Input_files', 'charged.sdf', 'params_format_test5.yaml', 20, 19, 0, [-252.7254], 1, 1, False, False, False, False), # test sdf
    ('Input_files', 'charged.smi', 'params_format_test6.yaml', 20, 19, 0, [-252.7254], 1, 1, False, False, False, False), # test smi
    ('Input_files', 'charged.mol', 'params_format_test7.yaml', 20, 19, 0, [-252.7254], 1, 1, False, False, False, False), # test mol
    ('Input_files', 'charged.mol2', 'params_format_test8.yaml', 20, 19, 0, [-252.7254], 1, 1, False, False, False, False), # test mol2
    ('Input_files', 'pentane.xyz', 'params_format_test9.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 1, 2, False, False, False, False), # test mol2
    ('Input_files', 'multicharge.smi', 'params_format_test10.yaml', 20, 19, 0, [-252.7254], 2, 1, False, False, False, False), # test smi with auto detecting of charges and +2 and -2 charges
    # tests that will check if the code crushes when using combinations of organic molecules and metal complexes
    ('Multiple', 'pentane_Pd_blank_lines.smi', 'params_comb_test1.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 1, False, False, False, False), # test pentane + Pd complex with blank lines
    ('Multiple', 'pentane_Pd.smi', 'params_comb_test2.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 1, False, False, False, False), # test pentane + Pd complex
    ('Multiple', 'pentane_Pd_template.smi', 'params_comb_test3.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 1, False, False, False, 'squareplanar'), # test pentane + Pd complex with template
    ('Multiple', 'pentane_Ag_Au.smi', 'params_comb_test4.yaml', 20, 17, 0, [-5.27175, -4.44184, -3.84858], 0, 1, False, False, False, False), # test pentane + 2 metals
    ('Multiple', 'pentane_Pd.smi', 'params_comb_test5.yaml', 432, 'nan', 4, [-5.27175,-4.44184,-3.84858,-1.57172], 0, 1, True, False, False, False), # test pentane + Pd + dihedral
])

def test_confgen_others(folder, smiles, params_file, n_confs_others, prefilter_confs_rdkit_others, filter_confs_rdkit, E_confs, charge, multiplicity_others, dihedral, xTB_ANI1, metal, template):
    # runs the program with the different tests
    cmd_others = ['python', '-m', 'pyconfort', '--varfile', params_file]

    test_init_rdkit_confs,test_prefilter_rdkit_confs,test_filter_rdkit_confs,round_confs_others,test_round_confs_others,test_charge,test_unique_confs,count,charge_com_others,multiplicity_com_others = conf_gen(path_others, precision_others, cmd_others, folder, smiles, E_confs, dihedral, xTB_ANI1, metal, template)

    # the assert statements are placed here, otherwise pytest doesn't explain the AssertionError
    # idx is used since I'm focusing on the pentane molecule at the end of the smiles list when using 2 smiles
    if folder == 'Multiple':
        idx = 1
    else:
        idx = 0

    assert str(n_confs_others) == str(test_init_rdkit_confs[idx])
    if prefilter_confs_rdkit_others != 'nan':
        assert str(prefilter_confs_rdkit_others) == str(test_prefilter_rdkit_confs[idx])
    assert str(filter_confs_rdkit) == str(test_filter_rdkit_confs[idx])

    assert str(round_confs_others) == str(test_round_confs_others)
    assert str(charge) == str(test_charge[idx])

    # make sure the COM files have the right charge and multiplicity
    assert str(charge_com_others) == str(charge)
    assert str(multiplicity_com_others) == str(multiplicity_others)
