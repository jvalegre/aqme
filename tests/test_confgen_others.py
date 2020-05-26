#!/usr/bin/env python

######################################################.
# 		  Main file for testing with pytest: 	     #
#     conformer generation of organic molecules      #
######################################################.

import os
import pytest
from definitions_testing import conf_gen, only_check

# saves the working directory
path_others = os.getcwd()
# decimal digits for comparing E
precision_others = 5

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1",
[
    # tests of input files with different formats
    ('Input_files', 'pentane.csv', 'params_format_test1.yaml', 20, 16, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 0, False, False), # test csv
    ('Input_files', 'pentane.cdx', 'params_format_test2.yaml', 20, 16, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 0, False, False), # test cdx
    ('Input_files', 'pentane.com', 'params_format_test3.yaml', 20, 16, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 1, False, False), # test com with charge 1
    ('Input_files', 'pentane.gjf', 'params_format_test4.yaml', 20, 16, 0, [-5.27175, -4.44184, -3.84858, -1.57172], -1, False, False), # test gjf with charge -1
    ('Input_files', 'pentane.sdf', 'params_format_test5.yaml', 20, 16, 0, [-5.27175, -4.44184, -3.84858, -1.57172], 0, False, False), # test sdf
    # tests that will check if the code crushes when using combinations of organic molecules and metal complexes
    ('Multiple', 'pentane_Pd_blank_lines.smi', 'params_comb_test1.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test pentane + Pd complex with blank lines
    ('Multiple', 'pentane_Pd.smi', 'params_comb_test2.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test pentane + Pd complex
    ('Multiple', 'pentane_Pd_template.smi', 'params_comb_test3.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test pentane + Pd complex with template
    ('Multiple', 'pentane_Ag_Au.smi', 'params_comb_test4.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test pentane + 2 metals
    ('Multiple', 'pentane_Pd.smi', 'params_comb_test5.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test pentane + Pd + dihedral
    # tests to check that gen and genecp work correctly
    ('Genecp', 'Pd_squareplanar.smi', 'params_genecp_test1.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test gen
    ('Genecp', 'Pd_squareplanar.smi', 'params_genecp_test2.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False, False), # test genecp
])

def test_confgen_others(folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1):
    # runs the program with the different tests
    cmd_others = ['python', '-m', 'pyconfort', '--varfile', params_file]

    conf_gen(path_others, precision_others, cmd_others, folder, smiles, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1)
