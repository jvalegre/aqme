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
@pytest.mark.parametrize("folder, smiles, params_file, charge_orca, multiplicity_orca, job_type",
[
    # tests for ORCA in QPREP
    ('Organic_molecules', 'pentane.smi', 'params_test28.yaml', 0, 1, 'QPREP'), # ref test for ORCA
    
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test6.yaml', 1, 1, 'QPREP'), # test with SUMM and degree = 30
])

def test_confgen_organic(folder, smiles, params_file, n_confs_organic, prefilter_confs_rdkit_organic_organic, filter_confs_rdkit_organic, E_confs, charge_organic, multiplicity_organic, dihedral, xTB_ANI):
    # runs the program with the different tests
    cmd_organic = ['python', '-m', 'pyconfort', '--varfile', params_file]

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
