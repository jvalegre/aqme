#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#               Miscellaneous tests                  #
######################################################.

import os
import subprocess
import pytest
from definitions_testing import misc_sdf_test, misc_com_test, misc_genecp_test, misc_freq_test, misc_lot_test, misc_nocom_test

# saves the working directory
path_misc = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("smiles, params_file, goal, goal_2, job_type",
[
    # tests of templates to check if the code create less files when there are repeated ligands
    ('Pd_squareplanar.smi', 'params_Pd_test1.yaml', 3, 'nan','count_sdf'), # test squareplanar template with gen
    ('Pd_squareplanar2.smi', 'params_Pd_test4.yaml', 2, 'nan','count_sdf'), # test with 2 same ligands
    ('Rh_squarepyramidal.smi', 'params_Rh_test1.yaml', 13, 'nan','count_sdf'), # test squarepyramidal with 2 same ligands
    ('Rh_squarepyramidal2.smi', 'params_Rh_test3.yaml', 9, 'nan','count_sdf'), # test squarepyramidal with 3 same ligands
    # test if the program removes the *_rdkit.sdf files when using SUMM
    ('pentane.smi', 'params_test15.yaml', 'None', 1, 'count_sdf'), # test SUMM
    # test the lowest_only option from SDF files and from the conformer generator workflow
    ('pentane.smi', 'params_test_misc1.yaml', 1, 'nan', 'count_com'), # test lowest_only with SUMM
    ('pentane.smi', 'params_test_misc2.yaml', 1, 'nan', 'count_com'), # test lowest_only with no SUMM
    ('pentane_n_lowest.sdf', 'params_test_misc3.yaml', 1, 'nan', 'count_com'), # test lowest_only with no SUMM
    # test the lowest_n option from SDF files and from the conformer generator workflow
    ('pentane.smi', 'params_test_misc4.yaml', 2, 'nan', 'count_com'), # test n_lowest with SUMM
    ('pentane.smi', 'params_test_misc5.yaml', 2, 'nan', 'count_com'), # test n_lowest with no SUMM
    ('pentane_n_lowest.sdf', 'params_test_misc6.yaml', 2, 'nan', 'count_com'), # test n_lowest with no SUMM
    # test with lowest_n = True and lowest_only = True (it should work with lowest_only)
    ('pentane.smi', 'params_test_misc12.yaml', 2, 'nan', 'count_com'), # test n_lowest with SUMM
    # test to create COM files with SDF files (i.e. from RDKit generated files)
    ('pentane_n_lowest.sdf', 'params_test_misc7.yaml', 4, 'nan', 'count_com'), # test create COM from SDF
    # test to see if .txt, .yaml and .yml to get genecp is working
#MISSING    ('pentane.smi', 'params_test_misc8.yaml', 1, 1, 'genecp_txt'), # test manual genecp
    # test to see if the input line, frequencies, chk and other options are working correctly
    ('pentane.smi', 'params_test_misc9.yaml', 0, 1, 'freq_maxcycles'), # test frequencies, maxcycles, chk, mem, nprocs, solvent, dispersion
    # test for using multiple levels of theory
    ('pentane.smi', 'params_test_misc10.yaml', 1, 1, 'multiple_lot'), # test multiple levels of theory and genecp
    # test for checking QCORR=False
    ('pentane.smi', 'params_test_misc11.yaml', 0, 1, 'no_com_files'), # test gauss_write=False
])

def test_confgen_misc(smiles, params_file, goal, goal_2, job_type):
    # runs the program with the different tests
    cmd_misc = ['python', '-m', 'aqme', '--varfile', params_file]
    os.chdir(path_misc+'/'+'Misc')
    subprocess.call(cmd_misc)

    if job_type == 'count_sdf':
        test_goal,test_goal_2 = misc_sdf_test(path_misc, smiles)
        # First goal: counts the number of *_rdkit.sdf files created
        # Second goal: counts the number of *_rdkit_rotated.sdf files created

    if job_type == 'count_com':
        test_goal = misc_com_test(path_misc, smiles)
        # First goal: calculates the number of COM files generated

    if job_type == 'genecp_txt':
        test_goal,test_goal_2 = misc_genecp_test(path_misc, smiles)
        # Finds gen for C from the C_genecp.txt file. If it finds it, test_goal is 1
        # Finds ecp for C from the C_genecp.txt file. If it finds it, test_goal_2 is 1

    if job_type == 'freq_maxcycles':
        test_goal,test_goal_2,test_goal_3,test_goal_4,test_goal_5,test_goal_6,test_goal_7 = misc_freq_test(path_misc, smiles)
        # Should not find freq in the keyword line of the COM file (test_goal should be 0)
        # Should find opt=(maxcycles=250) in the keyword line of the COM file (test_goal should be 1)

    if job_type == 'multiple_lot':
        test_goal,test_goal_2 = misc_lot_test(path_misc, smiles)

    if job_type == 'no_com_files':
        test_goal,test_goal_2 = misc_nocom_test(path_misc, smiles)

    assert test_goal == goal
    if goal_2 != 'nan':
        assert test_goal_2 == goal_2
    if job_type == 'freq_maxcycles':
        # checks for CHK, memory, nprocs, solvent and empirical dispersion
        assert test_goal_3 == 1
        assert test_goal_4 == 1
        assert test_goal_5 == 1
        assert test_goal_6 == 1
        assert test_goal_7 == 1
