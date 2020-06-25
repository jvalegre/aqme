#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#          template number of files tests            #
######################################################.

import os
import glob
import subprocess
import pytest
from definitions_testing import calc_genecp

# saves the working directory
path_misc = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("smiles, params_file, goal, goal_2, job_type",
[
    # multiple tests for the single-point option
    ('pentane.smi', 'params_sp_test1.yaml', 'nan', 'nan', 'sp_charge'),
    # check that analysis creates a CSV file with a summary of the results
    ('nan', 'nan', 'nan', 'nan', 'csv_analysis'),
    # check that the duplicate part moves the duplicates into a different folder
    ('nan', 'nan', 'nan', 'nan', 'dup'),
    # check that the duplicate part creates a CSV file with a summary of the results
    ('nan', 'nan', 'nan', 'nan', 'dup_csv'),
])

def test_confgen_misc(smiles, params_file, goal, goal_2, job_type):
    os.chdir(path_misc+'/Misc_2')
    # runs the program with the different tests
    if params_file == 'params_sp_test1.yaml':
        cmd_misc_2 = ['python', '-m', 'pyconfort', '--varfile', params_file]
        subprocess.call(cmd_misc_2)

    if job_type == 'sp_charge':
        os.chdir(path_misc+'/Misc_2/finished/single_point_input_files/wb97xd-def2svp')
        # this part will fail if the suffix is not working correctly
        file = 'CH4_Normal_termination_spc.com'
        atom = 'C' # not used really
        _,NBO,_,_,charge_com,multiplicity_com  = calc_genecp(file, atom)

        # Check that it adds the last line for NBO ($NBO $END)
        assert NBO == 1
        # Check charge and multiplicity explicitely added by the user
        assert charge_com == 3
        assert multiplicity_com == 5

    # this part is linked to the test_analysis_dup_sp.py part
    if job_type == 'csv_analysis':
        os.chdir(path_misc+'/Analysis')
        df = pd.read_csv('Duplicates Data.csv')
        assert df['finished'][0] == 3
        assert df['imaginary_frequencies'][0] == 1
        assert df['atomic_basis_error'][0] == 1
        assert df['SCF_error'][0] == 1
        assert df['unknown_error'][0] == 1
        assert df['failed_unfinished'][0] == 1

    if job_type == 'dup':
        # goes to the folder with the duplicates and try to find them
        os.chdir(path_misc+'/Misc_2/duplicates')
        assert 'CH4_Normal_termination_format.OUT' in glob.glob('*')

    if job_type == 'dup_csv':
        df = pd.read_csv('Duplicates Data.csv')
        assert df['finished'][0] == 1
        assert df['duplicates'][0] == 1
