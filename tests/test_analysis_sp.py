#!/usr/bin/env python

######################################################.
# 	        Testing analysis with pytest 	         #
######################################################.

import os
import pytest
from definitions_testing import analysis,single_point

# saves the working directory
path_analysis_dup_sp = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, file, params_file, type_of_job",
[
    # tests of the analysis part (I use file as the output LOG files)
    ('Analysis', 'CH4_Normal_termination.log', 'params_analysis_test.yaml', 'analysis'), # test normal termination
    ('Analysis', 'Basis_set_error1.LOG', 'params_analysis_test.yaml', 'analysis'), # test incompatibilities with gen/genecp
    ('Analysis', 'Basis_set_error2.LOG', 'params_analysis_test.yaml', 'analysis'), # test incompatibilities with gen/genecp
    ('Analysis', 'MeOH_Error_termination.LOG', 'params_analysis_test.yaml', 'analysis'), # test error terminations
    ('Analysis', 'Imag_freq.log', 'params_analysis_test.yaml', 'analysis'), # test imaginary frequencies
    ('Analysis', 'MeOH_SCF_error.out', 'params_analysis_test.yaml', 'analysis'), # test SCF errors
    ('Analysis', 'MeOH_Unfinished.OUT', 'params_analysis_test.yaml', 'analysis'), # test unfinished calculations
    # tests for single points
    ('Single_point', 'CH4_freq.log', 'params_sp_test.yaml', 'single_point'), # test single-point generation
    ('Single_point', 'Pd_SP.LOG', 'params_sp_test.yaml', 'single_point'), # test single-point generation with genecp
])

def test_analysis_dup_sp(folder, file, params_file, type_of_job):
    # runs the program with the different tests
    cmd_pyconfort = ['python', '-m', 'pyconfort', '--varfile', params_file]

    if type_of_job == 'analysis':
        analysis(path_analysis_dup_sp, cmd_pyconfort, folder, file)

    elif type_of_job == 'single_point':
        count,NBO,pop,opt = single_point(path_analysis_dup_sp, cmd_pyconfort, folder, file)

        if file == 'Pd_SP.LOG':
            assert count == 2 # finds genecp for Pd
        elif file == 'CH4_freq.log':
            assert count == 0 # does not find genecp part

        assert NBO == 1 # finds final line for sp
        assert pop == 1 # finds input line for sp
        assert opt == 0 # it does not find standard opt option
