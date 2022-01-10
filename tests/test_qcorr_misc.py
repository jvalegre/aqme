#!/usr/bin/env python

######################################################.
# 	   Miscellaneous QCORR tests with pytest 	     #
######################################################.

import os
import subprocess
import glob
import pytest
from definitions_testing import analysis,single_point

# saves the working directory
path_qcorr = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, params_file",
[
    # miscellaneous tests of the analysis part
    ('Analysis_misc', 'params_QCORR_misc_test.yaml'), # test QCORR and SP
])

def test_analysis_dup_sp(folder, params_file):
    # runs the program with the different tests
    cmd_qcorr_misc = ['python', '-m', 'aqme', '--varfile', params_file]

    os.chdir(path_qcorr+'/'+folder)
    subprocess.call(cmd_qcorr_misc)

    # check that the 2 failing calcs generated COM files
    os.chdir(path_qcorr+'/'+folder+'/input_files/run_2')
    assert len(glob.glob('*.com')) == 2

    outlines = com_lines('MeOH_Error_termination.com')

    assert outlines[2].find('# wb97xd/genecp  empiricaldispersion=GD3BJ opt=(maxcycles=500) scrf=(SMD,solvent=Chloroform)') > -1
    assert outlines[18].find('S   1   1.00') > -1
    assert outlines[37].find('C-ECP     4     78') > -1

    # check that the 1 normally terminating calc generated a COM file
    os.chdir(path_qcorr+'/'+folder+'/success/G16-SP_input_files/wb97xd-def2svp')
    assert len(glob.glob('*.com')) == 1

    outlines2 = com_lines('CH4_Normal_termination_SPC.com')

    assert outlines2[2].find('# wb97xd/genecp  empiricaldispersion=GD3BJ scrf=(SMD,solvent=Acetone)') > -1
    assert outlines2[17].find('S   1   1.00') > -1
    assert outlines2[36].find('C-ECP     4     78') > -1

    # check that the final destination of the files is correct
    os.chdir(path_analysis_dup_sp+'/'+folder+'/success/output_files')
    assert len(glob.glob('*.log')) == 1
    os.chdir(path_analysis_dup_sp+'/'+folder+'/failed/run_1/imag_freq')
    assert len(glob.glob('*.log')) == 1
    os.chdir(path_analysis_dup_sp+'/'+folder+'/failed/run_1/error/unknown_error')
    assert len(glob.glob('*.LOG')) == 1
