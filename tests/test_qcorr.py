#!/usr/bin/env python

######################################################.
# 	        Testing QCORR with pytest 	             #
######################################################.

import os
import glob
import pytest
from definitions_testing import analysis,single_point

# saves the working directory
path_analysis_dup_sp = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, file, params_file, type_of_job",
[
    # tests of the analysis part (I use file as the output LOG files)
    ('Analysis', 'CH4_Normal_termination.log', 'params_QCORR_test.yaml', 'analysis'), # test normal termination
    ('Analysis', 'CH4_z_Duplicate.LOG', 'params_QCORR_test.yaml', 'analysis'), # test duplicates
    ('Analysis', 'Basis_set_error1.LOG', 'params_QCORR_test.yaml', 'analysis'), # test incompatibilities with gen/genecp
    ('Analysis', 'Basis_set_error2.LOG', 'params_QCORR_test.yaml', 'analysis'), # test incompatibilities with gen/genecp
    ('Analysis', 'MeOH_Error_termination.LOG', 'params_QCORR_test.yaml', 'analysis'), # test error terminations
    ('Analysis', 'Imag_freq.log', 'params_QCORR_test.yaml', 'analysis'), # test imaginary frequencies
    ('Analysis', 'Imag_freq_low.log', 'params_QCORR_test.yaml', 'analysis'), # test cut-off for imag freqs
    ('Analysis', 'MeOH_SCF_error.out', 'params_QCORR_test.yaml', 'analysis'), # test SCF errors
    ('Analysis', 'MeOH_Unfinished.OUT', 'params_QCORR_test.yaml', 'analysis'), # test unfinished calculations
    ('Analysis', 'csv', 'params_QCORR_test.yaml', 'analysis'), # test final csv file with results
    ('Analysis', 'dat', 'params_QCORR_test.yaml', 'analysis'), # test final csv file with results
    ('Analysis2', 'Imag_freq.log', 'params_QCORR_test2.yaml', 'analysis'), # test amplitude for moving imag freqs
    # tests for single points
    ('Single_point', 'CH4_Normal_termination.log', 'params_QCORR_test.yaml', 'single_point'), # test single-point generation
])

def test_analysis_dup_sp(folder, file, params_file, type_of_job):
    # runs the program with the different tests
    cmd_pyconfort = ['python', '-m', 'pyconfort', '--varfile', params_file]

    if type_of_job == 'analysis':
        if file != 'csv' and file != 'dat':
            outlines = analysis(path_analysis_dup_sp, cmd_pyconfort, folder, file)

            if file == 'Basis_set_error1.LOG' or file == 'Basis_set_error2.LOG':
                os.chdir(path_analysis_dup_sp+'/'+folder+'/input_files/run_2')
                assert file.split('.')[0]+'.com' not in glob.glob('*.*')

            # full test of the generated com file (only coordinates for the next files)
            elif file == 'MeOH_Error_termination.LOG':
                line = '# wb97xd/genecp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=500) scrf=(SMD,solvent=Chloroform)'
                line_number = 2
                assert outlines[line_number].find(line) > -1

                line2 = '0 1'
                line_number2 = 6
                assert outlines[line_number2].find(line2) > -1

                line3 = 'H  -1.14928800  -0.80105100  -0.00024300'
                line_number3 = 12
                assert outlines[line_number3].find(line3) > -1

                line4 = 'H O 0'
                line_number4 = 14
                assert outlines[line_number4].find(line4) > -1

                line5 = 'def2svp'
                line_number5 = 15
                assert outlines[line_number5].find(line5) > -1

                line6 = 'C 0'
                line_number6 = 17
                assert outlines[line_number6].find(line6) > -1

                line7 = 'LANL2DZ'
                line_number7 = 18
                assert outlines[line_number7].find(line7) > -1

                line8 = 'C 0'
                line_number8 = 21
                assert outlines[line_number8].find(line8) > -1

            elif file == 'Imag_freq.log':
                if folder == 'Analysis':
                    line = 'H  -0.56133100   0.63933100  -0.67133100'
                    line_number = 10
                    assert outlines[line_number].find(line) > -1
                elif folder == 'Analysis2':
                    line = 'H  -0.35133100   0.66333100  -0.79133100'
                    line_number = 10
                    assert outlines[line_number].find(line) > -1

            elif file == 'MeOH_SCF_error.out':
                line = 'H  -1.04798700   0.80281000  -0.68030200'
                line_number = 10
                assert outlines[line_number].find(line) > -1

            elif file == 'MeOH_Unfinished.OUT':
                line = 'H  -1.04779100   0.87481300  -0.58663200'
                line_number = 10
                assert outlines[line_number].find(line) > -1

        else:
            df_QCORR,dat_files = analysis(path_analysis_dup_sp, cmd_pyconfort, folder, file)

        if file == 'csv':
            assert df_QCORR['Total files'][0] == 9
            assert df_QCORR['Normal termination'][0] == 2
            assert df_QCORR['Imaginary frequencies'][0] == 1
            assert df_QCORR['SCF error'][0] == 1
            assert df_QCORR['Basis set error'][0] == 2
            assert df_QCORR['Other errors'][0] == 1
            assert df_QCORR['Unfinished'][0] == 1
            assert df_QCORR['Duplicates'][0] == 1

        elif file == 'dat':
            assert 'pyCONFORT-QCORR-run_1_output.dat' in dat_files

        # check that the files are in their corresponding folders (just once)
        if file == 'CH4_Normal_termination.log':
            # all the files are moved from the original directory
            os.chdir(path_analysis_dup_sp+'/'+folder)
            assert len(glob.glob('*.LOG')) == 0
            assert len(glob.glob('*.log')) == 0
            assert len(glob.glob('*.out')) == 0
            assert len(glob.glob('*.OUT')) == 0
            # new paths for the files
            os.chdir(path_analysis_dup_sp+'/'+folder+'/duplicates/run_1')
            assert len(glob.glob('*.LOG')) == 1
            os.chdir(path_analysis_dup_sp+'/'+folder+'/success/output_files')
            assert len(glob.glob('*.log')) == 2
            os.chdir(path_analysis_dup_sp+'/'+folder+'/failed/run_1/unfinished')
            assert len(glob.glob('*.OUT')) == 1
            os.chdir(path_analysis_dup_sp+'/'+folder+'/failed/run_1/imag_freq')
            assert len(glob.glob('*.log')) == 1
            os.chdir(path_analysis_dup_sp+'/'+folder+'/failed/run_1/error/basis_set_error')
            assert len(glob.glob('*.LOG')) == 2
            os.chdir(path_analysis_dup_sp+'/'+folder+'/failed/run_1/error/scf_error')
            assert len(glob.glob('*.out')) == 1
            os.chdir(path_analysis_dup_sp+'/'+folder+'/failed/run_1/error/unknown_error')
            assert len(glob.glob('*.LOG')) == 1

    elif type_of_job == 'single_point':
        # test there are only 2 functional/basis set combinations (instead of 4)
        os.chdir(path_analysis_dup_sp+'/Analysis/success/G16-SP_input_files/')
        assert len(glob.glob('*')) == 2
        assert 'b3lyp-321g' in glob.glob('*')
        assert 'wb97xd-def2svp' in glob.glob('*')

        # test the input files are generate correctly
        os.chdir(path_analysis_dup_sp+'/Analysis/success/G16-SP_input_files/b3lyp-321g')
        files_sp = glob.glob('*.com')
        assert len(files_sp) == 2
        # test for the suffix option
        assert 'CH4_Normal_termination_SPC.com' in files_sp

        # test the content of the input files
        outlines = single_point(path_analysis_dup_sp, 'Analysis', file)

        line = '# b3lyp/genecp nmr = giao'
        line_number = 2
        assert outlines[line_number].find(line) > -1

        line2 = '5 3'
        line_number2 = 6
        assert outlines[line_number2].find(line2) > -1

        line3 = 'H   0.63133100  -0.63133100  -0.63133100'
        line_number3 = 11
        assert outlines[line_number3].find(line3) > -1

        line4 = 'C 0'
        line_number4 = 15
        assert outlines[line_number4].find(line4) > -1

        line5 = '321g'
        line_number5 = 16
        assert outlines[line_number5].find(line5) > -1

        line6 = 'H 0'
        line_number6 = 18
        assert outlines[line_number6].find(line6) > -1

        line7 = 'LANL2TZ'
        line_number7 = 19
        assert outlines[line_number7].find(line7) > -1

        line8 = 'H 0'
        line_number8 = 22
        assert outlines[line_number8].find(line8) > -1

        # extra test for final line in SP
        line9 = '-1'
        line_number9 = 13
        assert outlines[line_number9].find(line9) > -1
