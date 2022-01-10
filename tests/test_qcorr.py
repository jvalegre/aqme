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
@pytest.mark.parametrize("folder, file, command_line, type_of_job",
[
    # tests of the analysis part (I use file as the output LOG files)
    ('Analysis', 'CH4_Normal_termination.log', 'python -m aqme --qcorr gaussian', 'analysis'), # test normal termination
    ('Analysis', 'CH4_z_Duplicate.LOG', None, 'analysis'), # test duplicates
    ('Analysis', 'Basis_set_error1.LOG', None, 'analysis'), # test incompatibilities with gen/genecp
    ('Analysis', 'Basis_set_error2.LOG', None, 'analysis'), # test incompatibilities with gen/genecp
    ('Analysis', 'MeOH_Error_termination.LOG', None, 'analysis'), # test error terminations
    ('Analysis', 'Imag_freq.log', None, 'analysis'), # test imaginary frequencies
    ('Analysis', 'Imag_freq_low.log', None, 'analysis'), # test cut-off for imag freqs
    ('Analysis', 'MeOH_SCF_error.out', None, 'analysis'), # test SCF errors
    ('Analysis', 'MeOH_Unfinished.OUT', None, 'analysis'), # test unfinished calculations
    errortype = 'before_E_calculation'
    ('Analysis', 'TS_CH3HCH3.log', None, 'analysis'), # test for Normally terminated TSs
    error TS
    imag freq TS
    0 freqs in TS (check if the errortype is different than when you have 2 imag freqs)
    check ts_opt
    spin cont
    isomerization with com/gjf
    isomerization with csv
    isomeriz with csv for TSs
    bondi for isomeriz
    rcov for isomeriz
    tell if the imag freq from a TS is right based on displacement
    ('Analysis', 'CH4_Fail_freq_only.log', None, 'analysis'), # test for Normal terminated OPT and unfinished freq
    ('Analysis', 'csv', None, 'analysis'), # test final csv file with results
    ('Analysis', 'dat', None, 'analysis'), # test final csv file with results
    ('Analysis2', 'Imag_freq.log', 'params_QCORR_test2.yaml', 'analysis'), # test amplitude for moving imag freqs
    s**2 threshold
    e for duplicates threshold
    how do we check the folder create the right RUN numbers?
    # tests for single points
    ('Single_point', 'CH4_Normal_termination.log', 'params_QCORR_test.yaml', 'single_point'), # test single-point generation
    nocheck
    overwrites charge/mult
])

def test_analysis_dup_sp(folder, file, command_line, type_of_job):
    # runs the program with the different tests
    cmd_aqme = [command_line]

    if type_of_job == 'analysis':
        if file != 'csv' and file != 'dat':
            outlines = analysis(path_analysis_dup_sp, cmd_aqme, folder, file)

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
            df_QCORR,dat_files = analysis(path_analysis_dup_sp, cmd_aqme, folder, file)

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
        LoTs_SP = ['b3lyp-321g','wb97xd-def2svp']
        for LoT_SP in LoTs_SP:
            os.chdir(path_analysis_dup_sp+'/Analysis/success/G16-SP_input_files/'+LoT_SP)
            files_sp = glob.glob('*.com')
            # only normal terminations generate COM files for SP
            assert len(files_sp) == 2
            # test for the suffix option
            assert 'CH4_Normal_termination_SPC.com' in files_sp

            # test the content of the input files
            outlines = single_point(path_analysis_dup_sp, 'Analysis', file)

            if LoT_SP == 'b3lyp-321g':
                line = '# b3lyp/genecp nmr = giao'
            else:
                line = '# wb97xd/genecp nmr = giao'
            line_number = 2
            assert outlines[line_number].find(line) > -1

            line2 = '5 3'
            line_number2 = 6
            assert outlines[line_number2].find(line2) > -1

            line3 = 'H   0.63133100  -0.63133100  -0.63133100'
            line_number3 = 11
            assert outlines[line_number3].find(line3) > -1

            line4 = 'C 0'
            line_number4 = 13
            assert outlines[line_number4].find(line4) > -1

            if LoT_SP == 'b3lyp-321g':
                line5 = '321g'
            else:
                line5 = 'def2svp'
            line_number5 = 14
            assert outlines[line_number5].find(line5) > -1

            line6 = 'H 0'
            line_number6 = 16
            assert outlines[line_number6].find(line6) > -1

            if LoT_SP == 'b3lyp-321g':
                line7 = 'LANL2TZ'
            else:
                line7 = 'LANL2DZ'
            line_number7 = 17
            assert outlines[line_number7].find(line7) > -1

            line8 = 'H 0'
            line_number8 = 20
            assert outlines[line_number8].find(line8) > -1

            if LoT_SP == 'b3lyp-321g':
                line9 = 'LANL2TZ'
            else:
                line9 = 'LANL2DZ'
            line_number9 = 21
            assert outlines[line_number9].find(line9) > -1

            # extra test for final line in SP
            line10 = '-1'
            line_number10 = 23
            assert outlines[line_number10].find(line10) > -1
