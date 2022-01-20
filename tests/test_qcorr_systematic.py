#!/usr/bin/env python

######################################################.
# 	        Testing QCORR with pytest 	             #
######################################################.

import os
import glob
import subprocess
import shutil
import pandas as pd
import pytest
from definitions_testing import analysis,single_point

# saves the working directory
path_analysis_syst_initial = os.getcwd()+'/Analysis_syst'
path_analysis_syst = path_analysis_syst_initial+'/QCALC/G16'

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("file, type_of_job",
[
    # tests of the analysis part (I use file as the output LOG files)
    ('CH4_Normal_termination_syst.log', 'analysis'), # test normal termination
    ('CH4_z_Duplicate_syst.LOG', 'analysis'), # test duplicates
    ('Basis_set_error1_syst.LOG', 'analysis'), # test incompatibilities with gen/genecp
    ('Basis_set_error2_syst.LOG', 'analysis'), # test incompatibilities with gen/genecp
    ('MeOH_Error_termination_syst.LOG', 'analysis'), # test error terminations
    ('Imag_freq_syst.log', 'analysis'), # test imaginary frequencies
    ('MeOH_SCF_error_syst.out', 'analysis'), # test SCF errors
    ('MeOH_Unfinished_syst.OUT', 'analysis'), # test unfinished calculations
    ('csv', 'analysis'), # test final csv file with results
    ('dat', 'analysis'), # test final csv file with results
    # tests for single points
    ('CH4_Normal_termination_syst.log', 'single_point'), # test single-point generation
])

def test_analysis_dup_sp(file, type_of_job):
    # runs the program with the different tests
    cmd_aqme = ['python', '-m', 'aqme', '--varfile', 'params_QCORR_syst_test.yaml']
    LoTs = ['b3lyp-6-31G(d)','wb97xd-def2svp']

    if type_of_job == 'analysis':
        if file == 'CH4_Normal_termination_syst.log':
            # first run of aqme
            os.chdir(path_analysis_syst_initial)
            subprocess.call(cmd_aqme)

            # moving files from the first run to simulate that the COM files from input files were performed
            for LoT in LoTs:
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/run_1/error/basis_set_error')
                shutil.copy('Basis_set_error1_syst.LOG', path_analysis_syst+'/'+LoT+'/input_files/run_2')
                shutil.copy('Basis_set_error2_syst.LOG', path_analysis_syst+'/'+LoT+'/input_files/run_2')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/run_1/error/scf_error')
                shutil.copy('MeOH_SCF_error_syst.out', path_analysis_syst+'/'+LoT+'/input_files/run_2')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/run_1/error/unknown_error')
                shutil.copy('MeOH_Error_termination_syst.LOG', path_analysis_syst+'/'+LoT+'/input_files/run_2')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/run_1/imag_freq')
                shutil.copy('Imag_freq_syst.log', path_analysis_syst+'/'+LoT+'/input_files/run_2')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/run_1/unfinished')
                shutil.copy('MeOH_Unfinished_syst.OUT', path_analysis_syst+'/'+LoT+'/input_files/run_2')
                os.chdir(path_analysis_syst+'/'+LoT+'/success/output_files')
                shutil.copy('CH4_Normal_termination_syst.log', path_analysis_syst+'/'+LoT+'/input_files/run_2')
                os.chdir(path_analysis_syst+'/'+LoT+'/duplicates/run_1')
                shutil.copy('CH4_z_Duplicate_syst.LOG', path_analysis_syst+'/'+LoT+'/input_files/run_2')

            # second run of aqme
            os.chdir(path_analysis_syst_initial)
            subprocess.call(cmd_aqme)

        for LoT in LoTs:
            files_com = ['MeOH_Error_termination_syst.com','Imag_freq_syst.com','MeOH_SCF_error_syst.com','MeOH_Unfinished_syst.com']
            # check that the com files were generated correctly in runs 1 through 3
            os.chdir(path_analysis_syst+'/'+LoT+'/input_files/run_1')
            assert len(glob.glob('*.com')) == 8
            run_cycles_com = ['run_2','run_3']
            for run_cycle in run_cycles_com:
                os.chdir(path_analysis_syst+'/'+LoT+'/input_files/'+run_cycle)
                assert len(glob.glob('*.com')) == 4

                for file_com in files_com:
                    outfile = open(file_com,"r")
                    outlines = outfile.readlines()

                    line = '%mem=96GB'
                    line_number = 0
                    assert outlines[line_number].find(line) > -1

                    line = '%nprocshared=24'
                    line_number = 1
                    assert outlines[line_number].find(line) > -1

                    line = '0 1'
                    line_number = 6
                    assert outlines[line_number].find(line) > -1

                    if file_com == 'MeOH_Error_termination_syst.com' and file == 'MeOH_Error_termination_syst.LOG':
                        if LoT == 'b3lyp-6-31G(d)':
                            line2 = '# b3lyp/genecp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=500) scrf=(SMD,solvent=Chloroform)'

                        else:
                            line2 = '# wb97xd/genecp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=500) scrf=(SMD,solvent=Chloroform)'
                        line_number2 = 2
                        assert outlines[line_number2].find(line2) > -1

                        line3 = 'H  -1.14928800  -0.80105100  -0.00024300'
                        line_number3 = 12
                        assert outlines[line_number3].find(line3) > -1

                        line4 = 'H O 0'
                        line_number4 = 14
                        assert outlines[line_number4].find(line4) > -1

                        if LoT == 'b3lyp-6-31G(d)':
                            line5 = '6-31G(d)'
                        else:
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

                        line8 = 'LANL2DZ'
                        line_number8 = 22
                        assert outlines[line_number8].find(line8) > -1

                    elif file_com == 'Imag_freq_syst.com' and file == 'Imag_freq_syst.log':
                        line2 = 'H  -0.56133100   0.63933100  -0.67133100'
                        line_number2 = 10
                        assert outlines[line_number2].find(line2) > -1

                    elif file_com == 'MeOH_SCF_error_syst.com' and file == 'MeOH_SCF_error_syst.out':
                        if LoT == 'b3lyp-6-31G(d)':
                            line2 = '# b3lyp/genecp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=500) scrf=(SMD,solvent=Chloroform) scf=qc'

                        else:
                            line2 = '# wb97xd/genecp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=500) scrf=(SMD,solvent=Chloroform) scf=qc'

                        line_number2 = 2
                        assert outlines[line_number2].find(line2) > -1

                        line3 = 'H  -1.04798700   0.80281000  -0.68030200'
                        line_number3 = 10
                        assert outlines[line_number3].find(line3) > -1

                    elif file_com == 'MeOH_Unfinished_syst.com' and file == 'MeOH_Unfinished_syst.OUT':
                        line2 = 'H  -1.04779100   0.87481300  -0.58663200'
                        line_number2 = 10
                        assert outlines[line_number2].find(line2) > -1

                    outfile.close()

            # check that the log files are placed into the right place for the different runs
            run_cycles = ['run_1','run_2']
            for run_cycle in run_cycles:
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/'+run_cycle+'/error/basis_set_error')
                assert len(glob.glob('*.*')) == 2
                assert 'Basis_set_error1_syst.LOG' in glob.glob('*.*')
                assert 'Basis_set_error2_syst.LOG' in glob.glob('*.*')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/'+run_cycle+'/error/scf_error')
                assert len(glob.glob('*.*')) == 1
                assert 'MeOH_SCF_error_syst.out' in glob.glob('*.*')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/'+run_cycle+'/error/unknown_error')
                assert len(glob.glob('*.*')) == 1
                assert 'MeOH_Error_termination_syst.LOG' in glob.glob('*.*')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/'+run_cycle+'/imag_freq')
                assert len(glob.glob('*.*')) == 1
                assert 'Imag_freq_syst.log' in glob.glob('*.*')
                os.chdir(path_analysis_syst+'/'+LoT+'/failed/'+run_cycle+'/unfinished')
                assert len(glob.glob('*.*')) == 1
                assert 'MeOH_Unfinished_syst.OUT' in glob.glob('*.*')
                os.chdir(path_analysis_syst+'/'+LoT+'/success/output_files')
                assert len(glob.glob('*.*')) == 1
                assert 'CH4_Normal_termination_syst.log' in glob.glob('*.*')
                os.chdir(path_analysis_syst+'/'+LoT+'/duplicates/'+run_cycle)
                dup_files = []
                valid_output_formats = ['LOG','OUT','log','out']
                for potential_dup in glob.glob('*.*'):
                    if potential_dup not in dup_files:
                        if potential_dup.split('.')[1] in valid_output_formats:
                            dup_files.append(potential_dup)
                assert len(dup_files) == 1
                assert 'CH4_z_Duplicate_syst.LOG' in glob.glob('*.*')

            if file == 'csv':
                # Retrieving the generated CSV file
                os.chdir(path_analysis_syst+'/'+LoT+'/csv_files')

                run_cycles = ['run_1','run_2']
                for run_cycle in run_cycles:
                    df_QCORR = pd.read_csv('Analysis-Data-QCORR-'+run_cycle+'.csv')
                    assert df_QCORR['Total files'][0] == 8
                    assert df_QCORR['Normal termination'][0] == 1
                    assert df_QCORR['Imaginary frequencies'][0] == 1
                    assert df_QCORR['SCF error'][0] == 1
                    assert df_QCORR['Basis set error'][0] == 2
                    assert df_QCORR['Other errors'][0] == 1
                    assert df_QCORR['Unfinished'][0] == 1
                    assert df_QCORR['Duplicates'][0] == 1

            elif file == 'dat':
                os.chdir(path_analysis_syst+'/'+LoT+'/dat_files')
                run_cycles = ['run_1','run_2']
                for run_cycle in run_cycles:
                    assert 'aqme-QCORR-'+run_cycle+'_output.dat' in glob.glob('*.*')

    elif type_of_job == 'single_point':
        for LoT in LoTs:
            # test there are only 2 functional/basis set combinations (instead of 4)
            os.chdir(path_analysis_syst+'/'+LoT+'/success/G16-SP_input_files/')
            assert len(glob.glob('*')) == 2
            assert 'b3lyp-321g' in glob.glob('*')
            assert 'wb97xd-def2svp' in glob.glob('*')

            # test the input files are generated correctly
            LoTs_SP = ['b3lyp-321g','wb97xd-def2svp']
            for LoT_SP in LoTs_SP:
                os.chdir(path_analysis_syst+'/'+LoT+'/success/G16-SP_input_files/'+LoT_SP)
                # only normal terminations generate COM files for SP
                assert len(glob.glob('*.com')) == 1
                # test for the suffix option
                assert 'CH4_Normal_termination_syst_SPC.com' in glob.glob('*.com')

                # test the content of the input files
                outfile = open('CH4_Normal_termination_syst_SPC.com',"r")
                outlines = outfile.readlines()

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

                outfile.close()
