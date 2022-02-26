#!/usr/bin/env python

######################################################.
# 	        Testing QCORR with pytest 	             #
######################################################.

import os
from os import path
import glob
import pytest
import shutil
import subprocess
from pathlib import Path

# saves the working directory
path_main = os.getcwd()
path_qcorr = os.getcwd()+'/Example_workflows/QCORR_processing_QM_outputs'

# QCORR tests
@pytest.mark.parametrize("test_type, init_folder, file, command_line, target_folder",
[
    # QCORR analysis tests with standard options
    ('analysis', 'QCORR_1', 'CH4.log', 'python -m aqme --qcorr ', 'successful_QM_outputs'), # test successful termination
    ('analysis', 'QCORR_1', 'MeOH_G09.log', None, 'successful_QM_outputs'), # test successful termination
    ('analysis', 'QCORR_1', 'z_CH4_duplicate.log', None, 'unsuccessful_QM_outputs/run_1/duplicates/'), # test duplicates
    ('analysis', 'QCORR_1', 'Basis_set_error1.LOG', None, 'unsuccessful_QM_outputs/run_1/error/basis_set_error/'), # test incompatibilities with basis sets
    ('analysis', 'QCORR_1', 'Basis_set_error2.LOG', None, 'unsuccessful_QM_outputs/run_1/error/basis_set_error/'), # test incompatibilities with basis sets
    ('analysis', 'QCORR_1', 'MeOH_G09_FAIL.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error/'), # test error terminations
    ('analysis', 'QCORR_1', 'CH2OH2_unfinished.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error/'), # test unfinished calculations
    ('analysis', 'QCORR_1', 'Imag_freq.log', None, 'unsuccessful_QM_outputs/run_1/extra_imag_freq/'), # test imaginary frequencies
    ('analysis', 'QCORR_1', 'imag_freq_no_opt.log', None, 'unsuccessful_QM_outputs/run_1/extra_imag_freq/'), # test imaginary frequencies without OPT
    ('analysis', 'QCORR_1', 'MeOH_SCF_error.log', None, 'unsuccessful_QM_outputs/run_1/error/scf_error/'), # test SCF errors
    ('analysis', 'QCORR_1', 'CH4_before_E.log', None, 'unsuccessful_QM_outputs/run_1/error/no_data/'), # test calcs that finish before any coords are printed
    ('analysis', 'QCORR_1', 'freq_conv_YYYN.log', None, 'unsuccessful_QM_outputs/run_1/freq_no_conv/'), # test calcs with freq calcs that did not converge after OPT
    ('analysis', 'QCORR_1', 'bpinene_spin_contamin.log', None, 'unsuccessful_QM_outputs/run_1/spin_contaminated/'), # test calcs with spin contamination
    ('analysis', 'QCORR_1', 'CH4_Fail_freq_only.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error/'), # test for Normal terminated OPT and unfinished freq
    ('analysis', 'QCORR_1', 'TS_CH3HCH3.log', None, 'successful_QM_outputs'), # test successful termination in TSs
    ('analysis', 'QCORR_1', 'TS_CH3HCH3_unfinished.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error/'), # test unfinished TSs
    ('analysis', 'QCORR_1', 'TS_CH3HCH3_imag_freq.log', None, 'unsuccessful_QM_outputs/run_1/extra_imag_freq/'), # test imaginary frequencies for TSs
    ('analysis', 'QCORR_1', 'TS_CH3HCH3_no_imag_freq.log', None, 'unsuccessful_QM_outputs/run_1/ts_no_imag_freq/'), # test imaginary frequencies for TSs
    ('analysis', 'QCORR_1', 'CH4_SP.log', None, 'successful_QM_outputs/SP_calcs/'), # test for single-point calcs
    ('analysis', 'QCORR_1', 'H_freq.log', None, 'successful_QM_outputs'), # test successful termination with 1 atom
    ('analysis', 'QCORR_1', 'H_SP.log', None, 'successful_QM_outputs/SP_calcs/'), # test for single-point calcs with 1 atom
    # CO2_linear_3freqs_FAIL.log
    # CO2_linear_4freqs.log
    ('analysis', 'QCORR_1', 'json', None, 'successful_QM_outputs/json_files'), # test for correct creation of json files (all the successful terminations only)
    ('analysis', 'QCORR_1', 'fullcheck', None, 'successful_QM_outputs/json_files'), # test for correct fullcheck option
    ('analysis', 'QCORR_1', 'csv', None, 'dat_files'), # test final csv file with results
    ('analysis', 'QCORR_1', 'dat', None, 'dat_files'), # test final dat file with results
    # isomerization with com/gjf
    # isomerization with csv (ongoing)
    # isomeriz with csv for TSs (ongoing)
    # check that the initial_QM_inputs folder doesnt get created unless there are COM files (and check that it DOES get created when there are com files)
    # move 1 error and 1 normal file to the main folder, and run qcorr, make sure the file is placed in succes and run_2/error
    (None, None, 'restore_original', None, None), # test single-point generation
    # SAME FOR COMMAND LINE, only general (i.e. mount of files in each folder in 1 test)
    # (None, None, 'restore_original', None, None), # test single-point generation
    # SAME FOR YAML file, only general (i.e. mount of files in each folder in 1 test)
    # (None, None, 'restore_original', None, None), # test single-point generation
    # bondi for isomeriz
    # rcov for isomeriz
    # check ts_opt
    # tell if the imag freq from a TS is right based on displacement (ongoing)
    # s**2 threshold
    # ('analysis', 'QCORR_1', 'Imag_freq_low.log', None), # test cut-off for imag freqs
    # e for duplicates threshold
    # json2input tests
    # ('single_point', 'QCORR_1', 'CH4.log', None), # test single-point generation
    # overwrites charge/mult
])

def test_analysis_dup_sp(test_type, init_folder, file, command_line, target_folder):
    # copy the test folders
    if not path.exists(f'{path_main}/Example_workflows_original'):
        shutil.copytree(f'{path_main}/Example_workflows', f'{path_main}/Example_workflows_original')

    # runs the program with the different tests
    w_dir_main=f'{path_qcorr}/{init_folder}'

    if test_type == 'analysis':
        if command_line is not None:
            cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_main, '--qm_files', '*.log']
            if init_folder == 'QCORR_1':
                cmd_aqme = cmd_aqme + ['--author', 'Juan V. Alegre-Requena']

            elif init_folder == 'QCORR_2':
                input_folder=f'{w_dir_main}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs'
                cmd_aqme = cmd_aqme + ['--isom', 'com', '--isom_inputs', input_folder]
            os.chdir(w_dir_main)
            subprocess.run(cmd_aqme)

        if file.split('.')[-1].lower() == 'log':
            assert path.exists(f'{w_dir_main}/{target_folder}/{file}')

        elif file == 'json':
            os.chdir(f'{w_dir_main}/{target_folder}')
            json_files = glob.glob('*.json')
            target_files = ['CH4.json', 'CO2_linear_4freqs.json', 'H_freq.json', 'MeOH_G09.json', 'TS_CH3HCH3.json']
            assert len(json_files) == 5
            assert sorted(json_files) == sorted(target_files)

        elif file == 'fullcheck':
            target_fullcheck = ['-- Full check analysis --\n']
            target_fullcheck.append('x  Different program used in the calculations:\n')
            target_fullcheck.append('     * Gaussian 09, Revision A.02 in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('     * Gaussian 16, Revision C.01 in:\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('       - MeOH_G09\n')
            target_fullcheck.append('x  Different grid used in the calculations:\n')
            target_fullcheck.append('     * sg1 in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('     * ultrafine in:\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('     * fine in:\n')
            target_fullcheck.append('       - MeOH_G09\n')
            target_fullcheck.append('x  Different lot used in the calculations:\n')
            target_fullcheck.append('     * M062X/3-21G in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('     * M062X/def2TZVP in:\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('     * PBE1PBE/CC-pVTZ in:\n')
            target_fullcheck.append('       - MeOH_G09\n')
            target_fullcheck.append('     * B3LYP/3-21G in:\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('o  Same dispersion (none) used in all the calculations\n')
            target_fullcheck.append('x  Different solvation used in the calculations:\n')
            target_fullcheck.append('     * gas_phase in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('     * scrf=(solvent=chloroform,pcm) in:\n')
            target_fullcheck.append('       - MeOH_G09')

            outfile = open(f'{w_dir_main}/{target_folder}/--QCORR_Fullcheck_Analysis--.dat', "r")
            outlines = outfile.readlines()
            outfile.close()
            assert target_fullcheck == outlines

        elif file == 'dat':
            target_dat = ['o  Analyzing output files in C:\\Users\\juanv\\OneDrive - Colostate\\Software Programas\\aqme\\Example_workflows\\QCORR_processing_QM_outputs\\QCORR_1\n']
            target_dat.append('\n')
            target_dat.append('Basis_set_error1.LOG: Termination = normal, Error type = atomicbasiserror\n')
            target_dat.append('Basis_set_error2.LOG: Termination = normal, Error type = atomicbasiserror\n')
            target_dat.append('bpinene_spin_contamin.log: Termination = normal, Error type = spin_contaminated\n')

            outfile = open(f'{w_dir_main}/{target_folder}/QCORR-run_1.dat', "r")
            outlines = outfile.readlines()
            outfile.close()

            for i,line in enumerate(target_dat):
                if i == 0:
                    assert outlines[i].find('o  Analyzing output files in') > -1
                else:
                    assert line == outlines[i]


    # leave the folders as they were initially
    if file == 'restore_original':
        os.chdir(path_main)
        shutil.rmtree(f'{path_main}/Example_workflows')
        filepath = Path(f'{path_main}/Example_workflows_original')
        filepath.rename(f'{path_main}/Example_workflows')


#     check if the files get moved to the target_folder
#     check that the COM files are fixed in the way they're supposed to be fixed in each case (or not
#     created in tohers such as spin cont)
#     check CSV gets the right nunmbers of each type
#     check the json2input works correctly with all the examples from the ntoebook

# ANTeS
#     if type_of_job == 'analysis':
#         if file != 'csv' and file != 'dat':
#             outlines = analysis(path_main, cmd_aqme, folder, file)

#             if file == 'Basis_set_error1.LOG' or file == 'Basis_set_error2.LOG':
#                 os.chdir(path_main+'/'+folder+'/input_files/run_2')
#                 assert file.split('.')[0]+'.com' not in glob.glob('*.*')

#             # full test of the generated com file (only coordinates for the next files)
#             elif file == 'MeOH_Error_termination.LOG':
#                 line = '# wb97xd/genecp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=500) scrf=(SMD,solvent=Chloroform)'
#                 line_number = 2
#                 assert outlines[line_number].find(line) > -1

#                 line2 = '0 1'
#                 line_number2 = 6
#                 assert outlines[line_number2].find(line2) > -1

#                 line3 = 'H  -1.14928800  -0.80105100  -0.00024300'
#                 line_number3 = 12
#                 assert outlines[line_number3].find(line3) > -1

#                 line4 = 'H O 0'
#                 line_number4 = 14
#                 assert outlines[line_number4].find(line4) > -1

#                 line5 = 'def2svp'
#                 line_number5 = 15
#                 assert outlines[line_number5].find(line5) > -1

#                 line6 = 'C 0'
#                 line_number6 = 17
#                 assert outlines[line_number6].find(line6) > -1

#                 line7 = 'LANL2DZ'
#                 line_number7 = 18
#                 assert outlines[line_number7].find(line7) > -1

#                 line8 = 'C 0'
#                 line_number8 = 21
#                 assert outlines[line_number8].find(line8) > -1

#             elif file == 'Imag_freq.log':
#                 if folder == 'analysis':
#                     line = 'H  -0.56133100   0.63933100  -0.67133100'
#                     line_number = 10
#                     assert outlines[line_number].find(line) > -1
#                 elif folder == 'analysis2':
#                     line = 'H  -0.35133100   0.66333100  -0.79133100'
#                     line_number = 10
#                     assert outlines[line_number].find(line) > -1

#             elif file == 'MeOH_SCF_error.out':
#                 line = 'H  -1.04798700   0.80281000  -0.68030200'
#                 line_number = 10
#                 assert outlines[line_number].find(line) > -1

#             elif file == 'MeOH_Unfinished.OUT':
#                 line = 'H  -1.04779100   0.87481300  -0.58663200'
#                 line_number = 10
#                 assert outlines[line_number].find(line) > -1

#         else:
#             df_QCORR,dat_files = analysis(path_main, cmd_aqme, folder, file)

#         if file == 'csv':
#             assert df_QCORR['Total files'][0] == 9
#             assert df_QCORR['Normal termination'][0] == 2
#             assert df_QCORR['Imaginary frequencies'][0] == 1
#             assert df_QCORR['SCF error'][0] == 1
#             assert df_QCORR['Basis set error'][0] == 2
#             assert df_QCORR['Other errors'][0] == 1
#             assert df_QCORR['Unfinished'][0] == 1
#             assert df_QCORR['Duplicates'][0] == 1

#         elif file == 'dat':
#             assert 'aqme-QCORR-run_1_output.dat' in dat_files

#         # check that the files are in their corresponding folders (just once)
#         if file == 'CH4.log':
#             # all the files are moved from the original directory
#             os.chdir(path_main+'/'+folder)
#             assert len(glob.glob('*.LOG')) == 0
#             assert len(glob.glob('*.log')) == 0
#             assert len(glob.glob('*.out')) == 0
#             assert len(glob.glob('*.OUT')) == 0
#             # new paths for the files
#             os.chdir(path_main+'/'+folder+'/duplicates/run_1')
#             assert len(glob.glob('*.LOG')) == 1
#             os.chdir(path_main+'/'+folder+'/success/output_files')
#             assert len(glob.glob('*.log')) == 2
#             os.chdir(path_main+'/'+folder+'/failed/run_1/unfinished')
#             assert len(glob.glob('*.OUT')) == 1
#             os.chdir(path_main+'/'+folder+'/failed/run_1/imag_freq')
#             assert len(glob.glob('*.log')) == 1
#             os.chdir(path_main+'/'+folder+'/failed/run_1/error/basis_set_error')
#             assert len(glob.glob('*.LOG')) == 2
#             os.chdir(path_main+'/'+folder+'/failed/run_1/error/scf_error')
#             assert len(glob.glob('*.out')) == 1
#             os.chdir(path_main+'/'+folder+'/failed/run_1/error/unknown_error')
#             assert len(glob.glob('*.LOG')) == 1

#     elif type_of_job == 'single_point':
#         # test there are only 2 functional/basis set combinations (instead of 4)
#         os.chdir(path_main+'/analysis/success/G16-SP_input_files/')
#         assert len(glob.glob('*')) == 2
#         assert 'b3lyp-321g' in glob.glob('*')
#         assert 'wb97xd-def2svp' in glob.glob('*')

#         # test the input files are generate correctly
#         LoTs_SP = ['b3lyp-321g','wb97xd-def2svp']
#         for LoT_SP in LoTs_SP:
#             os.chdir(path_main+'/analysis/success/G16-SP_input_files/'+LoT_SP)
#             files_sp = glob.glob('*.com')
#             # only normal terminations generate COM files for SP
#             assert len(files_sp) == 2
#             # test for the suffix option
#             assert 'CH4_SPC.com' in files_sp

#             # test the content of the input files
#             outlines = single_point(path_main, 'analysis', file)

#             if LoT_SP == 'b3lyp-321g':
#                 line = '# b3lyp/genecp nmr = giao'
#             else:
#                 line = '# wb97xd/genecp nmr = giao'
#             line_number = 2
#             assert outlines[line_number].find(line) > -1

#             line2 = '5 3'
#             line_number2 = 6
#             assert outlines[line_number2].find(line2) > -1

#             line3 = 'H   0.63133100  -0.63133100  -0.63133100'
#             line_number3 = 11
#             assert outlines[line_number3].find(line3) > -1

#             line4 = 'C 0'
#             line_number4 = 13
#             assert outlines[line_number4].find(line4) > -1

#             if LoT_SP == 'b3lyp-321g':
#                 line5 = '321g'
#             else:
#                 line5 = 'def2svp'
#             line_number5 = 14
#             assert outlines[line_number5].find(line5) > -1

#             line6 = 'H 0'
#             line_number6 = 16
#             assert outlines[line_number6].find(line6) > -1

#             if LoT_SP == 'b3lyp-321g':
#                 line7 = 'LANL2TZ'
#             else:
#                 line7 = 'LANL2DZ'
#             line_number7 = 17
#             assert outlines[line_number7].find(line7) > -1

#             line8 = 'H 0'
#             line_number8 = 20
#             assert outlines[line_number8].find(line8) > -1

#             if LoT_SP == 'b3lyp-321g':
#                 line9 = 'LANL2TZ'
#             else:
#                 line9 = 'LANL2DZ'
#             line_number9 = 21
#             assert outlines[line_number9].find(line9) > -1

#             # extra test for final line in SP
#             line10 = '-1'
#             line_number10 = 23
#             assert outlines[line_number10].find(line10) > -1
