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
import pandas as pd

# saves the working directory
path_main = os.getcwd()
path_qcorr = os.getcwd()+'/Example_workflows/QCORR_processing_QM_outputs'

# QCORR tests
@pytest.mark.parametrize("init_folder, file, command_line, target_folder, restore_folder",
[
    # QCORR analysis tests with standard options
    ('QCORR_1', 'CH4.log', 'run_QCORR', 'successful_QM_outputs', False), # test successful termination
    ('QCORR_1', 'MeOH_G09.log', None, 'successful_QM_outputs', False), # test successful termination
    ('QCORR_1', 'z_CH4_duplicate.log', None, 'unsuccessful_QM_outputs/run_1/duplicates', False), # test duplicates
    ('QCORR_1', 'Basis_set_error1.log', None, 'unsuccessful_QM_outputs/run_1/error/basis_set_error', False), # test incompatibilities with basis sets
    ('QCORR_1', 'Basis_set_error2.log', None, 'unsuccessful_QM_outputs/run_1/error/basis_set_error', False), # test incompatibilities with basis sets
    ('QCORR_1', 'MeOH_G09_FAIL.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error', False), # test error terminations
    ('QCORR_1', 'CH2OH2_unfinished.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error', False), # test unfinished calculations
    ('QCORR_1', 'Imag_freq.log', None, 'unsuccessful_QM_outputs/run_1/extra_imag_freq', False), # test imaginary frequencies
    ('QCORR_1', 'imag_freq_no_opt.log', None, 'unsuccessful_QM_outputs/run_1/extra_imag_freq', False), # test imaginary frequencies without OPT
    ('QCORR_1', 'MeOH_SCF_error.log', None, 'unsuccessful_QM_outputs/run_1/error/scf_error', False), # test SCF errors
    ('QCORR_1', 'CH4_before_E.log', None, 'unsuccessful_QM_outputs/run_1/error/no_data', False), # test calcs that finish before any coords are printed
    ('QCORR_1', 'freq_conv_YYNN.log', None, 'unsuccessful_QM_outputs/run_1/freq_no_conv', False), # test calcs with freq calcs that did not converge after OPT
    ('QCORR_1', 'freq_ok_YYNN.log', None, 'successful_QM_outputs', False), # test calcs with freq calcs that did not converge after OPT
    ('QCORR_1', 'TS_CH3HCH3_no_conv_freq.log', None, 'unsuccessful_QM_outputs/run_1/freq_no_conv', False), # test calcs with freq calcs that did not converge after OPT in TSs
    ('QCORR_1', 'bpinene_spin_contamin.log', None, 'unsuccessful_QM_outputs/run_1/spin_contaminated', False), # test calcs with spin contamination
    ('QCORR_1', 'CH4_T1_SP_spin_contamin.log', None, 'unsuccessful_QM_outputs/run_1/spin_contaminated', False), # test calcs with spin contamination
    ('QCORR_1', 'CH4_Fail_freq_only.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error', False), # test for Normal terminated OPT and unfinished freq
    ('QCORR_1', 'TS_CH3HCH3.log', None, 'successful_QM_outputs', False), # test successful termination in TSs
    ('QCORR_1', 'TS_CH3HCH3_unfinished.log', None, 'unsuccessful_QM_outputs/run_1/error/not_specified_error', False), # test unfinished TSs
    ('QCORR_1', 'TS_CH3HCH3_imag_freq.log', None, 'unsuccessful_QM_outputs/run_1/extra_imag_freq', False), # test imaginary frequencies for TSs
    ('QCORR_1', 'TS_CH3HCH3_no_imag_freq.log', None, 'unsuccessful_QM_outputs/run_1/ts_no_imag_freq', False), # test imaginary frequencies for TSs
    ('QCORR_1', 'CH4_SP.log', None, 'successful_QM_outputs/SP_calcs', False), # test for single-point calcs
    ('QCORR_1', 'MeOH_NMR.log', None, 'successful_QM_outputs/SP_calcs', False), # test for single-point calcs 
    ('QCORR_1', 'H_freq.log', None, 'successful_QM_outputs', False), # test successful termination with 1 atom
    ('QCORR_1', 'H_SP.log', None, 'successful_QM_outputs/SP_calcs', False), # test for single-point calcs with 1 atom
    ('QCORR_1', 'CO2_linear_3freqs_FAIL.log', None, 'unsuccessful_QM_outputs/run_1/linear_mol_wrong', False), # test for linear mols with wrong number of freqs
    ('QCORR_1', 'CO2_linear_4freqs.log', None, 'successful_QM_outputs', False), # test successful termination for linear mols
    ('QCORR_1', 'json', None, 'successful_QM_outputs/json_files', False), # test for correct creation of json files (all the successful terminations only)
    ('QCORR_1', 'fullcheck', None, 'successful_QM_outputs/json_files', False), # test for correct fullcheck option
    ('QCORR_1', 'csv', None, None, False), # test final csv file with results
    ('QCORR_1', 'dat', None, None, False), # test final dat file with results
    ('QCORR_1', 'check_init', None, None, False), # test that the folder for initial QM inputs is not created
    # Test isomerization filter
    ('QCORR_2', 'CH4.log', 'run_QCORR', 'successful_QM_outputs', False), # test successful termination when using the isomerization filter
    ('QCORR_2', 'CH2OH2_isomerized.log', None, 'unsuccessful_QM_outputs/run_1/isomerization', False), # test isomerized calcs
    ('QCORR_2', 'check_fixed', None, None, False), # test that the fixed_QM_input folder is not generated when isomerized
    # Test if the unsuccessful runs are created in the same folder (files are contained inside the unsuccessful run_1 folder)
    ('QCORR_5', 'CH4.log', 'run_QCORR', 'successful_QM_outputs', False), # test that the fixed_QM_input folder is not generated when isomerized
    ('QCORR_5', 'CH2OH2_isomerized.log', None, 'unsuccessful_QM_outputs', False), # test that the fixed_QM_input folder is not generated when isomerized
    (None, None, None, None, True), # reset the initial folder to start another set of tests
    # Test if QCORR works using parameters from a YAML file
    ('QCORR_2b', 'CH4.log', 'run_QCORR', 'successful_QM_outputs', False), # test successful termination when using the isomerization filter
    ('QCORR_2b', 'CH2OH2_isomerized.log', None, 'unsuccessful_QM_outputs/run_1/isomerization', False), # test isomerized calcs
    # Test if the unsuccessful runs are created in the same folder (files are contained inside the parent QCORR folder)
    ('QCORR_5b', 'CH4.log', 'run_QCORR', 'successful_QM_outputs', False), # test that the fixed_QM_input folder is not generated when isomerized
    ('QCORR_5b', 'CH2OH2_isomerized.log', 'run_QCORR', 'unsuccessful_QM_outputs', False), # test that the fixed_QM_input folder is not generated when isomerized
    (None, None, None, None, True), # reset the initial folder to start another set of tests
    # Test if the vdwfrac and covfrac options work in QCORR
    ('QCORR_2c', 'CH4.log', 'run_QCORR', 'successful_QM_outputs', False), # test successful termination when using the isomerization filter
    ('QCORR_2c', 'CH2OH2_isomerized.log', None, 'successful_QM_outputs', False), # test that the vdwfrac and covfrac are disabled now with very low values
    # Test other options from QCORR (deactivate the freq_conv filter, S**2 threshold, energy threshold for duplicates, imag freqs cut-offs)
    ('QCORR_1b', 'freq_conv_YYNN.log', None, 'successful_QM_outputs', False), # test to disable freq calcs that did not converge after OPT
    ('QCORR_1b', 'bpinene_spin_contamin.log', None, 'successful_QM_outputs', False), # test to change filter for spin contamination
    ('QCORR_1b', 'z_CH4_duplicate.log', None, 'successful_QM_outputs', False), # test to change energy threshold to consider duplicates
    ('QCORR_1b', 'imag_freq_no_opt.log', None, 'successful_QM_outputs', False), # test to change cut-off to consider imaginary frequencies
    ('QCORR_1b', 'Imag_freq.log', None, 'unsuccessful_QM_outputs/run_1/extra_imag_freq', False), # test to change amplitude for displacing imaginary frequencies
    (None, None, None, None, True), # reset the initial folder to start another set of tests
    # add genECP test
    # isomerization with csv (ongoing)
    # isomeriz with csv for TSs (ongoing)
    # tell if the imag freq from a TS is right based on displacement (ongoing)
])

def test_QCORR_analysis(init_folder, file, command_line, target_folder, restore_folder):
    # copy the test folders
    if not path.exists(f'{path_main}/Example_workflows_original'):
        shutil.copytree(f'{path_main}/Example_workflows', f'{path_main}/Example_workflows_original')

    # runs the program with the different tests
    w_dir_main=f'{path_qcorr}/{init_folder}'
    cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_main, '--files', '*.log', '--freq_conv', 'opt=(calcfc,maxstep=5)']

    if init_folder == 'QCORR_1':
        if command_line is not None:
            subprocess.run(cmd_aqme)   

        if file.split('.')[-1].lower() == 'log':
            # ensure the output file moves to the right folder
            assert path.exists(f'{w_dir_main}/{target_folder}/{file}')

            # ensure that the com files are generated correctly
            if file.split('.')[0] in ['CH4', 'MeOH_G09', 'z_CH4_duplicate', 'Basis_set_error1', 
                                      'Basis_set_error2', 'CH4_before_E', 'bpinene_spin_contamin',
                                      'TS_CH3HCH3', 'TS_CH3HCH3_no_imag_freq', 'CH4_SP', 'H_freq',
                                      'H_SP', 'MeOH_NMR', 'CO2_linear_4freqs', 'freq_ok_YYNN',
                                      'CH4_T1_SP_spin_contamin']:
                assert not path.exists(f'{w_dir_main}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/{file.split(".")[0]}.com')

            else:
                assert path.exists(f'{w_dir_main}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/{file.split(".")[0]}.com')

                # ensure that QCORR applies the correct structural distortions to the errored calcs
                if file.split('.')[0] == 'MeOH_G09_FAIL':
                    line_2 = '# opt=calcfc freq=noraman cc-pvtz scrf=(solvent=chloroform,pcm) pbe1pbe g09defaults'
                    line_6 = '0 1'
                    line_8 = 'H   1.13330900   0.94774000  -0.00000300'
                    line_10 = 'H   0.97796200  -0.55747000   0.87365300'

                elif file.split('.')[0] == 'CH2OH2_unfinished':
                    line_2 = '# opt freq 3-21g m062x'
                    line_6 = '0 1'
                    line_8 = 'H   0.09630700   1.26666100  -0.75061200'
                    line_10 = 'O  -1.21496800  -0.19993700  -0.10923000'

                elif file.split('.')[0] == 'TS_CH3HCH3_no_conv_freq':
                    line_2 = '# opt=(ts,noeigen,calcfc,maxstep=5) freq b3lyp/3-21g'
                    line_6 = '0 2'
                    line_8 = 'H  -0.91171200   0.52637700  -1.63085200'
                    line_10 = 'H   0.00000000  -1.05275500  -1.63085200'

                elif file.split('.')[0] == 'CH4_Fail_freq_only':
                    line_2 = '# b3lyp/3-21G freq=noraman'
                    line_6 = '0 1'
                    line_8 = 'H   0.63133100   0.63133100   0.63133100'
                    line_10 = 'H  -0.63133100   0.63133100  -0.63133100'

                elif file.split('.')[0] == 'CO2_linear_3freqs_FAIL':
                    line_2 = '# opt=maxcycles=100 freq=noraman b3lyp 6-31G symmetry=(PG=Cinfv)'
                    line_6 = '0 1'
                    line_8 = 'O   0.00000000   1.18790900  -0.00032000'
                    line_10 = None

                elif file.split('.')[0] == 'TS_CH3HCH3_unfinished':
                    line_2 = '# opt=(calcfc,ts,noeigen) freq b3lyp/3-21g'
                    line_6 = '0 2'
                    line_8 = 'H  -0.90909900   0.52486800  -1.63442100'
                    line_10 = 'H   0.00000000  -1.04973700  -1.63442100'

                elif file.split('.')[0] == 'imag_freq_no_opt':
                    line_2 = '# M062X/Def2TZVP freq=noraman opt'
                    line_6 = '0 1'
                    line_8 = 'C  -0.90757400   0.00709700  -0.00994500'
                    line_10 = 'C   1.19952600   1.19528800   0.00698400'


                elif file.split('.')[0] == 'TS_CH3HCH3_imag_freq':
                    line_2 = '# opt=(calcfc,ts,noeigen) freq b3lyp/3-21g'
                    line_6 = '0 2'
                    line_8 = 'H  -0.87171200   0.59637700  -1.63085200'
                    line_10 = 'H  -0.08200000  -1.05275500  -1.63085200'

                elif file.split('.')[0] == 'Imag_freq':
                    line_2 = '# opt freq 3-21g m062x'
                    line_6 = '0 1'
                    line_8 = 'H   0.38503600  -0.39992100  -0.94851000'
                    line_10 = 'H   0.20952000   1.07184200  -0.02121100'

                elif file.split('.')[0] == 'MeOH_SCF_error':
                    line_2 = '# opt freq=noraman b3lyp/3-21g scf=xqc'
                    line_6 = '0 1'
                    line_8 = 'H  -1.08038200   0.99705900  -0.00806000'
                    line_10 = 'H  -1.06929400  -0.53149000   0.89474400'

                elif file.split('.')[0] == 'freq_conv_YYNN':
                    line_2 = '# M062X/Def2TZVP freq=noraman opt=(calcfc,maxstep=5)'
                    line_6 = '0 1'
                    line_8 = 'C  -0.90757400   0.00709700  -0.00594500'
                    line_10 = 'C   1.19952600   1.19528800   0.00098400'

                outfile = open(f'{w_dir_main}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/{file.split(".")[0]}.com', "r")
                outlines = outfile.readlines()
                outfile.close()

                assert outlines[2].strip() == line_2
                assert outlines[6].strip() == line_6
                assert outlines[8].strip() == line_8
                if line_10 is not None:
                    assert outlines[10].strip() == line_10

        elif file == 'json':
            os.chdir(f'{w_dir_main}/{target_folder}')
            json_files = glob.glob('*.json')
            target_files = ['CH4.json', 'CO2_linear_4freqs.json', 'freq_ok_YYNN.json', 'H_freq.json', 'MeOH_G09.json', 'TS_CH3HCH3.json']
            assert len(json_files) == 6
            assert sorted(json_files) == sorted(target_files)

        elif file == 'fullcheck':
            target_fullcheck = ['-- Full check analysis --\n']
            target_fullcheck.append('x  Different program used in the calculations:\n')
            target_fullcheck.append('     * Gaussian 09, Revision A.02 in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('       - CO2_linear_4freqs\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('     * Gaussian 16, Revision C.01 in:\n')
            target_fullcheck.append('       - freq_ok_YYNN\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('       - MeOH_G09\n')
            target_fullcheck.append('x  Different grid_type used in the calculations:\n')
            target_fullcheck.append('     * sg1 in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('       - CO2_linear_4freqs\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('     * ultrafine in:\n')
            target_fullcheck.append('       - freq_ok_YYNN\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('     * fine in:\n')
            target_fullcheck.append('       - MeOH_G09\n')
            target_fullcheck.append('x  Different level_of_theory used in the calculations:\n')
            target_fullcheck.append('     * M062X/3-21G in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('     * B3LYP/6-31G in:\n')
            target_fullcheck.append('       - CO2_linear_4freqs\n')
            target_fullcheck.append('     * M062X/def2TZVP in:\n')
            target_fullcheck.append('       - freq_ok_YYNN\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('     * PBE1PBE/CC-pVTZ in:\n')
            target_fullcheck.append('       - MeOH_G09\n')
            target_fullcheck.append('     * B3LYP/3-21G in:\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('o  Same dispersion (none) used in all the calculations\n')
            target_fullcheck.append('x  Different solvation used in the calculations:\n')
            target_fullcheck.append('     * gas_phase in:\n')
            target_fullcheck.append('       - CH4\n')
            target_fullcheck.append('       - CO2_linear_4freqs\n')
            target_fullcheck.append('       - freq_ok_YYNN\n')
            target_fullcheck.append('       - H_freq\n')
            target_fullcheck.append('       - TS_CH3HCH3\n')
            target_fullcheck.append('     * scrf=(solvent=chloroform,pcm) in:\n')
            target_fullcheck.append('       - MeOH_G09')

            outfile = open(f'{w_dir_main}/{target_folder}/--QCORR_Fullcheck_Analysis--.dat', "r")
            outlines = outfile.readlines()
            outfile.close()
            assert target_fullcheck == outlines

        elif file == 'dat':
            target_dat = ['\n']
            target_dat.append(['Command line used in AQME: aqme --qcorr\n'])
            target_dat.append('o  Analyzing output files in\n')
            target_dat.append('\n')
            target_dat.append('Basis_set_error1.log: Termination = other, Error type = atomicbasiserror\n')
            target_dat.append('Basis_set_error2.log: Termination = other, Error type = atomicbasiserror\n')
            target_dat.append('bpinene_spin_contamin.log: Termination = normal, Error type = spin_contaminated\n')

            outfile = open(f'{w_dir_main}/QCORR-run_1.dat', "r")
            outlines = outfile.readlines()
            outfile.close()

            for i,line in enumerate(target_dat):
                if i == 1:
                    assert outlines[i].find('Command line used in AQME: aqme --qcorr') > -1
                elif i == 2:
                    assert outlines[i].find('o  Analyzing output files in') > -1
                else:
                    assert line == outlines[i]

        elif file == 'csv':
            qcorr_stats = pd.read_csv(f'{w_dir_main}/QCORR-run_1-stats.csv')
            assert qcorr_stats['Total files'][0] == 27
            assert qcorr_stats['Normal termination'][0] == 6
            assert qcorr_stats['Single-point calcs'][0] == 3
            assert qcorr_stats['Extra imag. freq.'][0] == 3
            assert qcorr_stats['TS with no imag. freq.'][0] == 1
            assert qcorr_stats['Freq not converged'][0] == 2
            assert qcorr_stats['Linear mol with wrong n of freqs'][0] == 1
            assert qcorr_stats['SCF error'][0] == 1
            assert qcorr_stats['No data'][0] == 1
            assert qcorr_stats['Basis set error'][0] == 2
            assert qcorr_stats['Other errors'][0] == 4
            assert qcorr_stats['Spin contamination'][0] == 2
            assert qcorr_stats['Duplicates'][0] == 1
        
        elif file == 'check_init':
            assert not path.exists(f'{w_dir_main}/initial_QM_inputs/')

    elif init_folder == 'QCORR_1b':
        w_dir_main=f'{path_qcorr}/QCORR_1'
        cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_main, '--files', file, '--freq_conv', 'opt=(calcfc,maxstep=5)']
        if file.split('.')[0] == 'freq_conv_YYNN':
            cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_main, '--files', file]

        elif file.split('.')[0] == 'bpinene_spin_contamin':
            cmd_aqme = cmd_aqme + ['--s2_threshold', '50']

        elif file.split('.')[0] == 'z_CH4_duplicate':
            cmd_aqme = cmd_aqme + ['--dup_threshold', '0.000000000001']

        elif file.split('.')[0] == 'imag_freq_no_opt':
            cmd_aqme = cmd_aqme + ['--ifreq_cutoff', '50']

        elif file.split('.')[0] == 'Imag_freq':
            cmd_aqme = cmd_aqme + ['--amplitude_ifreq', '-0.4']  
        
        subprocess.run(cmd_aqme)      

        # ensure the output file moves to the right folder
        assert path.exists(f'{w_dir_main}/{target_folder}/{file}')

        if file.split('.')[0] == 'Imag_freq':
            # ensure that QCORR applies the correct structural distortions to the errored calcs
            line_2 = '# opt freq 3-21g m062x'
            line_6 = '0 1'
            line_8 = 'H   0.60103600  -0.59192100  -0.73251000'
            line_10 = 'H   0.46752000   1.00584200   0.19478900'
                
            outfile = open(f'{w_dir_main}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/{file.split(".")[0]}.com', "r")
            outlines = outfile.readlines()
            outfile.close()

            assert outlines[2].strip() == line_2
            assert outlines[6].strip() == line_6
            assert outlines[8].strip() == line_8
            assert outlines[10].strip() == line_10

    elif init_folder == 'QCORR_1c':
        w_dir_main=f'{path_qcorr}/QCORR_1'
        cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_main, '--files', file, '--fullcheck', 'False']
        subprocess.run(cmd_aqme) 

        assert not path.exists(f'{w_dir_main}/{target_folder}/--QCORR_Fullcheck_Analysis--.dat')

    elif init_folder == 'QCORR_2':
        if command_line is not None:
            cmd_aqme = cmd_aqme + ['--isom_type', 'com', '--isom_inputs', w_dir_main]
            subprocess.run(cmd_aqme)
        
        # ensure the output file moves to the right folder, including the initial COM file
        if file.split('.')[-1].lower() == 'log':
            assert path.exists(f'{w_dir_main}/{target_folder}/{file}')
            assert path.exists(f'{w_dir_main}/initial_QM_inputs/{file.split(".")[0]}.com')

        # ensure that no com files were generated
        elif file == 'check_fixed':
            assert not path.exists(f'{w_dir_main}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/')

    elif init_folder == 'QCORR_2b':
        w_dir_main=f'{path_qcorr}/QCORR_2'
        if command_line is not None:
            param_file = f'{w_dir_main}/QCORR_params.yaml'
            cmd_aqme = ['python', '-m', 'aqme', '--w_dir_main', w_dir_main, '--isom_inputs', w_dir_main, '--varfile', param_file]
            subprocess.run(cmd_aqme)
        # ensure the output file moves to the right folder, including the initial COM file
        if file.split('.')[-1].lower() == 'log':
            assert path.exists(f'{w_dir_main}/{target_folder}/{file}')
            assert path.exists(f'{w_dir_main}/initial_QM_inputs/{file.split(".")[0]}.com')

    elif init_folder == 'QCORR_2c':
        w_dir_main=f'{path_qcorr}/QCORR_2'
        if command_line is not None:
            cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_main, '--files', '*.log', '--freq_conv', 'opt=(calcfc,maxstep=5)']
            cmd_aqme = cmd_aqme + ['--isom_type', 'com', '--isom_inputs', w_dir_main, '--vdwfrac', '0.01', '--covfrac', '0.01']
            subprocess.run(cmd_aqme)
        # ensure the output file moves to the right folder, including the initial COM file
        if file.split('.')[-1].lower() == 'log':
            assert path.exists(f'{w_dir_main}/{target_folder}/{file}')

    elif init_folder == 'QCORR_5':
        w_dir_QCORR_5=f'{path_qcorr}/{init_folder}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/'

        if command_line is not None:
            cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_QCORR_5, '--files', '*.log', '--freq_conv', 'opt=(calcfc,maxstep=5)']
            cmd_aqme = cmd_aqme + ['--isom_type', 'gjf', '--isom_inputs', w_dir_QCORR_5]
            subprocess.run(cmd_aqme)

        if file.split('.')[0] == 'CH4':
            assert path.exists(f'{w_dir_main}/{target_folder}/{file}')
            assert path.exists(f'{w_dir_main}/{target_folder}/json_files/{file.split(".")[0]}.json')
        else:
            assert path.exists(f'{w_dir_main}/{target_folder}/run_2/isomerization/{file}')
            assert not path.exists(f'{w_dir_main}/{target_folder}/run_1/fixed_QM_inputs/unsuccessful_QM_outputs')

    elif init_folder == 'QCORR_5b':
        w_dir_QCORR_5b=f'{path_qcorr}/QCORR_5'
        file_QCORR_5b = Path(f'{w_dir_QCORR_5b}/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/{file}')
        file_QCORR_5b.rename(f'{w_dir_QCORR_5b}/{file}')
        
        if command_line is not None:
            cmd_aqme = ['python', '-m', 'aqme', '--qcorr', '--w_dir_main', w_dir_QCORR_5b, '--files', '*.log', '--freq_conv', 'opt=(calcfc,maxstep=5)']
            cmd_aqme = cmd_aqme + ['--isom_type', 'gjf', '--isom_inputs', w_dir_QCORR_5b+'/unsuccessful_QM_outputs/run_1/fixed_QM_inputs/']
            subprocess.run(cmd_aqme)

        if file.split('.')[0] == 'CH4':
            assert path.exists(f'{w_dir_QCORR_5b}/{target_folder}/{file}')
            assert path.exists(f'{w_dir_QCORR_5b}/{target_folder}/json_files/{file.split(".")[0]}.json')
        else:
            assert path.exists(f'{w_dir_QCORR_5b}/{target_folder}/run_2/isomerization/{file}')
            assert not path.exists(f'{w_dir_QCORR_5b}/{target_folder}/run_1/fixed_QM_inputs/unsuccessful_QM_outputs')

    # leave the folders as they were initially to run a different batch of tests
    elif restore_folder:
        os.chdir(path_main)
        shutil.rmtree(f'{path_main}/Example_workflows')
        filepath = Path(f'{path_main}/Example_workflows_original')
        filepath.rename(f'{path_main}/Example_workflows')
