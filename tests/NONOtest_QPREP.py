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
@pytest.mark.parametrize("init_folder, target_folder, command_request, restore_folder",
[
    # QCORR analysis tests with standard options
    ('json_files', 'com_files', None, False), # test successful termination
    ('log_files', 'com_files', None, False), # test successful termination
    ('sdf_files', 'com_files', None, False), # test successful termination
    # INCLUDE test for XYZ to COM
    (None, None, None, True), # reset the initial folder to start another set of tests
    # with destination (destination=CC)
    # from YAML file (varfile=XX)
    # using ORCA (program='orca')
    # changing charge and mult (charge=3, mult=2)
    # adding NBO final line (qm_end=XX)
    # adding suffix (suffix=XX)
    # adding chk (chk=XX)
    # changing mem (mem='100GB')
    # changing nprocs (nprocs=48)
    # using genecp (gen_atoms=['C'], bs_gen='LANL2TZ', bs='LANL2DZ')
])

def test_QCORR_analysis(init_folder, target_folder, command_request, restore_folder):
    # copy the test folders
    if not path.exists(f'{path_main}/Example_workflows_original'):
        shutil.copytree(f'{path_main}/Example_workflows', f'{path_main}/Example_workflows_original')

    # runs the program with the different tests
    w_dir_main=f'{path_qcorr}/{init_folder}'
    destination=f'{path_qcorr}/{init_folder}/{target_folder}'

    if command_request is None:
        if init_folder == 'json_files':
            files = '*.json'
            files_assert = ['CH4.com', 'MeOH_NMR.com']
        elif init_folder == 'log_files':
            files = '*.log'
            files_assert = ['CH4.com', 'H_freq.com']
        elif init_folder == 'sdf_files':
            files = '*.sdf'
            files_assert = ['quinine_rdkit_conf_1.com', 'quinine_rdkit_conf_10.com']

        cmd_aqme = ['python', '-m', 'aqme', '--qprep', '--w_dir_main', w_dir_main, '--destination', destination,
        '--files', files, '--program', 'gaussian', '--qm_input', 'wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)']

        subprocess.run(cmd_aqme)  

        for file in files_assert:    
            assert path.exists(f'{destination}/{file.split(".")[0]}.com')

            outfile = open(f'{destination}/{file.split(".")[0]}.com', "r")
            outlines = outfile.readlines()
            outfile.close()

            line_2 = '# wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)'

            if file.split('.')[0] != 'H_freq':
                # ensure that QCORR applies the correct structural distortions to the errored calcs
                line_6 = '0 1'

                if file.split('.')[0] == 'CH4':               
                    line_8 = 'H   0.45703600  -0.46392100  -0.87651000'
                    line_10 = 'H   0.29552000   1.04984200   0.05078900'

                if file.split('.')[0] == 'MeOH_NMR':
                    line_8 = 'H  -1.085929   0.983797  -0.00000100'
                    line_10 = 'H  -1.02538  -0.54454  -0.88941'

                if file.split('.')[0] == 'quinine_rdkit_conf_1':
                    line_8 = 'O   2.93580000   2.55850000   2.17990000'
                    line_10 = 'C   2.26850000   1.23230000   0.38520000'

                if file.split('.')[0] == 'quinine_rdkit_conf_10':
                    line_8 = 'O   3.23820000   2.67370000  -1.56720000'
                    line_10 = 'C   2.35510000   0.65880000  -0.80190000'
                    assert len(glob.glob(f'{destination}/*.com')) == 10

                assert outlines[8].strip() == line_8
                assert outlines[10].strip() == line_10

            else:
                line_6 = '0 2'       
                line_7 = 'H   0.00000000   0.00000000   0.00000000'

                assert outlines[7].strip() == line_7

            assert outlines[2].strip() == line_2
            assert outlines[6].strip() == line_6


    # leave the folders as they were initially to run a different batch of tests
    if restore_folder:
        os.chdir(path_main)
        shutil.rmtree(f'{path_main}/Example_workflows')
        filepath = Path(f'{path_main}/Example_workflows_original')
        filepath.rename(f'{path_main}/Example_workflows')

