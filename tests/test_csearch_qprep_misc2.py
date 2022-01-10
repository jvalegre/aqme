#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#              Miscellaneous tests 2                 #
######################################################.

import os
import subprocess
import glob
import pandas as pd
import pytest
from definitions_testing import calc_energy,get_not_empty_files,com_lines,remove_data,calc_genecp

# saves the working directory
path_misc2 = os.getcwd()
# decimal digits for comparing E
precision_misc2 = 2

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, final_confs_misc2, E_confs, folders_in_CSEARCH, job_type",
[
    # tests with the options of Fullmonte
    ('Organic_molecules', 'pentane_fullm.smi', 'params_test_fullm1.yaml', 4, [-5.27175,-4.44184,-3.84858,-1.57172], 2, 'Fullmonte'), # checks Fullmonte energies to original RDKit energies
    ('Organic_molecules', 'octane_fullm.smi', 'params_test_fullm2.yaml', 31, [-5.88586, -5.06191, -5.01297, -5.01041, -4.46582, -4.45445, -4.26456, -4.23993, -4.21557, -4.20096, -3.90739, -3.89512, -3.74114, -3.70362, -3.65763, -3.51251, -3.48719, -3.30641, -3.30554, -3.2525, -3.20074, -3.082, -2.14724, -2.14109, -1.46293, -1.41364, -1.40318, -1.35993, -1.34075, -1.27378, -1.12246], 2, 'Fullmonte'), # ref to test options for Fullmonte
    ('Organic_molecules', 'octane_fullm.smi', 'params_test_fullm3.yaml', 22, [-5.88586, -5.06191, -5.01297, -5.01041, -4.46582, -4.45445, -4.26456, -4.23993, -4.21557, -4.20096, -4.18556, -4.04982, -3.90739, -3.89512, -3.74114, -3.70362, -3.65763, -3.51251, -3.30641, -3.20074, -2.13838, -1.35259], 2, 'Fullmonte'), # test ewin_sample_fullmonte
# MISSING - INDEX ERROR    ('Organic_molecules', 'octane_fullm.smi', 'params_test_fullm4.yaml', 22, [-5.27175,-4.44184,-3.84858,-1.57172], 'Fullmonte'), # test ewin_fullmonte
    ('Organic_molecules', 'octane_fullm.smi', 'params_test_fullm5.yaml', 3, [-3.30641, -0.52597, -0.38729], 2, 'Fullmonte'), # test nsteps_fullmonte
    ('Organic_molecules', 'octane_fullm.smi', 'params_test_fullm6.yaml', 17, [-4.26456, -4.23993, -4.23584, -4.04982, -3.90739, -3.65763, -3.51251, -3.30641, -3.20074, -3.082, -2.95114, -1.41364, -0.82222, -0.61234, -0.52597, -0.38729, -0.20171], 2, 'Fullmonte'), # test nrot_fullmonte
    ('Organic_molecules', 'octane_fullm.smi', 'params_test_fullm7.yaml',21, [-3.89512, -3.51251, -3.30641, -3.30554, -3.2525, -3.20074, -3.082, -2.95114, -2.71742, -1.40318, -1.3888, -1.36318, -0.89859, -0.89236, -0.86685, -0.72651, -0.61234, -0.3618, -0.25553, 0.10326, 0.27629], 2, 'Fullmonte'), # test ang_fullmonte
    # test if fullmonte works for metals with and without templates
    ('Metal_complexes', 'Ir_fullm.smi', 'params_test_metal_fullm1.yaml', 3, [5.75593, 5.75822, 7.80888], 2, 'Fullmonte'), # test using metals
    ('Metal_complexes', 'Pd_squareplanar_fullm.smi', 'params_test_metal_fullm2.yaml', 2, [10.58415, 10.92013], 2, 'Fullmonte'), # test using templates
    # test to check if the program removes the rdkit folder from the CSEARCH folder when using SUMM
    ('Organic_molecules', 'pentane.smi', 'params_test15.yaml', 'nan', 'nan', 2, 'SUMM'), # test SUMM folders
    # check if the program removes the I atom when using metals and templates
    ('Metal_complexes', 'Ir_Iatom.smi', 'params_test_Iatom1.yaml', 'nan', 'nan', 'nan', 'check_I_atom_rdkit'), # test with RDKit
    ('Metal_complexes', 'Ir_Iatom.smi', 'params_test_Iatom2.yaml', 'nan', 'nan', 'nan', 'check_I_atom_fullmonte'), # test with fullmonte
    ('Metal_complexes', 'Ir_Iatom.smi', 'params_test_Iatom3.yaml', 'nan', 'nan', 'nan', 'check_I_atom_summ'), # test with SUMM
    ('Metal_complexes', 'Pd_squareplanar_Iatom.smi', 'params_test_Iatom4.yaml', 'nan', 'nan', 'nan', 'check_I_atom_rdkit'), # test with RDKit
    ('Metal_complexes', 'Pd_squareplanar_Iatom.smi', 'params_test_Iatom5.yaml', 'nan', 'nan', 'nan', 'check_I_atom_fullmonte'), # test with fullmonte
    ('Metal_complexes', 'Pd_squareplanar_Iatom.smi', 'params_test_Iatom6.yaml', 'nan', 'nan', 'nan', 'check_I_atom_summ') # test with SUMM
])

def test_confgen_misc2(folder, smiles, params_file, final_confs_misc2, E_confs, folders_in_CSEARCH, job_type):
    # runs the program with the different tests
    cmd_misc2 = ['python', '-m', 'aqme', '--varfile', params_file]

    os.chdir(path_misc2+'/'+folder+'/'+smiles.split('.')[0])
    subprocess.call(cmd_misc2)

    if job_type == 'Fullmonte':

        # Retrieving the generated CSV file
        os.chdir(path_misc2+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/csv_files')
        df_output = pd.read_csv(smiles.split('.')[0]+'-CSEARCH-Data.csv')

        # check n of final conformers from the CSV file
        assert str(df_output['FullMonte-Unique-conformers'][0]) == str(final_confs_misc2)

        # check energies from the SDF file
        os.chdir(path_misc2+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/fullmonte')
        if smiles == 'Pd_squareplanar_fullm.smi':
            test_E_confs = calc_energy(smiles.split('.')[0]+'_1_fullmonte.sdf')
            # this test ensures that the numbering starts in 1 not in 0 when using metal templates
            assert smiles.split('.')[0]+'_0_fullmonte.sdf' not in glob.glob('*.sdf')
        else:
            test_E_confs = calc_energy(smiles.split('.')[0]+'_fullmonte.sdf')

        test_E_round_confs = [round(num, precision_misc2) for num in test_E_confs]
        E_round_confs = [round(num, precision_misc2) for num in E_confs]
        assert str(test_E_round_confs) == str(E_round_confs)

        # check number of com files
        com_files_folder = path_misc2+'/'+folder+'/'+smiles.split('.')[0]+'/QMCALC/G16/wb97xd-6-31g(d)'
        test_com_files = 0
        test_com_files = get_not_empty_files(com_files_folder,test_com_files,'com')
        if smiles == 'Pd_squareplanar_fullm.smi':
            assert str(test_com_files) == '6'
        else:
            assert str(test_com_files) == str(final_confs_misc2)

        # extra check to ensure that all the CSEARCH methods start from number 0
        assert smiles.split('.')[0]+'_0.com' not in glob.glob('*.com')

        # check opt, charge and multiplicity
        if smiles == 'Pd_squareplanar_fullm.smi':
            _,_,_,opt,charge_com,multiplicity_com = calc_genecp(smiles.split('.')[0]+'_1_1.com', [])
        else:
            _,_,_,opt,charge_com,multiplicity_com = calc_genecp(smiles.split('.')[0]+'_1.com', [])

        if smiles == 'Ir_fullm.smi':
            assert str(charge_com) == '1'
        else:
            assert str(charge_com) == '0'
        assert str(multiplicity_com) == '1'
        assert str(opt) == '1'

    if job_type == 'Fullmonte' or job_type == 'SUMM':
        # checking amount of folders (ideally, there should be only 2 folders, one for the CSEARCH
        # method and another for the CSV files)
        os.chdir(path_misc2+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH')
        assert str((len(glob.glob('*')))) == str(folders_in_CSEARCH)

        if job_type == 'SUMM':
            # extra check to ensure that all the CSEARCH methods start from number 1
            os.chdir(path_misc2+'/'+folder+'/'+smiles.split('.')[0]+'/QMCALC/G16/wb97xd-6-31g(d)')
            starts_in_0 = False
            for file in glob.glob('*.com'):
                if file.split('.')[0].split('_')[1] == '0' or file.split('.')[0].split('_')[2] == '0':
                    starts_in_0 = True
            assert starts_in_0 == False

    if job_type.split('_')[0] == 'check':
        # ensures that the metal is back after CSEARCH
        os.chdir(path_misc2+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/'+job_type.split('_')[3])
        method = job_type.split('_')[3]
        if smiles == 'Pd_squareplanar_Iatom.smi':
            sdf_file = smiles.split('.')[0]+'_1_'+method+'.sdf'
        else:
            sdf_file = smiles.split('.')[0]+'_'+method+'.sdf'


        outfile = open(sdf_file,"r")
        outlines = outfile.readlines()

        atom_to_find = smiles.split('.')[0].split('_')[0]
        if atom_to_find == 'Pd':
            line_number = 5
        elif atom_to_find == 'Ir':
            line_number = 10

        assert outlines[line_number].find(atom_to_find) > -1

        outfile.close()

        # ensures that the metal is back after QPREP
        os.chdir(path_misc2+'/'+folder+'/'+smiles.split('.')[0]+'/QMCALC/G16/wb97xd-6-31g(d)')

        com_name = glob.glob('*.com')[0]
        outlines2 = com_lines(com_name)

        if atom_to_find == 'Pd':
            line_number2 = 8
        elif atom_to_find == 'Ir':
            line_number2 = 13

        assert outlines2[line_number2].find(atom_to_find) > -1

    remove_data(path_misc2, folder, smiles)
