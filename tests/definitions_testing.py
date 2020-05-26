#!/usr/bin/env python

######################################################.
# 		  This file stores all the functions 	     #
# 	   	      used in the pytest analysis	    	 #
######################################################.

import os
import glob
import pandas as pd
import subprocess

def calc_energy(file):
    energy = []
    f = open(file,"r")
    readlines = f.readlines()
    for i,line in enumerate(readlines):
        if readlines[i].find('>  <Energy>') > -1:
            energy.append(float(readlines[i+1].split()[0]))
    f.close()

    return energy

def calc_genecp(file, atom):
    # count the amount of times Pd 0 is repeated: for gen only 1, for gen_ecp 2
    count,NBO,pop,opt = 0,0,0,0
    f = open(file,"r")
    readlines = f.readlines()
    for i, line in enumerate(readlines):
        if line.find('pop') > -1:
            pop += 1
        if line.find('opt') > -1:
            opt += 1
        if line.find(atom+' 0') > -1:
            count += 1
        if line.find('$NBO $END') > -1:
            NBO += 1
    f.close()

    return count,NBO,pop,opt

def conf_gen(path, precision, cmd_pyconfort, folder, smiles, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1):
    # open right folder and run the code
    os.chdir(path+'/'+folder+'/'+smiles.split('.')[0])
    subprocess.run(cmd_pyconfort)

    # Retrieving the generated CSV file
    df_output = pd.read_csv(smiles.split('.')[0]+'-Duplicates Data.csv')

    # tests for RDKit
    if not dihedral:
        test_init_rdkit_confs = df_output['RDKIT-Initial-samples']
        test_prefilter_rdkit_confs = df_output['RDKit-initial_energy_threshold']
        test_filter_rdkit_confs = df_output['RDKit-RMSD-and-energy-duplicates']

        assert str(n_confs) == str(test_init_rdkit_confs[0])
        assert str(prefilter_confs_rdkit) == str(test_prefilter_rdkit_confs[0])
        assert str(filter_confs_rdkit) == str(test_filter_rdkit_confs[0])

    else:
        test_init_rdkit_confs = df_output['RDKIT-Rotated-conformers']
        test_unique_confs = df_output['RDKIT-Rotated-Unique-conformers']

        # I use the filter_confs_rdkit variable to assert for unique confs in dihedral scan
        assert str(n_confs) == str(test_init_rdkit_confs[0])
        assert str(filter_confs_rdkit) == str(test_unique_confs[0])

    # read the energies of the conformers
    os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/rdkit_generated_sdf_files')

    if not dihedral:
        test_rdkit_E_confs = calc_energy(smiles.split('.')[0]+'_rdkit.sdf')
    else:
        test_rdkit_E_confs = calc_energy(smiles.split('.')[0]+'_rdkit_rotated.sdf')

    # test for energies
    try:
        test_round_confs = [round(num, precision) for num in test_rdkit_E_confs]
        round_confs = [round(num, precision) for num in E_confs]
    except:
        test_round_confs = 'nan'
        round_confs = 'nan'

    assert str(round_confs) == str(test_round_confs)

    # tests charge
    test_charge = df_output['Overall charge']

    assert str(charge) == str(test_charge[0])


def only_check(path, precision, cmd_pyconfort, folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1):
    # open right folder and run the code
    os.chdir(path+'/'+folder+'/'+smiles.split('.')[0])
    subprocess.run(cmd_pyconfort)
    if params_file.find('_genecp_') > -1:
        os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/generated_gaussian_files/wb97xd-def2svp')
        file = glob.glob('*.com')[0]
        count,NBO,pop,opt = calc_genecp(file, 'Pd')

        if params_file == 'params_genecp_test1.yaml': # for gen
            assert count == 1
        else: # for genecp
            assert count == 2

def find_coordinates(file,coordinates):
    coordinates_found = 0
    com_file = file.split('.')[0]+'.com'
	outfile = open(com_file,"r")
	outlines = outfile.readlines()
	for i,outline in enumerate(outlines):
		if outline.find(coordinates) > -1:
            coordinates_found = 1
            break
    return coordinates_found

def check_log_files(path, folder, file):
    if file == 'CH4_Normal_termination.log':
        os.chdir(path+'/'+folder+'/finished')
        assert file in glob.glob('*.*')
    elif file == 'Basis_set_error1.LOG' or file == 'Basis_set_error2.LOG':
        os.chdir(path+'/'+folder+'/failed_error/atomic_basis_error')
        assert file in glob.glob('*.*')
    elif file == 'MeOH_Error_termination.LOG':
        os.chdir(path+'/'+folder+'/failed_error/unknown_error')
        assert file in glob.glob('*.*')
    elif file == 'Imag_freq.log':
        os.chdir(path+'/'+folder+'/imaginary_frequencies')
        assert file in glob.glob('*.*')
    elif file == 'MeOH_SCF_error.LOG':
        os.chdir(path+'/'+folder+'/failed_error/SCF_error')
        assert file in glob.glob('*.*')
    elif file == 'MeOH_Unfinished.LOG':
        os.chdir(path+'/'+folder+'/failed_unfinished')
        assert file in glob.glob('*.*')

def check_com_files(path, folder, file):
    if file == 'Basis_set_error1.LOG' or file == 'Basis_set_error2.LOG':
        os.chdir(path+'/'+folder+'/new_gaussian_input_files/def2svp-wb97xd')
        assert file not in glob.glob('*.*')
    elif file == 'MeOH_Error_termination.LOG':
        coordinates = 'H  -1.14928800  -0.80105100  -0.00024300'
        coordinates_error_found = find_coordinates(file,coordinates)
        assert coordinates_error_found == 1
        com_input_line_error = '# wb97xd/def2svp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=100) scrf=(SMD,solvent=Chloroform)'
        input_found_error = find_coordinates(file,com_input_line_error)
        assert input_found_error == 1
    elif file == 'Imag_freq.log':
        coordinates = 'H  -0.56133100   0.63933100  -0.67133100'
        coordinates_imag_found = find_coordinates(file,coordinates)
        assert coordinates_imag_found == 1
        com_input_line_imag = '# wb97xd/def2svp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=100) scrf=(SMD,solvent=Chloroform)'
        input_found_imag = find_coordinates(file,com_input_line_imag)
        assert input_found_imag == 1
    elif file == 'MeOH_SCF_error.LOG':
        coordinates = 'H  -1.04798700   0.80281000  -0.68030200'
        coordinates_scf_found = find_coordinates(file,coordinates)
        assert coordinates_scf_found == 1
        com_input_line_scf = '# wb97xd/def2svp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=100) scrf=(SMD,solvent=Chloroform) scf=qc'
        input_found_scf = find_coordinates(file,com_input_line_scf)
        assert input_found_scf == 1
    elif file == 'MeOH_Unfinished.LOG':
        coordinates = 'H  -1.04779100   0.87481300  -0.58663200'
        coordinates_unfinished_found = find_coordinates(file,coordinates)
        assert coordinates_unfinished_found == 1
        com_input_line_unfinished = '# wb97xd/def2svp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=100) scrf=(SMD,solvent=Chloroform)'
        input_found_unfinished = find_coordinates(file,com_input_line_unfinished)
        assert input_found_unfinished == 1

def analysis(path, cmd_pyconfort, folder, file):
    os.chdir(path+'/'+folder)
    # the code will move the files the first time, this 'if' avoids errors
    files = glob.glob('*.*')
    if len(files) > 0:
        subprocess.run(cmd_pyconfort)
    # make sure the LOG files are in the right folders after analysis
    check_log_files(path, folder, file)
    # make sure the generated COM files have the right level of theory and geometries
    os.chdir(path+'/'+folder+'/new_gaussian_input_files/def2svp-wb97xd')
    check_com_files(path, folder, file)

def single_point(path, cmd_pyconfort, folder, file):
    os.chdir(path+'/'+folder)
    files = glob.glob('*.*')
    if len(files) > 0:
        subprocess.run(cmd_pyconfort)
    os.chdir(path+'/'+folder+'/finished/single_point_input_files/wb97xd-def2svp')
    assert len(glob.glob('*.*')) == 2

    if file == 'Pd_SP.LOG':
        count,NBO,pop,opt = calc_genecp(file.split('.')[0]+'.com', 'Pd')
        assert count == 2 # finds genecp for Pd
        assert NBO == 1 # finds final line for sp
        assert pop == 1 # finds input line for sp
        assert opt == 0 # it does not find standard opt option

    elif file == 'CH4_freq.log':
        count,NBO,pop,opt = calc_genecp(file.split('.')[0]+'.com', 'C H')
        assert count == 0 # does not find genecp part
        assert NBO == 1 # finds final line for sp
        assert pop == 1 # finds input line for sp
        assert opt == 0 # it does not find standard opt option
