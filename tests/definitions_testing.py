#!/usr/bin/env python

######################################################.
# 		  This file stores all the functions 	     #
# 	   	      used in the pytest analysis	    	 #
######################################################.

import os
import glob
import shutil
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
    count,NBO,pop,opt,charge_com,multiplicity_com = 0,0,0,0,None,None
    f = open(file,"r")
    readlines = f.readlines()
    for i, line in enumerate(readlines):
        line_parts = line.split()
        if line.find('pop') > -1:
            pop += 1
        if line.find('opt') > -1:
            opt += 1
        if line.find('# wb97xd') > -1:
            charge_com = readlines[i+4].split()[0]
            multiplicity_com = readlines[i+4].split()[1]
        elif line.find('$NBO $END') > -1:
            NBO += 1
        for _,atom_ind in enumerate(atom):
            if atom_ind in line_parts and '0' in line_parts:
                count += 1
                break
    f.close()

    return count,NBO,pop,opt,charge_com,multiplicity_com

def remove_data(path, folder, smiles):
    # remove all the data created by the job
    if not smiles:
        os.chdir(path+'/'+folder)
    else:
        os.chdir(path+'/'+folder+'/'+smiles.split('.')[0])
    all_data = glob.glob('*')
    discard_ext = ['sdf','csv','dat']
    for _,file in enumerate(all_data):
        if len(file.split('.')) == 1:
            shutil.rmtree(file, ignore_errors=True)
        elif file != 'charged.csv' or file != 'charged.sdf':
            if file.split('.')[1] in discard_ext:
                os.remove(file)

def rdkit_tests(df_output,dihedral,xTB_ANI1):
    if not dihedral:
        if not xTB_ANI1:
            test_init_rdkit_confs = df_output['RDKIT-Initial-samples']
            test_prefilter_rdkit_confs = df_output['RDKit-initial_energy_threshold']
            test_filter_rdkit_confs = df_output['RDKit-RMSD-and-energy-duplicates']
            test_unique_confs = 'nan'
        elif xTB_ANI1 == 'xTB':
            test_init_rdkit_confs = df_output['xTB-Initial-samples']
            test_prefilter_rdkit_confs = df_output['xTB-initial_energy_threshold']
            test_filter_rdkit_confs = df_output['xTB-RMSD-and-energy-duplicates']
            test_unique_confs = 'nan'
        elif xTB_ANI1 == 'AN1ccx':
            test_init_rdkit_confs = df_output['AN1ccx-Initial-samples']
            test_prefilter_rdkit_confs = df_output['AN1ccx-initial_energy_threshold']
            test_filter_rdkit_confs = df_output['AN1ccx-RMSD-and-energy-duplicates']
            test_unique_confs = 'nan'
    else:
        test_init_rdkit_confs = df_output['RDKIT-Rotated-conformers']
        test_prefilter_rdkit_confs = 'nan'
        test_filter_rdkit_confs = 'nan'
        test_unique_confs = df_output['RDKIT-Rotated-Unique-conformers']

    return test_init_rdkit_confs, test_prefilter_rdkit_confs, test_filter_rdkit_confs, test_unique_confs

def conf_gen(path, precision, cmd_pyconfort, folder, smiles, E_confs, dihedral, xTB_ANI1, metal, template):
    # open right folder and run the code
    os.chdir(path+'/'+folder+'/'+smiles.split('.')[0])
    subprocess.call(cmd_pyconfort)

    # Retrieving the generated CSV file
    df_output = pd.read_csv(smiles.split('.')[0]+'-Duplicates Data.csv')

    # tests for RDKit
    test_init_rdkit_confs, test_prefilter_rdkit_confs, test_filter_rdkit_confs, test_unique_confs = rdkit_tests(df_output,dihedral,xTB_ANI1)

    # file_smi is a variable used for finding SDF and COM files
    if folder == 'Multiple':
        file_smi = 'pentane'
    else:
        file_smi = smiles.split('.')[0]

    # read the energies of the conformers
    try:
        if not xTB_ANI1:
            os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/rdkit_generated_sdf_files')
            # this is to tests smi files with multiple smiles
            if folder == 'Multiple':
                if not dihedral:
                    test_rdkit_E_confs = calc_energy(file_smi+'_rdkit.sdf')
                else:
                    test_rdkit_E_confs = calc_energy(file_smi+'_rdkit_rotated.sdf')
            else:
                if template == 'squareplanar' or template == 'squarepyramidal':
                    if not dihedral:
                        test_rdkit_E_confs = calc_energy(file_smi+'_0_rdkit.sdf')
                    else:
                        test_rdkit_E_confs = calc_energy(file_smi+'_0_rdkit_rotated.sdf')
                else:
                    if not dihedral:
                        test_rdkit_E_confs = calc_energy(file_smi+'_rdkit.sdf')
                    else:
                        test_rdkit_E_confs = calc_energy(file_smi+'_rdkit_rotated.sdf')
        elif xTB_ANI1 == 'xTB':
            os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/xtb_minimised_generated_sdf_files')
            test_rdkit_E_confs = calc_energy(file_smi+'_xtb.sdf')
        elif xTB_ANI1 == 'AN1ccx':
            os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/ani1ccx_minimised_generated_sdf_files')
            test_rdkit_E_confs = calc_energy(file_smi+'_ani.sdf')

        # test for energies
        try:
            test_round_confs = [round(num, precision) for num in test_rdkit_E_confs]
            round_confs = [round(num, precision) for num in E_confs]

        except TypeError:
            test_round_confs = 'nan'
            round_confs = 'nan'

        # tests charge
        test_charge = df_output['Overall charge']

        os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/generated_gaussian_files/wb97xd-def2svp')

        file_gen = glob.glob(file_smi+'*.com')[0]
        if metal != False:
            count,_,_,_,charge_com,multiplicity_com = calc_genecp(file_gen, metal)
        else:
            count,_,_,_,charge_com,multiplicity_com = calc_genecp(file_gen, ['C'])

    # this is included for the tests where no confs are generated
    except (FileNotFoundError,IndexError):
        test_init_rdkit_confs, test_filter_rdkit_confs, test_prefilter_rdkit_confs, test_unique_confs = 'nan','nan','nan','nan'
        test_round_confs, round_confs, test_charge = 'nan','nan','nan'
        count,charge_com,multiplicity_com = 'nan','nan','nan'

    remove_data(path, folder, smiles)

    return test_init_rdkit_confs,test_prefilter_rdkit_confs,test_filter_rdkit_confs,round_confs,test_round_confs,test_charge,test_unique_confs,count,charge_com,multiplicity_com

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
    elif file == 'MeOH_SCF_error.out':
        os.chdir(path+'/'+folder+'/failed_error/SCF_error')
        assert file in glob.glob('*.*')
    elif file == 'MeOH_Unfinished.OUT':
        os.chdir(path+'/'+folder+'/failed_unfinished')
        assert file in glob.glob('*.*')

def check_com_files(path, folder, file):
    if file == 'Basis_set_error1.LOG' or file == 'Basis_set_error2.LOG':
        os.chdir(path+'/'+folder+'/new_gaussian_input_files/wb97xd-def2svp')
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
    elif file == 'MeOH_SCF_error.out':
        coordinates = 'H  -1.04798700   0.80281000  -0.68030200'
        coordinates_scf_found = find_coordinates(file,coordinates)
        assert coordinates_scf_found == 1
        com_input_line_scf = '# wb97xd/def2svp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=100) scrf=(SMD,solvent=Chloroform) scf=qc'
        input_found_scf = find_coordinates(file,com_input_line_scf)
        assert input_found_scf == 1
    elif file == 'MeOH_Unfinished.OUT':
        coordinates = 'H  -1.04779100   0.87481300  -0.58663200'
        coordinates_unfinished_found = find_coordinates(file,coordinates)
        assert coordinates_unfinished_found == 1
        com_input_line_unfinished = '# wb97xd/def2svp freq=noraman empiricaldispersion=GD3BJ opt=(calcfc,maxcycles=100) scrf=(SMD,solvent=Chloroform)'
        input_found_unfinished = find_coordinates(file,com_input_line_unfinished)
        assert input_found_unfinished == 1

def analysis(path, cmd_pyconfort, folder, file):
    os.chdir(path+'/'+folder)
    # the code will move the files the first time, this 'if' avoids errors
    files = glob.glob('*.log')
    if len(files) > 0:
        subprocess.call(cmd_pyconfort)
    # make sure the LOG files are in the right folders after analysis
    check_log_files(path, folder, file)
    # make sure the generated COM files have the right level of theory and geometries
    os.chdir(path+'/'+folder+'/new_gaussian_input_files/')
    check_com_files(path, folder, file)

def single_point(path, cmd_pyconfort, folder, file):
    os.chdir(path+'/'+folder)
    files = glob.glob('*.*')
    if len(files) > 0:
        subprocess.call(cmd_pyconfort)
    os.chdir(path+'/'+folder+'/finished/single_point_input_files/wb97xd-def2svp')
    assert len(glob.glob('*.*')) == 2

    if file == 'Pd_SP.LOG':
        count,NBO,pop,opt,_,_ = calc_genecp(file.split('.')[0]+'.com', ['Pd'])

    elif file == 'CH4_freq.log':
        count,NBO,pop,opt,_,_ = calc_genecp(file.split('.')[0]+'.com', ['C H'])

    return count,NBO,pop,opt

def conf_gen_exp_rules(path, folder, precision, cmd_exp_rules, smiles, E_confs_no_rules, E_confs_rules):
    # open right folder and run the code
    os.chdir(path+'/'+folder)
    if smiles == 'Ir_1':
        subprocess.call(cmd_exp_rules)

    # read the energies of the conformers with and without the exp_rules filter
    os.chdir(path+'/'+folder+'/rdkit_generated_sdf_files')

    test_E_confs_no_rules = calc_energy(smiles+'_rdkit.sdf')
    test_E_confs_rules = calc_energy(smiles+'_rdkit_filter_exp_rules.sdf')

    # test the energies and number of conformers with and without filtering
    test_round_E_confs_no_rules = [round(num, precision) for num in test_E_confs_no_rules]
    round_E_confs_no_rules = [round(num, precision) for num in E_confs_no_rules]

    test_round_E_confs_rules = [round(num, precision) for num in test_E_confs_rules]
    round_E_confs_rules = [round(num, precision) for num in E_confs_rules]

    # tests charge and genecp
    os.chdir(path+'/'+folder+'/generated_gaussian_files/wb97xd-def2svp')
    file_exp_rules = glob.glob(smiles+'*')[0]
    test_com_files = len(glob.glob(smiles+'*'))
    _,_,_,_,test_charge,_ = calc_genecp(file_exp_rules, ['Ir'])

    remove_data(path, folder, smiles=False)

    return round_E_confs_no_rules,round_E_confs_rules,test_round_E_confs_no_rules,test_round_E_confs_rules,test_charge,test_com_files
