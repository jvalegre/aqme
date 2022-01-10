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
        if line.find('>  <Energy>') > -1:
            energy_value = readlines[i+1].split()[0]
            energy.append(float(energy_value))
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
    # remove all the data created by the job except those files included in exceptions
    if not smiles:
        os.chdir(path+'/'+folder)
    else:
        os.chdir(path+'/'+folder+'/'+smiles.split('.')[0])
    all_data = glob.glob('*')
    discard_ext = ['sdf','csv','dat']
    exceptions = ['charged.csv','charged.sdf','pentane_n_lowest.sdf','Ir_4.sdf','Ir_4.csv','Ir_1_charge_0_passes_geom_rules.sdf','Ir_2_charge_1_passes_geom_rules.sdf','Ir_3_charge_1_fails_geom_rules.sdf']
    for _,file in enumerate(all_data):
        if len(file.split('.')) == 1:
            shutil.rmtree(file, ignore_errors=True)
        elif file not in exceptions:
            if file.split('.')[1] in discard_ext or file == 'cdx.smi':
                os.remove(file)

def rdkit_tests(df_output,dihedral,xTB_ANI,cmd_aqme):
    if not dihedral:
        if not xTB_ANI:
            test_init_rdkit_confs = df_output['RDKit-Initial-samples']
            test_prefilter_rdkit_confs = df_output['RDKit-initial_energy_threshold']
            test_filter_rdkit_confs = df_output['RDKit-RMSD-and-energy-duplicates']
            test_unique_confs = 'nan'
        elif xTB_ANI == 'xTB':
            test_init_rdkit_confs = df_output['xTB-Initial-samples']
            test_prefilter_rdkit_confs = df_output['xTB-initial_energy_threshold']
            test_filter_rdkit_confs = df_output['xTB-RMSD-and-energy-duplicates']
            test_unique_confs = 'nan'
        elif xTB_ANI == 'ANI':
            test_init_rdkit_confs = df_output['ANI-Initial-samples']
            test_prefilter_rdkit_confs = df_output['ANI-initial_energy_threshold']
            test_filter_rdkit_confs = df_output['ANI-RMSD-and-energy-duplicates']
            test_unique_confs = 'nan'
    else:
        if cmd_aqme[4] == 'params_Cu_test2.yaml':
            test_unique_confs = df_output['summ-conformers']
            test_init_rdkit_confs = df_output['RDKit-Initial-samples']
        else:
            test_init_rdkit_confs = df_output['summ-conformers']
            test_unique_confs = df_output['summ-Unique-conformers']
        test_prefilter_rdkit_confs = 'nan'
        test_filter_rdkit_confs = 'nan'


    return test_init_rdkit_confs, test_prefilter_rdkit_confs, test_filter_rdkit_confs, test_unique_confs

def conf_gen(path, precision, cmd_aqme, folder, smiles, E_confs, dihedral, xTB_ANI, metal, template):
    # open right folder and run the code
    os.chdir(path+'/'+folder+'/'+smiles.split('.')[0])
    subprocess.call(cmd_aqme)

    # Retrieving the generated CSV file
    os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/csv_files')

    df_output = pd.read_csv(smiles.split('.')[0]+'-CSEARCH-Data.csv')


    # tests for RDKit
    test_init_rdkit_confs, test_prefilter_rdkit_confs, test_filter_rdkit_confs, test_unique_confs = rdkit_tests(df_output,dihedral,xTB_ANI,cmd_aqme)

    # file_smi is a variable used for finding SDF and COM files
    if folder == 'Multiple':
        file_smi = 'pentane'
    else:
        file_smi = smiles.split('.')[0]

    # read the energies of the conformers
    try:
        if not xTB_ANI:
            if not dihedral:
                os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/rdkit')
            else:
                os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/summ')
            # this is to tests smi files with multiple smiles
            if folder == 'Multiple':
                if not dihedral:
                    test_rdkit_E_confs = calc_energy(file_smi+'_rdkit.sdf')
                else:
                    test_rdkit_E_confs = calc_energy(file_smi+'_summ.sdf')
            else:
                if template == 'squareplanar' or template == 'squarepyramidal':
                    if not dihedral:
                        test_rdkit_E_confs = calc_energy(file_smi+'_0_rdkit.sdf')
                    else:
                        test_rdkit_E_confs = calc_energy(file_smi+'_0_summ.sdf')
                else:
                    if not dihedral:
                        if smiles.split('.')[1] == 'sdf' or smiles.split('.')[1] == 'cdx':
                            test_rdkit_E_confs = calc_energy(file_smi+'_0_rdkit.sdf')
                        else:
                            test_rdkit_E_confs = calc_energy(file_smi+'_rdkit.sdf')
                    else:
                        if smiles.split('.')[1] == 'sdf' or smiles.split('.')[1] == 'cdx':
                            test_rdkit_E_confs = calc_energy(file_smi+'_0_summ.sdf')
                        else:
                            test_rdkit_E_confs = calc_energy(file_smi+'_summ.sdf')
        elif xTB_ANI == 'xTB':
            os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/xtb')
            test_rdkit_E_confs = calc_energy(file_smi+'_xtb.sdf')
        elif xTB_ANI == 'ANI':
            os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/CSEARCH/ani')
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

        os.chdir(path+'/'+folder+'/'+smiles.split('.')[0]+'/QMCALC/G16/wb97xd-6-31g(d)')

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
    for _,outline in enumerate(outlines):
        if outline.find(coordinates) > -1:
            coordinates_found = 1
            break
    return coordinates_found

def com_lines(file):
    try:
        com_file = file.split('.')[0]+'.com'
        outfile = open(com_file,"r")
        outlines = outfile.readlines()
    except:
        outlines = None

    return outlines

def analysis(path, cmd_aqme, folder, file):
    os.chdir(path+'/'+folder)
    df_QCORR,dat_files = [],[]

    # the code will move the files the first time, this 'if' statement avoids errors
    files = glob.glob('*.log')
    if len(files) > 0:
        subprocess.call(cmd_aqme)

    if file != 'csv' and file != 'dat':
        # copy the lines from the generated COM files
        os.chdir(path+'/'+folder+'/input_files/run_2')
        outlines = com_lines(file)

        return outlines

    elif file == 'csv':
        os.chdir(path+'/'+folder+'/csv_files')
        df_QCORR = pd.read_csv('Analysis-Data-QCORR-run_1.csv')
        return df_QCORR, dat_files

    elif file == 'dat':
        os.chdir(path+'/'+folder+'/dat_files')
        dat_files = glob.glob('*.dat')
        return df_QCORR, dat_files

def single_point(path_analysis_dup_sp, folder, file):
    # copy the lines from the generated COM files
    os.chdir(path_analysis_dup_sp+'/'+folder+'/success/G16-SP_input_files/b3lyp-321g')
    outlines = com_lines(file.split('.')[0]+'_SPC.com')

    return outlines

def Ir_geom_rules(path, folder, precision, cmd_geom_rules, smiles, E_confs_no_rules, E_confs_rules):
    # open right folder and run the code
    os.chdir(path+'/'+folder)
    if folder == 'Ir_geom_rules' and smiles == 'Ir_1':
        subprocess.call(cmd_geom_rules)
    elif folder == 'Ir_geom_rules2':
        subprocess.call(cmd_geom_rules)

    # read the energies of the conformers with and without the geom_rules filter
    os.chdir(path+'/'+folder+'/CSEARCH/rdkit')
    test_E_confs_no_rules = calc_energy(smiles+'_rdkit.sdf')

    os.chdir(path+'/'+folder+'/CSEARCH/rdkit/filter_geom_rules')
    test_E_confs_rules = calc_energy(smiles+'_rdkit_filter_geom_rules.sdf')

    # test the energies and number of conformers with and without filtering
    test_round_E_confs_no_rules = [round(num, precision) for num in test_E_confs_no_rules]
    round_E_confs_no_rules = [round(num, precision) for num in E_confs_no_rules]

    test_round_E_confs_rules = [round(num, precision) for num in test_E_confs_rules]
    round_E_confs_rules = [round(num, precision) for num in E_confs_rules]

    # tests charge and genecp
    os.chdir(path+'/'+folder+'/QMCALC/G16/wb97xd-6-31g(d)')
    file_geom_rules = glob.glob(smiles+'*')[0]
    test_com_files = len(glob.glob(smiles+'*'))
    _,_,_,_,test_charge,_ = calc_genecp(file_geom_rules, ['Ir'])

    if smiles == 'Ir_7' or folder == 'Ir_geom_rules2':
        remove_data(path, folder, smiles=False)

    return round_E_confs_no_rules,round_E_confs_rules,test_round_E_confs_no_rules,test_round_E_confs_rules,test_charge,test_com_files

def get_not_empty_files(folder_for_files,variable,format):
    os.chdir(folder_for_files)
    for file in glob.glob('*.'+format):
        if os.stat(file).st_size > 0:
            variable += 1
    return variable

def Pd_geom_rules(path, folder, precision_geom_rules, cmd_geom_rules, smiles):

    test_sdf_created,test_sdf_final,test_com_files = 0,0,0

    # run the code
    os.chdir(path+'/'+folder)
    subprocess.call(cmd_geom_rules)

    # check the amount of rdkit sdf files before and after geom_rules with size > 0
    rdkit_sdf_folder = path+'/'+folder+'/CSEARCH/rdkit'
    test_sdf_created = get_not_empty_files(rdkit_sdf_folder,test_sdf_created,'sdf')

    geom_rules_sdf_folder = path+'/'+folder+'/CSEARCH/rdkit/filter_geom_rules'
    test_sdf_final = get_not_empty_files(geom_rules_sdf_folder,test_sdf_final,'sdf')

    # check the amount of com files created
    com_files_folder = path+'/'+folder+'/QMCALC/G16/wb97xd-6-31g(d)'
    test_com_files = get_not_empty_files(com_files_folder,test_com_files,'com')

    remove_data(path, folder, smiles=False)

    return test_sdf_created,test_sdf_final,test_com_files

def misc_sdf_test(path_misc, smiles):

    try:
        # run the code and get to the folder where the SDF files are created
        os.chdir(path_misc+'/'+'Misc/CSEARCH/rdkit')
        # find how many SDF files were generated
        test_goal = len(glob.glob(smiles.split('.')[0]+'*_rdkit.sdf'))
    except:
        test_goal = 'None'

    try:
        os.chdir(path_misc+'/'+'Misc/CSEARCH/summ')
        test_goal_2 = len(glob.glob(smiles.split('.')[0]+'*_summ.sdf'))
    except:
        test_goal_2 = 0

    remove_data(path_misc, 'Misc', smiles=False)

    return test_goal, test_goal_2

def misc_com_test(path_misc, smiles):
    # get the amount of COM files created
    os.chdir(path_misc+'/'+'Misc/QMCALC/G16/wb97xd-6-31g(d)')
    test_goal = len(glob.glob(smiles.split('.')[0]+'*.com'))

    remove_data(path_misc, 'Misc', smiles=False)

    return test_goal

def misc_genecp_test(path_misc, smiles):
    os.chdir(path_misc+'/'+'Misc/QMCALC/G16/wb97xd-6-31g(d)')
    com_gen_files = glob.glob(smiles.split('.')[0]+'*.com')
    # I take just the first COM file (the others should be the same)
    file_gen = com_gen_files[0]

    # this function will find the lines from the manually inputed genecp
    gen_input = '      0.0745000              1.0000000'
    ecp_input = 'C-ECP     4     78'
    gen_found = find_coordinates(file_gen,gen_input)
    ecp_found = find_coordinates(file_gen,ecp_input)

    remove_data(path_misc, 'Misc', smiles=False)

    return gen_found,ecp_found

def misc_freq_test(path_misc, smiles):
    os.chdir(path_misc+'/'+'Misc/QMCALC/G16/wb97xd-6-31g(d)')
    com_freq_files = glob.glob(smiles.split('.')[0]+'*.com')
    # I take just the first COM file (the others should be the same)
    file_freq = com_freq_files[0]

    freq_input = 'freq'
    max_input = 'opt=(maxcycles=250)'
    chk_input = '%chk=pentane'
    mem_input = '%mem=20GB'
    nprocs_input = '%nprocshared=15'
    solvent_input = 'scrf=(SMD,solvent=hexane)'
    dispersion_input = 'empiricaldispersion=GD3'

    freq_found = find_coordinates(file_freq,freq_input)
    max_found = find_coordinates(file_freq,max_input)
    chk_found = find_coordinates(file_freq,chk_input)
    mem_found = find_coordinates(file_freq,mem_input)
    nprocs_found = find_coordinates(file_freq,nprocs_input)
    solvent_found = find_coordinates(file_freq,solvent_input)
    dispersion_found = find_coordinates(file_freq,dispersion_input)

    remove_data(path_misc, 'Misc', smiles=False)

    return freq_found,max_found,chk_found,mem_found,nprocs_found,solvent_found,dispersion_found

def misc_lot_test(path_misc, smiles):
    lots = ['wb97xd-def2svp', 'b3lyp-6-31G','m062x-def2tzvp']
    genecp_lots = ['LANL2DZ','midix','LANL2TZ']
    for i,lot in enumerate(lots):
        try:
            os.chdir(path_misc+'/'+'Misc/QMCALC/G16/'+lot)
            com_lot_files = glob.glob(smiles.split('.')[0]+'*.com')
            if len(com_lot_files) > 0:
                # I take just the first COM file (the others should be the same)
                file_lot = com_lot_files[0]
                files_found = 1
            else:
                # a two here means that the folder are created with no files inside
                files_found = 2
                break

        except FileNotFoundError:
            files_found = 0
            break
        if files_found == 1:
            gen_lot_found = find_coordinates(file_lot,genecp_lots[i])
            if gen_lot_found == 0:
                break
        else:
            # a two here means that the problem comes from the previous part
            gen_lot_found = 2
            break

    remove_data(path_misc, 'Misc', smiles=False)

    return files_found,gen_lot_found

def misc_nocom_test(path_misc, smiles):
    # finds the folder where COM files are generated normally
    try:
        os.chdir(path_misc+'/'+'Misc/QMCALC/G16/wb97xd-6-31g(d)')
        test_goal = 1
    except:
        test_goal = 0
        # run the code and get to the folder where the SDF files are created
        os.chdir(path_misc+'/'+'Misc/CSEARCH/rdkit')

    # find if the code generates SDF files
    if len(glob.glob(smiles.split('.')[0]+'_rdkit.sdf')) > 0:
        test_goal_2 = 1
    else:
        test_goal_2 = 0

    remove_data(path_misc, 'Misc', smiles=False)

    return test_goal,test_goal_2
