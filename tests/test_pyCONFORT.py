#!/usr/bin/env python
# -*- coding: utf-8 -*-

# RDKit imports
import os
import pytest
import pandas as pd
import subprocess

# xTB and ANI1 imports
import ase
import ase.optimize
import torch, torchani
os.environ['KMP_DUPLICATE_LIB_OK']='True'
device = torch.device('cpu')
model = torchani.models.ANI1ccx()
from ase.units import Hartree
from xtb import GFN2

# saves the working directory
path = os.getcwd()
# decimal digits for comparing E
precision = 5

def read_energies(file): # parses the energies from sdf files - then used to filter conformers
	energies = []
	f = open(file,"r")
	readlines = f.readlines()
	for i in range(len(readlines)):
		if readlines[i].find('>  <Energy>') > -1:
			energies.append(float(readlines[i+1].split()[0]))
	f.close()
	return energies

# tests for the DBGEN module with organic molecules and metal complexes
@pytest.mark.parametrize("smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, xTB_ANI1",
[
    ('pentane.smi', 'params_test1.yaml', 240, 236, 0, [-5.27175,-4.44184,-3.84858,-1.57172],0, False), # test sample = 'auto', auto_sample = 20
    ('pentane.smi', 'params_test2.yaml', 20, 16, 0, [-5.27175, -4.44184, -3.84858, -1.57172],0, False), # test sample = 20
    ('pentane.smi', 'params_test3.yaml', 20, 2, 13, [-5.27175, -4.44184, -4.44184, -3.84858, -1.57172],0, False), # test initial_energy_threshold = 1E-10
    ('pentane.smi', 'params_test4.yaml', 20, 10, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858, -1.57172],0, False), # test energy_threshold = 1E-15
    ('pentane.smi', 'params_test5.yaml', 20, 10, 0, [-5.27175, -5.27175, -5.27175, -5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858, -1.57172],0, False), # test rms_threshold = 1E-15
    ('pentane.smi', 'params_test6.yaml', 20, 2, 11, [-5.27175, -4.44184, -4.44184, -4.44184, -4.44184, -3.84858, -1.57172],0, False),
    ('pentane.smi', 'params_test7.yaml', 60, 56, 0, [-5.27175, -4.44184, -3.84858, -1.57172],0, False), # test sample = 'auto', auto_sample = 5
    ('pentane.smi', 'params_test8.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False), # test num_rot_bonds = 1
    ('pentane.smi', 'params_test9.yaml', 'nan', 'nan', 'nan', 'nan', 'nan', False), # test max_MolWt = 1
    ('pentane.smi', 'params_test10.yaml', 20, 16, 0, [2.52059, 3.68961, 4.94318, 6.51778],0, False), # test ff = 'UFF'
    ('pentane.smi', 'params_test11.yaml', 20, 0, 8, [-5.26093, -4.41687, -4.39313, -4.10961, -3.93585, -2.95568, -2.43353, -2.03709, -1.51856, -1.45757, -0.22202, 0.46406],0, False), # test opt_steps_RDKit = 40
    ('pentane.smi', 'params_test12.yaml', 20, 16, 0, [-5.27175,-4.44184,-3.84858,-1.57172],0, True), # test xTB = True
    ('pentane.smi', 'params_test13.yaml', 20, 16, 0, [-5.27175,-4.44184,-3.84858,-1.57172],0, True), # test ANI1ccx = True
    ('Ir_hexacoord.smi', 'params_Ir_test1.yaml', 1440, 1434, 0, [-5.26093, -4.41687, -4.39313, -4.10961, -3.93585, -2.95568, -2.43353, -2.03709, -1.51856, -1.45757, -0.22202, 0.46406],1, False), # test 1 single metal
    ('Ir_Ag_complex.smi', 'params_Ag_Au_test1.yaml', 2, 0, 0, [-5.26093, -4.41687, -4.39313, -4.10961, -3.93585, -2.95568, -2.43353, -2.03709, -1.51856, -1.45757, -0.22202, 0.46406],1, False), # test 2 metals
    ('Pd_squareplanar.smi', 'params_Pd_test1.yaml', 120, 111, 0, [-5.26093, -4.41687, -4.39313, -4.10961, -3.93585, -2.95568, -2.43353, -2.03709, -1.51856, -1.45757, -0.22202, 0.46406],1, False), # test squareplanar template
])

def test_confgen(smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, xTB_ANI1):
	# gets into the directory for testing SMILES
	os.chdir(path+'/'+smiles.split('.')[0])

	# Conformer generation using different parameters. It creates a CSV file
	subprocess.run(['python', '-m', 'DBGEN', '--varfile', params_file])

	# Retrieving the generated CSV file
	df_output = pd.read_csv(smiles.split('.')[0]+'-Duplicates Data.csv')

	file = params_file.split('.')[0]

	# tests for RDKit
	test_init_rdkit_confs = df_output['RDKIT-Initial-samples']
	test_prefilter_rdkit_confs = df_output['RDKit-energy-duplicates']
	test_filter_rdkit_confs = df_output['RDKit-RMSD-and-energy-duplicates']

	assert str(n_confs) == str(test_init_rdkit_confs[0])
	assert str(prefilter_confs_rdkit) == str(test_prefilter_rdkit_confs[0])
	assert str(filter_confs_rdkit) == str(test_filter_rdkit_confs[0])

	# read the energies of the conformers
	os.chdir(path+'/'+smiles.split('.')[0]+'/RDKit_generated_SDF_files')
	test_rdkit_E_confs = read_energies(smiles.split('.')[0]+'_rdkit.sdf')

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

# MISSING CHECKS:
# START FROM COM, SDF, XYZ
# CHECK COM FILES generation
# CHECK THAT THE METAL PART IS WORKING
# CHECK THAT MULTIPLE METALS ARE COMPATIBLE
# CHECK THAT THE TEMPLATE PART IS WORKING
# CHECK THAT THE EXPEIRMENTAL RULES IS WORKING
# CHECK THAT ANALYSIS IS WORKING FROM ALL THE DIFFERENT FOLDERS (USE SMALL TEST LOG FILES WITH ALL THE TYPE OF ERRORS)
# CHECK THAT AUTOPREP IS WORKING
# CHECK THAT GENECP AND GEN IS WORKING
