#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
import pandas as pd
from conftest import datapath
import subprocess

# The target value is gonna be the number of conformers (n_conf). The other
# parameters are gonna be variables used by DBGEN
@pytest.mark.parametrize("smiles, params_file, n_conf, E_confs, n_confs_xtb, E_confs_xtb",
[
    # Pentane example
    ('pentane.smi', 'params_test1.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test2.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test3.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test4.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test5.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test6.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test7.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test8.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
    # ('pentane.smi', 'params_test9.yaml', 2, [1,2,3,5], 2, [1,2,3,5]),
])

# PARAMETERS TESTED FROM PARAMS_TEST*.YAML FILES:
# sample, rms_threshold, energy_threshold, initial_energy_threshold, auto_sample, ff, ewin, xtb, dihedralscan,

def test_confgen(smiles, params_file, n_conf, E_confs, n_confs_xtb, E_confs_xtb):
    # saves the working directory
    path = os.getcwd()

    # Conformer generation using different parameters
    subprocess.call('python', '-m', 'DBGEN', '--varfile', params_file)

    df_output = pd.from_csv(params_file.SPLIT('.')[0])
    file = params_file.SPLIT('.')[0]

    # tests for RDKit
    os.chdir(path+'\\RDKit_generated_SDF_files')
    sdf_file_rdkit = Chem.SDMolSupplier(file+'_rdkit.sdf')
    test_init_rdkit_confs = df_output['rdKit']
    test_prefilter_rdkit_confs = df_output['rdKit']
    test_filter_rdkit_confs = df_output['rdKit']
    test_rdkit_E_confs = read_energies(file,log)

    assert n_confs == test_n_confs
    assert E_confs == test_E_confs

    # tests for xtb
    if xtb:
        sdf_file_xtb = smiles.split('.')[0]+'_xtb.sdf'
        sdf_file_xtb = Chem.SDMolSupplier(file_xtb)
        test_E_confs_xtb = read_energies(file_xtb,log)

        assert n_confs_xtb == test_n_confs_xtb
        assert E_confs_xtb == test_E_confs_xtb

# MISSING CHECKS:
# CHECK THAT THE AUTO FUNCTION IS WORKING
# CHECK THAT CSV GETS THE CORRECT INFO OF NUMBERS AND TIME!
# CHECK COM FILES generation
# CHECK THAT N OF DUPLICATES FOR THE DIFFERENT FILTERS IS ALWAYS THE samE
# CHECK THAT THE METAL PART IS WORKING
# CHECK THAT MULTIPLE METALS ARE COMPATIBLE
# CHECK THAT THE TEMPLATE PART IS WORKING
# CHECK THAT THE EXPEIRMENTAL RULES IS WORKING
# CHECK THAT ANALYSIS IS WORKING FROM ALL THE DIFFERENT FOLDERS (USE SMALL TEST LOG FILES WITH ALL THE TYPE OF ERRORS)
# CHECK THAT AUTOPREP IS WORKING
# CHECK THAT GENECP AND GEN IS WORKING
