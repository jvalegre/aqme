#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#     Conformer generation of organic molecules      #
######################################################.

import os
import glob
import subprocess
import pytest
from definitions_testing import remove_data

# saves the working directory
path_orca = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, charge_orca, multiplicity_orca, param, job_type",
[
    # tests for ORCA in QPREP
    ('Organic_molecules', 'pentane.smi', 'params_test28.yaml', 0, 1, False, 'QPREP'), # ref test for ORCA (all default options)
    ('Organic_molecules', 'pentane.smi', 'params_test29.yaml', 0, 1, 'nprocs', 'QPREP'), # test for nprocs in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test30.yaml', 0, 1, 'mem', 'QPREP'), # test for mem in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test31.yaml', 0, 1, 'mdci_orca', 'QPREP'), # test for mdci_orca in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test32.yaml', 0, 1, 'print_mini_orca', 'QPREP'), # test for print_mini_orca in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test33.yaml', 0, 1, 'set_input_line', 'QPREP'), # test for set_input_line in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test34.yaml', 0, 1, 'solvent_CPCM', 'QPREP'), # test for solvent_CPCM in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test35.yaml', 0, 1, 'solvent_SMD', 'QPREP'), # test for solvent_SMD in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test36.yaml', 0, 1, 'cpcm_input', 'QPREP'), # test for cpcm_input in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test37.yaml', 0, 1, 'orca_scf_iters', 'QPREP'), # test for orca_scf_iters in ORCA
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test6.yaml', 1, 1, 'QPREP'), # test with metals and auxiliary basis sets in ORCA
])

def test_confgen_organic(folder, smiles, params_file, n_confs_organic, prefilter_confs_rdkit_organic_organic, filter_confs_rdkit_organic, E_confs, charge_organic, multiplicity_organic, dihedral, xTB_ANI):
    # runs the program with the different tests
    cmd_orca = ['python', '-m', 'pyconfort', '--varfile', params_file]

    os.chdir(path_orca+'/'+folder)
    subprocess.call(cmd_orca)

    os.chdir(path_orca+'/'+folder+'/QMCALC/ORCA/wb97xd-6-31g(d)')
    inp_file = glob.glob('*.inp')[0]

    # check that the INP files contain the right parameters
    outlines = com_lines(inp_file)

    if job_type == False:
        assert outlines[2].find('%maxcore 96000')
        assert outlines[4].find('%pal nprocs 24 end')
        assert outlines[5].find('! cc-pVTZ/C DLPNO-CCSD(T)')
        assert outlines[6].find('%scf maxiter 500')
        assert outlines[9].find('printlevel mini')
        assert outlines[16].find('Dipole False')
        assert outlines[18].find('* xyz 0 1')
        assert outlines[36].find('*')

    elif job_type == 'nprocs':
        assert outlines[4].find('%pal nprocs 12 end')

    elif job_type == 'mem':
        assert outlines[2].find('%maxcore 4000')

    elif job_type == 'mdci_orca':
        assert outlines[8].find('% mdci')
        assert outlines[9].find('Density None')
        assert outlines[10].find('end')

    elif job_type == 'print_mini_orca':
        assert outlines[8].find('* xyz 0 1')

    elif job_type == 'set_input_line':
        assert outlines[5].find('! cc-pVTZ/C DLPNO-CCSD(T) NormalPNO TightSCF Extrapolate(2/3,cc)')

    elif job_type == 'solvent_CPCM':
        assert outlines[6].find('! CPCM(THF)')
        assert outlines[7].find('%scf maxiter 500')

    elif job_type == 'solvent_SMD':
        assert outlines[6].find('%cpcm')
        assert outlines[7].find('smd true')
        assert outlines[8].find('SMDsolvent "acetone"')
        assert outlines[9].find('end')

    elif job_type == 'cpcm_input':
        assert outlines[6].find('! CPCM(THF)')
        assert outlines[7].find('%cpcm')
        assert outlines[8].find('surfacetype gepol_ses_gaussian')
        assert outlines[9].find('end')

    elif job_type == 'orca_scf_iters':
        assert outlines[6].find('%scf maxiter 100')

    remove_data(path_orca, folder, smiles)
