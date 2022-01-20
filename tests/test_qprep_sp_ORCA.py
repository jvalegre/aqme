#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#     Conformer generation of organic molecules      #
######################################################.

import os
import glob
import subprocess
import pytest
from definitions_testing import com_lines,remove_data

# saves the working directory
path_orca = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, job_type",
[
    # tests for ORCA in QPREP
    ('Organic_molecules', 'pentane.smi', 'params_test28.yaml', False), # ref test for ORCA (all default options)
    ('Organic_molecules', 'pentane.smi', 'params_test29.yaml', 'nprocs'), # test for nprocs in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test30.yaml', 'mem'), # test for mem in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test31.yaml', 'mdci_orca'), # test for mdci_orca in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test32.yaml', 'print_mini_orca'), # test for print_mini_orca in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test33.yaml', 'qm_input'), # test for qm_input in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test34.yaml', 'solvent_CPCM'), # test for solvent_CPCM in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test35.yaml', 'solvent_SMD'), # test for solvent_SMD in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test36.yaml', 'cpcm_input'), # test for cpcm_input in ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test37.yaml', 'orca_scf_iters'), # test for orca_scf_iters in ORCA
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test6.yaml', 'aux'), # test with metals and auxiliary basis sets in ORCA
    ('Analysis_ORCA', 'NaN', 'params_QCORR_ORCA_test.yaml', 'sp'), # test with the SP option with ORCA
])

def test_confgen_organic(folder, smiles, params_file, job_type):
    # runs the program with the different tests
    cmd_orca = ['python', '-m', 'aqme', '--varfile', params_file]

    if folder != 'Analysis_ORCA':

        os.chdir(path_orca+'/'+folder+'/'+smiles.split('.')[0])
        subprocess.call(cmd_orca)

        os.chdir(path_orca+'/'+folder+'/'+smiles.split('.')[0]+'/QCALC/ORCA/DLPNO-CCSD(T)-cc-pVTZ')
        inp_file = glob.glob('*.inp')[0]

        # check that the INP files contain the right parameters
        outfile = open(inp_file,"r")
        outlines = outfile.readlines()

        if job_type == False:
            assert outlines[2].find('%maxcore 96000') > -1
            assert outlines[4].find('%pal nprocs 24 end') > -1
            assert outlines[5].find('! cc-pVTZ/C DLPNO-CCSD(T)') > -1
            assert outlines[6].find('%scf maxiter 500') > -1
            assert outlines[9].find('printlevel mini') > -1
            assert outlines[16].find('Dipole False') > -1
            assert outlines[18].find('* xyz 0 1') > -1
            assert outlines[36].find('*') > -1

        elif job_type == 'nprocs':
            assert outlines[4].find('%pal nprocs 12 end') > -1

        elif job_type == 'mem':
            assert outlines[2].find('%maxcore 4000') > -1

        elif job_type == 'mdci_orca':
            assert outlines[8].find('% mdci') > -1
            assert outlines[9].find('Density None') > -1
            assert outlines[10].find('end') > -1

        elif job_type == 'print_mini_orca':
            assert outlines[8].find('* xyz 0 1') > -1

        elif job_type == 'qm_input':
            assert outlines[5].find('! cc-pVTZ/C DLPNO-CCSD(T) NormalPNO TightSCF Extrapolate(2/3,cc)') > -1

        elif job_type == 'solvent_CPCM':
            assert outlines[6].find('! CPCM(THF)') > -1
            assert outlines[7].find('%scf maxiter 500') > -1

        elif job_type == 'solvent_SMD':
            assert outlines[6].find('%cpcm') > -1
            assert outlines[7].find('smd true') > -1
            assert outlines[8].find('SMDsolvent "acetone"') > -1
            assert outlines[9].find('end') > -1

        elif job_type == 'cpcm_input':
            assert outlines[6].find('! CPCM(acetone)') > -1
            assert outlines[7].find('%cpcm') > -1
            assert outlines[8].find('surfacetype gepol_ses_gaussian') > -1
            assert outlines[9].find('end') > -1

        elif job_type == 'orca_scf_iters':
            assert outlines[6].find('%scf maxiter 100') > -1

        elif job_type == 'aux':
            assert outlines[6].find('%basis') > -1
            assert outlines[7].find('NewGTO 77 "def2-TZVPP" end') > -1
            assert outlines[8].find('NewAuxCGTO 77 "def2-TZVPP/C" end') > -1
            assert outlines[9].find('end') > -1

        outfile.close()

        remove_data(path_orca, folder, smiles)

    # test for SP with ORCA
    else:
        os.chdir(path_orca+'/'+folder)
        subprocess.call(cmd_orca)

        os.chdir(path_orca+'/'+folder+'/success/ORCA-SP_input_files/wb97xd-6-31g(d)')
        inp_file = 'CH4_Normal_termination_SPC.inp'

        # check that the INP files contain the right parameters
        outfile = open(inp_file,"r")
        outlines = outfile.readlines()

        assert len(glob.glob('*.inp')) == 1
        assert outlines[2].find('%maxcore 96000') > -1
        assert outlines[4].find('%pal nprocs 24 end') > -1
        assert outlines[5].find('! 6-31g(d) wb97xd') > -1
        assert outlines[6].find('%scf maxiter 500') > -1
        assert outlines[9].find('printlevel mini') > -1
        assert outlines[16].find('Dipole False') > -1
        assert outlines[18].find('* xyz 0 1') > -1
        assert outlines[24].find('*') > -1

        outfile.close()
