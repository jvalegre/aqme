#!/usr/bin/env python

######################################################.
# 	        Testing QCORR with pytest 	             #
######################################################.

import os
import glob
import subprocess
import pytest
from definitions_testing import remove_data

# saves the working directory
path_QPREP_syst = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, LoT",
[
    # tests of the analysis part (I use file as the output LOG files)
    ('Organic_molecules', 'pentane.smi', 'params_test38.yaml', 'b3lyp-6-31g(d)'), # tests for multiple levels of theory
    ('Organic_molecules', 'pentane.smi', 'params_test38.yaml', 'm062x-LANL2DZ'), # tests for multiple levels of theory
    ('Organic_molecules', 'pentane.smi', 'params_test38.yaml', 'wb97xd-def2tzvpp'), # tests for multiple levels of theory
])

def test_analysis_dup_sp(folder, smiles, params_file, LoT):

    cmd_aqme = ['python', '-m', 'aqme', '--varfile', params_file]

    if LoT == 'b3lyp-6-31g(d)':
        # run aqme
        os.chdir(path_QPREP_syst+'/'+folder+'/'+smiles.split('.')[0])
        subprocess.call(cmd_aqme)

    os.chdir(path_QPREP_syst+'/'+folder+'/'+smiles.split('.')[0]+'/QMCALC/G16/'+LoT)
    file_com = glob.glob('*.com')[0]

    outfile = open(file_com,"r")
    outlines = outfile.readlines()

    line = '%mem=96GB'
    line_number = 0
    assert outlines[line_number].find(line) > -1

    line = '%nprocshared=24'
    line_number = 1
    assert outlines[line_number].find(line) > -1

    line = '0 1'
    line_number = 6
    assert outlines[line_number].find(line) > -1

    if LoT == 'b3lyp-6-31g(d)':
        line2 = 'b3lyp/genecp freq=noraman opt=(maxcycles=100)'
    elif LoT == 'm062x-LANL2DZ':
        line2 = 'm062x/genecp freq=noraman opt=(maxcycles=100)'
    elif LoT == 'wb97xd-def2tzvpp':
        line2 = 'wb97xd/genecp freq=noraman opt=(maxcycles=100)'
    line_number2 = 2
    assert outlines[line_number2].find(line2) > -1

    line4 = 'H 0'
    line_number4 = 25
    assert outlines[line_number4].find(line4) > -1

    if LoT == 'b3lyp-6-31g(d)':
        line5 = '6-31g(d)'
    elif LoT == 'm062x-LANL2DZ':
        line5 = 'LANL2DZ'
    elif LoT == 'wb97xd-def2tzvpp':
        line5 = 'def2tzvpp'
    line_number5 = 26
    assert outlines[line_number5].find(line5) > -1

    line6 = 'C 0'
    line_number6 = 28
    assert outlines[line_number6].find(line6) > -1

    if LoT == 'b3lyp-6-31g(d)':
        line7 = 'LANL2DZ'
    elif LoT == 'm062x-LANL2DZ':
        line7 = 'def2svp'
    elif LoT == 'wb97xd-def2tzvpp':
        line7 = 'LANL2TZ'
    line_number7 = 29
    assert outlines[line_number7].find(line7) > -1

    line8 = 'C 0'
    line_number8 = 32
    assert outlines[line_number8].find(line8) > -1

    if LoT == 'b3lyp-6-31g(d)':
        line9 = 'LANL2DZ'
    elif LoT == 'm062x-LANL2DZ':
        line9 = 'def2svp'
    elif LoT == 'wb97xd-def2tzvpp':
        line9 = 'LANL2TZ'
    line_number9 = 33
    assert outlines[line_number9].find(line9) > -1

    outfile.close()

    if LoT == 'wb97xd-def2tzvpp':
        remove_data(path_QPREP_syst, folder, smiles)
