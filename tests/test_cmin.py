#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#                   CMIN module                      #
######################################################.

import os
import glob
import pytest
from aqme.cmin import cmin
import rdkit
import shutil

# saves the working directory
w_dir_main = os.getcwd()
cmin_methods_dir = w_dir_main + "/tests/cmin_methods"
if not os.path.exists(cmin_methods_dir):
    os.mkdir(cmin_methods_dir)

# tests of basic ANI and xTB optimizations
@pytest.mark.parametrize(
    "path, program, sdf, output_nummols",
    [
        # tests for conformer generation with RDKit
        ("complete", "ani", "pentane_rdkit_methods.sdf", 4),
        ("complete", "xtb", "pentane_rdkit_methods.sdf", 4),
        ("partial", "ani", "tests/cmin_methods/pentane_rdkit_methods.sdf", 4), # test for partial path in the files option
        ("name", "ani", "pentane_rdkit_methods.sdf", 4), # test for direct name in the files option
    ],
)
def test_cmin_methods(
    path, program, sdf, output_nummols
):

    # runs the program with the different tests
    os.chdir(w_dir_main)
    if path == 'complete':
        cmin(program=program,files=f'{cmin_methods_dir}/{sdf}')
        os.chdir(cmin_methods_dir)
    elif path == 'partial':
        cmin(program=program,files=f'{sdf}')
        os.chdir(cmin_methods_dir)
        sdf = 'pentane_rdkit_methods.sdf'
    elif path == 'name':
        os.chdir(cmin_methods_dir) # first go to the folder with SDF
        cmin(program=program,files=f'{sdf}')


    file = f'{cmin_methods_dir}/CMIN/{sdf.split(".")[0]}_{program}.sdf'
    file2 = f'{cmin_methods_dir}/CMIN/{sdf.split(".")[0]}_{program}_all_confs.sdf'

    assert os.path.exists(file)
    assert os.path.exists(file2)

    mols = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
    mols_all = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
    # in this case, RDKit, ANI and xTB should lead to the same 4 conformers
    assert len(mols_all) == output_nummols
    assert len(mols) == output_nummols

    # check that the optimizations work (different geometry than initial)
    outfile = open(file, "r")
    outlines = outfile.readlines()
    outfile.close()
    coords = ['2.5236','0.0073','0.1777']
    for coord in coords:
        assert coord not in outlines[4]
    os.chdir(w_dir_main)

# tests for removing foler
@pytest.mark.parametrize(
    "folder_list, file_list",
    [
        (
            ["tests/cmin_methods/CMIN"],["tests/cmin_methods/CMIN*"]
        ),
    ],
)
def test_remove(folder_list, file_list):
    for i,folder in enumerate(folder_list):
        shutil.rmtree(w_dir_main + "/" + folder)
        for f in glob.glob(file_list[i]):
            os.remove(f)
    os.chdir(w_dir_main)
