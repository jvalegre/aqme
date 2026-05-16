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

# Mapeo de rutas relativas a la ubicación del archivo del test
tests_dir = os.path.dirname(os.path.abspath(__file__))
w_dir_main = os.path.dirname(tests_dir) # Raíz del repositorio (aqme)

cmin_methods_dir = os.path.join(tests_dir, "cmin_methods")
cmin_empty_dir = os.path.join(tests_dir, "cmin_empty")

if not os.path.exists(cmin_methods_dir):
    os.mkdir(cmin_methods_dir)

# tests of basic QME optimizations (xtb, mace, aimnet2)
@pytest.mark.parametrize(
    "path, program, sdf, output_nummols",
    [
        ("complete", "xtb", "pentane_rdkit_methods.sdf", 4),
        ("complete", "mace", "pentane_rdkit_methods.sdf", 4), 
        ("complete", "aimnet2", "pentane_rdkit_methods.sdf", 4),
        ("partial", "xtb", "tests/cmin_methods/pentane_rdkit_methods.sdf", 4), 
        ("name", "xtb", "pentane_rdkit_methods.sdf", 4), 
    ],
)
def test_cmin_methods(
    path, program, sdf, output_nummols
):
    # runs the program with the different tests
    os.chdir(cmin_methods_dir)
    if path == 'complete':
        file_path = os.path.normpath(f'{cmin_methods_dir}/{sdf}')
        cmin(program=program, files=file_path)
        os.chdir(w_dir_main)
    elif path == 'partial':
        os.chdir(w_dir_main)
        file_path = os.path.normpath(sdf)
        cmin(program=program, files=file_path)
        os.chdir(cmin_methods_dir)
        sdf = 'pentane_rdkit_methods.sdf'
    elif path == 'name':
        os.chdir(cmin_methods_dir) 
        cmin(program=program, files=sdf)

    file = os.path.normpath(f'{cmin_methods_dir}/CMIN/{sdf.split(".")[0]}.sdf')
    file2 = os.path.normpath(f'{cmin_methods_dir}/CMIN/All_confs/{sdf.split(".")[0]}_all_confs.sdf')

    assert os.path.exists(file)
    assert os.path.exists(file2)

    mols = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
    mols_all = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
    
    assert len(mols_all) == output_nummols
    assert len(mols) == output_nummols

    # check that the optimizations work (different geometry than initial)
    outfile = open(file, "r")
    outlines = outfile.readlines()
    outfile.close()
    coords = ['2.5236','0.0073','0.1777']
    for coord in coords:
        assert coord not in outlines[4]

    # Force the release of the RDKit file handle for Windows
    del mols
    del mols_all
    os.chdir(w_dir_main)

# Test to check the application of constraints (distance, angle, and dihedral)
# CMIN is executed applying the three types of constraints supported by QME FixInternals
def test_cmin_constraints():
    os.chdir(cmin_methods_dir)
    sdf = "pentane_rdkit_methods.sdf"
    file_path = os.path.normpath(f'{cmin_methods_dir}/{sdf}')

    cmin(
        program="xtb",
        files=file_path,
        constraints_dist=[[0, 1, 1.54]],
        constraints_angle=[[0, 1, 2, 109.5]],
        constraints_dihedral=[[0, 1, 2, 3, 180.0]]
    )
    os.chdir(w_dir_main)

    file = os.path.normpath(f'{cmin_methods_dir}/CMIN/{sdf.split(".")[0]}.sdf')
    file2 = os.path.normpath(f'{cmin_methods_dir}/CMIN/All_confs/{sdf.split(".")[0]}_all_confs.sdf')

    assert os.path.exists(file)
    assert os.path.exists(file2)

@pytest.mark.parametrize(
    "test",
    [
        ("empty_values"),
    ],
)
def test_empty_values(
    test
):
    # runs the program with the different tests
    os.chdir(cmin_empty_dir)
    cmin(program='xtb', files=f'{cmin_empty_dir}/*.sdf')
    os.chdir(w_dir_main)

    file = os.path.normpath(f'{cmin_empty_dir}/CMIN/a.sdf')
    file2 = os.path.normpath(f'{cmin_empty_dir}/CMIN/All_confs/a_all_confs.sdf')
    file3 = os.path.normpath(f'{cmin_empty_dir}/CMIN/b_fail.sdf')
    file4 = os.path.normpath(f'{cmin_empty_dir}/CMIN/All_confs/b_fail_all_confs.sdf')
    file5 = os.path.normpath(f'{cmin_empty_dir}/CMIN/c.sdf')
    file6 = os.path.normpath(f'{cmin_empty_dir}/CMIN/All_confs/c_all_confs.sdf')

    assert os.path.exists(file)
    assert os.path.exists(file2)
    assert os.path.exists(file3)
    assert os.path.getsize(file3) == 0
    assert os.path.exists(file4)
    assert os.path.getsize(file4) == 0
    assert os.path.exists(file5)
    assert os.path.exists(file6)


# tests for removing folder
@pytest.mark.parametrize(
    "folder_list, file_list",
    [
        (
            ["tests/cmin_methods/CMIN","tests/cmin_empty/CMIN"],["tests/cmin_methods/CMIN*","tests/cmin_empty/CMIN*"]
        ),
    ],
)
def test_remove(folder_list, file_list):
    for i,folder in enumerate(folder_list):
        target_folder = os.path.join(w_dir_main, folder)
        if os.path.exists(target_folder):
            shutil.rmtree(target_folder)
        
        parent_dir = os.path.dirname(target_folder)
        for f in glob.glob(os.path.join(parent_dir, "*.dat")):
            try:
                os.remove(f)
            except PermissionError:
                pass
    os.chdir(w_dir_main)