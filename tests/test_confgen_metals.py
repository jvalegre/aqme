#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#       Conformer generation of metal complexes      #
######################################################.

import os
import pytest
from definitions_testing import conf_gen

# saves the working directory
path_metals = os.getcwd()
# decimal digits for comparing E
precision_metals = 5

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1",
[
    # tests for conformer generation with RDKit, xTB and ANI1
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test1.yaml', 10, 1, 1, [3.82210, 3.82412, 4.13001, 4.13016, 4.13232, 4.18135, 7.54147, 7.58881], 1, False, False), # test single metal with genecp
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test2.yaml', 96, 'nan', 33, [1.53156, 1.55012, 1.5506, 1.55114, 1.55163, 1.57662], 1, True, False), # test with dihedral scan
    ('Metal_complexes', 'Ag_Au_complex.smi', 'params_Ag_Au_test1.yaml', 360, 352, 0, [10.45361, 11.58698, 20.01407, 20.66195, 20.7945, 22.27646, 23.33896, 23.39252], 0, False, False), # test 2 metals with genecp
    ('Metal_complexes', 'Ag_Au_complex_2.smi', 'params_Ag_Au_test1.yaml', 20, 0, 0, [37.95387, 41.58506, 41.8691, 42.4647, 42.84581, 46.75471, 48.4257, 48.67218, 49.49887, 51.41683, 52.68135, 53.25264, 54.42994, 54.5214, 56.52762, 59.4592, 59.75275, 61.4215, 63.41206, 69.54026], 1, False, False), # test 2 metals with genecp and charge
    ('Metal_complexes', 'Ag_Au_complex_2.smi', 'params_Ag_Au_test2.yaml', 240, 'nan', 146, [36.49354, 36.77305, 36.82544, 36.97313, 37.23573, 37.24436, 37.56838, 38.25121, 39.10735, 41.20973, 41.24798, 41.26245, 41.35228, 41.35546, 41.35927, 41.46424, 41.53602, 41.59332, 41.62411, 41.6405, 41.72955, 41.73488, 41.79588, 41.86305, 41.90294, 41.92186, 41.92265, 41.93435, 42.05995, 42.06945, 42.27756, 42.37207, 42.45762, 42.48821, 42.50402, 42.53456, 42.95388, 43.13758, 43.25449, 44.78602, 44.93579, 44.96228, 45.26116, 45.26598, 45.44017, 45.55071, 45.8679, 46.061, 47.5134, 47.90304, 48.0069, 48.14663, 48.16596, 48.23897, 48.28195, 48.32833, 48.35315, 48.38528, 48.42291, 48.46023, 48.54182, 48.91924, 48.99078, 49.08647, 49.13862, 49.29747, 49.35792, 49.7613, 50.33035, 50.35371, 50.47871, 50.48243, 50.52847, 50.78177, 50.89001, 50.91743, 51.16714, 51.29948, 51.40979, 51.42053, 51.46385, 51.76021, 51.9065, 52.01792, 52.09117, 52.15601, 52.17765, 52.244, 52.2868, 52.51926, 52.69236, 53.57123, 53.60618, 53.77328, 53.86299, 53.88996, 53.9276, 53.93788, 53.94374, 53.97461, 53.97969, 54.07632, 54.10034, 54.13162, 54.15401, 54.25457, 54.26089, 54.2831, 54.28863, 55.60037, 55.64299, 55.77447, 55.88946, 55.9445, 56.08696, 56.91175, 57.01257, 57.11954, 57.21402, 57.38656, 57.47613, 58.04783, 58.16629, 58.24595, 58.27262, 58.3785, 58.40188, 58.46164, 58.63853, 58.70071, 58.72827, 58.8047, 58.91022, 59.02796, 59.36356, 60.17173, 60.54619, 61.79955, 61.91091, 61.95232, 62.04226, 62.15073, 62.2454, 62.6222, 63.04568, 67.79518], 1, True, False), # test with dihedral scan
    ('Metal_complexes', 'Au_linear.smi', 'params_Au_test1.yaml', 20, 18, 0, [249.80053, 250.17904, 253.42573], -1, False, False), # test squarepyramidal with genecp
    ('Metal_complexes', 'Au_linear.smi', 'params_Au_test2.yaml', 3456, 'nan', 3, [249.80053, 250.17904, 253.42573], -1, False, False), # test squarepyramidal with genecp
    ('Metal_complexes', 'Pd_squareplanar.smi', 'params_Pd_test1.yaml', 360, 349, 0, [8.97676, 9.01872, 9.10379, 9.15897, 9.16366, 9.17272, 9.18237, 9.19448, 9.22382, 9.26467, 9.43331], 0, False, False), # test squareplanar template with gen
    ('Metal_complexes', 'Pd_squareplanar.smi', 'params_Pd_test2.yaml', 48, 'nan', 6, [8.97676, 9.01872, 9.15897, 9.17272, 9.18237, 9.19448], 0, True, False), # test nodihedral
    ('Metal_complexes', 'Rh_squarepyramidal.smi', 'params_Rh_test1.yaml', 10, 7, 0, [249.80053, 250.17904, 253.42573], 0, False, False), # test squarepyramidal with genecp
])

def test_confgen_metals(folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1):
    # runs the program with the different tests
    cmd_metals = ['python', '-m', 'pyconfort', '--varfile', params_file]

    conf_gen(path_metals, precision_metals, cmd_metals, folder, smiles, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, dihedral, xTB_ANI1)

# MISSING CHECKS:
# experimental rules for confgen SDF to COM
