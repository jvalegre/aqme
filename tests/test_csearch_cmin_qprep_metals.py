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
precision_metals = 2

# tests for individual metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, multiplicity, dihedral, xTB_ANI1, genecp, metal, template",
[
    # tests for conformer generation with RDKit, xTB and ANI1
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test1.yaml', 10, 5, 0, [5.75331, 5.75562, 5.75593, 5.77069, 5.78123], 1, 1, False, False, 'genecp', ['Ir'], False), # test single metal with genecp
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test2.yaml', 60, 'nan', 12, [5.71918, 5.75331, 5.7556, 5.75593, 5.75822, 5.77069, 5.78123, 5.80183, 7.78989, 7.80634, 7.80888, 7.83923], 1, 1, True, False, 'genecp', ['Ir'], False), # test with SUMM and degree = 30
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test3.yaml', 10, 5, 0, [5.75331, 5.75562, 5.75593, 5.77069, 5.78123], 1, 3, False, False, 'genecp', ['Ir'], False), # multiplicity 3
#MISSING    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test4.yaml', 5, 0, 0, [-29977.81073, -29977.32134, -29976.96919, -29976.86387, -29955.55021], 1, 1, False, 'xTB', 'genecp', ['Ir'], False), # test single metal with genecp and xTB
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test5.yaml', 15, 'nan', 12, [5.71918, 5.75331, 5.7556, 5.75593, 5.75822, 5.77069, 5.78123, 5.80183, 7.78989, 7.80634, 7.80888, 7.83923], 1, 1, True, False, 'genecp', ['Ir'], False), # test with SUMM and degree = 30
    ('Metal_complexes', 'Ag_Au_complex.smi', 'params_Ag_Au_test1.yaml', 360, 352, 0, [10.45361, 11.58698, 20.01407, 20.66195, 20.7945, 22.27646, 23.33896, 23.39252], 0, 1, False, False, 'genecp', ['Ag', 'Au'], False), # test 2 metals with genecp
    ('Metal_complexes', 'Ag_Au_complex_2.smi', 'params_Ag_Au_test2.yaml', 20, 0, 0, [39.43794, 41.28533, 43.29455, 44.53245, 47.28272, 48.32505, 48.8603, 50.19419, 52.50503, 54.24193, 55.38771, 56.12443, 56.66391, 58.54504, 58.75516, 60.24671, 62.73743, 63.6848, 67.87692, 71.8144], 1, 1, False, False, 'genecp', ['Ag', 'Au'], False), # test 2 metals with genecp and charge
#FIXED (STILL ANSWER Q) - why the E are not like in the above case if it's the same smi but from rdkit to summ?
    ('Metal_complexes', 'Ag_Au_complex_2.smi', 'params_Ag_Au_test3.yaml', 240, 'nan', 152, [37.45745, 37.47559, 37.52452, 37.62471, 37.71385, 37.7933, 37.84798, 38.06928, 38.21904, 38.55595, 38.97523, 40.50324, 41.11259, 41.14737, 41.22266, 41.22671, 41.47715, 41.52431, 41.6141, 41.97117, 41.9778, 42.2083, 42.24988, 42.36919, 42.48911, 42.67075, 43.16273, 43.82749, 44.04067, 44.10506, 44.20492, 44.31674, 44.33584, 45.55114, 47.01417, 47.21215, 47.31181, 47.43075, 47.79447, 48.01871, 48.28802, 48.32217, 48.47347, 48.4973, 48.55565, 48.76442, 48.82108, 48.95752, 49.0197, 49.16268, 49.37342, 49.51542, 49.66585, 49.86223, 49.93011, 49.96764, 50.22653, 50.38697, 50.42082, 50.46815, 50.59108, 50.7647, 50.79361, 50.81871, 50.96488, 50.99408, 51.10115, 51.21785, 51.25334, 51.26594, 51.37509, 51.53163, 51.88127, 52.03421, 52.5599, 52.58515, 52.71705, 52.7716, 52.89906, 52.93695, 52.99679, 53.16577, 53.55891, 53.82861, 53.92672, 54.0144, 54.21466, 54.2264, 54.31041, 54.41852, 54.58651, 54.74842, 54.85375, 54.89399, 54.90465, 54.96143, 54.9828, 55.04113, 55.05167, 55.06721, 55.10108, 55.17805, 55.28463, 55.30989, 55.43862, 55.48628, 55.54901, 55.69621, 55.94229, 55.96426, 56.39534, 56.43791, 56.47583, 56.5471, 56.54814, 56.55645, 56.64721, 56.64735, 57.13082, 57.8234, 57.92668, 57.92992, 57.95021, 58.06632, 58.07678, 58.11178, 58.21233, 58.613, 58.75516, 59.57906, 61.46019, 62.21297, 62.25169, 62.45727, 62.45955, 62.54085, 62.64427, 63.00202, 63.52375, 63.82314, 67.00769, 67.35636, 67.65801, 67.84382, 67.95804, 68.00704, 69.02229, 70.62388, 73.69538, 77.13289, 80.70287, 82.57872], 1, 1, True, False, 'genecp', ['Ag', 'Au'], False), # test with SUMM and degree = 30
    ('Metal_complexes', 'Ag_Au_complex.smi', 'params_Ag_Au_test4.yaml', 360, 352, 0, [10.45361, 11.58698, 20.01407, 20.66195, 20.7945, 22.27646, 23.33896, 23.39252], 0, 3, False, False, 'genecp', ['Ag', 'Au'], False), # test 2 metals with multiplicity = 3
    ('Metal_complexes', 'Au_linear.smi', 'params_Au_test1.yaml', 20, 18, 0, [20.46458, 21.36073], -1, 1, False, False, 'gen', ['Au'], 'linear'), # test linear
    ('Metal_complexes', 'Au_linear.smi', 'params_Au_test2.yaml', 3456, 'nan', 4, [20.46458, 20.7802, 21.36073, 21.88162], -1, 1, True, False, 'gen', ['Au'], 'linear'), # test linear with SUMM and degree = 30
    ('Metal_complexes', 'Au_linear.smi', 'params_Au_test3.yaml', 54, 'nan', 2, [20.46458, 21.36073], -1, 1, True, False, 'gen', ['Au'], 'linear'), # test linear with SUMM and degree = 120 (standard)
    ('Metal_complexes', 'Cu_trigonal.smi', 'params_Cu_test1.yaml', 20, 19, 0, [88.48245], -2, 1, False, False, 'gen', ['Cu'], 'trigonalplanar'), # test trigonalplanar
# FIXED - gen is not working, I atom is still there in the com file
    ('Metal_complexes', 'Cu_trigonal.smi', 'params_Cu_test2.yaml', 20, 'nan', 1, [88.48245], -2, 1, True, False, 'gen', ['Cu'], 'trigonalplanar'), # test with SUMM
    ('Metal_complexes', 'Pd_squareplanar.smi', 'params_Pd_test1.yaml', 360, 356, 3, [129.54222], 0, 1, False, False, 'gen', ['Pd'], 'squareplanar'), # test squareplanar template with gen
    ('Metal_complexes', 'Pd_squareplanar.smi', 'params_Pd_test2.yaml', 3, 'nan', 1, [129.54222], 0, 1, True, False, 'gen', ['Pd'], 'squareplanar'), # test squareplanar template with SUMM
    ('Metal_complexes', 'Pd_squareplanar2.smi', 'params_Pd_test3.yaml', 360, 343, 15, [116.11335, 117.15133], 0, 1, False, False, 'gen', ['Pd'], 'squareplanar'), # test with 2 same ligands
    ('Metal_complexes', 'Rh_squarepyramidal.smi', 'params_Rh_test1.yaml', 10, 6, 0, [250.30172, 250.67941, 253.9259, 254.21842], 0, 1, False, False, 'genecp', ['Rh'], 'squarepyramidal'), # test squarepyramidal with 2 same ligands
    ('Metal_complexes', 'Rh_squarepyramidal.smi', 'params_Rh_test2.yaml', 12, 'nan', 4, [250.30172, 250.67941, 253.9259, 254.21842], 0, 1, True, False, 'genecp', ['Rh'], 'squarepyramidal'), # test SUMM
    ('Metal_complexes', 'Rh_squarepyramidal2.smi', 'params_Rh_test3.yaml', 10, 5, 2, [225.29346, 225.51244, 225.7681], 1, 1, False, False, 'genecp', ['Rh'], 'squarepyramidal'), # test squarepyramidal with 3 same ligands
    # checks that the program tolerates atoms in the metal option with different hybridizations
    ('Metal_complexes', 'P_hybrid.smi', 'params_hybrid_test1.yaml', 10, 8, 0, ['NOT ASSIGNED YET!'], 0, 1, False, False, False, False, False), # test squarepyramidal with 3 same ligands
])

def test_confgen_metals(folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, multiplicity, dihedral, xTB_ANI1, genecp, metal, template):
    # runs the program with the different tests
    cmd_metals = ['python', '-m', 'aqme', '--varfile', params_file]

    test_init_rdkit_confs,test_prefilter_rdkit_confs,test_filter_rdkit_confs,round_confs,test_round_confs,test_charge,test_unique_confs,count,charge_com,multiplicity_com = conf_gen(path_metals, precision_metals, cmd_metals, folder, smiles, E_confs, dihedral, xTB_ANI1, metal, template)

    # the assert statements are placed here, otherwise pytest doesn't explain the AssertionError
    # dihedral (SUMM) vs no dihedral scans
    if not dihedral:
        assert str(n_confs) == str(test_init_rdkit_confs[0])
        assert str(prefilter_confs_rdkit) == str(test_prefilter_rdkit_confs[0])
        assert str(filter_confs_rdkit) == str(test_filter_rdkit_confs[0])
    else:
        assert str(n_confs) == str(test_init_rdkit_confs[0])
        # I use the filter_confs_rdkit variable to assert for unique confs in SUMM
        assert str(filter_confs_rdkit) == str(test_unique_confs[0])

    assert str(round_confs) == str(test_round_confs)
    assert str(charge) == str(test_charge[0])

    # make sure gen and genecp are written correctly
    if genecp == 'gen':
        assert count == 1
    elif genecp == 'genecp':
        assert count == 2

    # make sure the COM files have the right charge and multiplicity
    assert str(charge_com) == str(charge)
    assert str(multiplicity_com) == str(multiplicity)
