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

# tests for individual metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, multiplicity, dihedral, xTB_ANI1, genecp, metal, template",
[
    # tests for conformer generation with RDKit, xTB and ANI1
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test1.yaml', 10, 5, 0, [5.75331, 5.75562, 5.75593, 5.77069, 5.78123], 1, 1, False, False, 'genecp', ['Ir'], False), # test single metal with genecp
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test2.yaml', 60, 'nan', 12, [5.71918, 5.75331, 5.75562, 5.75593, 5.75822, 5.77069, 5.78123, 5.80183, 7.78989, 7.80634, 7.80888, 7.83923], 1, 1, True, False, 'genecp', ['Ir'], False), # test with dihedral scan
    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test3.yaml', 10, 5, 0, [5.75331, 5.75562, 5.75593, 5.77069, 5.78123], 1, 3, False, False, 'genecp', ['Ir'], False), # test single metal with genecp
    ('Metal_complexes', 'Ag_Au_complex.smi', 'params_Ag_Au_test1.yaml', 360, 352, 0, [10.45361, 11.58698, 20.01407, 20.66195, 20.7945, 22.27646, 23.33896, 23.39252], 0, 1, False, False, 'genecp', ['Ag', 'Au'], False), # test 2 metals with genecp
    ('Metal_complexes', 'Ag_Au_complex_2.smi', 'params_Ag_Au_test2.yaml', 20, 0, 0, [39.43794, 41.28533, 43.29455, 44.53245, 47.28272, 48.32505, 48.8603, 50.19419, 52.50503, 54.24193, 55.38771, 56.12443, 56.66391, 58.54504, 58.75516, 60.24671, 62.73743, 63.6848, 67.87692, 71.8144], 1, 1, False, False, 'genecp', ['Ag', 'Au'], False), # test 2 metals with genecp and charge
    ('Metal_complexes', 'Ag_Au_complex_2.smi', 'params_Ag_Au_test3.yaml', 240, 'nan', 142, [37.4232, 37.42871, 37.44378, 37.65729, 37.71947, 37.79716, 38.0208, 38.40607, 39.21487, 40.41989, 40.64656, 40.82472, 40.97965, 41.04102, 41.08276, 41.16519, 41.19906, 41.47685, 41.50339, 41.51278, 41.8783, 42.25166, 42.34784, 42.40702, 42.66898, 43.13526, 43.5535, 43.60214, 43.80758, 43.8721, 43.91527, 44.34929, 44.91371, 46.85054, 47.0147, 47.31223, 47.42957, 47.74682, 47.85505, 48.05115, 48.30384, 48.47486, 48.49919, 48.68608, 48.98184, 49.00307, 49.0089, 49.12682, 49.50258, 49.90334, 49.95503, 49.98737, 50.16542, 50.35621, 50.44718, 50.60955, 50.63001, 50.6705, 50.7351, 50.76217, 50.76414, 50.80697, 50.99666, 51.00213, 51.01071, 51.02795, 51.1664, 51.34206, 51.39781, 51.53344, 51.6767, 51.68941, 51.73872, 51.81989, 51.97521, 52.02365, 52.06313, 52.79718, 53.0295, 53.15055, 53.85663, 54.14556, 54.18232, 54.24734, 54.31579, 54.5353, 54.54111, 54.63228, 54.8586, 54.92746, 55.01887, 55.06957, 55.07984, 55.14975, 55.18397, 55.35605, 55.4423, 55.60383, 55.67481, 55.90979, 56.11947, 56.21611, 56.32081, 56.40114, 56.53302, 56.60429, 56.69997, 56.70791, 56.84571, 57.08747, 57.28229, 57.37708, 57.84527, 57.90821, 57.99961, 58.03689, 58.08491, 58.08754, 58.09288, 58.75516, 58.90323, 59.33155, 59.73856, 60.82994, 61.17009, 61.43969, 61.50301, 62.18022, 62.30909, 62.38567, 62.45852, 62.48041, 62.63065, 62.94094, 63.01722, 63.56445, 63.61779, 66.02703, 66.5325, 68.24772, 69.15283, 69.50527], 1, 1, True, False, 'genecp', ['Ag', 'Au'], False), # test with dihedral scan
    ('Metal_complexes', 'Ag_Au_complex.smi', 'params_Ag_Au_test4.yaml', 360, 352, 0, [10.45361, 11.58698, 20.01407, 20.66195, 20.7945, 22.27646, 23.33896, 23.39252], 0, 3, False, False, 'genecp', ['Ag', 'Au'], False), # test 2 metals with multiplicity = 3
    ('Metal_complexes', 'Au_linear.smi', 'params_Au_test1.yaml', 20, 18, 0, [20.46458, 21.36073], -1, 1, False, False, 'gen', ['Au'], 'linear'), # test linear
    ('Metal_complexes', 'Au_linear.smi', 'params_Au_test2.yaml', 3456, 'nan', 4, [-0.05391, 0.07215, 0.64203, 0.64215], -1, 1, True, False, 'gen', ['Au'], 'linear'), # test linear with dihedral
    ('Metal_complexes', 'Cu_trigonal.smi', 'params_Cu_test1.yaml', 20, 19, 0, [88.48245], -2, 1, False, False, 'gen', ['Cu'], 'trigonalplanar'), # test trigonalplanar
    ('Metal_complexes', 'Cu_trigonal.smi', 'params_Cu_test2.yaml', 20, 'nan', 1, [88.48245], -2, 1, True, False, 'gen', ['Cu'], 'trigonalplanar'), # test trigonalplanar
    ('Metal_complexes', 'Pd_squareplanar.smi', 'params_Pd_test1.yaml', 360, 356, 0, [129.54222, 129.54379, 129.78549, 129.7858], 0, 1, False, False, 'gen', ['Pd'], 'squareplanar'), # test squareplanar template with gen
    ('Metal_complexes', 'Pd_squareplanar.smi', 'params_Pd_test2.yaml', 24, 'nan', 2, [10.58415, 10.92013], 0, 1, True, False, 'gen', ['Pd'], 'squareplanar'), # test nodihedral
    ('Metal_complexes', 'Pd_squareplanar.smi', 'params_Pd_test3.yaml', 360, 356, 0, [129.54222, 129.54379, 129.78549, 129.7858], 0, 3, False, False, 'gen', ['Pd'], 'squareplanar'), # test squareplanar template with gen
    ('Metal_complexes', 'Pd_squareplanar2.smi', 'params_Pd_test4.yaml', 360, 343, 11, [116.11335, 116.12354, 116.12364, 116.1243, 116.12595, 117.15133], 0, 1, False, False, 'gen', ['Pd'], 'squareplanar'), # test with 2 same ligands
    ('Metal_complexes', 'Rh_squarepyramidal.smi', 'params_Rh_test1.yaml', 10, 6, 0, [250.30172, 250.67941, 253.9259, 254.21842], 0, 1, False, False, 'genecp', ['Rh'], 'squarepyramidal'), # test squarepyramidal with 2 same ligands
    ('Metal_complexes', 'Rh_squarepyramidal.smi', 'params_Rh_test2.yaml', 48, 'nan', 3, [84.40719, 85.59864, 87.87227], 0, 1, True, False, 'genecp', ['Rh'], 'squarepyramidal'), # test squarepyramidal with 2 same ligands
    ('Metal_complexes', 'Rh_squarepyramidal2.smi', 'params_Rh_test3.yaml', 10, 5, 1, [225.29346, 225.51244, 225.59172, 225.92328], 1, 1, False, False, 'genecp', ['Rh'], 'squarepyramidal'), # test squarepyramidal with 3 same ligands
])

def test_confgen_metals(folder, smiles, params_file, n_confs, prefilter_confs_rdkit, filter_confs_rdkit, E_confs, charge, multiplicity, dihedral, xTB_ANI1, genecp, metal, template):
    # runs the program with the different tests
    cmd_metals = ['python', '-m', 'pyconfort', '--varfile', params_file]

    test_init_rdkit_confs,test_prefilter_rdkit_confs,test_filter_rdkit_confs,round_confs,test_round_confs,test_charge,test_unique_confs,count,charge_com,multiplicity_com = conf_gen(path_metals, precision_metals, cmd_metals, folder, smiles, E_confs, dihedral, xTB_ANI1, metal, template)

    # the assert statements are placed here, otherwise pytest doesn't explain the AssertionError
    # dihedral vs no dihedral scans
    if not dihedral:
        assert str(n_confs) == str(test_init_rdkit_confs[0])
        assert str(prefilter_confs_rdkit) == str(test_prefilter_rdkit_confs[0])
        assert str(filter_confs_rdkit) == str(test_filter_rdkit_confs[0])
    else:
        assert str(n_confs) == str(test_init_rdkit_confs[0])
        # I use the filter_confs_rdkit variable to assert for unique confs in dihedral scan
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

# MISSING:
# xTB with metals
