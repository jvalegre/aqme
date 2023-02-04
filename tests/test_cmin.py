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

# # saves the working directory
# w_dir_main = os.getcwd()
# cmin_methods_dir = w_dir_main + "/tests/cmin_methods"
# if not os.path.exists(cmin_methods_dir):
#     os.mkdir(cmin_methods_dir)
# cmin_xtb_dir = w_dir_main + "/tests/cmin_xtb"
# if not os.path.exists(cmin_xtb_dir):
#     os.mkdir(cmin_xtb_dir)

# tests of basic ANI and xTB optimizations
# @pytest.mark.parametrize(
#     "program, sdf, output_nummols",
#     [
#         # tests for conformer generation with RDKit
#         ("ani", "pentane_rdkit_methods.sdf", 4),
#         ("xtb", "pentane_rdkit_methods.sdf", 4),
#     ],
# )
# def test_cmin_methods(
#     program, sdf, output_nummols
# ):

#     # runs the program with the different tests
#     base_folder = f'CMIN_{program}'
#     cmin(
#         w_dir_main=cmin_methods_dir,
#         program=program,
#         destination=f'{cmin_methods_dir}/{base_folder}',
#         files=f'{cmin_methods_dir}/{sdf}'
#     )

#     file = f'{cmin_methods_dir}/{base_folder}/{sdf.split(".")[0]}_{program}.sdf'
#     file2 = f'{cmin_methods_dir}/{base_folder}/{sdf.split(".")[0]}_{program}_all_confs.sdf'

#     assert os.path.exists(file)
#     assert os.path.exists(file2)

#     mols = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
#     assert len(mols) == output_nummols

# # tests for PATHs
# @pytest.mark.parametrize(
#     "program, sdf, output_nummols",
#     [
#         # tests for conformer generation with RDKit
#         ("destination", 4),
#         ("working_dir", 4),
#     ],
# )
# def test_cmin_methods(
#     program, sdf, output_nummols
# ):

#     # runs the program with the different tests
#     base_folder = f'CMIN_{program}'
#     cmin(
#         w_dir_main=cmin_methods_dir,
#         program=program,
#         destination=f'{cmin_methods_dir}/{base_folder}',
#         files=f'{cmin_methods_dir}/{sdf}'
#     )

#     file = f'{cmin_methods_dir}/{base_folder}/{sdf.split(".")[0]}_{program}.sdf'
#     file2 = f'{cmin_methods_dir}/{base_folder}/{sdf.split(".")[0]}_{program}_all_confs.sdf'

#     assert os.path.exists(file)
#     assert os.path.exists(file2)

#     mols = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
#     assert len(mols) == output_nummols


# # tests for parameters of cmin paramters
# @pytest.mark.parametrize(
#     "program, sdf, metal_complex,metal_atoms,metal_oxi,complex_type, charge, mult, ewin_cmin, initial_energy_threshold, energy_threshold,rms_threshold, output_nummols",
#     [
#         # tests for conformer generation with RDKit
#         (
#             "xtb",
#             "pentane_rdkit.sdf",
#             False,
#             None,
#             None,
#             None,
#             0,
#             1,
#             5,
#             0.003,
#             0.4,
#             0.5,
#             4,
#         ),
#         (
#             "xtb",
#             "Pd_complex_0_rdkit.sdf",
#             True,
#             ["Pd"],
#             [2],
#             "squareplanar",
#             0,
#             1,
#             8,
#             0.004,
#             0.3,
#             0.2,
#             2,
#         ),
#     ],
# )
# def test_cmin_xtb_parameters(
#     program,
#     sdf,
#     metal_complex,
#     metal_atoms,
#     metal_oxi,
#     complex_type,
#     charge,
#     mult,
#     ewin_cmin,
#     initial_energy_threshold,
#     energy_threshold,
#     rms_threshold,
#     output_nummols,
# ):
#     os.chdir(cmin_xtb_dir)
#     # runs the program with the different tests
#         cmin(
#             w_dir_main=cmin_xtb_dir,
#             program=program,
#             files=sdf,
#             charge=charge,
#             mult=mult,
#             ewin_cmin=ewin_cmin,
#             initial_energy_threshold=initial_energy_threshold,
#             energy_threshold=energy_threshold,
#             rms_threshold=rms_threshold,
#             xtb_accuracy=xtb_accuracy,
#             xtb_electronic_temperature=xtb_electronic_temperature,
#             xtb_max_iterations=xtb_max_iterations,
#         )

#     # tests here
#     file = str("CMIN/" + sdf.split(".")[0] + "_" + program + ".sdf")
#     file2 = str(
#         "CMIN/" + sdf.split(".")[0] + "_" + program + "_all_confs.sdf"
#     )
#     assert os.path.exists(file)
#     assert os.path.exists(file2)

#     mols = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
#     assert len(mols) == output_nummols

#     assert int(mols[0].GetProp("Real charge")) == charge
#     assert int(mols[0].GetProp("Mult")) == mult

#     os.chdir(w_dir_main)


# tests for removing foler and creation of CMIN DAT and CSV files
# @pytest.mark.parametrize(
#     "folder_list, file_list",
#     [
#         ["/tests/cmin_methods/CMIN_xtb","/tests/cmin_methods/CMIN_ani"],
#         ["/tests/cmin_methods/CMIN_xtb*","/tests/cmin_methods/CMIN_ani*"]

#     ],
# )
# def test_remove(folder_list, file_list):
#     for i,folder in enumerate(folder_list):
#         shutil.rmtree(w_dir_main + "/" + folder)
#         for f in glob.glob(file_list[i]):
#             os.remove(f)