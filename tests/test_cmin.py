#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#                   CMIN module                      #
######################################################.

import os
import pytest
from aqme.cmin import cmin
import rdkit
import shutil
import glob

# saves the working directory
w_dir_main = os.getcwd()
cmin_methods_dir = w_dir_main + "/tests/cmin_methods"
if not os.path.exists(cmin_methods_dir):
    os.mkdir(cmin_methods_dir)
cmin_xtb_dir = w_dir_main + "/tests/cmin_xtb"
if not os.path.exists(cmin_xtb_dir):
    os.mkdir(cmin_xtb_dir)

# tests for parameters of csearch random initialzation
@pytest.mark.parametrize(
    "program, sdf, ani_method, xtb_method, opt_steps, opt_fmax, output_nummols",
    [
        # tests for conformer generation with RDKit
        ("ani", "pentane_rdkit_methods.sdf", "ANI1ccx", None, 100, 0.08, 4),
        ("xtb", "pentane_rdkit_methods.sdf", None, "GFN2-xTB", 400, 0.03, 4),
    ],
)
def test_cmin_methods(
    program, sdf, ani_method, xtb_method, opt_steps, opt_fmax, output_nummols
):
    os.chdir(cmin_methods_dir)
    # runs the program with the different tests
    if program == "ani":
        cmin(
            w_dir_main=cmin_methods_dir,
            ani_method=ani_method,
            program=program,
            files=sdf,
            opt_steps=opt_steps,
            opt_fmax=opt_fmax,
        )
    elif program == "xtb":
        cmin(
            w_dir_main=cmin_methods_dir,
            xtb_method=xtb_method,
            program=program,
            files=sdf,
            opt_steps=opt_steps,
            opt_fmax=opt_fmax,
        )

    # tests here
    file = str("CMIN/" + program + "/" + sdf.split(".")[0] + "_" + program + ".sdf")
    file2 = str(
        "CMIN/" + program + "/" + sdf.split(".")[0] + "_" + program + "_all_confs.sdf"
    )
    assert os.path.exists(file)
    assert os.path.exists(file2)

    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    os.chdir(w_dir_main)


# # tests for parameters of cmin paramters
# @pytest.mark.parametrize(
#     "program, sdf, metal_complex,metal_atoms,metal_oxi,complex_type, charge, mult, xtb_solvent, ewin_cmin, initial_energy_threshold, energy_threshold,rms_threshold,xtb_accuracy,xtb_electronic_temperature, xtb_max_iterations, output_nummols",
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
#             "Chloroform",
#             5,
#             0.003,
#             0.4,
#             0.5,
#             0.01,
#             298,
#             500,
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
#             "none",
#             8,
#             0.004,
#             0.3,
#             0.2,
#             0.02,
#             300,
#             400,
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
#     xtb_solvent,
#     ewin_cmin,
#     initial_energy_threshold,
#     energy_threshold,
#     rms_threshold,
#     xtb_accuracy,
#     xtb_electronic_temperature,
#     xtb_max_iterations,
#     output_nummols,
# ):
#     os.chdir(cmin_xtb_dir)
#     # runs the program with the different tests
#     if not metal_complex:
#         cmin(
#             w_dir_main=cmin_xtb_dir,
#             program=program,
#             files=sdf,
#             charge=charge,
#             mult=mult,
#             xtb_solvent=xtb_solvent,
#             ewin_cmin=ewin_cmin,
#             initial_energy_threshold=initial_energy_threshold,
#             energy_threshold=energy_threshold,
#             rms_threshold=rms_threshold,
#             xtb_accuracy=xtb_accuracy,
#             xtb_electronic_temperature=xtb_electronic_temperature,
#             xtb_max_iterations=xtb_max_iterations,
#         )
#     else:
#         cmin(
#             w_dir_main=cmin_xtb_dir,
#             program=program,
#             files=sdf,
#             metal_complex=metal_complex,
#             metal_atoms=metal_atoms,
#             metal_oxi=metal_oxi,
#             complex_type=complex_type,
#             charge=charge,
#             mult=mult,
#             xtb_solvent=xtb_solvent,
#             ewin_cmin=ewin_cmin,
#             initial_energy_threshold=initial_energy_threshold,
#             energy_threshold=energy_threshold,
#             rms_threshold=rms_threshold,
#             xtb_accuracy=xtb_accuracy,
#             xtb_electronic_temperature=xtb_electronic_temperature,
#             xtb_max_iterations=xtb_max_iterations,
#         )

#     # tests here
#     file = str("CMIN/" + program + "/" + sdf.split(".")[0] + "_" + program + ".sdf")
#     file2 = str(
#         "CMIN/" + program + "/" + sdf.split(".")[0] + "_" + program + "_all_confs.sdf"
#     )
#     assert os.path.exists(file)
#     assert os.path.exists(file2)

#     mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
#     assert len(mols) == output_nummols

#     assert int(mols[0].GetProp("Real charge")) == charge
#     assert int(mols[0].GetProp("Mult")) == mult

#     os.chdir(w_dir_main)


# tests for removing foler
@pytest.mark.parametrize(
    "folder",
    [
        # tests for conformer generation with RDKit
        ("remove")
    ],
)
def test_remove(folder):
    os.chdir(w_dir_main)
    if os.path.exists(f"{w_dir_main}/tests/cmin_methods/CMIN"):
        files_remove = [f"{w_dir_main}/tests/cmin_methods/AQME_data.dat",f"{w_dir_main}/tests/cmin_methods/ase.opt"]
        files_remove.append(f"{w_dir_main}/tests/cmin_methods/ANI1_opt.traj")
        shutil.rmtree(f"{w_dir_main}/tests/cmin_methods/CMIN")
        for f in glob.glob(f"{w_dir_main}/tests/cmin_methods/CMIN*")+files_remove:
            os.remove(f)
    if os.path.exists(f"{w_dir_main}/tests/cmin_xtb/CMIN"):
        files_remove = [f"{w_dir_main}/tests/cmin_xtb/AQME_data.dat",f"{w_dir_main}/tests/cmin_xtb/ase.opt"]
        shutil.rmtree(f"{w_dir_main}/tests/cmin_xtb/CMIN")
        for f in glob.glob(f"{w_dir_main}/tests/cmin_xtb/CMIN*")+files_remove:
            os.remove(f)