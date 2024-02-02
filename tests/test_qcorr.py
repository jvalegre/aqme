#!/usr/bin/env python

######################################################.
# 	        Testing QCORR with pytest 	             #
######################################################.

import os
from os import path
import glob
import pytest
import shutil
import subprocess
from pathlib import Path
import pandas as pd

# saves the working directory
path_main = os.getcwd()
path_qcorr = os.getcwd() + "/Example_workflows/QCORR_processing_QM_outputs"

# QCORR tests
@pytest.mark.parametrize(
    "init_folder, file, command_line, target_folder, restore_folder",
    [
        # QCORR analysis tests with standard options
        (
            "QCORR_1",
            "CH4.log",
            "run_QCORR",
            "success",  # there is a mix up between this and duplicate either one gets taken
            False,
        ),  # test successful termination
        (
            "QCORR_1",
            "MeOH_G09.log",
            None,
            "success",
            False,
        ),  # test successful termination
        (
            "QCORR_1",
            "z_CH4_duplicate.log",
            None,
            "failed/run_1/duplicates",
            False,
        ),  # test duplicates
        (
            "QCORR_1",
            "Basis_set_error1.log",
            None,
            "failed/run_1/error/basis_set_error",
            False,
        ),  # test incompatibilities with basis sets
        (
            "QCORR_1",
            "Basis_set_error2.log",
            None,
            "failed/run_1/error/basis_set_error",
            False,
        ),  # test incompatibilities with basis sets
        (
            "QCORR_1",
            "MeOH_G09_FAIL.log",
            None,
            "failed/run_1/error/not_specified_error",
            False,
        ),  # test error terminations
        (
            "QCORR_1",
            "CH2OH2_unfinished.log",
            None,
            "failed/run_1/error/not_specified_error",
            False,
        ),  # test unfinished calculations
        (
            "QCORR_1",
            "Imag_freq.log",
            None,
            "failed/run_1/extra_imag_freq",
            False,
        ),  # test imaginary frequencies
        (
            "QCORR_1",
            "imag_freq_no_opt.log",
            None,
            "failed/run_1/extra_imag_freq",
            False,
        ),  # test imaginary frequencies without OPT
        (
            "QCORR_1",
            "MeOH_SCF_error.log",
            None,
            "failed/run_1/error/scf_error",
            False,
        ),  # test SCF errors
        (
            "QCORR_1",
            "CH4_before_E.log",
            None,
            "failed/run_1/error/no_data",
            False,
        ),  # test calcs that finish before any coords are printed
        (
            "QCORR_1",
            "freq_conv_YYNN.log",
            None,
            "failed/run_1/freq_no_conv",
            False,
        ),  # test calcs with freq calcs that did not converge after OPT
        (
            "QCORR_1",
            "freq_ok_YYNN.log",
            None,
            "success",
            False,
        ),  # test calcs with freq calcs that did not converge after OPT
        (
            "QCORR_1",
            "TS_CH3HCH3_no_conv_freq.log",
            None,
            "failed/run_1/freq_no_conv",
            False,
        ),  # test calcs with freq calcs that did not converge after OPT in TSs
        (
            "QCORR_1",
            "bpinene_spin_contamin.log",
            None,
            "failed/run_1/spin_contaminated",
            False,
        ),  # test calcs with spin contamination
        (
            "QCORR_1",
            "CH4_T1_SP_spin_contamin.log",
            None,
            "failed/run_1/spin_contaminated",
            False,
        ),  # test calcs with spin contamination
        (
            "QCORR_1",
            "CH4_Fail_freq_only.log",
            None,
            "failed/run_1/error/not_specified_error",
            False,
        ),  # test for Normal terminated OPT and unfinished freq
        (
            "QCORR_1",
            "TS_CH3HCH3.log",
            None,
            "success",
            False,
        ),  # test successful termination in TSs
        (
            "QCORR_1",
            "TS_CH3HCH3_unfinished.log",
            None,
            "failed/run_1/error/not_specified_error",
            False,
        ),  # test unfinished TSs
        (
            "QCORR_1",
            "TS_CH3HCH3_imag_freq.log",
            None,
            "failed/run_1/extra_imag_freq",
            False,
        ),  # test imaginary frequencies for TSs
        (
            "QCORR_1",
            "TS_CH3HCH3_no_imag_freq.log",
            None,
            "failed/run_1/ts_no_imag_freq",
            False,
        ),  # test imaginary frequencies for TSs
        (
            "QCORR_1",
            "CH4_SP.log",
            None,
            "success/SP_calcs",
            False,
        ),  # test for single-point calcs
        (
            "QCORR_1",
            "MeOH_NMR.log",
            None,
            "success/SP_calcs",
            False,
        ),  # test for single-point calcs
        (
            "QCORR_1",
            "H_freq.log",
            None,
            "success",
            False,
        ),  # test successful termination with 1 atom
        (
            "QCORR_1",
            "H_SP.log",
            None,
            "success/SP_calcs",
            False,
        ),  # test for single-point calcs with 1 atom
        (
            "QCORR_1",
            "CO2_linear_3freqs_FAIL.log",
            None,
            "failed/run_1/linear_mol_wrong",
            False,
        ),  # test for linear mols with wrong number of freqs
        (
            "QCORR_1",
            "CO2_linear_4freqs.log",
            None,
            "success",
            False,
        ),  # test successful termination for linear mols
        (
            "QCORR_1",
            "nosymm.log",
            None,
            "failed/run_1/error/not_specified_error",
            False,
        ),  # test successful termination for linear mols
        (
            "QCORR_1",
            "json",
            None,
            "success/json_files",
            False,
        ),  # test for correct creation of json files (all the successful terminations only)
        (
            "QCORR_1",
            "fullcheck",
            None,
            "success/json_files",
            False,
        ),  # test for correct fullcheck option
        ("QCORR_1", "csv", None, None, False),  # test final csv file with results
        ("QCORR_1", "dat", None, None, False),  # test final dat file with results
        (
            "QCORR_1",
            "check_init",
            None,
            None,
            False,
        ),  # test that the folder for initial QM inputs is not created
        # Test isomerization filter
        (
            "QCORR_2",
            "CH4.log",
            "run_QCORR",
            "success",
            False,
        ),  # test successful termination when using the isomerization filter
        (
            "QCORR_2",
            "CH2OH2_isomerized.log",
            None,
            "failed/run_1/isomerization",
            False,
        ),  # test isomerized calcs
        (
            "QCORR_2",
            "check_fixed",
            None,
            None,
            False,
        ),  # test that the fixed_QM_input folder is not generated when isomerized
        # Test if the unsuccessful runs are created in the same folder (files are contained inside the unsuccessful run_1 folder)
        (
            "QCORR_5",
            "CH4.log",
            "run_QCORR",
            "success",
            False,
        ),  # test that the fixed_QM_input folder is not generated when isomerized
        (
            "QCORR_5",
            "CH2OH2_isomerized.log",
            None,
            "failed",
            False,
        ),  # test that the fixed_QM_input folder is not generated when isomerized
        (
            None,
            None,
            None,
            None,
            True,
        ),  # reset the initial folder to start another set of tests
        # Test if QCORR works using parameters from a YAML file
        (
            "QCORR_2b",
            "CH4.log",
            "run_QCORR",
            "success",
            False,
        ),  # test successful termination when using the isomerization filter
        (
            "QCORR_2b",
            "CH2OH2_isomerized.log",
            None,
            "failed/run_1/isomerization",
            False,
        ),  # test isomerized calcs
        # Test if the unsuccessful runs are created in the same folder (files are contained inside the parent QCORR folder)
        (
            "QCORR_5b",
            "CH4.log",
            "run_QCORR",
            "success",
            False,
        ),  # test that the fixed_QM_input folder is not generated when isomerized
        (
            "QCORR_5b",
            "CH2OH2_isomerized.log",
            "run_QCORR",
            "failed",
            False,
        ),  # test that the fixed_QM_input folder is not generated when isomerized
        (
            None,
            None,
            None,
            None,
            True,
        ),  # reset the initial folder to start another set of tests
        # Test if the vdwfrac and covfrac options work in QCORR
        (
            "QCORR_2c",
            "CH4.log",
            "run_QCORR",
            "success",
            False,
        ),  # test successful termination when using the isomerization filter
        (
            "QCORR_2c",
            "CH2OH2_isomerized.log",
            None,
            "success",
            False,
        ),  # test that the vdwfrac and covfrac are disabled now with very low values
        # Test other options from QCORR (deactivate the freq_conv filter, S**2 threshold, energy threshold for duplicates, imag freqs cut-offs)
        (
            "QCORR_1b",
            "freq_conv_YYNN.log",
            None,
            "success",
            False,
        ),  # test to disable freq calcs that did not converge after OPT
        (
            "QCORR_1b",
            "bpinene_spin_contamin.log",
            None,
            "success",
            False,
        ),  # test to change filter for spin contamination
        (
            "QCORR_1b",
            "z_CH4_duplicate.log",
            None,
            "success",
            False,
        ),  # test to change energy threshold to consider duplicates
        (
            "QCORR_1b",
            "imag_freq_no_opt.log",
            None,
            "success",
            False,
        ),  # test to change cut-off to consider imaginary frequencies
        (
            "QCORR_1b",
            "Imag_freq.log",
            None,
            "failed/run_1/extra_imag_freq",
            False,
        ),  # test to change amplitude for displacing imaginary frequencies
        (
            None,
            None,
            None,
            None,
            True,
        ),  # reset the initial folder to start another set of tests
        (
            "QCORR_1b",
            "Imag_freq_no_corr.log",
            None,
            "failed/run_1/extra_imag_freq",
            False,
        ),  # test the im_freq_input option
        (
            None,
            None,
            None,
            None,
            True,
        ),  # reset the initial folder to start another set of tests
        # QCORR analysis with no w_dir_main
        (
            "QCORR_1d",
            "CH4.log",
            "parent_run",
            "success",
            False,
        ),  # test successful termination with no w_dir_main
        # QCORR analysis with no w_dir_main
        (
            "QCORR_6",
            "CH4_dup.log",
            "run_QCORR",
            "failed/run_1/duplicates",
            False,
        ),  # test successful termination with no w_dir_main
        (
            "QCORR_7",
            "orca_TS_success.out",
            "run_QCORR",
            "success",
            False,
        ),  # test QCORR with ORCA optimizations
        (
            "QCORR_7",
            "orca_imag_freq.out",
            None,
            "failed/run_1/extra_imag_freq",
            False,
        ),  # test QCORR with ORCA optimizations
        (
            None,
            None,
            None,
            None,
            True,
        ),  # reset the initial folder to start another set of tests
        # add genECP test
        # isomerization with csv (ongoing)
        # isomeriz with csv for TSs (ongoing)
        # tell if the imag freq from a TS is right based on displacement (ongoing)
    ],
)
def test_QCORR_analysis(init_folder, file, command_line, target_folder, restore_folder):

    # start from main folder
    os.chdir(path_main)

    # copy the test folders
    if not path.exists(f"{path_main}/Example_workflows_original"):
        shutil.copytree(
            f"{path_main}/Example_workflows", f"{path_main}/Example_workflows_original"
        )

    # runs the program with the different tests
    w_dir_main = f"{path_qcorr}/{init_folder}"
    cmd_aqme = [
        "python",
        "-m",
        "aqme",
        "--qcorr",
        "--files",
        f"{w_dir_main}/*.log",
        "--freq_conv",
        "opt=(calcfc,maxstep=5)",
    ]

    if init_folder == "QCORR_1":
        if command_line is not None:
            subprocess.run(cmd_aqme)

        if file.split(".")[-1].lower() == "log":
            # ensure the output file moves to the right folder
            assert path.exists(f"{w_dir_main}/{target_folder}/{file}")

            # ensure that the com files are generated correctly
            if file.split(".")[0] in [
                "CH4",
                "MeOH_G09",
                "z_CH4_duplicate",
                "Basis_set_error1",
                "Basis_set_error2",
                "CH4_before_E",
                "bpinene_spin_contamin",
                "TS_CH3HCH3",
                "TS_CH3HCH3_no_imag_freq",
                "CH4_SP",
                "H_freq",
                "H_SP",
                "MeOH_NMR",
                "CO2_linear_4freqs",
                "freq_ok_YYNN",
                "CH4_T1_SP_spin_contamin",
            ]:
                assert not path.exists(
                    f'{w_dir_main}/failed/run_1/fixed_QM_inputs/{file.split(".")[0]}.com'
                )

            else:
                assert path.exists(
                    f'{w_dir_main}/failed/run_1/fixed_QM_inputs/{file.split(".")[0]}.com'
                )

                # ensure that QCORR applies the correct structural distortions to the errored calcs
                if file.split(".")[0] == "MeOH_G09_FAIL":
                    line_2 = "# opt=calcfc freq=noraman cc-pvtz scrf=(solvent=chloroform,pcm) pbe1pbe g09defaults"
                    line_6 = "0 1"
                    line_8 = "H   1.13330900   0.94774000  -0.00000300"
                    line_10 = "H   0.97796200  -0.55747000   0.87365300"

                elif file.split(".")[0] == "CH2OH2_unfinished":
                    line_2 = "# opt freq 3-21g m062x"
                    line_6 = "0 1"
                    line_8 = "H   0.09630700   1.26666100  -0.75061200"
                    line_10 = "O  -1.21496800  -0.19993700  -0.10923000"

                elif file.split(".")[0] == "TS_CH3HCH3_no_conv_freq":
                    line_2 = "# opt=(ts,noeigen,calcfc,maxstep=5) freq b3lyp/3-21g"
                    line_6 = "0 2"
                    line_8 = "H  -0.91171200   0.52637700  -1.63085200"
                    line_10 = "H   0.00000000  -1.05275500  -1.63085200"

                elif file.split(".")[0] == "CH4_Fail_freq_only":
                    line_2 = "# b3lyp/3-21G freq=noraman"
                    line_6 = "0 1"
                    line_8 = "H   0.63133100   0.63133100   0.63133100"
                    line_10 = "H  -0.63133100   0.63133100  -0.63133100"

                elif file.split(".")[0] == "CO2_linear_3freqs_FAIL":
                    line_2 = "# opt=maxcycles=100 freq=noraman b3lyp 6-31G symmetry=(PG=Cinfv)"
                    line_6 = "0 1"
                    line_8 = "O   0.00000000   1.18790900  -0.00032000"
                    line_10 = None

                elif file.split(".")[0] == "TS_CH3HCH3_unfinished":
                    line_2 = "# opt=(calcfc,ts,noeigen) freq b3lyp/3-21g"
                    line_6 = "0 2"
                    line_8 = "H  -0.90909900   0.52486800  -1.63442100"
                    line_10 = "H   0.00000000  -1.04973700  -1.63442100"

                elif file.split(".")[0] == "imag_freq_no_opt":
                    line_2 = "# M062X/Def2TZVP freq=noraman opt=(calcfc,maxstep=5)"
                    line_6 = "0 1"
                    line_8 = "C  -0.90757400   0.00709700  -0.00994500"
                    line_10 = "C   1.19952600   1.19528800   0.00698400"

                elif file.split(".")[0] == "TS_CH3HCH3_imag_freq":
                    line_2 = "# opt=(ts,noeigen,calcfc,maxstep=5) freq b3lyp/3-21g"
                    line_6 = "0 2"
                    line_8 = "H  -0.87171200   0.59637700  -1.63085200"
                    line_10 = "H  -0.08200000  -1.05275500  -1.63085200"

                elif file.split(".")[0] == "Imag_freq":
                    line_2 = "# opt=(calcfc,maxstep=5) freq 3-21g m062x"
                    line_6 = "0 1"
                    line_8 = "H   0.38503600  -0.39992100  -0.94851000"
                    line_10 = "H   0.20952000   1.07184200  -0.02121100"

                elif file.split(".")[0] == "MeOH_SCF_error":
                    line_2 = "# opt freq=noraman b3lyp/3-21g scf=xqc"
                    line_6 = "0 1"
                    line_8 = "H  -1.08038200   0.99705900  -0.00806000"
                    line_10 = "H  -1.06929400  -0.53149000   0.89474400"

                elif file.split(".")[0] == "freq_conv_YYNN":
                    line_2 = "# M062X/Def2TZVP freq=noraman opt=(calcfc,maxstep=5)"
                    line_6 = "0 1"
                    line_8 = "C  -0.90757400   0.00709700  -0.00594500"
                    line_10 = "C   1.19952600   1.19528800   0.00098400"

                elif file.split(".")[0] == "nosymm":
                    line_2 = "# m062x def2svp nosymm int=(ultrafine) scrf=(smd,solvent=tetrahydrofuran) opt freq=(noraman)"
                    line_6 = "0 1"
                    line_8 = "C   1.54100700  -0.43780100  -1.49642700"
                    line_10 = "N  -0.37309700  -1.37343200  -0.29066900"

                outfile = open(
                    f'{w_dir_main}/failed/run_1/fixed_QM_inputs/{file.split(".")[0]}.com',
                    "r",
                )
                outlines = outfile.readlines()
                outfile.close()

                assert outlines[2].strip() == line_2
                assert outlines[6].strip() == line_6
                assert outlines[8].strip() == line_8
                if line_10 is not None:
                    assert outlines[10].strip() == line_10

        elif file == "json":
            os.chdir(f"{w_dir_main}/{target_folder}")
            json_files = glob.glob("*.json")
            assert len(json_files) == 6

        elif file == "fullcheck":
            target_fullcheck = ["-- Full check analysis --\n"]
            target_fullcheck.append("x  Different program used in the calculations:\n")
            target_fullcheck.append("     * Gaussian 09, Revision A.02 in:\n")
            target_fullcheck.append("       - z_CH4_duplicate\n")
            target_fullcheck.append("       - CO2_linear_4freqs\n")
            target_fullcheck.append("       - TS_CH3HCH3\n")
            target_fullcheck.append("     * Gaussian 16, Revision C.01 in:\n")
            target_fullcheck.append("       - freq_ok_YYNN\n")
            target_fullcheck.append("       - H_freq\n")
            target_fullcheck.append("       - MeOH_G09\n")
            target_fullcheck.append(
                "x  Different grid_type used in the calculations:\n"
            )
            target_fullcheck.append("     * sg1 in:\n")
            target_fullcheck.append("       - z_CH4_duplicate\n")
            target_fullcheck.append("       - CO2_linear_4freqs\n")
            target_fullcheck.append("       - TS_CH3HCH3\n")
            target_fullcheck.append("     * ultrafine in:\n")
            target_fullcheck.append("       - freq_ok_YYNN\n")
            target_fullcheck.append("       - H_freq\n")
            target_fullcheck.append("     * fine in:\n")
            target_fullcheck.append("       - MeOH_G09\n")
            target_fullcheck.append(
                "x  Different level_of_theory used in the calculations:\n"
            )
            target_fullcheck.append("     * M062X/3-21G in:\n")
            target_fullcheck.append("       - z_CH4_duplicate\n")
            target_fullcheck.append("     * B3LYP/6-31G in:\n")
            target_fullcheck.append("       - CO2_linear_4freqs\n")
            target_fullcheck.append("     * M062X/def2TZVP in:\n")
            target_fullcheck.append("       - freq_ok_YYNN\n")
            target_fullcheck.append("       - H_freq\n")
            target_fullcheck.append("     * PBE1PBE/CC-pVTZ in:\n")
            target_fullcheck.append("       - MeOH_G09\n")
            target_fullcheck.append("     * B3LYP/3-21G in:\n")
            target_fullcheck.append("       - TS_CH3HCH3\n")
            target_fullcheck.append(
                "o  Same dispersion (none) used in all the calculations\n"
            )
            target_fullcheck.append(
                "x  Different solvation used in the calculations:\n"
            )
            target_fullcheck.append("     * gas_phase in:\n")
            target_fullcheck.append("       - z_CH4_duplicate\n")
            target_fullcheck.append("       - CO2_linear_4freqs\n")
            target_fullcheck.append("       - freq_ok_YYNN\n")
            target_fullcheck.append("       - H_freq\n")
            target_fullcheck.append("       - TS_CH3HCH3\n")
            target_fullcheck.append("     * scrf=(solvent=chloroform,pcm) in:\n")
            target_fullcheck.append("       - MeOH_G09")

            outfile = open(
                f"{w_dir_main}/{target_folder}/--QCORR_Fullcheck_Analysis--.dat", "r"
            )
            outlines = outfile.readlines()
            outfile.close()
            assert len(outlines) == len(target_fullcheck) + 1
            # for line in target_fullcheck:
            #     if line in outlines:
            #         pass
            #     else:
            #         assert False

        elif file == "dat":
            outfile = open(f"{path_main}/QCORR-run_1.dat", "r")
            outlines = outfile.readlines()
            outfile.close()

            assert "AQME v" in outlines[0]
            assert "Citation: AQME v" in outlines[1]
            assert "Command line used in AQME: python -m aqme --qcorr" in outlines[3]
            assert "o  Analyzing output files in" in outlines[5]
            assert (
                "Basis_set_error1.log: Termination = other, Error type = atomicbasiserror\n"
                in outlines
            )
            assert (
                "Basis_set_error2.log: Termination = other, Error type = atomicbasiserror\n"
                in outlines
            )
            assert (
                "bpinene_spin_contamin.log: Termination = normal, Error type = spin_contaminated\n"
                in outlines
            )
            assert "Time QCORR:" in outlines[-2]

        elif file == "csv":
            qcorr_stats = pd.read_csv(f"{path_main}/QCORR-run_1-stats.csv")
            assert qcorr_stats["Total files"][0] == 29
            assert qcorr_stats["Normal termination"][0] == 6
            assert qcorr_stats["Single-point calcs"][0] == 3
            assert qcorr_stats["Extra imag. freq."][0] == 4
            assert qcorr_stats["TS with no imag. freq."][0] == 1
            assert qcorr_stats["Freq not converged"][0] == 2
            assert qcorr_stats["Linear mol with wrong n of freqs"][0] == 1
            assert qcorr_stats["SCF error"][0] == 1
            assert qcorr_stats["No data"][0] == 1
            assert qcorr_stats["Basis set error"][0] == 2
            assert qcorr_stats["Other errors"][0] == 5
            assert qcorr_stats["Spin contamination"][0] == 2
            assert qcorr_stats["Duplicates"][0] == 1

        elif file == "check_init":
            assert not path.exists(f"{w_dir_main}/inputs/")

    elif init_folder == "QCORR_1b":
        w_dir_main = f"{path_qcorr}/QCORR_1"
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qcorr",
            "--files",
            f"{w_dir_main}/{file}",
            "--freq_conv",
            "opt=(calcfc,maxstep=5)",
        ]
        if file.split(".")[0] == "freq_conv_YYNN":
            cmd_aqme = [
                "python",
                "-m",
                "aqme",
                "--qcorr",
                "--files",
                f"{w_dir_main}/{file}",
            ]

        elif file.split(".")[0] == "bpinene_spin_contamin":
            cmd_aqme = cmd_aqme + ["--s2_threshold", "50"]

        elif file.split(".")[0] == "z_CH4_duplicate":
            cmd_aqme = cmd_aqme + ["--dup_threshold", "0.000000000001"]

        elif file.split(".")[0] == "imag_freq_no_opt":
            cmd_aqme = cmd_aqme + ["--ifreq_cutoff", "50"]

        elif file.split(".")[0] == "Imag_freq":
            cmd_aqme = cmd_aqme + ["--amplitude_ifreq", "-0.4"]

        elif file.split(".")[0] == "Imag_freq_no_corr":
            cmd_aqme = cmd_aqme + ["--im_freq_input", "None"]

        subprocess.run(cmd_aqme)

        # ensure the output file moves to the right folder
        assert path.exists(f"{w_dir_main}/{target_folder}/{file}")

        if file.split(".")[0] in ["Imag_freq","Imag_freq_no_corr"]:
            # ensure that QCORR applies the correct structural distortions to the errored calcs
            if file.split(".")[0] in "Imag_freq":
                line_2 = "# opt=(calcfc,maxstep=5) freq 3-21g m062x"
                line_6 = "0 1"
                line_8 = "H   0.60103600  -0.59192100  -0.73251000"
                line_10 = "H   0.46752000   1.00584200   0.19478900"
            elif file.split(".")[0] in "Imag_freq_no_corr":
                line_2 = "# opt freq 3-21g m062x"
                line_6 = "0 1"

            outfile = open(
                f'{w_dir_main}/failed/run_1/fixed_QM_inputs/{file.split(".")[0]}.com',
                "r",
            )
            outlines = outfile.readlines()
            outfile.close()

            if file.split(".")[0] in "Imag_freq":
                assert outlines[2].strip() == line_2
                assert outlines[6].strip() == line_6
                assert outlines[8].strip() == line_8
                assert outlines[10].strip() == line_10
            elif file.split(".")[0] in "Imag_freq_no_corr":
                assert outlines[2].strip() == line_2
                assert outlines[6].strip() == line_6

    elif init_folder == "QCORR_1c":
        w_dir_main = f"{path_qcorr}/QCORR_1"
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qcorr",
            "--files",
            f"{w_dir_main}/{file}",
            "--fullcheck",
            "False",
        ]
        subprocess.run(cmd_aqme)

        assert not path.exists(
            f"{w_dir_main}/{target_folder}/--QCORR_Fullcheck_Analysis--.dat"
        )

    elif init_folder == "QCORR_2":
        if command_line is not None:
            cmd_aqme = cmd_aqme + ["--isom_type", "com", "--isom_inputs", w_dir_main]
            subprocess.run(cmd_aqme)

        # ensure the output file moves to the right folder, including the initial COM file
        if file.split(".")[-1].lower() == "log":
            assert path.exists(f"{w_dir_main}/{target_folder}/{file}")
            assert path.exists(f'{w_dir_main}/inputs/{file.split(".")[0]}.com')

        # ensure that no com files were generated
        elif file == "check_fixed":
            assert not path.exists(f"{w_dir_main}/failed/run_1/fixed_QM_inputs/")

    elif init_folder == "QCORR_2b":
        w_dir_main = f"{path_qcorr}/QCORR_2"
        if command_line is not None:
            param_file = f"{w_dir_main}/QCORR_params.yaml"
            cmd_aqme = [
                "python",
                "-m",
                "aqme",
                "--isom_inputs",
                w_dir_main,
                "--varfile",
                param_file,
            ]
            subprocess.run(cmd_aqme)
        # ensure the output file moves to the right folder, including the initial COM file
        if file.split(".")[-1].lower() == "log":
            assert path.exists(f"{w_dir_main}/{target_folder}/{file}")
            assert path.exists(f'{w_dir_main}/inputs/{file.split(".")[0]}.com')

    elif init_folder == "QCORR_2c":
        w_dir_main = f"{path_qcorr}/QCORR_2"
        if command_line is not None:
            cmd_aqme = [
                "python",
                "-m",
                "aqme",
                "--qcorr",
                "--files",
                f"{w_dir_main}/*.log",
                "--freq_conv",
                "opt=(calcfc,maxstep=5)",
            ]
            cmd_aqme = cmd_aqme + [
                "--isom_type",
                "com",
                "--isom_inputs",
                w_dir_main,
                "--vdwfrac",
                "0.01",
                "--covfrac",
                "0.01",
            ]
            subprocess.run(cmd_aqme)

        # ensure the output file moves to the right folder, including the initial COM file
        if file.split(".")[-1].lower() == "log":
            assert path.exists(f"{w_dir_main}/{target_folder}/{file}")

    elif init_folder == "QCORR_5":
        w_dir_QCORR_5 = f"{path_qcorr}/{init_folder}/failed/run_1/fixed_QM_inputs/"

        if command_line is not None:
            cmd_aqme = [
                "python",
                "-m",
                "aqme",
                "--qcorr",
                "--files",
                f"{w_dir_QCORR_5}/*.log",
                "--freq_conv",
                "opt=(calcfc,maxstep=5)",
            ]
            cmd_aqme = cmd_aqme + ["--isom_type", "gjf", "--isom_inputs", w_dir_QCORR_5]
            subprocess.run(cmd_aqme)

        if file.split(".")[0] == "CH4":
            assert path.exists(f"{w_dir_main}/{target_folder}/{file}")
            assert path.exists(
                f'{w_dir_main}/{target_folder}/json_files/{file.split(".")[0]}.json'
            )
        else:
            assert path.exists(
                f"{w_dir_main}/{target_folder}/run_2/isomerization/{file}"
            )
            assert not path.exists(
                f"{w_dir_main}/{target_folder}/run_1/fixed_QM_inputs/failed"
            )

    elif init_folder == "QCORR_5b":
        w_dir_QCORR_5b = f"{path_qcorr}/QCORR_5"

        if command_line is not None:
            cmd_aqme = [
                "python",
                "-m",
                "aqme",
                "--qcorr",
                "--files",
                f"{w_dir_QCORR_5b}/failed/run_1/fixed_QM_inputs/*.log",
                "--freq_conv",
                "opt=(calcfc,maxstep=5)",
            ]
            cmd_aqme = cmd_aqme + [
                "--isom_type",
                "gjf",
                "--isom_inputs",
                w_dir_QCORR_5b + "/failed/run_1/fixed_QM_inputs",
            ]
            subprocess.run(cmd_aqme)

        # the energy of the json file was changed to avoid the duplicate filter
        if file.split(".")[0] == "CH4":
            assert path.exists(f"{w_dir_QCORR_5b}/{target_folder}/{file}")
            assert path.exists(
                f'{w_dir_QCORR_5b}/{target_folder}/json_files/{file.split(".")[0]}.json'
            )
        else:
            assert path.exists(
                f"{w_dir_QCORR_5b}/{target_folder}/run_2/isomerization/{file}"
            )
            assert not path.exists(
                f"{w_dir_QCORR_5b}/{target_folder}/run_1/fixed_QM_inputs/failed"
            )

    elif init_folder == "QCORR_1d":
        # go to the parent folder and run the program
        w_dir_main = f"{path_qcorr}/QCORR_1"
        os.chdir(w_dir_main)
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qcorr",
            "--files",
            "CH4.log",
            "--freq_conv",
            "opt=(calcfc,maxstep=5)",
        ]
        subprocess.run(cmd_aqme)

        # ensure the output file moves to the right folder
        assert path.exists(f"{w_dir_main}/{target_folder}/{file}")

        # ensure that the com file is not generated
        assert not path.exists(
            f'{w_dir_main}/failed/run_1/fixed_QM_inputs/{file.split(".")[0]}.com'
        )

    elif init_folder == "QCORR_6":
        subprocess.run(cmd_aqme)

        # ensure the output file moves to the right folder
        assert path.exists(f"{w_dir_main}/{target_folder}/{file}")

    elif init_folder == "QCORR_7":
        w_dir_main = f"{path_qcorr}/QCORR_7"
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qcorr",
            "--files",
            f"{w_dir_main}/{file}",
        ]
        subprocess.run(cmd_aqme)

        # ensure the output file moves to the right folder
        assert path.exists(f"{w_dir_main}/{target_folder}/{file}")

        if file == "orca_imag_freq.out":
            assert path.exists(
                f'{w_dir_main}/failed/run_1/fixed_QM_inputs/{file.split(".")[0]}.inp'
            )

            # ensure that QCORR applies the correct structural distortions to the errored calcs
            outfile = open(f'{w_dir_main}/failed/run_1/fixed_QM_inputs/{file.split(".")[0]}.inp',"r",)
            outlines = outfile.readlines()
            outfile.close()

            assert 'MaxStep 0.05' in outlines[9]
            assert 'C  -4.32349420   1.79733020  -0.42830100' in outlines[12]


    # leave the folders as they were initially to run a different batch of tests
    elif restore_folder:
        os.chdir(path_main)
        shutil.rmtree(f"{path_main}/Example_workflows")
        filepath = Path(f"{path_main}/Example_workflows_original")
        filepath.rename(f"{path_main}/Example_workflows")
        # remove DAT and CSV files generated by QCORR
        dat_files = glob.glob("*.dat")
        for dat_file in dat_files:
            if "QCORR" in dat_file:
                os.remove(dat_file)
        dat_files = glob.glob("*.csv")
        for dat_file in dat_files:
            if "QCORR" in dat_file:
                os.remove(dat_file)

# two tests 1) check that frozen atoms are in success json file
# 2) check that frozen flags are included in generate .com for failed 