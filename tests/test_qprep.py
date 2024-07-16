#!/usr/bin/env python

######################################################.
# 	        Testing QPREP with pytest 	             #
######################################################.

import os
from os import path
import glob
import pytest
import shutil
import subprocess
from pathlib import Path

# saves the working directory
path_main = os.getcwd()
path_qprep = path_main + "/Example_workflows/QPREP_generating_input_files"

# QPREP tests
@pytest.mark.parametrize(
    "test_type, init_folder, target_folder, restore_folder",
    [
        # QPREP analysis tests with standard options
        ("com_gen", "json_files", "com_files", False),  # test json inputs
        ("com_gen", "log_files", "com_files", False),  # test log inputs
        ("com_gen", "sdf_files", "com_files", False),  # test sdf inputs
        ("com_gen", "xyz_files", "com_files", False),  # test xyz inputs
        ("com_gen", "pdb_files", "com_files", False),  # test pdb inputs
        (
            "charge_mult",
            "json_files",
            "com_files",
            False,
        ),  # test json inputs with defined charge and mult
        (
            "charge_mult",
            "log_files",
            "com_files",
            False,
        ),  # test log inputs with defined charge and mult
        ("com_gen", "p_print", "com_files", False),  # test #p in qm_input (replaces previous files)
        (
            "charge_mult",
            "sdf_files",
            "com_files",
            False,
        ),  # test sdf inputs with defined charge and mult
        (
            "charge_mult",
            "xyz_files",
            "com_files",
            False,
        ),  # test xyz inputs with defined charge and mult
        # ('charge_mult', 'pdb_files', 'com_files', False), # test pdb inputs with defined charge and mult, ONGOING
        # creating ORCA input files
        ("orca", "log_files", "orca_files", False),  # test ORCA input files
        # changing memory, nprocs, chk, charge, mult and suffix
        ("input_params", "json_files", "params_files", False),  # test multiple params
        # changing chk_path
        ("chk_path", "json_files", "params_files", False),  # test chk_path
        # using genecp and final lines
        ("final_line", "json_files", "gen_final_files", False),  # test final line
        (
            "genecp_and_final",
            "json_files",
            "gen_final_files",
            False,
        ),  # test genecp and final line
        ("gen", "json_files", "gen_final_files", False),  # test gen
        # from YAML file (varfile=XX)
        ("yaml", "json_files", "yaml_files", False),  # test for yaml files
        # calling AQME from the parent folder where the files are located (omitting the w_dir_main keyword)
        ("parent", "json_files", "parent_run", False),  # test for yaml files
        # calling AQME with no destination keyword (creates a folder called QCALC by default)
        ("no_dest", "json_files", "QCALC", False),  # test for yaml files
        # check if the DAT file is generated
        ("dat_file", "json_files", "com_files", False),  # test for DAT file
        (
            None,
            None,
            None,
            True,
        ),  # reset the initial folder to start another set of tests
    ],
)
def test_QPREP_analysis(test_type, init_folder, target_folder, restore_folder):
    # copy the test folders
    if not path.exists(f"{path_main}/Example_workflows_original"):
        shutil.copytree(
            f"{path_main}/Example_workflows", f"{path_main}/Example_workflows_original"
        )

    # runs the program with the different tests
    if init_folder == 'p_print':
        qm_input = "p wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)"
        line_2 = "#p wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)"
        init_folder = 'sdf_files'

    else:
        qm_input = "wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)"
        line_2 = "# wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)"

    w_dir_main = f"{path_qprep}/{init_folder}"
    destination = f"{path_qprep}/{init_folder}/{target_folder}"

    if test_type == "com_gen":
        if init_folder == "json_files":
            files = "*.json"
            files_assert = ["CH4.com", "MeOH_NMR.com"]
        elif init_folder == "log_files":
            files = "*.log"
            files_assert = ["CH4.com", "H_freq.com"]
        elif init_folder == "sdf_files":
            files = "*.sdf"
            files_assert = ["quinine_rdkit_conf_1.com", "quinine_rdkit_conf_10.com"]
        elif init_folder == "xyz_files":
            files = "*.xyz"
            files_assert = ["Int-I_conf_1.com", "Int-I_conf_3.com"]
        elif init_folder == "pdb_files":
            files = "*.pdb"
            files_assert = ["7ac8_chainsEF_conf_1.com"]

        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            f"{w_dir_main}/{files}",
            "--program",
            "gaussian",
            "--qm_input",
            qm_input,
        ]
        subprocess.run(cmd_aqme)

        if init_folder in ["xyz_files", "pdb_files"]:
            # make sure all the generated SDF files are removed
            assert len(glob.glob(f"{w_dir_main}/*.sdf")) == 0

        for file in files_assert:
            assert path.exists(f'{destination}/{file.split(".")[0]}.com')

            outfile = open(f'{destination}/{file.split(".")[0]}.com', "r")
            outlines = outfile.readlines()
            outfile.close()
 
            if file.split(".")[0] == "H_freq":
                line_6 = "0 2"
                line_7 = "H   0.00000000   0.00000000   0.00000000"

                assert outlines[7].strip() == line_7

            else:
                # ensure that QPREP applies the correct structural distortions to the errored calcs
                line_6 = "0 1"

                if file.split(".")[0] == "CH4":
                    line_8 = "H   0.45703600  -0.46392100  -0.87651000"
                    line_10 = "H   0.29552000   1.04984200   0.05078900"

                elif file.split(".")[0] == "MeOH_NMR":
                    line_8 = "H  -1.08592900   0.98379700  -0.00000100"
                    line_10 = "H  -1.02538000  -0.54454000  -0.88941000"

                elif file.split(".")[0] == "quinine_rdkit_conf_1":
                    line_8 = "O   2.93580000   2.55850000   2.17990000"
                    line_10 = "C   2.26850000   1.23230000   0.38520000"

                elif file.split(".")[0] == "quinine_rdkit_conf_10":
                    line_8 = "O   3.23820000   2.67370000  -1.56720000"
                    line_10 = "C   2.35510000   0.65880000  -0.80190000"
                    assert (
                        len(glob.glob(f"{destination}/quinine_rdkit_conf_*.com")) == 10
                    )

                elif file.split(".")[0] == "Int-I_conf_1":
                    line_8 = "C   3.74440000  -1.75670000   0.37800000"
                    line_10 = "H   3.38180000  -1.72450000   1.40790000"

                elif file.split(".")[0] == "Int-I_conf_3":
                    line_8 = "C   3.59790000  -2.21510000   0.62110000"
                    line_10 = "H   2.68030000  -2.43910000   1.16970000"
                    assert len(glob.glob(f"{destination}/Int-I_conf_*.com")) == 3

                elif file.split(".")[0] == "7ac8_chainsEF_conf_1":
                    line_8 = "C   5.61200000 -38.15900000   9.29600000"
                    line_10 = "O   5.88600000 -40.54700000   9.19200000"

                assert outlines[8].strip() == line_8
                assert outlines[10].strip() == line_10
        
            assert outlines[2].strip() == line_2
            if init_folder == "pdb_files":
                line_6 = "-2 1"
                assert outlines[6].strip() == line_6
            else:
                assert outlines[6].strip() == line_6

        

    elif test_type == "charge_mult":
        if init_folder == "json_files":
            files_assert = ["MeOH_NMR_charged.com"]
        elif init_folder == "log_files":
            files_assert = ["CH4_charged.com"]
        elif init_folder == "sdf_files":
            files_assert = [
                "quinine_rdkit_charged_conf_1.com",
                "quinine_rdkit_charged_conf_10.com",
            ]
        elif init_folder == "xyz_files":
            files_assert = ["Int-I_charged_conf_1.com", "Int-I_charged_conf_3.com"]

        for file in files_assert:
            outfile = open(f'{destination}/{file.split(".")[0]}.com', "r")
            outlines = outfile.readlines()
            outfile.close()

            if file.split(".")[0] == "MeOH_NMR_charged":
                line_6 = "2 3"
            elif file.split(".")[0] == "CH4_charged":
                line_6 = "2 3"
            elif file.split(".")[0] == "quinine_rdkit_charged_conf_1":
                line_6 = "2 3"
            elif file.split(".")[0] == "quinine_rdkit_charged_conf_10":
                line_6 = "4 5"
            elif file.split(".")[0] == "Int-I_charged_conf_1":
                line_6 = "2 3"
            elif file.split(".")[0] == "Int-I_charged_conf_3":
                line_6 = "4 5"

            assert outlines[6].strip() == line_6

    elif test_type == "orca":
        ORCA_SP = "Extrapolate(2/3,cc) def2/J cc-pVTZ/C DLPNO-CCSD(T) NormalPNO TightSCF RIJCOSX GridX7\n"
        ORCA_SP += "%cpcm\n"
        ORCA_SP += "smd true\n"
        ORCA_SP += 'SMDsolvent "CH2Cl2"\n'
        ORCA_SP += "end"

        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            f"{w_dir_main}/CH4.log",
            "--program",
            "orca",
            "--qm_input",
            ORCA_SP,
        ]
        subprocess.run(cmd_aqme)

        outfile = open(f"{destination}/CH4.inp", "r")
        outlines = outfile.readlines()
        outfile.close()

        line_0 = "# CH4"
        line_1 = "%maxcore 16000"
        line_2 = "%pal nprocs 8 end"
        line_3 = "! Extrapolate(2/3,cc) def2/J cc-pVTZ/C DLPNO-CCSD(T) NormalPNO TightSCF RIJCOSX GridX7"
        line_4 = "%cpcm"
        line_5 = "smd true"
        line_6 = 'SMDsolvent "CH2Cl2"'
        line_7 = "end"
        line_8 = "* xyz 0 1"
        line_9 = "C  -0.00037100   0.00001500   0.00006600"
        line_10 = "H   0.45703600  -0.46392100  -0.87651000"

        assert len(glob.glob(f"{destination}/*.inp")) == 1
        assert (
            glob.glob(f"{destination}/*.inp")[0].replace("\\", "/").split("/")[-1]
            == "CH4.inp"
        )
        assert outlines[0].strip() == line_0
        assert outlines[1].strip() == line_1
        assert outlines[2].strip() == line_2
        assert outlines[3].strip() == line_3
        assert outlines[4].strip() == line_4
        assert outlines[5].strip() == line_5
        assert outlines[6].strip() == line_6
        assert outlines[7].strip() == line_7
        assert outlines[8].strip() == line_8
        assert outlines[9].strip() == line_9
        assert outlines[10].strip() == line_10

    elif test_type == "input_params":

        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            f"{w_dir_main}/CH4.json",
            "--program",
            "gaussian",
            "--qm_input",
            "wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)",
            "--charge",
            "2",
            "--mult",
            "3",
            "--suffix",
            "params",
            "--chk",
            "--mem",
            "100GB",
            "--nprocs",
            "32",
        ]
        subprocess.run(cmd_aqme)

        outfile = open(f"{destination}/CH4_params.com", "r")
        outlines = outfile.readlines()
        outfile.close()

        line_0 = "%chk=CH4_params.chk"
        line_1 = "%nprocshared=32"
        line_2 = "%mem=100GB"
        line_5 = "CH4_params"
        line_7 = "2 3"

        assert outlines[0].strip() == line_0
        assert outlines[1].strip() == line_1
        assert outlines[2].strip() == line_2
        assert outlines[5].strip() == line_5
        assert outlines[7].strip() == line_7

    elif test_type == "chk_path":
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            f"{w_dir_main}/CH4.json",
            "--program",
            "gaussian",
            "--qm_input",
            "wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)",
            "--suffix",
            "chk_path",
            "--chk",
            "--chk_path",
            "test/PATH/CH4_chk_path.chk"
        ]
        subprocess.run(cmd_aqme)

        outfile = open(f"{destination}/CH4_chk_path.com", "r")
        outlines = outfile.readlines()
        outfile.close()

        line_0 = "%chk=test/PATH/CH4_chk_path.chk"
        assert outlines[0].strip() == line_0

    elif test_type == "final_line":

        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            f"{w_dir_main}/CH4.json",
            "--program",
            "gaussian",
            "--qm_input",
            "pop=(nbo6read,savenbos) wb97xd/def2svp",
            "--qm_end",
            "$nbo bndidx $end",
            "--suffix",
            "final",
        ]
        subprocess.run(cmd_aqme)

        outfile = open(f"{destination}/CH4_final.com", "r")
        outlines = outfile.readlines()
        outfile.close()

        assert outlines[2].strip() == "# pop=(nbo6read,savenbos) wb97xd/def2svp"
        assert outlines[-2].strip() == "$nbo bndidx $end"
        assert outlines[-1] == "\n"

    elif test_type == "genecp_and_final":
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            f"{w_dir_main}/CH4.json",
            "--program",
            "gaussian",
            "--qm_input",
            "pop=(nbo6read,savenbos) wb97xd/genecp",
            "--qm_end",
            "$nbo bndidx $end",
            "--suffix",
            "genecp_and_final",
            "--gen_atoms",
            "['C']",
            "--bs_gen",
            "LANL2TZ",
            "--bs_nogen",
            "LANL2DZ",
        ]
        subprocess.run(cmd_aqme)

        outfile = open(f"{destination}/CH4_genecp_and_final.com", "r")
        outlines = outfile.readlines()
        outfile.close()

        assert outlines[2].strip() == "# pop=(nbo6read,savenbos) wb97xd/genecp"
        assert outlines[13].strip() == "$nbo bndidx $end"
        assert outlines[14].strip() == ""
        assert outlines[15].strip() == "H 0"
        assert outlines[16].strip() == "LANL2DZ"
        assert outlines[17].strip() == "****"
        assert outlines[18].strip() == "C 0"
        assert outlines[19].strip() == "LANL2TZ"
        assert outlines[20].strip() == "****"
        assert outlines[21].strip() == ""
        assert outlines[22].strip() == "C 0"
        assert outlines[23].strip() == "LANL2TZ"
        assert outlines[24] == "\n"

    elif test_type == "gen":
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            f"{w_dir_main}/CH4.json",
            "--program",
            "gaussian",
            "--qm_input",
            "wb97xd/gen",
            "--suffix",
            "gen",
            "--gen_atoms",
            "['C']",
            "--bs_gen",
            "LANL2TZ",
            "--bs_nogen",
            "LANL2DZ",
        ]
        subprocess.run(cmd_aqme)

        outfile = open(f"{destination}/CH4_gen.com", "r")
        outlines = outfile.readlines()
        outfile.close()

        assert outlines[2].strip() == "# wb97xd/gen"
        assert outlines[13].strip() == "H 0"
        assert outlines[14].strip() == "LANL2DZ"
        assert outlines[15].strip() == "****"
        assert outlines[16].strip() == "C 0"
        assert outlines[17].strip() == "LANL2TZ"
        assert outlines[18].strip() == "****"
        assert outlines[19] == "\n"

    elif test_type == "yaml":
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--destination",
            destination,
            "--varfile",
            f"{w_dir_main}/pytest_varfile.yaml",
        ]
        subprocess.run(cmd_aqme)

        outfile = open(f"{destination}/CH4_varfile.com", "r")
        outlines = outfile.readlines()
        outfile.close()

        line_2 = "# wb97xd/lanl2dz"
        line_6 = "2 3"
        line_8 = "H   0.45703600  -0.46392100  -0.87651000"
        line_10 = "H   0.29552000   1.04984200   0.05078900"

        assert outlines[2].strip() == line_2
        assert outlines[6].strip() == line_6
        assert outlines[8].strip() == line_8
        assert outlines[10].strip() == line_10

    elif test_type == "parent":
        os.chdir(w_dir_main)
        file = "CH4.json"
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--destination",
            destination,
            "--files",
            file,
            "--program",
            "gaussian",
            "--qm_input",
            "wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)",
        ]
        subprocess.run(cmd_aqme)

        assert path.exists(f'{destination}/{file.split(".")[0]}.com')

        outfile = open(f'{destination}/{file.split(".")[0]}.com', "r")
        outlines = outfile.readlines()
        outfile.close()

        line_2 = "# wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)"
        line_6 = "0 1"
        line_8 = "H   0.45703600  -0.46392100  -0.87651000"
        line_10 = "H   0.29552000   1.04984200   0.05078900"

        assert outlines[8].strip() == line_8
        assert outlines[10].strip() == line_10
        assert outlines[2].strip() == line_2
        assert outlines[6].strip() == line_6

        # remove the DAT file from the QPREP run
        dat_files = glob.glob("*.dat")
        for dat_file in dat_files:
            if "QPREP" in dat_file:
                os.remove(dat_file)
        os.chdir(path_main)

    elif test_type == "no_dest":
        file = "CH4.json"
        cmd_aqme = [
            "python",
            "-m",
            "aqme",
            "--qprep",
            "--files",
            f"{w_dir_main}/{file}",
            "--program",
            "gaussian",
            "--qm_input",
            "wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)",
        ]
        subprocess.run(cmd_aqme)

        assert path.exists(f'{os.getcwd()}/QCALC/{file.split(".")[0]}.com')

        outfile = open(f'{os.getcwd()}/QCALC/{file.split(".")[0]}.com', "r")
        outlines = outfile.readlines()
        outfile.close()

        line_2 = "# wb97xd/lanl2dz scrf=(smd,solvent=acetonitrile)"
        line_6 = "0 1"
        line_8 = "H   0.45703600  -0.46392100  -0.87651000"
        line_10 = "H   0.29552000   1.04984200   0.05078900"

        assert outlines[8].strip() == line_8
        assert outlines[10].strip() == line_10
        assert outlines[2].strip() == line_2
        assert outlines[6].strip() == line_6

    elif test_type == "dat_file":
        outfile = open(f"{path_main}/QPREP_data.dat", "r")
        outlines = outfile.readlines()
        outfile.close()

        assert "AQME v" in outlines[0]
        assert "Citation: AQME v" in outlines[1]
        assert "Time QPREP:" in outlines[-2]

    # leave the folders as they were initially to run a different batch of tests
    if restore_folder:
        os.chdir(path_main)
        shutil.rmtree(f"{path_main}/Example_workflows")
        filepath = Path(f"{path_main}/Example_workflows_original")
        filepath.rename(f"{path_main}/Example_workflows")
        # remove dat files generated by QPREP
        dat_files = glob.glob("*.dat")
        for dat_file in dat_files:
            if "QPREP" in dat_file:
                os.remove(dat_file)
