#!/usr/bin/env python

#########################################################################################.
#########################################################################################
###                                                                                   ###
###  AQME is a tool that allows to carry out automated:                               ###
###  (1) Conformational searches and creation of COM files using RDKit, xTB and ANI   ###
###  (2) LOG file processing (detects imaginary freqs and error terminations          ###
###      and creates new COM files)                                                   ###
###  (3) Use LOG files to create new COM files with new keywords (i.e. single-point   ###
###      corrections after geometry optimization)                                     ###
###                                                                                   ###
#########################################################################################
###                                                                                   ###
###  Version: v1.0.1, Release date: 22-May-2020                                       ###
###                                                                                   ###
#########################################################################################
###                                                                                   ###
###  Authors: Shree Sowndarya S. V., Juan V. Alegre Requena, Robert S. Paton          ###
###                                                                                   ###
###  Please, report any bugs or suggestions to:                                       ###
###  svss@colostate.edu or juanvi89@hotmail.com                                       ###
###                                                                                   ###
#########################################################################################
#########################################################################################.

import os
import time
import json
import pandas as pd
from pathlib import Path
from aqme.argument_parser import parser_args

from aqme.mainf import (
    csearch_main,
    geom_rules_main,
    qprep_main,
    move_sdf_main,
    graph_main,
    geom_par_main,
    nmr_main,
    energy_main,
    load_from_yaml,
    dbstep_par_main,
    nics_par_main,
    cclib_main,
    cmin_main,
)
from aqme.utils import Logger, get_filenames
from aqme.qcorr import qcorr, check_for_final_folder
from aqme.qprep import qprep


def main():
    # working directory and arguments
    w_dir_initial = os.getcwd()
    base_path = Path(w_dir_initial)

    args = parser_args()

    log_overall = Logger("aqme", args.output_name)
    # if needed to load from a yaml file
    args, log = load_from_yaml(args, log_overall)

    name = args.input.split(".")[0]

    # CSEARCH AND CMIN
    if args.CSEARCH in [
        "rdkit",
        "summ",
        "fullmonte",
        "crest",
    ]:  # RAUL: Is there any other posibilities or just None?
        start_time_overall = time.time()
        csearch_dup_data = csearch_main(w_dir_initial, args, log_overall)
        if args.time:
            elapsed_time = round(time.time() - start_time_overall, 2)
            log_overall.write(
                f"\n All molecules execution time CSEARCH: {elapsed_time} seconds"
            )
        os.chdir(w_dir_initial)
        if args.CMIN is None:
            csearch_csv_folder = base_path.joinpath("CSEARCH/csv_files")
            csearch_csv_folder.mkdir(exist_ok=True)
            csearch_csv_file = csearch_csv_folder.joinpath(f"{name}-CSEARCH-Data.csv")
            csearch_dup_data.to_csv(csearch_csv_file, index=False)

    # Separating CMIN
    if args.CSEARCH != None and args.CMIN in ["xtb", "ani"]:
        cmin_dup_data = cmin_main(w_dir_initial, args, log_overall, csearch_dup_data)
        if args.time:
            elapsed_time = round(time.time() - start_time_overall, 2)
            log_overall.write(
                f"\n All molecules execution time CMIN: {elapsed_time} seconds"
            )
        os.chdir(w_dir_initial)
        cmin_csv_folder = base_path.joinpath("CMIN/csv_files")
        cmin_csv_folder.mkdir(exist_ok=True)
        cmin_csv_file = cmin_csv_folder.joinpath(f"{name}-CMIN-Data.csv")
        cmin_dup_data.to_csv(cmin_csv_file, index=False)

    # applying rules to discard certain conformers based on rules that the user define
    if len(args.geom_rules) >= 1:
        geom_rules_active = True
        if args.qcorr == "gaussian":
            geom_rules_active = False
        geom_rules_main(args, log_overall, geom_rules_active)
        os.chdir(w_dir_initial)

    # QPREP
    if args.QPREP == "gaussian" or args.QPREP == "orca" or args.QPREP == "turbomole":
        qprep_main(w_dir_initial, args, log_overall)
        os.chdir(w_dir_initial)

    # if args.CSEARCH in ['rdkit', 'summ', 'fullmonte','crest'] or args.QPREP is not None:
    #     # moving files after compute and/or write_gauss
    #     move_sdf_main(args)
    #     os.chdir(w_dir_initial)


    #QCORR
    if args.qcorr:
        log_overall.write("\no  Writing analysis of output files in respective folders\n")
        qm_files = get_filenames('output',None)
        df_qcorr_success,df_qcorr_error = pd.DataFrame(),pd.DataFrame()

        # detects round of optimizations
        round_num = check_for_final_folder(w_dir_initial)

        # runs the qcorr module, which organizes the files and returns all the necessary information
        # as a json file
        if args.qcorr_json == "":
            df_qcorr_success, df_qcorr_error = qcorr(
                qm_files=qm_files, round_num=round_num
            )
        else:
            with open(f"{args.qcorr_json}_success") as f:
                df_qcorr_success = json.load(f)

        if not df_qcorr_error.empty():
            for i, file_name in enumerate(df_qcorr_error["file_name"]):
                # creates input files to fix imaginary freqs and not normal terminations
                if df_qcorr_error["termination"][i] != "normal" or df_qcorr_error[
                    "errortype"
                ][i] not in [None, "isomerization"]:
                    dir_for_analysis = (
                        w_dir_initial
                        + "/fixed_QM_input_files/run_"
                        + str(round_num + 1)
                    )
                    if not os.path.isdir(dir_for_analysis):
                        os.makedirs(dir_for_analysis)
                    log.write(
                        f"-> Creating new gaussian input file for {file_name} in {dir_for_analysis}"
                    )
                    if not args.nocom:
                        qprep(
                            destination=dir_for_analysis,
                            molecule=file_name,
                            charge=df_qcorr_error["charge"][i],
                            mult=df_qcorr_error["mult"][i],
                            atom_types=df_qcorr_error["atom_types"][i],
                            gen_atoms=args.gen_atoms,
                            bs_gen=args.genecp_bs,
                            bs=args.bs,
                            program=args.program,
                            cartesians=df_qcorr_error["cartesians"][i],
                            qm_keywords=df_qcorr_error["keywords_line"][i],
                            qm_end=args.qm_input_end_sp,
                        )

        if not df_qcorr_success.empty() and args.sp != None:
            # creates input files for single-point energy calcs
            for i, file_name in enumerate(df_qcorr_success["file_name"]):
                if (
                    df_qcorr_success["termination"][i] == "normal"
                    and df_qcorr_success["errortype"][i] == None
                ):
                    dir_for_sp = w_dir_initial + "/single_point_input_files"
                    if not os.path.isdir(dir_for_sp):
                        os.makedirs(dir_for_sp)
                    if args.charge_sp != "None":
                        charge = args.charge_sp
                    if args.mult_sp != "None":
                        mult = args.mult_sp
                    if args.suffix != "None":
                        file_name = file_name + "_" + args.suffix
                    log.write(
                        f"-> Creating new gaussian input file for {file_name} in {dir_for_sp}"
                    )
                    qprep(
                        destination=dir_for_analysis,
                        molecule=file_name,
                        charge=df_qcorr_success["charge"][i],
                        mult=df_qcorr_success["mult"][i],
                        atom_types=df_qcorr_success["atom_types"][i],
                        gen_atoms=args.gen_atoms,
                        bs_gen=args.gen_bs_sp,
                        bs=args.bs_sp,
                        program=args.program_sp,
                        cartesians=df_qcorr_success["cartesians"][i],
                        qm_keywords=args.qm_input_sp,
                        qm_end=args.qm_input_end_sp,
                    )

    # QPRED
    if args.QPRED == "nmr":
        nmr_main(args, log_overall, w_dir_initial)
    if args.QPRED == "energy":
        energy_main(args, log_overall, w_dir_initial)
    if args.QPRED == "dbstep":
        dbstep_par_main(args, log_overall, w_dir_initial)
    if args.QPRED == "nics":
        nics_par_main(args, log_overall, w_dir_initial)
    if args.QPRED == "cclib-json":
        cclib_main(args, log_overall, w_dir_initial)
    os.chdir(w_dir_initial)

    # QSTAT
    if args.QSTAT == "descp":
        geom_par_main(args, log_overall, w_dir_initial)
    if args.QSTAT == "graph":
        graph_main(args, log_overall, w_dir_initial)
    os.chdir(w_dir_initial)

    log_overall.finalize()

    out_data_file = Path(f"aqme_output.dat")
    if out_data_file.exists():
        out_data_file.replace(f"aqme_{args.output_name}.dat")


if __name__ == "__main__":
    main()
