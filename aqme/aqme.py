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
from pathlib import Path
from aqme.argument_parser import parser_args

from aqme.mainf import (
    csearch_main,
    geom_rules_main,
    qprep_main,
    load_from_yaml,
    cmin_main,
)
from aqme.utils import Logger
from aqme.qcorr import qcorr, json2input
from aqme.qdescp import qdescp


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

    # QCORR
    if args.qcorr:
        qcorr(
            qm_files=args.qm_files,
            w_dir_main=args.w_dir_main,
            dup_threshold=args.dup_threshold,
            mem=args.mem,
            nprocs=args.nprocs,
            chk=args.chk,
            qm_input=args.qm_input,
            s2_threshold=args.s2_threshold,
            isom=args.isom,
            isom_inputs=args.isom_inputs,
            vdwfrac=args.vdwfrac,
            covfrac=args.covfrac,
            bs_gen=args.bs_gen,
            bs=args.bs,
            gen_atoms=args.gen_atoms,
            qm_end=args.qm_end,
            amplitude_ifreq=args.amplitude_ifreq,
            freq_conv=args.freq_conv,
            ifreq_cutoff=args.ifreq_cutoff,
            fullcheck=args.fullcheck,
            author=args.author,
            program=args.program,
            varfile=None,
        )

    if args.json2input:
        json2input(
            json_files=args.json_files,
            w_dir_main=args.w_dir_main,
            destination=args.destination,
            suffix=args.suffix,
            charge=args.charge,
            mult=args.mult,
            mem=args.mem,
            nprocs=args.nprocs,
            chk=args.chk,
            qm_input=args.qm_input,
            bs_gen=args.bs_gen,
            bs=args.bs,
            gen_atoms=args.gen_atoms,
            qm_end=args.qm_end,
            program=args.program,
        )

    # qdescp
    if args.qdescp in ["geometricdescp", "nmr", "dbstep", "nbo"]:
        qdescp(
            w_dir_main=args.w_dir_main,
            destination=args.destination,
            qm_files=args.qm_files,
            json_files=args.json_files,
            task=args.qdescp,
            varfile=None,
        )

    # if args.QPRED == "nmr":
    #     nmr_main(args, log_overall, w_dir_initial)
    # if args.QPRED == "energy":
    #     energy_main(args, log_overall, w_dir_initial)
    # if args.QPRED == "dbstep":
    #     dbstep_par_main(args, log_overall, w_dir_initial)
    # if args.QPRED == "nics":
    #     nics_par_main(args, log_overall, w_dir_initial)
    # if args.QPRED == "cclib-json":
    #     cclib_main(args, log_overall, w_dir_initial)
    # os.chdir(w_dir_initial)
    #
    # # QSTAT
    # if args.QSTAT == "descp":
    #     geom_par_main(args, log_overall, w_dir_initial)
    # if args.QSTAT == "graph":
    #     graph_main(args, log_overall, w_dir_initial)
    # os.chdir(w_dir_initial)
    #
    log_overall.finalize()

    out_data_file = Path(f"aqme_output.dat")
    if out_data_file.exists():
        out_data_file.replace(f"aqme_{args.output_name}.dat")


if __name__ == "__main__":
    main()
