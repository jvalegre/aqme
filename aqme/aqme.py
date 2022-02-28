#!/usr/bin/env python

###########################################################################################.
###########################################################################################
###                                                                                     ###
###  AQME is a tool that allows to carry out automated:                                 ###
###  (CSEARCH) Conformational searches and creation of COM files using RDKit and CREST  ###
###  (CMIN) Geometry refinement of initial conformers with xTB and ANI                  ###
###  (QCORR) Out put file processing from QM calculations and automated issue fixing,   ###
###  including imaginary freqs, spin contamination, isomerization issues and            ###
###  error terminations, among others                                                   ###
###  (QPREP) Use QM output (.log or .out) and json files to create new COM files with   ###
###  new keywords (i.e. for single-point corrections after geometry optimization)       ###                               ###
###                                                                                     ###
###########################################################################################
###                                                                                     ###
###  Version: v0.2, Release date: 28-Feb-2022                                           ###
###                                                                                     ###
###########################################################################################
###                                                                                     ###
###  Authors: Shree Sowndarya S. V., Juan V. Alegre Requena                             ###
###                                                                                     ###
###  Please, report any bugs or suggestions to:                                         ###
###  svss@colostate.edu or juanvi89@hotmail.com                                         ###
###                                                                                     ###
###########################################################################################
###########################################################################################.

import os
import time
from pathlib import Path
from aqme.argument_parser import parser_args
from aqme.mainf import (
    csearch_main,
    geom_rules_main,
    qprep_main,
    cmin_main)
from aqme.utils import (Logger, 
    load_from_yaml)
from aqme.qcorr import qcorr


def main():

    args = parser_args()   

    # working directory and arguments
    w_dir_main = Path(args.w_dir_main)

    log_overall = Logger("aqme", args.output_name)

    # if needed to load from a yaml file
    args,_ = load_from_yaml(args, log_overall)

    name = args.input.split(".")[0]

    # CSEARCH AND CMIN
    if args.CSEARCH in [
        "rdkit",
        "summ",
        "fullmonte",
        "crest",
    ]:
        start_time_overall = time.time()
        csearch_dup_data = csearch_main(w_dir_main, args, log_overall)
        if args.time:
            elapsed_time = round(time.time() - start_time_overall, 2)
            log_overall.write(
                f"\n All molecules execution time CSEARCH: {elapsed_time} seconds"
            )
        os.chdir(w_dir_main)
        if args.CMIN is None:
            csearch_csv_folder = w_dir_main.joinpath("CSEARCH/csv_files")
            csearch_csv_folder.mkdir(exist_ok=True)
            csearch_csv_file = csearch_csv_folder.joinpath(f"{name}-CSEARCH-Data.csv")
            csearch_dup_data.to_csv(csearch_csv_file, index=False)

    # Separating CMIN
    if args.CSEARCH != None and args.CMIN in ["xtb", "ani"]:
        cmin_dup_data = cmin_main(w_dir_main, args, log_overall, csearch_dup_data)
        if args.time:
            elapsed_time = round(time.time() - start_time_overall, 2)
            log_overall.write(
                f"\n All molecules execution time CMIN: {elapsed_time} seconds"
            )
        os.chdir(w_dir_main)
        cmin_csv_folder = w_dir_main.joinpath("CMIN/csv_files")
        cmin_csv_folder.mkdir(exist_ok=True)
        cmin_csv_file = cmin_csv_folder.joinpath(f"{name}-CMIN-Data.csv")
        cmin_dup_data.to_csv(cmin_csv_file, index=False)

    # applying rules to discard certain conformers based on rules that the user define
    if len(args.geom_rules) >= 1:
        geom_rules_active = True
        if args.qcorr:
            geom_rules_active = False
        geom_rules_main(args, log_overall, geom_rules_active)
        os.chdir(w_dir_main)

    # QPREP
    if args.qprep:
        qprep_main(w_dir_main, args, log_overall)

    # QCORR
    if args.qcorr:
        qcorr(
            files=args.files,
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
            program=args.program,
            varfile=None
        )

    # # qdescp
    # if args.qdescp in ["geometricdescp", "nmr", "dbstep", "nbo"]:
    #     qdescp(
    #         w_dir_main=args.w_dir_main,
    #         destination=args.destination,
    #         files=args.files,
    #         json_files=args.json_files,
    #         task=args.qdescp,
    #         varfile=None,
    #     )

    # if args.QPRED == "nmr":
    #     nmr_main(args, log_overall, w_dir_main)
    # if args.QPRED == "energy":
    #     energy_main(args, log_overall, w_dir_main)
    # if args.QPRED == "dbstep":
    #     dbstep_par_main(args, log_overall, w_dir_main)
    # if args.QPRED == "nics":
    #     nics_par_main(args, log_overall, w_dir_main)
    # if args.QPRED == "cclib-json":
    #     cclib_main(args, log_overall, w_dir_main)
    # os.chdir(w_dir_main)
    #
    # # QSTAT
    # if args.QSTAT == "descp":
    #     geom_par_main(args, log_overall, w_dir_main)
    # if args.QSTAT == "graph":
    #     graph_main(args, log_overall, w_dir_main)
    # os.chdir(w_dir_main)
    #
    log_overall.finalize()

    out_data_file = Path("aqme_output.dat")
    if out_data_file.exists():
        out_data_file.replace(f"aqme_{args.output_name}.dat")


if __name__ == "__main__":
    main()
