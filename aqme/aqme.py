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
from aqme.mainf import (
    csearch_main,
    qprep,
    cmin_main)
from aqme.utils import (Logger,
    command_line_args)
from aqme.qcorr import qcorr


def main():
    '''
    Main function of AQME, acts as the starting point when the program is run through a terminal
    '''

    # load user-defined arguments from command line
    args = command_line_args()
    
    # working directory and arguments
    w_dir_main = Path(args.w_dir_main)

    log_overall = Logger("aqme", args.output_name)
    start_time_overall = time.time()

    name = args.input.split(".")[0]

    # CSEARCH AND CMIN
    if args.csearch in [
        "rdkit",
        "summ",
        "fullmonte",
        "crest",
    ]:
        csearch_dup_data = csearch_main(w_dir_main, args, log_overall)
        os.chdir(w_dir_main)
        elapsed_time = round(time.time() - start_time_overall, 2)
        log_overall.write(f"\n Time CSEARCH: {elapsed_time} seconds")
        if args.cmin is None:
            csearch_csv_folder = w_dir_main.joinpath("CSEARCH/csv_files")
            csearch_csv_folder.mkdir(exist_ok=True)
            csearch_csv_file = csearch_csv_folder.joinpath(f"{name}-CSEARCH-Data.csv")
            csearch_dup_data.to_csv(csearch_csv_file, index=False)

    # Separating CMIN
    if args.csearch != None and args.cmin in ["xtb", "ani"]:
        cmin_dup_data = cmin_main(w_dir_main, args, log_overall, csearch_dup_data)
        os.chdir(w_dir_main)
        cmin_csv_folder = w_dir_main.joinpath("CMIN/csv_files")
        cmin_csv_folder.mkdir(exist_ok=True)
        cmin_csv_file = cmin_csv_folder.joinpath(f"{name}-CMIN-Data.csv")
        cmin_dup_data.to_csv(cmin_csv_file, index=False)
        elapsed_time = round(time.time() - start_time_overall, 2)
        log_overall.write(f"\n Time CMIN: {elapsed_time} seconds")

    # QPREP
    if args.qprep:
        qprep(files=args.files,
            atom_types=args.atom_types,
            cartesians=args.cartesians,
            w_dir_main=args.w_dir_main,
            destination=args.destination,
            varfile=args.varfile,
            program=args.program,
            qm_input=args.qm_input,
            qm_end=args.qm_end,
            charge=args.charge,
            mult=args.mult,
            suffix=args.suffix,
            chk=args.chk,
            mem=args.mem,
            nprocs=args.nprocs,
            gen_atoms=args.gen_atoms,
            bs_gen=args.bs_gen,
            bs=args.bs)
        elapsed_time = round(time.time() - start_time_overall, 2)
        log_overall.write(f"\n Time QPREP: {elapsed_time} seconds")

    # QCORR
    if args.qcorr:
        qcorr(files=args.files,
            w_dir_main=args.w_dir_main,
            fullcheck=args.fullcheck,
            varfile=args.varfile,
            ifreq_cutoff=args.ifreq_cutoff,
            amplitude_ifreq=args.amplitude_ifreq,
            freq_conv=args.freq_conv,
            s2_threshold=args.s2_threshold,
            dup_threshold=args.dup_threshold,
            isom=args.isom,
            isom_inputs=args.isom_inputs,
            vdwfrac=args.vdwfrac,
            covfrac=args.covfrac,
            program=args.program,
            mem=args.mem,
            nprocs=args.nprocs,
            qm_input=args.qm_input,
            qm_end=args.qm_end,
            chk=args.chk,
            gen_atoms=args.gen_atoms,
            bs_gen=args.bs_gen,
            bs=args.bs)
        elapsed_time = round(time.time() - start_time_overall, 2)
        log_overall.write(f"\n Time QCORR: {elapsed_time} seconds")

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

    # if args.qpred == "nmr":
    #     nmr_main(args, log_overall, w_dir_main)
    # if args.qpred == "energy":
    #     energy_main(args, log_overall, w_dir_main)
    # if args.qpred == "dbstep":
    #     dbstep_par_main(args, log_overall, w_dir_main)
    # if args.qpred == "nics":
    #     nics_par_main(args, log_overall, w_dir_main)
    # if args.qpred == "cclib-json":
    #     cclib_main(args, log_overall, w_dir_main)
    # os.chdir(w_dir_main)
    #
    # # qstat
    # if args.qstat == "descp":
    #     geom_par_main(args, log_overall, w_dir_main)
    # if args.qstat == "graph":
    #     graph_main(args, log_overall, w_dir_main)
    # os.chdir(w_dir_main)
    #
    log_overall.finalize()

    out_data_file = Path("aqme_output.dat")
    if out_data_file.exists():
        out_data_file.replace(f"aqme_{args.output_name}.dat")


if __name__ == "__main__":
    main()
