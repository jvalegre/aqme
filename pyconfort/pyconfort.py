#!/usr/bin/env python

#########################################################################################.
#########################################################################################
###                                                                                   ###
###  pyCONFORT is a tool that allows to carry out automated:                          ###
###  (1) Conformational searches and creation of COM files using RDKit, xTB and ANI1  ###
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
from pyconfort.argument_parser import parser_args

from pyconfort.mainf import (csearch_main, exp_rules_main, qprep_main, 
                             move_sdf_main, qcorr_gaussian_main,dup_main,
                             graph_main,geom_par_main,nmr_main,energy_main,
                             load_from_yaml,creation_of_ana_csv,dbstep_par_main,
                             nics_par_main,cclib_main,cmin_main)
from pyconfort.utils import Logger

def main():
    # working directory and arguments
    w_dir_initial = os.getcwd()
    base_path = Path(w_dir_initial)

    args = parser_args()

    log_overall = Logger("pyCONFORT", args.output_name)
    #if needed to load from a yaml file
    load_from_yaml(args,log_overall)

    name = args.input.split('.')[0]

    # IS THIS NECESSARY?!
    #setting variable if needed
    for _ in args.basis_set:
        if len(args.basis_set) != len(args.basis_set_genecp_atoms):
            args.basis_set_genecp_atoms.append('')
    for _ in args.basis_set_sp:
        if len(args.basis_set_sp) != len(args.basis_set_genecp_atoms_sp):
            args.basis_set_genecp_atoms_sp.append('')

    #CSEARCH AND CMIN
    if args.CSEARCH in ['rdkit', 'summ', 'fullmonte']: # RAUL: Is there any other posibilities or just None?
        start_time_overall = time.time()
        csearch_dup_data = csearch_main(w_dir_initial,args,log_overall)
        if args.time:
            elapsed_time = round(time.time() - start_time_overall,2) 
            log_overall.write(f"\n All molecules execution time CSEARCH: {elapsed_time} seconds")
        os.chdir(w_dir_initial)
        if args.CMIN is None:
            csearch_csv_folder = base_path.joinpath('CSEARCH/csv_files')
            csearch_csv_folder.mkdir(exist_ok=True)
            csearch_csv_file = csearch_csv_folder.joinpath(f'{name}-CSEARCH-Data.csv')
            csearch_dup_data.to_csv(csearch_csv_file,index=False)

    #Separating CMIN
    if args.CSEARCH != None and args.CMIN in ['xtb','ani']:
        cmin_dup_data = cmin_main(w_dir_initial,args,log_overall,csearch_dup_data)
        if args.time:
            elapsed_time = round(time.time() - start_time_overall,2) 
            log_overall.write(f"\n All molecules execution time CMIN: {elapsed_time} seconds")
        os.chdir(w_dir_initial)
        cmin_csv_folder = base_path.joinpath('/CMIN/csv_files')
        cmin_csv_folder.mkdir(exist_ok=True)
        cmin_csv_file = cmin_csv_folder.joinpath(f'{name}-CMIN-Data.csv')
        cmin_dup_data.to_csv(cmin_csv_file,index=False)


    #applying rules to discard certain conformers based on rules that the user define
    if len(args.exp_rules) >= 1:
        exp_rules_active = True
        if args.QCORR=='gaussian':
            exp_rules_active = False
        exp_rules_main(args,log_overall,exp_rules_active)
        os.chdir(w_dir_initial)

    #QPREP
    if args.QPREP=='gaussian' or args.QPREP=='orca' or args.QPREP=='turbomole':
        qprep_main(w_dir_initial,args,log_overall)
        os.chdir(w_dir_initial)

    if args.CSEARCH in ['rdkit', 'summ', 'fullmonte'] or args.QPREP is not None:
        # moving files after compute and/or write_gauss
        move_sdf_main(args)
        os.chdir(w_dir_initial)

    #QCORR
    if args.QCORR=='gaussian':
        log_overall.write("\no  Writing analysis of output files in respective folders\n")
        # main part of the duplicate function
        if args.dup:
            try:
                import goodvibes
            except (ModuleNotFoundError,AttributeError):
                log_overall.write("\nx  GoodVibes is not installed as a module (pip or conda), the duplicate option will be disabled in QCORR\n")
            else:
                # If the AttributeError can happen here do not accept the change
                duplicates = dup_main(args, log_overall, w_dir_initial) 
                os.chdir(w_dir_initial)
        else:
            duplicates = False

        # main part of the output file analyzer for errors/imag freqs
        qcorr_gaussian_main(duplicates,w_dir_initial,args,log_overall)
        os.chdir(w_dir_initial)

    #QPRED
    if args.QPRED=='nmr':
        nmr_main(args,log_overall,w_dir_initial)
    if args.QPRED=='energy':
        energy_main(args,log_overall,w_dir_initial)
    if args.QPRED=='dbstep':
        dbstep_par_main(args,log_overall,w_dir_initial)
    if args.QPRED=='nics':
        nics_par_main(args,log_overall,w_dir_initial)
    if args.QPRED=='cclib-json':
        cclib_main(args,log_overall,w_dir_initial)
    os.chdir(w_dir_initial)

    #QSTAT
    if args.QSTAT=='descp':
        geom_par_main(args,log_overall,w_dir_initial)
    if args.QSTAT=='graph':
        graph_main(args,log_overall,w_dir_initial)
    os.chdir(w_dir_initial)

    log_overall.finalize()

    out_data_file = Path(f'pyCONFORT_output.dat')
    if out_data_file.exists():
        out_data_file.replace(f'pyCONFORT_{args.output_name}.dat')


if __name__ == "__main__":
    main()
