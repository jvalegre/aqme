#!/usr/bin/env python

import os
import sys
import time
from pathlib import Path

###########################################################################################.
###########################################################################################
###                                                                                     ###
###  AQME is a tool that allows to carry out automated:                                 ###
###  (CSEARCH) Conformational searches and creation of COM files using RDKit and CREST  ###
###  (CMIN) Geometry refinement of initial conformers with xTB and ANI                  ###
###  (QCORR) Out put file processing from QM calculations and automated issue fixing,   ###
###  including imaginary freqs, spin contamination, isomerization issues and            ###
###  error terminations, among others                                                   ###
###  (QPREP) Use QM outputs, XYZ, SDF, PDB, JSON and other 3D formats to create input   ###
###  files for multiple QM programs                                                     ###
###  (QDESCP) Generate xTB molecular descriptors, including Boltzmann averaged values,  ###
###  to use in machine learning models                                                  ###
###                                                                                     ###
###########################################################################################
###                                                                                     ###
###  Authors: Juan V. Alegre Requena, Shree Sowndarya S. V., Brenda Manzanilla          ###
###                                                                                     ###
###  Please, report any bugs or suggestions to:                                         ###
###  jv.alegre@csic.es                                                                  ###
###                                                                                     ###
###########################################################################################
###########################################################################################.

from aqme.csearch import csearch
from aqme.cmin import cmin
from aqme.qprep import qprep
from aqme.utils import command_line_args
from aqme.qcorr import qcorr
from aqme.qdescp import qdescp
from aqme.utils import Logger, aqme_ref, aqme_version, time_run, _format_command_line_for_log


def main():
    """
    Main function of AQME, acts as the starting point when the program is run through a terminal
    """

    # load user-defined arguments from command line
    args = command_line_args()
    args.command_line = True

    if not args.csearch and not args.cmin and not args.qprep and not args.qcorr and not args.qdescp and not args.milo:
        print('x  No module was specified in the command line! (i.e. --csearch for conformer generation). If you did specify a module, check that you are using quotation marks when using options (i.e. --files "*.sdf").\n')

    # CSEARCH
    if args.csearch:
        csearch(**vars(args))

    # CMIN
    if args.cmin:
        cmin(**vars(args))

    # QPREP
    if args.qprep:
        qprep(**vars(args))

    # QCORR
    if args.qcorr:
        qcorr(**vars(args))

    # QDESCP
    if args.qdescp:
        qdescp(**vars(args))

    # MILO
    if args.milo:
        from aqme.Anat_Milo.core import run_anat_milo_workflow
        milo_log = Logger(Path(os.getcwd()) / "MILO", "data", verbose=True)
        milo_log.write(f"AQME v {aqme_version} {time_run} \nCitation: {aqme_ref}\n")
        milo_log.write(f"Command line used in AQME: python -m aqme {_format_command_line_for_log(sys.argv[1:])}\n")

        start_time = time.time()
        try:
            if args.outputs and args.yaml_file:
                run_anat_milo_workflow(
                    yaml_file=args.yaml_file,
                    outputs_dir=args.outputs,
                    log=milo_log,
                )
            elif args.outputs:
                run_anat_milo_workflow(
                    outputs_dir=args.outputs,
                    log=milo_log,
                )
            elif args.yaml_file and args.input:
                run_anat_milo_workflow(
                    yaml_file=args.yaml_file,
                    data_dir=args.input,
                    log=milo_log,
                )
            else:
                milo_log.write("x  MILO needs either --outputs, or --yaml_file together with --input.")
                raise SystemExit(1)
        except Exception as exc:
            milo_log.write(f"x  MILO failed: {exc}")
            raise SystemExit(1) from exc
        finally:
            elapsed = round(time.time() - start_time, 2)
            milo_log.write(f"\nTime MILO: {elapsed} seconds\n")
            milo_log.finalize()

if __name__ == "__main__":
    main()
import os
import time
from pathlib import Path
