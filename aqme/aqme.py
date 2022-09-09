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
###  (QPREP) Use QM outputs, XYZ, SDF, PDB, JSON and other 3D formats to create input   ###
###  files for multiple QM programs                                                     ###
###  (QDESCP) Generate xTB molecular descriptors, including Boltzmann averaged values,  ###
###  to use in machine learning models                                                  ###
###                                                                                     ###
###########################################################################################
###                                                                                     ###
###  Authors: Shree Sowndarya S. V., Juan V. Alegre Requena                             ###
###                                                                                     ###
###  Please, report any bugs or suggestions to:                                         ###
###  svss@colostate.edu or jvalegre@unizar.es                                           ###
###                                                                                     ###
###########################################################################################
###########################################################################################.

from aqme.csearch import csearch
from aqme.cmin import cmin
from aqme.qprep import qprep
from aqme.utils import command_line_args
from aqme.qcorr import qcorr
from aqme.qdescp import qdescp


def main():
    """
    Main function of AQME, acts as the starting point when the program is run through a terminal
    """

    # load user-defined arguments from command line
    args = command_line_args()
    args.command_line = True

    # CSEARCH
    if args.csearch:
        csearch(
            input=args.input,
            command_line=args.command_line,
            smi=args.smi,
            name=args.name,
            w_dir_main=args.w_dir_main,
            destination=args.destination,
            charge=args.charge,
            mult=args.mult,
            sample=args.sample,
            max_workers=args.max_workers,
            metal_complex=args.metal_complex,
            metal_atoms=args.metal_atoms,
            metal_idx=args.metal_idx,
            complex_coord=args.complex_coord,
            metal_sym=args.metal_sym,
            complex_type=args.complex_type,
            opt_steps_rdkit=args.opt_steps_rdkit,
            heavyonly=args.heavyonly,
            max_matches_rmsd=args.max_matches_rmsd,
            max_mol_wt=args.max_mol_wt,
            ewin_csearch=args.ewin_csearch,
            initial_energy_threshold=args.initial_energy_threshold,
            energy_threshold=args.energy_threshold,
            rms_threshold=args.rms_threshold,
            auto_sample=args.auto_sample,
            ff=args.ff,
            degree=args.degree,
            verbose=args.verbose,
            output=args.output,
            seed=args.seed,
            max_torsions=args.max_torsions,
            varfile=args.varfile,
            program=args.program,
            constraints_atoms=args.constraints_atoms,
            constraints_dist=args.constraints_dist,
            constraints_angle=args.constraints_angle,
            constraints_dihedral=args.constraints_dihedral,
            prefix=args.prefix,
            stacksize=args.stacksize,
            ewin_fullmonte=args.ewin_fullmonte,
            ewin_sample_fullmonte=args.ewin_sample_fullmonte,
            nsteps_fullmonte=args.nsteps_fullmonte,
            nrot_fullmonte=args.nrot_fullmonte,
            ang_fullmonte=args.ang_fullmonte,
            crest_keywords=args.crest_keywords,
            angle_off=args.angle_off,
            nprocs=args.nprocs,
            cregen=args.cregen,
            cregen_keywords=args.cregen_keywords,
            crest_force=args.crest_force,
        )

    # CMIN
    if args.cmin:
        cmin(
            files=args.files,
            command_line=args.command_line,
            w_dir_main=args.w_dir_main,
            destination=args.destination,
            verbose=args.verbose,
            program=args.program,
            xtb_method=args.xtb_method,
            xtb_solvent=args.xtb_solvent,
            output=args.output,
            charge=args.charge,
            mult=args.mult,
            metal_complex=args.metal_complex,
            metal_atoms=args.metal_atoms,
            metal_oxi=args.metal_oxi,
            complex_type=args.complex_type,
            ewin_cmin=args.ewin_cmin,
            initial_energy_threshold=args.initial_energy_threshold,
            energy_threshold=args.energy_threshold,
            rms_threshold=args.rms_threshold,
            xtb_accuracy=args.xtb_accuracy,
            xtb_electronic_temperature=args.xtb_electronic_temperature,
            xtb_max_iterations=args.xtb_max_iterations,
            opt_steps=args.opt_steps,
            opt_fmax=args.opt_fmax,
            ani_method=args.ani_method,
            stacksize=args.stacksize,
        )

    # QPREP
    if args.qprep:
        qprep(
            files=args.files,
            command_line=args.command_line,
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
            bs_nogen=args.bs_nogen,
            verbose=args.verbose,
        )

    # QCORR
    if args.qcorr:
        qcorr(
            files=args.files,
            command_line=args.command_line,
            w_dir_main=args.w_dir_main,
            fullcheck=args.fullcheck,
            varfile=args.varfile,
            ifreq_cutoff=args.ifreq_cutoff,
            amplitude_ifreq=args.amplitude_ifreq,
            freq_conv=args.freq_conv,
            s2_threshold=args.s2_threshold,
            dup_threshold=args.dup_threshold,
            isom_type=args.isom_type,
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
            bs_nogen=args.bs_nogen,
        )

    # qdescp
    # required to do CSEARCH or have SDF files of coformers
    if args.qdescp:
        qdescp(
            w_dir_main=args.w_dir_main,
            destination=args.destination,
            files=args.files,
            charge=args.charge,
            mult=args.mult,
            program=args.program,
            qdescp_temp=args.qdescp_temp,
            qdescp_acc=args.qdescp_acc,
            qdescp_solvent=args.qdescp_solvent,
            boltz=args.boltz,
            nmr_atoms=args.nmr_atoms,
            nmr_slope=args.nmr_slope,
            nmr_intercept=args.nmr_intercept,
            nmr_experim=args.nmr_experim,
        )


if __name__ == "__main__":
    main()
