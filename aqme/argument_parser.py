#####################################################.
#      This file contains the argument parser         #
#####################################################.

import argparse
import os

def parser_args():

    parser = argparse.ArgumentParser(
        description="Arguments for the different modules of AQME. The arguments can be changed using keywords in the command line or a yaml file."
    )

    # necessary input details
    parser.add_argument(
        "--varfile", 
        default=None, 
        help="Parameters in YAML format"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="File containing molecular structure(s)",
        default=" ",
    )
    parser.add_argument(
        "--output_name",
        action="store",
        default="output",
        help='Change output filename to AQME-"output".dat',
        type=str,
    )

    parser.add_argument(
        "--path",
        help="Path for analysis/boltzmann factor/combining files where the gaussian folder created is present",
        default="",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False, help="verbose output"
    )
    parser.add_argument(
        "--output",
        default=".sdf",
        help="The extension of the SDF files written",
    )

    # EXPAND QPREP TO GAUSSIAN, ORCA; TURBOMOLE!
    # work the script has to do
    parser.add_argument(
        "--CSEARCH",
        action="store",
        default=None,
        help="Perform conformational analysis with or without dihedrals",
        choices=["rdkit", "summ", "fullmonte", "crest"],
    )
    parser.add_argument(
        "--CMIN",
        action="store",
        default=None,
        help="Perform minimization after conformational analysis",
        choices=["xtb", "ani"],
    )
    parser.add_argument(
        "--QPREP",
        action="store",
        default=None,
        help="Create input files for QM calculations",
        choices=["gaussian", "orca"],
    )
    parser.add_argument(
        "--qcorr",
        action="store_true",
        default=False,
        help="Post-processing of output files from QM calculations",
    )
    parser.add_argument(
        "--QSTAT",
        action="store",
        default=None,
        help="Generate parameters for different conformers",
        choices=["graph", "descp"],
    )
    parser.add_argument(
        "--QPRED",
        action="store",
        default=None,
        help="Perform predictions for different conformers",
        choices=["nmr", "energy", "dbstep", "nics", "cclib-json"],
    )

    # arguments for TMBUILD
    parser.add_argument(
        "--metal_complex",
        action="store_true",
        default=False,
        help="Request metal complex with coord. no. 4, 5 or 6",
    )
    parser.add_argument(
        "--metal",
        help="Specify metallic element",
        default=[],
        type=str
    )
    parser.add_argument(
        "--mult",
        help="Multiplicity of metal complex or organic complexes",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--complex_coord",
        help="Coord. no. of metal complex (automatically updates)",
        default=[],
        type=int,
    )
    parser.add_argument(
        "--complex_type",
        help="Force geometry of the metal complex (options: linear, trigonalplanar, squareplanar, squarepyramidal)",
        default="",
        type=str,
    )
    parser.add_argument(
        "--m_oxi", 
        help="Metal oxidation state",
        default=[],
        type=int
    )
    parser.add_argument(
        "--metal_idx",
        help="Metal index (automatically updates)",
        default=[],
        type=int,
    )
    parser.add_argument(
        "--charge",
        help="Charge of metal complex (automatically updates)",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--charge_default",
        help="Charge default to be considered",
        default="auto",
    )
    parser.add_argument(
        "--metal_sym",
        help="Symbols of metals to be considered from list (automatically updates)",
        default=[],
        type=str,
    )

    # argumets for CSEARCH and CMIN
    parser.add_argument(
        "--ewin_cmin",
        action="store",
        default=5.0,
        help="energy window to print conformers for minimization using xTB or ANI (kcal/mol)",
        type=float,
    )
    parser.add_argument(
        "--ewin_csearch",
        action="store",
        default=5.0,
        help="energy window to print conformers for RDKit (kcal/mol)",
        type=float,
    )
    parser.add_argument(
        "--opt_fmax",
        action="store",
        default=0.05,
        help="fmax value used in xTB and AN1 optimizations",
        type=float,
    )
    parser.add_argument(
        "--opt_steps",
        action="store",
        default=1000,
        help="max cycles used in xTB and AN1 optimizations",
        type=int,
    )
    parser.add_argument(
        "--opt_steps_RDKit",
        action="store",
        default=1000,
        help="max cycles used in RDKit optimizations",
        type=int,
    )
    parser.add_argument(
        "--time",
        "-t",
        action="store_true",
        default=True,
        help="request program runtime",
    )
    parser.add_argument(
        "--heavyonly",
        help="only consider torsion angles involving heavy (non H) elements (default=True)",
        default=True,
    )
    parser.add_argument(
        "-d",
        "--degree",
        type=float,
        help="Amount, in degrees, to enumerate torsions by (default 120.0)",
        default=120.0,
    )
    parser.add_argument(
        "--max_torsions",
        type=int,
        help="Skip any molecules with more than this many torsions (default 20)",
        default=20,
    )
    parser.add_argument(
        "--sample",
        help="number of conformers to sample to get non-torsional differences (default 100)",
        default="auto",
    )
    parser.add_argument(
        "--auto_sample",
        help="final factor to multiply in the auto mode for the sample option (default 20)",
        default=20,
        type=int,
    )
    parser.add_argument(
        "--ff",
        help="force field (MMFF or UFF)",
        default="MMFF",
    )
    parser.add_argument(
        "--seed",
        help="random seed (default 062609)",
        default="062609",
        type=int,
    )
    parser.add_argument(
        "--rms_threshold",
        help="cutoff for considering sampled conformers the same (default 0.25)",
        default=0.25,
        type=float,
    )
    parser.add_argument(
        "--max_matches_RMSD",
        help="iteration cutoff for considering  matches in sampled conformers the same (default 1000)",
        default=1000,
        type=int,
    )
    parser.add_argument(
        "--energy_threshold",
        action="store",
        default=0.25,
        help="energy difference between unique conformers (default 0.25)",
    )
    parser.add_argument(
        "--initial_energy_threshold",
        action="store",
        default=0.0001,
        help="energy difference between unique conformers for the first filter of only E (default 0.0001)",
    )
    parser.add_argument(
        "--max_MolWt",
        help="Max. molecular weight of molecule",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "--ani_method",
        help="Specify ANI method used (i.e. ANI1x, ANI1ccx, ANI2x)",
        default="ANI2x",
        type=str,
    )
    parser.add_argument(
        "--STACKSIZE", help="Stack size available for xTB calculations", default="1G"
    )
    parser.add_argument(
        "--xtb_method",
        help="Specify xTB method used",
        default="GFN2-xTB",
        type=str,
    )
    parser.add_argument(
        "--xtb_solvent",
        help="Specify GBSA solvent used",
        default="none",
        type=str,
    )
    parser.add_argument(
        "--xtb_accuracy",
        help="Numerical accuracy of the xTB calculation",
        action="store",
        default=1.0,
    )
    parser.add_argument(
        "--xtb_electronic_temperature",
        help="Electronic temperature for TB methods",
        action="store",
        default=300.0,
    )
    parser.add_argument(
        "--xtb_max_iterations",
        help="Numerical accuracy of the xTB calculation",
        action="store",
        default=250,
    )
    parser.add_argument(
        "--cpus",
        action="store",
        default=12,
        help="Maximum number of threads to parallelize on while running CSEARCH and CMIN",
        type=int,
    )

    # arguments for FULLMONTE
    parser.add_argument(
        "--ewin_sample_fullmonte",
        action="store",
        default=2.0,
        help="energy window to consider conformers for sampling in FULLMONTE (default 2 kcal/mol)",
        type=float,
    )
    parser.add_argument(
        "--ewin_fullmonte",
        action="store",
        default=5.0,
        help="energy window to consider conformers for FULLMONTE (default 5 kcal/mol)",
        type=float,
    )
    parser.add_argument(
        "--nsteps_fullmonte",
        action="store",
        default=100,
        help="Number of steps to consider for FULLMONTE (default 100)",
        type=int,
    )
    parser.add_argument(
        "--nrot_fullmonte",
        action="store",
        default=3,
        help="Number of diherals to rotate for FULLMONTE (default 3) ",
        type=int,
    )
    parser.add_argument(
        "--ang_fullmonte",
        action="store",
        default=30,
        help="Angle to rotate each diheral of for FULLMONTE (default 30)",
        type=float,
    )

    # arguments for CREST
    parser.add_argument(
        "--cregen",
        action="store_true",
        help="Do cregen after crest",
        default=False,
    )
    parser.add_argument(
        "--cregen_ethr",
        action="store",
        default=0.2,
        help="Energy thershold for CREGEN after crest",
    )
    parser.add_argument(
        "--cregen_rthr",
        action="store",
        default=0.125,
        help="RMS thershold for CREGEN fter crest",
    )
    parser.add_argument(
        "--cregen_bthr",
        action="store",
        default=0.01,
        help="Rotational constant thershold for CREGEN after crest",
    )
    parser.add_argument(
        "--cregen_ewin",
        action="store",
        default=6,
        help="Energy window for CREGEN after crest",
    )

    # arguments for QPREP
    parser.add_argument(
        "--program",
        help="Program required to create the new input files",
        default="gaussian",
        type=str,
    )
    parser.add_argument(
        "--nprocs",
        help="Number of processors used in the QM calculations",
        default=2,
        type=int,
    )
    parser.add_argument(
        "--mem",
        help="Memory for the QM calculations (i) Gaussian: total memory; (ii) ORCA: memory per processor",
        default="4GB",
        type=str,
    )
    parser.add_argument(
        "--mol",
        help="Mol object to prepare the input QM file",
        default=None,
    )
    parser.add_argument(
        "--destination",
        help="Destination to create the new input files",
        default=None,
    )
    parser.add_argument(
        "--qm_input",
        help="Keywords line for new input files",
        default="",
    )
    parser.add_argument(
        "--ts_input",
        help="OPT options in Gaussian in QCORR for input files of TSs (disable this option with None)",
        default="opt=(calcfc,noeigen,ts,maxstep=5)",
    )
    parser.add_argument(
        "--qm_end",
        help="Final line in the new input files",
        default="",
    )
    parser.add_argument(
        "--gen_atoms",
        help="Atoms included in the gen(ECP) basis set",
        default=[],
    )
    parser.add_argument(
        "--bs",
        default='',
        help="Basis set used for non gen(ECP) atoms in gen(ECP) calculations",
        type=str,
    )
    parser.add_argument(
        "--bs_gen",
        default='',
        help="Basis set used for gen(ECP) atoms",
        type=str,
    )
    parser.add_argument(
        "--lowest_only",
        action="store_true",
        default=False,
        help="Lowest conformer to write in Gaussian",
    )
    parser.add_argument(
        "--lowest_n",
        action="store_true",
        default=False,
        help="Lowest Number of conformers to write in Gaussian",
    )
    parser.add_argument(
        "--energy_threshold_for_gaussian",
        help="Cut-off for considering sampled conformers in Gaussian inputs",
        default="100.0",
        type=float,
    )
    parser.add_argument(
        "--chk",
        action="store_true",
        default=False,
        help="Include the chk input line in new input files for Gaussian calculations",
    )

    # QPREP with Turbomole
    parser.add_argument(
        "--tmfunctional",
        help="Turbomole functionals",
        type=str,
        nargs="*",
        default=[
            "TPSS",
        ],
    )
    parser.add_argument(
        "--tmfunctionalsp",
        help="Turbomole functionals for SP",
        type=str,
        nargs="*",
        default=[
            "TPSS",
        ],
    )
    parser.add_argument(
        "--tmbasis", help="Turbomole basis set for all atoms for the optimization."
    )
    parser.add_argument(
        "--tmbasisfile", help="Turbomole basis set file for the optimization."
    )
    parser.add_argument(
        "--tmbasissp", help="Turbomole basis set for all atoms for SP", nargs="*"
    )
    parser.add_argument(
        "--tmbasisfilesp", help="Turbomole basis set files for SP", nargs="*"
    )
    parser.add_argument(
        "--tmgrid",
        help="Turbomole integration grid",
        default="m4",
        choices=["1", "2", "3", "4", "5", "6", "7", "m3", "m4", "m5"],
    )
    parser.add_argument(
        "--tmdispersion",
        help="Turbomole dispersion",
        default="off",
        choices=["off", "on", "bj", "d4"],
    )
    parser.add_argument(
        "--tmepsilon", help="Turbomole dielectric constant for COSMO", default="gas"
    )
    parser.add_argument(
        "--tmcavity", help="Turbomole solvent cavity for COSMO-RS", default="none"
    )
    parser.add_argument(
        "--tmricore", help="Turbomole max RIcore memory in MB", default=200, type=int
    )
    parser.add_argument(
        "--tmmaxcore", help="Turbomole max per core memory in MB", default=200, type=int
    )

    # other options for QPREP
    parser.add_argument(
        "--com_from_xyz",
        action="store_true",
        default=False,
        help="Create input files for Gaussian from an xyz file",
    )

    # arguments for QCORR (besides related functions from QPREP)
    # analysis of files
    parser.add_argument(
        "--w_dir_main",
        action="store",
        default=os.getcwd(),
        help="Working directory",
    )   

    parser.add_argument(
        "--qm_files",
        action="store",
        default=[],
        help="Filenames of QM output files to analyze",
    )   
    parser.add_argument(
        "--dup",
        action="store",
        default=True,
        help="Remove duplicates after DFT optimization",
    )
    parser.add_argument(
        "--dup_threshold",
        action="store",
        default=0.0001,
        help="Energy (in hartree) used as the energy difference in E, H and G to detect duplicates",
        type=float,
    )
    parser.add_argument(
        "--amplitude_ifreq",
        action="store",
        default=0.2,
        help="Amplitude used to displace the imaginary frequencies to fix",
        type=float,
    )
    parser.add_argument(
        "--ifreq_cutoff",
        action="store",
        default=0.0,
        help="Cut off for to consider whether a frequency is imaginary (absolute of the specified value is used)",
        type=float,
    )
    parser.add_argument(
        "--freq_conv",
        action="store",
        default='opt=(calcfc,maxstep=5)',
        help="If a string is defined, it will remove calculations that converged during optimization but did not convergence in the subsequent frequency calculation. Options: opt sections as strings i.e. (opt=(calcfc,maxstep=5)). If readfc is specified in the string, the chk option must be included as well. Turn this option off by using freq_conv=False.",
    )
    parser.add_argument(
        "--s2_threshold",
        action="store",
        default=10.0,
        help="Cut off for spin contamination during analysis in per cent of the expected value (i.e. multiplicity 3 has an the expected <S**2> of 2.0, if s2_threshold = 10 the <S**2> value is allowed to be 2.0 +- 0.2). Set s2_threshold = 0 to deactivate this option.",
        type=float,
    )
    parser.add_argument(
        "--isom",
        action="store",
        default=False,
        help="Check for isomerization from the initial input file to the resulting QM output files in QCORR. It requires the extension of the initial input files (i.e. isom='com') and the folder of the input files must be added in the isom_inputs option",
    )
    parser.add_argument(
        "--isom_inputs",
        action="store",
        default=os.getcwd(),
        help="Folder containing the initial input files to check for isomerization",
        type=str,
    )
    parser.add_argument(
        "--vdwfrac",
        action="store",
        help="Fraction of the summed VDW radii that constitutes a bond between two atoms in the isomerization filter",
        default=0.50,
        type=float,
    )
    parser.add_argument(
        "--covfrac",
        action="store",
        help="Fraction of the summed covalent radii that constitutes a bond between two atoms in the isomerization filter",
        default=1.10,
        type=float,
    )
    parser.add_argument(
        "--fullcheck",
        action="store",
        default=True,
        help="Perform an analysis to detect whether the calculations were done homogeneously (i.e. same level of theory, solvent, grid size, etc)",
    )
    parser.add_argument(
        "--author",
        action="store",
        default="",
        help="Author of the calculations",
        type=str,
    )

    # writing input files from json format
    parser.add_argument(
        "--json2input",
        action="store_true",
        help="Create QM single point input files from json files",
        default=False,
    )
    parser.add_argument(
        "--json_files",
        help="Filenames of json files to analyze",
        default=[],
    )
    parser.add_argument(
        "--suffix",
        help="Suffix for the new input files",
        default='',
    )

    # argumets for QSTAT
    parser.add_argument(
        "--rot_dihedral",
        action="store_true",
        default=False,
        help="Turn on for tracking the geometric parameters for the rotatable dihedrals (Need not specify anything in the dihedral list)",
    )
    parser.add_argument(
        "--dihedral",
        help="Specify the atom indexes to track dihedrals for different conformes only for specific dihedrals (For all rotatable dihedrals turn rot_dihedral to True)",
        default=[],
        type=str,
        nargs=4,
        action="append",
    )
    parser.add_argument(
        "--bond",
        help="Specify the atom indexes to track bond lengths for different conformers",
        default=[],
        type=str,
        nargs=2,
        action="append",
    )
    parser.add_argument(
        "--angle",
        help="Specify the atom indexes to track angles for different conformers",
        default=[],
        type=str,
        nargs=3,
        action="append",
    )
    parser.add_argument(
        "--geom_par_name",
        action="store",
        default="descp",
        help="Change the prefix for the descriptors obtained",
    )
    parser.add_argument(
        "--dbstep_cen_lig_file",
        help="Center for DBSTEP steric paramters in a txt ( FORMAT : name, center, ligand)",
        action="store",
        default="No file passed",
    )

    # arguments for nmr
    parser.add_argument(
        "--nmr_exp",
        default="fromsdf",
        help="From where the experimental NMR details will be obtained",
    )
    parser.add_argument(
        "--nmr_online",
        action="store_true",
        default=False,
        help="Turn to true for checking NMR scaling factors from ChesHire Database",
    )
    parser.add_argument(
        "--nmr_aos",
        help="Specify the type of atomic basis used for nmr calculation (default = giao) ",
        default="giao",
        type=str,
    )
    parser.add_argument(
        "--nmr_nucleus",
        help="Specify the nucleus for nmr analysis default (['C','H'])",
        default=["C", "H"],
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--nmr_slope",
        help="Specify the slope for each nucleus for nmr analysis default([1.0673,1.0759])",
        default=[1.0673, 1.0759],
        type=float,
        nargs="*",
    )
    parser.add_argument(
        "--nmr_intercept",
        help="Specify the intercept for each nucleus for nmr analysis default([-15.191,-2.2094])",
        default=[-15.191, -2.2094],
        type=float,
        nargs="*",
    )
    parser.add_argument(
        "--nmr_tms_ref",
        help="Specify the reference for TMS for each nucleus for nmr analysis default([191.79,31.39])",
        default=[191.79, 31.39],
        type=float,
        nargs="*",
    )

    # arguments for NICS
    parser.add_argument(
        "--nics_range",
        help="Range to calculate NICS along a given axis",
        default=4,
    )
    parser.add_argument(
        "--nics_number",
        help="Step size to calculate NICS along a given axis",
        default=16,
    )
    parser.add_argument(
        "--nics_atoms_file",
        help="NICS atoms in a txt ( FORMAT : name, atom1, atom2, atom3..)",
        action="store",
        default="No file passed",
    )

    # arguments for cclib

    # submission of Gaussian files
    parser.add_argument(
        "--qsub",
        action="store_true",
        default=False,
        help="Submit Gaussian files when they are created",
    )
    parser.add_argument(
        "--qsub_ana",
        action="store_true",
        default=False,
        help="Submit Gaussian files after analysis",
    )
    parser.add_argument(
        "--submission_command",
        help="Queueing system that the submission is done on",
        default="qsub_summit",
        type=str,
    )

    # apply exp rules
    parser.add_argument(
        "--geom_rules",
        default=[],
        help="Discarding rules applied to filter-off conformers (based on experimental observation for example). Format: i) Automatic rules: ['Ir_bidentate_x3'], ii) manual rules: ['ATOM1-ATOM2-ATOM3, ANGLE'] (i.e. ['C-Pd-C, 180'])",
    )
    parser.add_argument(
        "--angle_off",
        type=float,
        help="Deviation to discard in geom_rules (i.e. 180 +- 30 degrees)",
        default=30,
    )

    ##### further additions #####
    # NCI complex
    parser.add_argument(
        "--nci_complex",
        action="store_true",
        default=False,
        help="Request NCI complexes conformational search",
    )
    # TS complex
    parser.add_argument(
        "--ts_complex",
        action="store_true",
        default=False,
        help="Request TS complexes conformational search",
    )
    parser.add_argument(
        "--cbonds",
        action="store",
        default=0.5,
        help="cbonds fro NCI complexes",
    )
    parser.add_argument(
        "--prefix",
        help="Prefix for naming files",
        default="None",
        type=str,
    )

    args = parser.parse_args()

    return args


### part for using the options in a script or jupyter notebook
class options_add:
    pass


def set_options(kwargs):
    # set default options and options provided
    options = options_add()
    # dictionary containing default values for options
    var_dict = {
        "varfile": None,
        "input": " ",
        "output_name": "output",
        "path": "",
        "verbose": False,
        "output": ".sdf",
        "CSEARCH": None,
        "CMIN": None,
        "QPREP": None,
        "QCORR": None,
        "QSTAT": None,
        "QPRED": None,
        "metal_complex": False,
        "metal": [],
        "mult": 1,
        "complex_coord": [],
        "complex_type": "",
        "m_oxi": [],
        "metal_idx": [],
        "charge": [],
        "charge_default": "auto",
        "metal_sym": [],
        "ewin_cmin": 5.0,
        "ewin_csearch": 5.0,
        "opt_fmax": 0.05,
        "opt_steps": 1000,
        "opt_steps_RDKit": 1000,
        "time": True,
        "heavyonly": True,
        "degree": 120.0,
        "max_torsions": 20,
        "sample": "auto",
        "auto_sample": 20,
        "ff": "MMFF",
        "seed": 62609,
        "rms_threshold": 0.25,
        "max_matches_RMSD": 1000,
        "energy_threshold": 0.25,
        "initial_energy_threshold": 0.0001,
        "max_MolWt": 10000,
        "ani_method": "ANI2x",
        "STACKSIZE": "1G",
        "xtb_method": "GFN2-xTB",
        "xtb_solvent": "none",
        "xtb_accuracy": 1.0,
        "xtb_electronic_temperature": 300.0,
        "xtb_max_iterations": 250,
        "cpus": 12,
        "ewin_sample_fullmonte": 2.0,
        "ewin_fullmonte": 5.0,
        "nsteps_fullmonte": 100,
        "nrot_fullmonte": 3,
        "ang_fullmonte": 30,
        "cregen": False,
        "cregen_ethr": 0.2,
        "cregen_rthr": 0.125,
        "cregen_bthr": 0.01,
        "cregen_ewin": 6,
        "nprocs": 24,
        "mem": "96GB",
        "genecp_bs": None,
        "bs": None,
        "aux_atoms_orca": [],
        "aux_genecp_bs": [],
        "aux_fit_gen_atoms": [],
        "cpcm_input": "None",
        "orca_scf_iters": 500,
        "mdci_orca": "None",
        "print_mini_orca": True,
        "qm_input": "",
        "ts_input": "opt=(calcfc,noeigen,ts,maxstep=5)",
        "qm_end": "",
        "gen_atoms": None,
        "lowest_only": False,
        "lowest_n": False,
        "energy_threshold_for_gaussian": 100.0,
        "chk": False,
        "tmfunctional": ["TPSS"],
        "tmfunctionalsp": ["TPSS"],
        "tmbasis": None,
        "tmbasisfile": None,
        "tmbasissp": None,
        "tmbasisfilesp": None,
        "tmgrid": "m4",
        "tmdispersion": "off",
        "tmepsilon": "gas",
        "tmcavity": "none",
        "tmricore": 200,
        "tmmaxcore": 200,
        "com_from_xyz": False,
        "dup": True,
        "dup_threshold": 0.0001,
        "check_geom": False,
        "length_criteria": 1.4,
        "amplitude_ifreq": 0.2,
        "ifreq_cutoff": 0.0,
        "s2_threshold": 10.0,
        "sp": None,
        "nics": False,
        "charge_sp": None,
        "mult_sp": None,
        "qm_input_sp": "",
        "qm_end_sp": "",
        "gen_bs_sp": None,
        "bs_sp": None,
        "suffix": "None",
        "aux_atoms_orca_sp": [],
        "aux_gen_bs_sp": [],
        "aux_fit_gen_atoms_sp": [],
        "cpcm_input_sp": "None",
        "orca_scf_iters_sp": 500,
        "mdci_orca_sp": "None",
        "print_mini_orca_sp": True,
        "rot_dihedral": False,
        "dihedral": [],
        "bond": [],
        "angle": [],
        "geom_par_name": "descp",
        "dbstep_cen_lig_file": "No file passed",
        "nmr_exp": "fromsdf",
        "nmr_online": False,
        "nmr_aos": "giao",
        "nmr_nucleus": ["C", "H"],
        "nmr_slope": [1.0673, 1.0759],
        "nmr_intercept": [-15.191, -2.2094],
        "nmr_tms_ref": [191.79, 31.39],
        "nics_range": 4,
        "nics_number": 16,
        "nics_atoms_file": "No file passed",
        "qsub": False,
        "qsub_ana": False,
        "submission_command": "qsub_summit",
        "geom_rules": [],
        "angle_off": 30,
        "nci_complex": False,
        "ts_complex": False,
        "cbonds": 0.5,
        "prefix": "None",
        "constraints_dist": None,
        "constraints_angle": None,
        "constraints_dihedral": None,
        "isom": None,
        "nocheck": False,
        "nocom": False,
        "program": "gaussian",
        "program_sp": "gaussian",
        "qcorr_json": "",
        "bs_gen" : None
    }

    for key in var_dict:
        vars(options)[key] = var_dict[key]
    for key in kwargs:
        if key in var_dict:
            vars(options)[key] = kwargs[key]
        else:
            print(
                "Warning! Option: [",
                key,
                ":",
                kwargs[key],
                "] provided but no option exists, try -h to see available options.",
            )

    return options
