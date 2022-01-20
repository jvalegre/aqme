#####################################################.
#      This file contains the argument parser         #
#####################################################.

import argparse


def parser_args():

    parser = argparse.ArgumentParser(
        description="Generate conformers depending on type of optimization (change parameters in the params yaml file)."
    )

    # necessary input details
    parser.add_argument(
        "--varfile", dest="varfile", default=None, help="Parameters in YAML format"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="File containing molecular structure(s)",
        dest="input",
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
        dest="path",
        default="",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False, help="verbose output"
    )
    parser.add_argument(
        "--output",
        dest="output",
        default=".sdf",
        metavar="output",
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
        action="store",
        default=None,
        help="Fix the output files from QM calculations",
        choices=["gaussian"],
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
        "--metal", help="Specify metallic element", default=[], dest="metal", type=str
    )
    parser.add_argument(
        "--mult",
        help="Multiplicity of metal complex or organic complexes",
        default="1",
        dest="mult",
        type=int,
    )
    parser.add_argument(
        "--complex_coord",
        help="Coord. no. of metal complex (automatically updates)",
        default=[],
        dest="complex_coord",
        type=int,
    )
    parser.add_argument(
        "--complex_type",
        help="Force geometry of the metal complex (options: linear, trigonalplanar, squareplanar, squarepyramidal)",
        default="",
        dest="complex_type",
        type=str,
    )
    parser.add_argument(
        "--m_oxi", help="Metal oxidation state", default=[], dest="m_oxi", type=int
    )
    parser.add_argument(
        "--metal_idx",
        help="Metal index (automatically updates)",
        default=[],
        dest="metal_idx",
        type=int,
    )
    parser.add_argument(
        "--charge",
        help="Charge of metal complex (automatically updates)",
        default=[],
        dest="charge",
        type=int,
    )
    parser.add_argument(
        "--charge_default",
        help="Charge default to be considered",
        default="auto",
        dest="charge_default",
    )
    parser.add_argument(
        "--metal_sym",
        help="Symbols of metals to be considered from list (automatically updates)",
        default=[],
        dest="metal_sym",
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
        metavar="heavyonly",
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
        metavar="sample",
    )
    parser.add_argument(
        "--auto_sample",
        help="final factor to multiply in the auto mode for the sample option (default 20)",
        default=20,
        type=int,
        metavar="auto_sample",
    )
    parser.add_argument(
        "--ff", help="force field (MMFF or UFF)", default="MMFF", metavar="ff"
    )
    parser.add_argument(
        "--seed",
        help="random seed (default 062609)",
        default="062609",
        type=int,
        metavar="s",
    )
    parser.add_argument(
        "--rms_threshold",
        help="cutoff for considering sampled conformers the same (default 0.25)",
        default=0.25,
        type=float,
        metavar="R",
    )
    parser.add_argument(
        "--max_matches_RMSD",
        help="iteration cutoff for considering  matches in sampled conformers the same (default 1000)",
        default=1000,
        type=int,
        metavar="max_matches_RMSD",
    )
    parser.add_argument(
        "--energy_threshold",
        dest="energy_threshold",
        action="store",
        default=0.25,
        help="energy difference between unique conformers (default 0.25)",
    )
    parser.add_argument(
        "--initial_energy_threshold",
        dest="initial_energy_threshold",
        action="store",
        default=0.0001,
        help="energy difference between unique conformers for the first filter of only E (default 0.0001)",
    )
    parser.add_argument(
        "--max_MolWt",
        help="Max. molecular weight of molecule",
        default=10000,
        type=int,
        metavar="max_MolWt",
    )
    parser.add_argument(
        "--ani_method",
        help="Specify ANI method used (i.e. ANI1x, ANI1ccx, ANI2x)",
        default="ANI2x",
        dest="ani_method",
        type=str,
    )
    parser.add_argument(
        "--STACKSIZE", help="Stack size available for xTB calculations", default="1G"
    )
    parser.add_argument(
        "--xtb_method",
        help="Specify xTB method used",
        default="GFN2-xTB",
        dest="xtb_method",
        type=str,
    )
    parser.add_argument(
        "--xtb_solvent",
        help="Specify GBSA solvent used",
        default="none",
        dest="xtb_solvent",
        type=str,
    )
    parser.add_argument(
        "--xtb_accuracy",
        help="Numerical accuracy of the xTB calculation",
        action="store",
        default=1.0,
        dest="xtb_accuracy",
    )
    parser.add_argument(
        "--xtb_electronic_temperature",
        help="Electronic temperature for TB methods",
        action="store",
        default=300.0,
        dest="xtb_electronic_temperature",
    )
    parser.add_argument(
        "--xtb_max_iterations",
        help="Numerical accuracy of the xTB calculation",
        action="store",
        default=250,
        dest="xtb_max_iterations",
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
        dest="cregen",
        action="store_true",
        help="Do cregen after crest",
        default=False,
    )
    parser.add_argument(
        "--cregen_ethr",
        dest="cregen_ethr",
        action="store",
        default=0.2,
        help="Energy thershold for CREGEN after crest",
    )
    parser.add_argument(
        "--cregen_rthr",
        dest="cregen_rthr",
        action="store",
        default=0.125,
        help="RMS thershold for CREGEN fter crest",
    )
    parser.add_argument(
        "--cregen_bthr",
        dest="cregen_bthr",
        action="store",
        default=0.01,
        help="Rotational constant thershold for CREGEN after crest",
    )
    parser.add_argument(
        "--cregen_ewin",
        dest="cregen_ewin",
        action="store",
        default=6,
        help="Energy window for CREGEN after crest",
    )

    # arguments for QPREP
    parser.add_argument(
        "--program",
        help="Target program for the creation of input files in QPREP",
        default="gaussian",
        type=str,
        dest="program",
    )
    parser.add_argument(
        "--program_sp",
        help="Target program for the creation of input files in QPREP for single-points in QCORR",
        default="gaussian",
        type=str,
        dest="program_sp",
    )
    parser.add_argument(
        "--nprocs",
        help="Number of processors for the DFT calculations",
        default=24,
        type=int,
        dest="nprocs",
    )
    parser.add_argument(
        "--mem",
        help="Memory for the DFT calculations (i) Gaussian: total memory; (ii) ORCA: memory per processor",
        default="96GB",
        type=str,
        dest="mem",
    )
    parser.add_argument(
        "--aux_atoms_orca",
        default=[],
        help="List of atoms included in the aux part when using multiple basis sets in ORCA",
        dest="aux_atoms_orca",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--aux_genecp_bs",
        default=[],
        help="Auxiliary basis set for genecp/gen in ORCA",
        dest="aux_genecp_bs",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--aux_fit_gen_atoms",
        default=[],
        help="Fitting for the auxiliary basis set in ORCA (i.e. ['def2-TZVPP/C'])",
        dest="aux_fit_gen_atoms",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--cpcm_input",
        default="None",
        help="Additional lines for ORCA input files in the cpcm section. Format: ['LINE1','LINE2',etc]",
        dest="cpcm_input",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--orca_scf_iters",
        default=500,
        help="Number of SCF iterations in ORCA",
        dest="orca_scf_iters",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--mdci_orca",
        default="None",
        help="mdci section in ORCA. Format: ['LINE1','LINE2',etc]",
        dest="mdci_orca",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--print_mini_orca",
        action="store_true",
        default=True,
        help="Option to print 'mini' (reduced outputs) in ORCA",
    )
    parser.add_argument(
        "--qm_input",
        help="(i) keywords used in Gaussian input files (overiding opt and freq) or (ii) additional keywords for the ORCA input line",
        default="",
        dest="qm_input",
    )
    parser.add_argument(
        "--ts_input",
        help="OPT options in Gaussian in QCORR for input files of TSs (disable this option with None)",
        default="opt=(calcfc,noeigen,ts,maxstep=5)",
        dest="ts_input",
    )
    parser.add_argument(
        "--qm_input_end",
        help="Last input line for Gaussian",
        default="",
        dest="qm_input_end",
    )
    parser.add_argument(
        "--gen_atoms",
        help="Gen(ECP) atoms for Gaussian",
        default=None,
        dest="gen_atoms",
    )
    parser.add_argument(
        "--bs",
        default=None,
        help="Basis set for regular atoms during Gen(ECP) calculations",
        dest="bs",
        type=str,
    )
    parser.add_argument(
        "--bs_gen",
        default=None,
        help="Basis set for Gen(ECP) atoms",
        dest="bs_gen",
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
        dest="energy_threshold_for_gaussian",
    )
    parser.add_argument(
        "--chk",
        action="store_true",
        default=False,
        help="Create .chk files for Gaussian",
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
        "--dup",
        action="store",
        default=True,
        help="Remove duplicates after DFT optimization",
    )
    parser.add_argument(
        "--dup_threshold",
        action="store",
        default=0.0001,
        help="Energy difference (in Hartree) for E+ZPE, H and G considered in the duplicate filter",
        type=float,
    )
    parser.add_argument(
        "--amplitude_ifreq",
        action="store",
        default=0.2,
        help="Amplitude used to displace the imaginary frequencies to fix during analysis",
        type=float,
    )
    parser.add_argument(
        "--ifreq_cutoff",
        action="store",
        default=0.0,
        help="Cut off for imaginary frequencies during analysis",
        type=float,
    )
    parser.add_argument(
        "--s2_threshold",
        action="store",
        default=10.0,
        help="Cut off for spin contamination during analysis in \%\ of the expected value (i.e. multiplicity 3 has an the expected <S**2> of 2.0, if s2_threshold = 10 the <S**2> value is allowed to be 2.0 +- 0.2). Set s2_threshold = 0 to deactivate this option.",
        type=float,
    )
    parser.add_argument(
        "--isom",
        action="store",
        default=None,
        help="Checks that geometries mantain the same connectivity after xTB and DFT optimization",
        type=str,
    )
    parser.add_argument(
        "--vdwfrac",
        action="store",
        help="What fraction of the summed VDW radii constitutes a bond between two atoms in the isomerization filter?",
        default=0.50,
        type=float,
    )
    parser.add_argument(
        "--covfrac",
        action="store",
        help="What fraction of the summed covalent radii constitutes a bond between two atoms in the isomerization filter?",
        default=1.10,
        type=float,
    )
    parser.add_argument(
        "--nocheck",
        action="store_true",
        default=False,
        help="Skips analysis for QM output files (i.e. it generates inputs for all the files not only the normally terminated files",
    )
    parser.add_argument(
        "--nocom",
        action="store_true",
        default=False,
        help="Skips generation of QM input files to fix errors",
    )
    parser.add_argument(
        "--qcorr_json",
        action="store",
        default="",
        help="Starts QCORR from a previously generated json file (skips analysis)",
        type=str,
    )

    # writing single point files
    parser.add_argument(
        "--sp",
        help="Create QM single point input files",
        default=None,
        dest="sp",
        type=str,
    )
    parser.add_argument(
        "--nics", action="store_true", default=False, help="Create input files for NICS"
    )
    parser.add_argument(
        "--charge_sp",
        help="The charge for single point calculation in Gaussian and ORCA",
        default=None,
        metavar="charge_sp",
    )
    parser.add_argument(
        "--mult_sp",
        help="The multiplicity for single point calculation in Gaussian and ORCA",
        default=None,
        metavar="mult_sp",
    )
    parser.add_argument(
        "--qm_input_sp",
        help="Extra keywords for single-point calculations in Gaussian and ORCA",
        default="",
        dest="qm_input_sp",
    )
    parser.add_argument(
        "--qm_end",
        help="Last input line for single point calculations in Gaussian",
        default="",
        dest="qm_end",
    )
    parser.add_argument(
        "--bs_sp",
        default=None,
        help="Regular basis set(s) for single point calculations in Gaussian when Gen(ECP) is used",
        dest="bs_sp",
    )
    parser.add_argument(
        "--gen_bs_sp",
        default=None,
        help="Gen(ECP) basis set(s) for single point calculations in Gaussian",
        dest="gen_bs_sp",
    )
    parser.add_argument(
        "--suffix",
        help="The suffix for single point calculation in Gaussian and ORCA",
        default="None",
        type=str,
        metavar="suffix",
    )
    parser.add_argument(
        "--aux_atoms_orca_sp",
        default=[],
        help="List of atoms included in the aux part when using multiple basis sets in ORCA single-point calculations",
        dest="aux_atoms_orca_sp",
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--aux_gen_bs_sp",
        default=[],
        help="Auxiliary basis set for genecp/gen in ORCA single-point calculations",
        dest="aux_gen_bs_sp",
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--aux_fit_gen_atoms_sp",
        default=[],
        help="Fitting for the auxiliary basis set in ORCA single-point calculations (i.e. ['def2-TZVPP/C'])",
        dest="aux_fit_gen_atoms_sp",
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--cpcm_input_sp",
        default="None",
        help="Additional lines for ORCA single-point calculations in the cpcm section",
        dest="cpcm_input_sp",
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--orca_scf_iters_sp",
        default=500,
        help="Number of SCF iterations in ORCA single-point calculations",
        dest="orca_scf_iters_sp",
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--mdci_orca_sp",
        default="None",
        help="mdci section in ORCA single-point calculations",
        dest="mdci_orca_sp",
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--print_mini_orca_sp",
        action="store_true",
        default=True,
        help="Option to print 'mini' (reduced outputs) in ORCA single-point calculations",
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
        dest="dihedral",
        type=str,
        nargs=4,
        action="append",
    )
    parser.add_argument(
        "--bond",
        help="Specify the atom indexes to track bond lengths for different conformers",
        default=[],
        dest="bond",
        type=str,
        nargs=2,
        action="append",
    )
    parser.add_argument(
        "--angle",
        help="Specify the atom indexes to track angles for different conformers",
        default=[],
        dest="angle",
        type=str,
        nargs=3,
        action="append",
    )
    parser.add_argument(
        "--geom_par_name",
        action="store",
        dest="geom_par_name",
        default="descp",
        help="Change the prefix for the descriptors obtained",
    )
    parser.add_argument(
        "--dbstep_cen_lig_file",
        help="Center for DBSTEP steric paramters in a txt ( FORMAT : name, center, ligand)",
        action="store",
        default="No file passed",
        dest="dbstep_cen_lig_file",
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
        dest="nmr_nucleus",
        type=str,
        nargs="*",
    )
    parser.add_argument(
        "--nmr_slope",
        help="Specify the slope for each nucleus for nmr analysis default([1.0673,1.0759])",
        default=[1.0673, 1.0759],
        dest="nmr_slope",
        type=float,
        nargs="*",
    )
    parser.add_argument(
        "--nmr_intercept",
        help="Specify the intercept for each nucleus for nmr analysis default([-15.191,-2.2094])",
        default=[-15.191, -2.2094],
        dest="nmr_intercept",
        type=float,
        nargs="*",
    )
    parser.add_argument(
        "--nmr_tms_ref",
        help="Specify the reference for TMS for each nucleus for nmr analysis default([191.79,31.39])",
        default=[191.79, 31.39],
        dest="nmr_tms_ref",
        type=float,
        nargs="*",
    )

    # arguments for NICS
    parser.add_argument(
        "--nics_range",
        help="Range to calculate NICS along a given axis",
        default=4,
        dest="nics_range",
    )
    parser.add_argument(
        "--nics_number",
        help="Step size to calculate NICS along a given axis",
        default=16,
        dest="nics_number",
    )
    parser.add_argument(
        "--nics_atoms_file",
        help="NICS atoms in a txt ( FORMAT : name, atom1, atom2, atom3..)",
        action="store",
        default="No file passed",
        dest="nics_atoms_file",
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
        metavar="submission_command",
        type=str,
    )

    # apply exp rules
    parser.add_argument(
        "--geom_rules",
        dest="geom_rules",
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
        metavar="prefix",
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
        "qm_input_end_sp": "",
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
