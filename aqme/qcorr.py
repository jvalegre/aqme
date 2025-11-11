r"""
Parameters
----------

   files : list of str, default=''
      Filenames of QM output files to analyze. If *.log (or other strings that 
      are not lists such as *.out) are specified, the program will look for all 
      the log files in the working directory through glob.glob(*.log)
   w_dir_main : str, default=os.getcwd()
      Working directory
   fullcheck : bool, default=True
      Perform an analysis to detect whether the calculations were done 
      homogeneously (i.e. same level of theory, solvent, grid size, etc)
   varfile : str, default=None
      Option to parse the variables using a yaml file (specify the filename)
   ifreq_cutoff : float, default=0.0
      Cut off for to consider whether a frequency is imaginary (absolute of the 
      specified value is used)
   amplitude_ifreq : float, default=0.2
      Amplitude used to displace the imaginary frequencies to fix
   freq_conv : str, default=None
      If a string is defined, it will remove calculations that converged during 
      optimization but did not convergence in the subsequent frequency 
      calculation. Options: opt keyword as string (i.e. 'opt=(calcfc,maxstep=5)'). 
      If readfc is specified in the string, the chk option must be included as well.
   im_freq_input : str, default='opt=(calcfc,maxstep=5)' (Gaussian), '\n%geom\nCalc_Hess true\nMaxStep 0.05\nend' (ORCA)
      When extra imaginery frequencies are detected by QCORR, it automatically adds
      hessian calcs before starting geometry optimizations. This option can be 
      disabled using im_freq_input=None.
   s2_threshold : float, default=10.0
      Cut off for spin contamination during analysis in % of the expected value 
      (i.e. multiplicity 3 has an the expected <S**2> of 2.0, 
      if s2_threshold = 10, the <S**2> value is allowed to be 2.0 +- 0.2). 
      Set s2_threshold = 0 to deactivate this option.
   dup_threshold : float, default=0.0001
      Energy (in hartree) used as the energy difference in E, H and G to detect 
      duplicates
   ro_threshold : float, default=0.1
      Rotational constant value used as the threshold to detect duplicates 
   isom_type : str, default=None
      Check for isomerization from the initial input file to the resulting 
      output files. It requires the extension of the initial input files 
      (i.e. isom_type='com' or 'gjf') and the folder of the input files must be 
      added in the isom_inputs option
   isom_inputs : str, default=os.getcwd()
      Folder containing the initial input files to check for isomerization
   vdwfrac : float, default=0.50
      Fraction of the summed VDW radii that constitutes a bond between two atoms 
      in the isomerization filter
   covfrac : float, default=1.10
      Fraction of the summed covalent radii that constitutes a bond between two 
      atoms in the isomerization filter
   nodup_check : bool, default=False
      If True, the duplicate filter is disabled
      
.. note::

   New input files are generated through the QPREP module and, therefore, all 
   QPREP arguments can be used when calling QCORR and will overwrite default 
   options. For example, if the user specifies qm_input='wb97xd/def2svp', 
   all the new input files generated to fix issues will contain this keywords 
   line. See examples in the 'Example_workflows' folder for more information.

"""
######################################################.
#        This file stores the QCORR class            #
######################################################.

import os
import sys
import glob
import time
import pandas as pd
import json
import cclib
import numpy as np
from pathlib import Path
from aqme.utils import (
    move_file,
    get_info_input,
    load_variables,
    read_file,
    cclib_atoms_coords,
    check_files,
    check_dependencies
)
from aqme.qcorr_utils import (
    detect_linear,
    check_isomerization,
    full_check,
    get_json_data,
    get_cclib_params
)
from aqme.qprep import qprep


class qcorr:
    """
    Class containing all the functions from the QCORR module.

    Parameters
    ----------
    kwargs : argument class
        Specify any arguments from the QCORR module (for a complete list of variables, visit the AQME documentation)
    """

    def __init__(self, **kwargs):
        self._load_and_validate_variables(kwargs)
        self._check_dependencies_and_files()
        self._validate_file_formats()
        self.qcorr_processing()
        # this is added to avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

    def _load_and_validate_variables(self, kwargs):
        """Load default and user-specified variables, set number of processors"""
        self.args = load_variables(kwargs, "qcorr")
        # set number of processors
        if self.args.nprocs is None:
            self.args.nprocs = 8

    def _check_dependencies_and_files(self):
        """Check whether dependencies are installed and retrieve files to run in QCORR"""
        _ = check_dependencies(self)
        _ = check_files(self, 'qcorr')

    def _validate_file_formats(self):
        """Validate that input file formats are compatible with QCORR"""
        file_ext = os.path.basename(self.args.files[0]).split('.')[-1].lower()
        if file_ext not in ['log', 'out', 'json']:
            self.args.log.write(f"\nx  The format used ({file_ext}) is not compatible with QCORR! Formats accepted: log, out, json")
            self.args.log.finalize()
            sys.exit()

    def qcorr_processing(self):
        """
        General function of the QCORR module that:

        1. Analyzes the QM output files and moves output files with normal termination and no extra imaginary frequencies to the same folder
        2. Generates input files to fix errors and extra imaginary frequencies
        3. Generates input files with new keywords line(s) from the normally terminated files from point 1 (i.e. single-point energy corrections)
        """
        start_time_overall = time.time()
        file_terms, duplicate_data = self._initialize_tracking_data()
        destination_success_json = None

        self.args.log.write(f"o  Analyzing output files in {self.args.w_dir_main}\n")
        os.chdir(self.args.w_dir_main)

        # analyze files
        for file in sorted(self.args.files):
            file_name = os.path.basename(Path(file)).split(".")[0]
            destination_success_json = self._process_single_file(
                file, file_name, file_terms, duplicate_data, destination_success_json
            )

        # write information about the QCORR analysis in a csv
        csv_qcorr = self.write_qcorr_csv(file_terms)

        # performs a full analysis to ensure that the calcs were run with the same parameters
        self._perform_fullcheck_analysis(csv_qcorr, destination_success_json)
        self._finalize_qcorr_analysis(start_time_overall)

    def _initialize_tracking_data(self):
        """Initialize data structures for tracking file terminations and duplicates"""
        file_terms = {
            "finished": 0,
            "sp_calcs": 0,
            "extra_imag_freq": 0,
            "ts_no_imag_freq": 0,
            "freq_no_conv": 0,
            "spin_contaminated": 0,
            "duplicate_calc": 0,
            "atom_error": 0,
            "scf_error": 0,
            "no_data": 0,
            "linear_mol_wrong": 0,
            "not_specified": 0,
            "geom_qcorr": 0,
            "isomerized": 0,
        }

        duplicate_data = {
            "File": [],
            "Energies": [],
            "Enthalpies": [],
            "Gibbs": [],
            "RO_constant": [],
        }

        return file_terms, duplicate_data

    def _process_single_file(self, file, file_name, file_terms, duplicate_data, destination_success_json):
        """Process a single QM output file through the QCORR workflow"""
        # get initial cclib data and termination/error types and discard calcs with no data
        termination, errortype, cclib_data, file = self.cclib_init(file, file_name)

        if errortype in ["no_data", "atomicbasiserror"]:
            file_terms, destination = self.organize_outputs(
                file, termination, errortype, file_terms
            )
            if os.path.exists(file_name + ".json"):
                destination_json = destination.joinpath("json_files/")
                move_file(destination_json, self.args.w_dir_main, file_name + ".json")
            if errortype == "atomicbasiserror":
                self.args.log.write(f"{os.path.basename(file)}: Termination = {termination}, Error type = {errortype}")
            return destination_success_json

        # check for duplicates and fix wrong number of freqs in normally terminated calculations
        if termination == "normal":
            (
                atom_types,
                cartesians,
                duplicate_data,
                errortype,
                cclib_data,
                dup_off,
            ) = self.analyze_normal(duplicate_data, errortype, cclib_data, file_name)
        # fix calcs that did not terminated normally
        elif termination != "normal":
            atom_types, cartesians, cclib_data = self.analyze_abnormal(
                errortype, cclib_data
            )

        # check for isomerization
        if self.args.isom_type is not None:
            errortype = self.analyze_isom(file, cartesians, atom_types, errortype)

        # move initial QM input files (if the files are placed in the same folder as the output files)
        self._handle_input_file_movement(file_name, cclib_data)

        # create input files through QPREP to fix the errors (some errors require user intervention)
        if errortype not in [
            "ts_no_imag_freq",
            "isomerization",
            "duplicate_calc",
            "spin_contaminated",
            "none",
            "sp_calc",
        ]:
            self.qcorr_fixing(cclib_data, file, atom_types, cartesians)

        # This part places the calculations and json files in different folders depending on the type of termination
        if errortype == "duplicate_calc":
            self.args.log.write(f"{os.path.basename(file)}: Termination = {termination}, Error type = {errortype}, Duplicate of = {dup_off}")
        else:
            self.args.log.write(f"{os.path.basename(file)}: Termination = {termination}, Error type = {errortype}")

        file_terms, destination = self.organize_outputs(
            file, termination, errortype, file_terms
        )

        destination_json = destination.joinpath("json_files/")
        move_file(destination_json, self.args.w_dir_main, file_name + ".json")

        if errortype == "none":
            destination_success_json = destination_json

        return destination_success_json

    def _handle_input_file_movement(self, file_name, cclib_data):
        """Move initial QM input files if they are in the same folder as output files"""
        if cclib_data["metadata"]["QM program"].lower().find("gaussian") > -1:
            input_suffix = "com"
        elif cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
            input_suffix = "inp"

        if (
            os.path.exists(f"{self.args.w_dir_main}/{file_name}.{input_suffix}")
            and self.args.round_num == 1
        ):
            move_file(
                self.args.w_dir_main.joinpath("inputs/"),
                self.args.w_dir_main,
                f"{file_name}.com",
            )

    def _perform_fullcheck_analysis(self, csv_qcorr, destination_success_json):
        """Perform full analysis to ensure calculations were run with the same parameters"""
        # Currently, fullcheck is not working with ORCA calcs - we'll check this from json files if needed
        if self.args.fullcheck == "False":
            self.args.fullcheck = False
        elif self.args.fullcheck == "True":
            self.args.fullcheck = True

        # Check if we should skip fullcheck for ORCA by examining json files
        if self.args.fullcheck and destination_success_json is not None:
            # Quick check if any json files exist and if they're ORCA
            json_files = glob.glob(f"{destination_success_json}/*.json")
            if len(json_files) > 0:
                try:
                    with open(json_files[0]) as f:
                        test_data = json.load(f)
                        if test_data.get("metadata", {}).get("QM program", "").lower().find("orca") > -1:
                            self.args.fullcheck = False
                except:
                    pass

        if self.args.fullcheck:
            no_normal_terms = False
            try:
                df_qcorr = pd.read_csv(csv_qcorr)
                if df_qcorr["Normal termination"][0] > 0:
                    json_files = glob.glob(f"{destination_success_json}/*.json")
                    full_check(
                        w_dir_main=destination_success_json,
                        destination_fullcheck=destination_success_json,
                        files=json_files,
                        log=self.args.log,
                    )
                else:
                    no_normal_terms = True
            except UnboundLocalError:
                no_normal_terms = True
            if no_normal_terms:
                self.args.log.write("\nx  No normal terminations with no errors to run the full check analysis")

    def _finalize_qcorr_analysis(self, start_time_overall):
        """Finalize QCORR analysis with timing and log"""
        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\n Time QCORR: {elapsed_time} seconds\n")
        self.args.log.finalize()

    # include geom filters (ongoing work)

    # 			if len(self.args.geom) >= 1:
    # 				passing_rules = True
    # 				valid_mol_gen = True
    # 				self.args.log.write("  ----- geom filter(s) will be applied to the output file -----\n")
    # 				try:
    # 					format_file = file.split('.')[1]
    # 					mol = output_to_mol(file,format_file)
    # 					print_error_geom=False
    # 					if ob_compat and rdkit_compat:
    # 						passing_rules = geom_output(mol,self.args,self.args.log,file,print_error_geom)
    # 						if not passing_rules:
    # 							errortype = 'fail_geom'
    # 					os.remove(file.split('.')[0]+'.mol')
    # 				except AttributeError:
    # 					valid_mol_gen = False
    # 					os.remove(file.split('.')[0]+'.mol')
    # 					self.args.log.write("The file could not be converted into a mol object, geom filter(s) will be disabled\n")

    def cclib_init(self, file, file_name):
        """
        Determine termination and error types (initial determination), create json files
        with cclib and load the data in the cclib json files
        """
        # cclib generation of json files with ccwrite
        termination, errortype, cclib_data, file = self.cclib_gen(file, file_name)

        if errortype == "no_data":
            return termination, errortype, None, file

        # calculations with 1 atom
        if cclib_data["natom"] == 1:
            termination, errortype, cclib_data = self._handle_single_atom_calc(cclib_data)
        # general errors
        elif "vibfreqs" not in cclib_data:
            termination, errortype, cclib_data = self._handle_missing_vibfreqs(
                cclib_data, file, errortype
            )

        # normal terminations - check spin contamination
        if "vibfreqs" in cclib_data or errortype == "sp_calc":
            errortype = self._check_spin_contamination(cclib_data, errortype)

        return termination, errortype, cclib_data, file

    def _handle_single_atom_calc(self, cclib_data):
        """Handle calculations with a single atom"""
        cclib_data["vibfreqs"] = {}
        termination = "normal"
        errortype = "none"

        if "scfenergies" not in cclib_data:
            termination = "other"
            errortype = "not_specified"
        elif "freeenergy" not in cclib_data:
            errortype = "sp_calc"

        return termination, errortype, cclib_data

    def _handle_missing_vibfreqs(self, cclib_data, file, errortype):
        """Handle calculations with missing vibfreqs (general errors)"""
        termination = "other"
        errortype = "not_specified"

        if "optdone" in cclib_data:
            # if the optimization finished, only a freq job is required
            if cclib_data["optdone"] in [True, 'true']:
                errortype = "no_freq"

        # use very short reversed loop to find basis set incompatibilities and SCF errors
        outlines = read_file(os.getcwd(), self.args.w_dir_main, file)
        termination, errortype, cclib_data = self._detect_basis_scf_errors(
            outlines, termination, errortype, cclib_data
        )

        return termination, errortype, cclib_data

    def _detect_basis_scf_errors(self, outlines, termination, errortype, cclib_data):
        """Detect basis set incompatibilities and SCF errors from output file"""
        for i in reversed(range(len(outlines) - 15, len(outlines))):
            if (
                outlines[i].find("Normal termination") > -1
                and errortype != "no_freq"
            ):
                termination = "normal"
                errortype = "sp_calc"
                cclib_data["metadata"][
                    "ground or transition state"
                ] = "SP calculation"
                break
            elif (
                outlines[i - 1].find("Atomic number out of range") > -1
                or outlines[i - 1].find("basis sets are only available") > -1
            ):
                errortype = "atomicbasiserror"
                break
            elif outlines[i].find("SCF Error") > -1:
                errortype = "SCFerror"
                break

        return termination, errortype, cclib_data

    def _check_spin_contamination(self, cclib_data, errortype):
        """Check for spin contamination using user-defined thresholds"""
        if "S2 after annihilation" not in cclib_data:
            return errortype

        unpaired_e = cclib_data["mult"] - 1
        # this first part accounts for singlet diradicals (threshold is 10% of the spin before annihilation)
        if unpaired_e == 0:
            if (
                float(cclib_data["S2 after annihilation"])
                > abs(float(self.args.s2_threshold) / 100)
                * cclib_data["S2 before annihilation"]
            ):
                errortype = "spin_contaminated"
        else:
            spin = unpaired_e * 0.5
            s2_expected_value = spin * (spin + 1)
            spin_diff = abs(
                float(cclib_data["S2 after annihilation"])
                - s2_expected_value
            )
            if (
                spin_diff
                > abs(float(self.args.s2_threshold) / 100) * s2_expected_value
            ):
                errortype = "spin_contaminated"

        return errortype 

    def analyze_normal(self, duplicate_data, errortype, cclib_data, file_name):
        """
        Analyze errors from normally terminated calculations
        """
        # retrieve previous successful results in case new calculations are duplicates
        duplicate_data = self._load_previous_success_data(duplicate_data, errortype)

        atom_types, cartesians = cclib_atoms_coords(cclib_data, -1)
        dup_off = None

        if errortype == "none":
            errortype, dup_off = self._check_for_duplicates(
                duplicate_data, cclib_data, errortype, file_name
            )

        if errortype == "none":
            errortype = self._analyze_imaginary_frequencies(
                duplicate_data, cclib_data, atom_types, file_name, cartesians
            )
            errortype = self._detect_freq_convergence_issues(cclib_data, errortype)

        if errortype in ["extra_imag_freq", "freq_no_conv", "linear_mol_wrong"]:
            cartesians, cclib_data = self._apply_frequency_fixes(
                errortype, cclib_data, cartesians
            )

        return atom_types, cartesians, duplicate_data, errortype, cclib_data, dup_off

    def _load_previous_success_data(self, duplicate_data, errortype):
        """Load previous successful results to check for duplicates"""
        if self.args.resume_qcorr:
            destination_json = self.args.w_dir_main.joinpath("../../../success/json_files/")
        else:
            destination_json = self.args.w_dir_main.joinpath("success/json_files/")

        if os.path.exists(destination_json):
            previous_success = glob.glob(f"{destination_json}/*.json")
            for previous_json in previous_success:
                with open(previous_json) as json_file:
                    cclib_data_json = json.load(json_file)
                E_json, H_json, G_json, ro_json, _ = get_cclib_params(cclib_data_json, errortype)
                duplicate_data["File"].append(os.path.basename(previous_json))
                duplicate_data["Energies"].append(E_json)
                duplicate_data["Enthalpies"].append(H_json)
                duplicate_data["Gibbs"].append(G_json)
                duplicate_data["RO_constant"].append(ro_json)

        return duplicate_data

    def _check_for_duplicates(self, duplicate_data, cclib_data, errortype, file_name):
        """Check if the current calculation is a duplicate of a previous one"""
        E_dup, H_dup, G_dup, ro_dup, errortype = get_cclib_params(cclib_data, errortype)
        dup_off = None

        if self.args.nodup_check == False:
            # detects if this calculation is a duplicate
            for i, _ in enumerate(duplicate_data["Energies"]):
                E_diff = abs(E_dup - duplicate_data["Energies"][i])
                H_diff = abs(H_dup - duplicate_data["Enthalpies"][i])
                G_diff = abs(G_dup - duplicate_data["Gibbs"][i])
                ro_diff = 0
                if (ro_dup is not None) and (
                    duplicate_data["RO_constant"][i] is not None
                ):
                    ro_diff = np.linalg.norm(
                        np.array(ro_dup) - np.array(duplicate_data["RO_constant"][i])
                    )
                if max([E_diff, H_diff, G_diff]) < abs(float(self.args.dup_threshold)):
                    if (ro_dup is not None) and (ro_diff < self.args.ro_threshold):
                        errortype = "duplicate_calc"
                        dup_off = duplicate_data["File"][i]
                        break

        # Store this calculation's data only if it's not a duplicate
        # (moved outside to be called from analyze_imaginary_frequencies)
        return errortype, dup_off

    def _analyze_imaginary_frequencies(self, duplicate_data, cclib_data, atom_types, file_name, cartesians):
        """Analyze imaginary frequencies and detect issues"""
        # First, get and store the current calculation's energies
        E_dup, H_dup, G_dup, ro_dup, errortype = get_cclib_params(cclib_data, "none")
        duplicate_data["File"].append(file_name)
        duplicate_data["Energies"].append(E_dup)
        duplicate_data["Enthalpies"].append(H_dup)
        duplicate_data["Gibbs"].append(G_dup)
        duplicate_data["RO_constant"].append(ro_dup)

        # Count imaginary frequencies
        initial_ifreqs = 0
        for freq in cclib_data["vibfreqs"]:
            if float(freq) < 0 and abs(float(freq)) > abs(
                float(self.args.ifreq_cutoff)
            ):
                initial_ifreqs += 1

        # exclude TS imag frequency
        if (
            cclib_data["metadata"]["ground or transition state"]
            == "transition_state"
        ):
            initial_ifreqs -= 1

        # gives new coordinates by displacing the normal mode(s) of the negative freq(s)
        if initial_ifreqs > 0:
            errortype = "extra_imag_freq"
        elif initial_ifreqs < 0:
            errortype = "ts_no_imag_freq"

        # Check for linear molecules
        if len(atom_types) in [3, 4]:
            errortype = detect_linear(errortype, atom_types, cclib_data)

        return errortype

    def _detect_freq_convergence_issues(self, cclib_data, errortype):
        """Detect frequency convergence issues (Gaussian only)"""
        if "gaussian" in cclib_data["metadata"]["QM program"].lower():
            if self.args.freq_conv is not None:
                if (
                    errortype == "none"
                    and cclib_data["opt times converged"] == 1
                ):
                    errortype = "freq_no_conv"
        return errortype

    def _apply_frequency_fixes(self, errortype, cclib_data, cartesians):
        """Apply fixes for frequency-related errors"""
        if errortype == "extra_imag_freq":
            cartesians = self.fix_imag_freqs(cclib_data, cartesians)

        # in case no previous OPT was done (only works if it's not a TS)
        opt_found = False
        for keyword in cclib_data["metadata"]["keywords line"].split():
            if keyword.lower().startswith("opt"):
                opt_found = True

        if not opt_found:
            cclib_data["metadata"]["keywords line"] += " opt"

        # adding the Hessian calculation before OPT increases the rate of success
        if errortype in ["freq_no_conv", "extra_imag_freq"]:
            cclib_data = self._add_hessian_to_keywords(errortype, cclib_data)
        elif errortype == "linear_mol_wrong":
            cclib_data["metadata"]["keywords line"] += " symmetry=(PG=Cinfv)"

        return cartesians, cclib_data

    def _add_hessian_to_keywords(self, errortype, cclib_data):
        """Add Hessian calculation to keywords line"""
        new_opt = None
        if errortype == "freq_no_conv" and self.args.freq_conv not in [None, 'None']:
            new_opt = self.args.freq_conv
        elif errortype == "extra_imag_freq" and self.args.im_freq_input not in [None, 'None']:
            new_opt = self.args.im_freq_input

        if cclib_data["metadata"]["QM program"].lower().find("gaussian") > -1:
            new_keywords_line = ""
            for keyword in cclib_data["metadata"]["keywords line"].split():
                if keyword.lower().startswith("opt"):
                    if new_opt is not None:
                        keyword = new_opt
                        if cclib_data["metadata"]["ground or transition state"] == "transition_state":
                            keyword = keyword.replace("=(", "=(ts,noeigen,")
                new_keywords_line += keyword
                new_keywords_line += " "
        elif cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
            new_keywords_line = cclib_data["metadata"]["keywords line"]
            if self.args.im_freq_input not in [None, 'None']:
                # change the default value
                if self.args.im_freq_input == 'opt=(calcfc,maxstep=5)':
                    self.args.im_freq_input = '\n%geom\nCalc_Hess true\nMaxStep 0.05\nend'
                new_keywords_line += self.args.im_freq_input

        cclib_data["metadata"]["keywords line"] = new_keywords_line
        return cclib_data

    def analyze_abnormal(self, errortype, cclib_data):
        """
        Analyze errors from calculations that did not finish normally
        """
        program = self._determine_qm_program(cclib_data)

        # for calcs with finished OPT but no freqs, adjust the keywords so only FREQ is calculated
        if errortype == "no_freq":
            atom_types, cartesians, cclib_data = self._fix_no_freq_error(cclib_data)
        else:
            # help to fix SCF convergence errors
            if errortype == "SCFerror":
                cclib_data = self._fix_scf_error(program, cclib_data)

            if errortype in ["not_specified", "SCFerror"]:
                atom_types, cartesians = self._find_best_geometry(cclib_data)

        return atom_types, cartesians, cclib_data

    def _determine_qm_program(self, cclib_data):
        """Determine which QM program was used"""
        if cclib_data["metadata"]["QM program"].lower().find("gaussian") > -1:
            return 'gaussian'
        elif cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
            return "orca"
        return None

    def _fix_no_freq_error(self, cclib_data):
        """Fix calculations that finished OPT but have no frequency calculation"""
        new_keywords_line = ""
        for keyword in cclib_data["metadata"]["keywords line"].split():
            if not keyword.lower().startswith("opt"):
                new_keywords_line += keyword
                new_keywords_line += " "
        cclib_data["metadata"]["keywords line"] = new_keywords_line
        atom_types, cartesians = cclib_atoms_coords(cclib_data, -1)
        return atom_types, cartesians, cclib_data

    def _fix_scf_error(self, program, cclib_data):
        """Apply SCF convergence fixes based on the QM program"""
        if program == 'gaussian':
            cclib_data = self._fix_scf_error_gaussian(cclib_data)
        elif program == 'orca':
            cclib_data = self._fix_scf_error_orca(cclib_data)
        return cclib_data

    def _fix_scf_error_gaussian(self, cclib_data):
        """Fix SCF convergence errors for Gaussian"""
        if (
            cclib_data["metadata"]["keywords line"].find(" scf=xqc") > -1
            or cclib_data["metadata"]["keywords line"].find(" scf=qc") > -1
        ):
            new_keywords_line = ""
            for keyword in cclib_data["metadata"]["keywords line"].split():
                if keyword == "scf=xqc":
                    keyword = "scf=qc"
                new_keywords_line += keyword
                new_keywords_line += " "
            cclib_data["metadata"]["keywords line"] = new_keywords_line
        else:
            cclib_data["metadata"]["keywords line"] += " scf=xqc"
        return cclib_data

    def _fix_scf_error_orca(self, cclib_data):
        """Fix SCF convergence errors for ORCA"""
        if 'SlowConv' not in cclib_data["metadata"]["keywords line"]:
            cclib_data["metadata"]["keywords line"] = 'SlowConv ' + cclib_data["metadata"]["keywords line"]
        return cclib_data

    def _find_best_geometry(self, cclib_data):
        """Find the geometry with the best (lowest) RMS forces"""
        if "geovalues" in cclib_data:
            RMS_forces = [row[1] for row in cclib_data["geovalues"]]
            # cclib uses None when the values are corrupted in the output files, replace None for a large number
            RMS_forces = [10000 if val is None else val for val in RMS_forces]
            min_RMS = RMS_forces.index(min(RMS_forces))
        else:
            # for optimizations that fail in the first step
            min_RMS = 0

        atom_types, cartesians = cclib_atoms_coords(cclib_data, min_RMS)
        return atom_types, cartesians


    def analyze_isom(self, file, cartesians, atom_types, errortype):
        """
        Check if the initial structure isomerized during QM geometry optimization
        """
        if not self._validate_isom_directory():
            os.chdir(self.args.initial_dir)
            self.args.log.finalize()
            sys.exit()

        try:
            atoms_com, coords_com = self._load_initial_structure(file)
            isom_data = self._prepare_isom_data(atoms_com, coords_com, cartesians, atom_types)
            isomerized = self._perform_isomerization_check(isom_data, file)
        except FileNotFoundError:
            self.args.log.write(f"x  No com file were found for {os.path.basename(file)}, the check_geom test will be disabled for this calculation")
            isomerized = False

        if isomerized:
            errortype = "isomerization"

        os.chdir(self.args.w_dir_main)
        return errortype

    def _validate_isom_directory(self):
        """Validate that the isomerization input directory exists"""
        try:
            os.chdir(self.args.isom_inputs)
            return True
        except FileNotFoundError:
            self.args.log.write("x  The PATH specified in isom_inputs doesn't exist!")
            return False

    def _load_initial_structure(self, file):
        """Load the initial structure from input file or CSV"""
        atoms_com, coords_com, atoms_and_coords = [], [], []
        init_csv = pd.DataFrame()

        if len(self.args.isom_type.split(".")) == 1:
            atoms_and_coords, _, _ = get_info_input(
                f'{os.path.basename(Path(file)).split(".")[0]}.{self.args.isom_type}'
            )
        elif self.args.isom_type.split(".")[1] != "csv":
            init_csv = pd.read_csv(self.args.isom_type)

        for line in atoms_and_coords:
            atoms_com.append(line.split()[0])
            coords_com.append(
                [
                    float(line.split()[1]),
                    float(line.split()[2]),
                    float(line.split()[3]),
                ]
            )

        return atoms_com, coords_com

    def _prepare_isom_data(self, atoms_com, coords_com, cartesians, atom_types):
        """Prepare data dictionary for isomerization check"""
        init_csv = pd.DataFrame()
        if len(self.args.isom_type.split(".")) > 1 and self.args.isom_type.split(".")[1] == "csv":
            init_csv = pd.read_csv(self.args.isom_type)

        isom_data = {
            "Coords input": coords_com,
            "Coords output": cartesians,
            "Atoms input": atoms_com,
            "Atoms output": atom_types,
            "VdW radii fraction": self.args.vdwfrac,
            "Covalent radii fraction": self.args.covfrac,
            "Initial csv": init_csv,
        }
        return isom_data

    def _perform_isomerization_check(self, isom_data, file):
        """Perform the actual isomerization check"""
        return check_isomerization(isom_data, file)

    def qcorr_fixing(self, cclib_data, file, atom_types, cartesians):
        """
        Create com files for resubmission with the suggested protocols to correct the errors
        """
        cclib_data = self._apply_user_parameters(cclib_data)
        destination_fix = self._determine_output_destination()
        program = self._determine_qm_program(cclib_data)

        if program in ["gaussian", "orca"]:
            if program == 'gaussian':
                self._validate_gen_basis_set(cclib_data)
            self._generate_qprep_input(
                destination_fix, cclib_data, file, atom_types, cartesians, program
            )
        else:
            self.args.log.write(f"x  Couldn't create an input file to fix {os.path.basename(file)} (compatible programs: Gaussian and ORCA)\n")

    def _apply_user_parameters(self, cclib_data):
        """Apply user-defined keywords line, memory and number of processors"""
        # user-defined keywords line, mem and nprocs overwrites previously used parameters
        if self.args.qm_input != "":
            cclib_data["metadata"]["keywords line"] = self.args.qm_input

        if self.args.mem != "16GB":
            cclib_data["metadata"]["memory"] = self.args.mem
        elif "memory" not in cclib_data["metadata"]:
            cclib_data["metadata"]["memory"] = "16GB"

        if self.args.nprocs != 8:
            cclib_data["metadata"]["processors"] = self.args.nprocs
        elif "processors" not in cclib_data["metadata"]:
            cclib_data["metadata"]["processors"] = 8

        return cclib_data

    def _validate_gen_basis_set(self, cclib_data):
        """Validate gen/genecp basis set specifications for Gaussian"""
        for keyword in cclib_data["metadata"]["keywords line"].split():
            for subkey in keyword.split('/'):
                if subkey in ['gen', 'genecp']:
                    if '' in [self.args.bs_gen, self.args.bs_nogen]:
                        self.args.log.write("x  WARNING! You are using gen(ECP) but you are not specifying two basis sets. Please, add them with the bs_gen and bs_nogen options.")
                        self.args.log.finalize()
                        sys.exit()
                    elif len(self.args.gen_atoms) == 0:
                        self.args.log.write("x  WARNING! You are using gen(ECP) but you are not specifying the atoms included for gen(ECP). Please, add them with the gen_atoms option.")
                        self.args.log.finalize()
                        sys.exit()

    def _determine_output_destination(self):
        """Determine destination folder for fixed QM input files"""
        if self.args.resume_qcorr:
            return Path(
                f"{self.args.w_dir_main}/../../run_{self.args.round_num}/fixed_QM_inputs"
            )
        else:
            return Path(
                f"{self.args.w_dir_main}/failed/run_{self.args.round_num}/fixed_QM_inputs"
            )

    def _generate_qprep_input(self, destination_fix, cclib_data, file, atom_types, cartesians, program):
        """Generate QM input file using qprep"""
        qprep(
            destination=destination_fix,
            w_dir_main=self.args.w_dir_main,
            files=os.path.basename(file),
            charge=cclib_data["charge"],
            mult=cclib_data["mult"],
            program=program,
            atom_types=atom_types,
            cartesians=cartesians,
            qm_input=cclib_data["metadata"]["keywords line"],
            mem=cclib_data["metadata"]["memory"],
            nprocs=cclib_data["metadata"]["processors"],
            chk=self.args.chk,
            qm_end=self.args.qm_end,
            bs_gen=self.args.bs_gen,
            bs_nogen=self.args.bs_nogen,
            gen_atoms=self.args.gen_atoms,
            create_dat=False,
        )

    def cclib_gen(self, file, file_name):
        """
        Create a json file with cclib and load a dictionary
        """

        termination, errortype = "normal", "none"
        cclib_data = {}

        if not os.path.exists(file):
            file = f'{self.args.initial_dir}/{file}'
        try:
                cclib_data = cclib.io.ccopen(file).parse()
                # Convert to dictionary
                cclib_data = cclib_data.__dict__.copy()

        except ValueError:
            termination = "other"
            errortype = "no_data"

        if "natom" not in cclib_data:
            termination = "other"
            errortype = "no_data"

        # add parameters that might be missing from cclib (depends on the version)
        if errortype != "no_data":
            cclib_data = get_json_data(self, file, cclib_data)

        # this is just a "dirty hack" until cclib is updated to be compatible for print mini in ORCA
        if hasattr(cclib_data, "metadata"):
            if cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
                if "final single point energy" in cclib_data["energy"]:
                    termination, errortype = "normal", "none"

        if errortype == "no_data":
            self.args.log.write(f"x  Potential cclib compatibility problem or no data found for file {file_name} (Termination = {termination}, Error type = {errortype})")

        return termination, errortype, cclib_data, file

    def fix_imag_freqs(self, cclib_data, cartesians):
        """
        Fixes undersired (extra) imaginary frequencies from QM calculations. This function multiplies the imaginary normal mode vectors by the selected amplitude (0.2 is the default amplitude in the pyQRC script from GitHub, user: bobbypaton).	By default, all the extra imaginary modes are used (i.e. in calculations with three	extra imaginary frequencies, all the three modes will be used to displace the atoms). This can be tuned with the --ifreq_cutoff option (i.e. only use freqs lower than -50 cm-1).

        Parameters
        ----------
        cclib_data : cclib object
            Variables parsed with cclib
        cartesians : list of lists
            List of lists containing the molecular coordinates as floats

        Returns
        -------
        cartesians : list of lists
            New set of cartesian coordinates generated after displacing the original coordinates along the normal modes of the corresponding imaginary frequencies
        """

        shift = []

        # could get rid of atomic units here, if zpe_rat definition is changed
        for mode, _ in enumerate(cclib_data["vibfreqs"]):
            # moves along all imaginary freqs (ignoring the TS imag freq, assumed to be the most negative freq)
            if (
                mode == 0
                and cclib_data["metadata"]["ground or transition state"]
                == "transition_state"
            ):
                shift.append(0.0)
            else:
                if cclib_data["vibfreqs"][mode] < 0.0:
                    shift.append(float(self.args.amplitude_ifreq))
                else:
                    shift.append(0.0)

            # The starting geometry is displaced along each normal mode according to the random shift
            for atom in range(0, cclib_data["natom"]):
                for coord in range(0, 3):
                    cartesians[atom][coord] = (
                        cartesians[atom][coord]
                        + cclib_data["vibdisps"][mode][atom][coord]
                        * shift[mode]
                    )

        return cartesians

    def organize_outputs(self, file, termination, errortype, file_terms):
        """
        1. Moves the QM output files to their corresponding folders after the analysis.
        2. Keeps track of the number of calculations with the different types of terminations and error types

        Parameters
        ----------
        file : str
            Output file
        termination : string
            Type of termination of the QM output file (i.e. normal, error, unfinished)
        errortype : string
            Type of error type of the QM output file (i.e. None, not_specified, extra_imag_freq, etc)
        file_terms : dict
            Keeps track of the number of calculations for each termination and error type

        Returns
        -------
        file_terms : dict
            Keeps track of the number of calculations for each termination and error type
        """
        destination = self._determine_destination_folder(termination, errortype)
        file_terms = self._update_file_statistics(errortype, file_terms)
        self._move_output_file(file, destination)

        return file_terms, destination

    def _determine_destination_folder(self, termination, errortype):
        """Determine the destination folder based on termination and error type"""
        if self.args.resume_qcorr:
            destination_error = self.args.w_dir_main.joinpath(
                f"../../run_{self.args.round_num}/"
            )
            destination_normal = self.args.w_dir_main.joinpath("../../../success/")
        else:
            destination_error = self.args.w_dir_main.joinpath(
                f"failed/run_{self.args.round_num}/"
            )
            destination_normal = self.args.w_dir_main.joinpath("success/")

        # Map error types to their destination folders
        error_destinations = {
            "none": (destination_normal, termination == "normal"),
            "sp_calc": (destination_normal.joinpath("SP_calcs/"), termination == "normal"),
            "extra_imag_freq": (destination_error.joinpath("extra_imag_freq/"), True),
            "ts_no_imag_freq": (destination_error.joinpath("ts_no_imag_freq/"), True),
            "spin_contaminated": (destination_error.joinpath("spin_contaminated/"), True),
            "duplicate_calc": (destination_error.joinpath("duplicates/"), True),
            "atomicbasiserror": (destination_error.joinpath("error/basis_set_error/"), True),
            "SCFerror": (destination_error.joinpath("error/scf_error/"), True),
            "no_data": (destination_error.joinpath("error/no_data/"), True),
            "fail_geom": (destination_error.joinpath("geom_filter/"), True),
            "isomerization": (destination_error.joinpath("isomerization/"), True),
            "freq_no_conv": (destination_error.joinpath("freq_no_conv/"), True),
            "linear_mol_wrong": (destination_error.joinpath("linear_mol_wrong/"), True),
        }

        dest_tuple = error_destinations.get(errortype)
        if dest_tuple and dest_tuple[1]:
            return dest_tuple[0]
        else:
            return destination_error.joinpath("error/not_specified_error/")

    def _update_file_statistics(self, errortype, file_terms):
        """Update file statistics counter based on error type"""
        error_counters = {
            "none": "finished",
            "sp_calc": "sp_calcs",
            "extra_imag_freq": "extra_imag_freq",
            "ts_no_imag_freq": "ts_no_imag_freq",
            "spin_contaminated": "spin_contaminated",
            "duplicate_calc": "duplicate_calc",
            "atomicbasiserror": "atom_error",
            "SCFerror": "scf_error",
            "no_data": "no_data",
            "fail_geom": "geom_qcorr",
            "isomerization": "isomerized",
            "freq_no_conv": "freq_no_conv",
            "linear_mol_wrong": "linear_mol_wrong",
        }

        counter_key = error_counters.get(errortype, "not_specified")
        file_terms[counter_key] += 1
        return file_terms

    def _move_output_file(self, file, destination):
        """Move output file to the determined destination"""
        move_file(destination, self.args.w_dir_main, os.path.basename(file))

    def write_qcorr_csv(self, file_terms):
        """
        Write information about the QCORR analysis in a csv
        """

        ana_data = pd.DataFrame()
        ana_data.at[0, "Total files"] = len(self.args.files)
        ana_data.at[0, "Normal termination"] = file_terms["finished"]
        ana_data.at[0, "Single-point calcs"] = file_terms["sp_calcs"]
        ana_data.at[0, "Extra imag. freq."] = file_terms["extra_imag_freq"]
        ana_data.at[0, "TS with no imag. freq."] = file_terms["ts_no_imag_freq"]
        ana_data.at[0, "Freq not converged"] = file_terms["freq_no_conv"]
        ana_data.at[0, "Linear mol with wrong n of freqs"] = file_terms[
            "linear_mol_wrong"
        ]
        ana_data.at[0, "SCF error"] = file_terms["scf_error"]
        ana_data.at[0, "No data"] = file_terms["no_data"]
        ana_data.at[0, "Basis set error"] = file_terms["atom_error"]
        ana_data.at[0, "Other errors"] = file_terms["not_specified"]
        if float(self.args.s2_threshold) > 0.0:
            ana_data.at[0, "Spin contamination"] = file_terms["spin_contaminated"]
        ana_data.at[0, "Duplicates"] = file_terms["duplicate_calc"]
        if len(self.args.geom) >= 1:
            ana_data.at[0, "geom filter"] = file_terms["geom_qcorr"]
        if self.args.isom_type is not None:
            ana_data.at[0, "Isomerization"] = file_terms["isomerized"]
        path_as_str = self.args.initial_dir.as_posix()
        csv_qcorr = path_as_str + f"/QCORR-run_{self.args.round_num}-stats.csv"
        if self.args.verbose:
            ana_data.to_csv(csv_qcorr, index=False)

        return csv_qcorr
