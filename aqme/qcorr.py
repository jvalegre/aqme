"""
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
   s2_threshold : float, default=10.0
      Cut off for spin contamination during analysis in % of the expected value 
      (i.e. multiplicity 3 has an the expected <S**2> of 2.0, 
      if s2_threshold = 10, the <S**2> value is allowed to be 2.0 +- 0.2). 
      Set s2_threshold = 0 to deactivate this option.
   dup_threshold : float, default=0.0001
      Energy (in hartree) used as the energy difference in E, H and G to detect 
      duplicates
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
import subprocess
import numpy as np

try:
    import cclib
except ModuleNotFoundError:
    print("x  cclib is not installed! You can install the program with 'conda install -c conda-forge cclib' or 'pip install cclib'")
    sys.exit()
from pathlib import Path
from aqme.utils import (
    move_file,
    QM_coords,
    get_info_input,
    load_variables,
    read_file,
    cclib_atoms_coords,
)
from aqme.qcorr_utils import (
    detect_linear,
    check_isomerization,
    full_check,
    get_json_data,
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

        # load default and user-specified variables
        self.args = load_variables(kwargs, "qcorr")

        if len(self.args.files) == 0:
            self.args.log.write('\nx  No files were found! Make sure you use quotation marks if you are using * (i.e. --files "*.log")')
            self.args.log.finalize()
            sys.exit()

        # QCORR analysis
        if self.args.files[0].split('.')[1].lower() not in ['log','out','json']:
            self.args.log.write(f"\nx  The format used ({self.args.files[0].split('.')[1].lower()}) is not compatible with QCORR! Formats accepted: log, out, json")
            self.args.log.finalize()
            sys.exit()

        self.qcorr_processing()

        # this is added to avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

    def qcorr_processing(self):
        """
        General function of the QCORR module that:

        1. Analyzes the QM output files and moves output files with normal termination and no extra imaginary frequencies to the same folder
        2. Generates input files to fix errors and extra imaginary frequencies
        3. Generates input files with new keywords line(s) from the normally terminated files from point 1 (i.e. single-point energy corrections)
        """

        start_time_overall = time.time()

        # generate some data
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
            "geom_rules_qcorr": 0,
            "isomerized": 0,
        }

        duplicate_data = {
            "File": [],
            "Energies": [],
            "Enthalpies": [],
            "Gibbs": [],
            "RO_constant": [],
        }

        self.args.log.write(f"o  Analyzing output files in {self.args.w_dir_main}\n")
        os.chdir(self.args.w_dir_main)
        # analyze files
        for file in sorted(self.args.files):
            # get initial cclib data and termination/error types and discard calcs with no data
            file_name = os.path.basename(file).split(".")[0]
            termination, errortype, cclib_data, outlines, file = self.cclib_init(
                file, file_name
            )
            if errortype in ["no_data", "atomicbasiserror"]:
                file_terms, _ = self.organize_outputs(
                    file, termination, errortype, file_terms
                )
                if errortype == "atomicbasiserror":
                    os.remove(file_name + ".json")
                    self.args.log.write(f"{os.path.basename(file)}: Termination = {termination}, Error type = {errortype}")
                continue

            # check for duplicates and fix wrong number of freqs in normally terminated calculations and
            elif termination == "normal":
                (
                    atom_types,
                    cartesians,
                    duplicate_data,
                    errortype,
                    cclib_data,
                    dup_off,
                ) = self.analyze_normal(
                    duplicate_data, errortype, cclib_data, file_name
                )

            # fix calcs that did not terminated normally

            elif termination != "normal":
                atom_types, cartesians, cclib_data = self.analyze_abnormal(
                    errortype, cclib_data, outlines
                )

            # check for isomerization
            if self.args.isom_type is not None:
                errortype = self.analyze_isom(file, cartesians, atom_types, errortype)

            # move initial QM input files (if the files are placed in the same folder as the output files)
            if (
                os.path.exists(f"{self.args.w_dir_main}/{file_name}.com")
                and self.args.round_num == 1
            ):
                move_file(
                    self.args.w_dir_main.joinpath("inputs/"),
                    self.args.w_dir_main,
                    f"{file_name}.com",
                )

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

            if errortype in ["none", "sp_calc"]:
                destination_json = destination.joinpath("json_files/")
                move_file(destination_json, self.args.w_dir_main, file_name + ".json")
            else:
                os.remove(file_name + ".json")

            # write information about the QCORR analysis in a csv
            csv_qcorr = self.write_qcorr_csv(file_terms)

        # performs a full analysis to ensure that the calcs were run with the same parameters
        if self.args.fullcheck == "False":
            self.args.fullcheck = False
        elif self.args.fullcheck == "True":
            self.args.fullcheck = True
        if self.args.fullcheck:
            no_normal_terms = False
            try:
                df_qcorr = pd.read_csv(csv_qcorr)
                if df_qcorr["Normal termination"][0] > 0:
                    json_files = glob.glob(f"{destination_json}/*.json")
                    full_check(
                        w_dir_main=destination_json,
                        destination_fullcheck=destination_json,
                        files=json_files,
                        log=self.args.log,
                    )
                else:
                    no_normal_terms = True
            except UnboundLocalError:
                no_normal_terms = True
            if no_normal_terms:
                self.args.log.write("\nx  No normal terminations with no errors to run the full check analysis")

        elapsed_time = round(time.time() - start_time_overall, 2)
        self.args.log.write(f"\n Time QCORR: {elapsed_time} seconds\n")
        self.args.log.finalize()

        # NOT needed as already created in initial_dir
        # move dat and csv file containing the QCORR information if this is a sequential QCORR analysis
        # if self.args.resume_qcorr:
        #     destination_data = self.args.w_dir_main.joinpath("../../../")
        #     move_file(
        #         destination_data,
        #         self.args.w_dir_main,
        #         f"QCORR-run_{self.args.round_num}.dat",
        #     )
        #     move_file(
        #         destination_data,
        #         self.args.w_dir_main,
        #         f"QCORR-run_{self.args.round_num}-stats.csv",
        #     )

    # include geom filters (ongoing work)

    # 			if len(self.args.geom_rules) >= 1:
    # 				passing_rules = True
    # 				valid_mol_gen = True
    # 				self.args.log.write("  ----- geom_rules filter(s) will be applied to the output file -----\n")
    # 				try:
    # 					format_file = file.split('.')[1]
    # 					mol = output_to_mol(file,format_file)
    # 					print_error_geom_rules=False
    # 					if ob_compat and rdkit_compat:
    # 						passing_rules = geom_rules_output(mol,self.args,self.args.log,file,print_error_geom_rules)
    # 						if not passing_rules:
    # 							errortype = 'fail_geom_rules'
    # 					os.remove(file.split('.')[0]+'.mol')
    # 				except AttributeError:
    # 					valid_mol_gen = False
    # 					os.remove(file.split('.')[0]+'.mol')
    # 					self.args.log.write("The file could not be converted into a mol object, geom_rules filter(s) will be disabled\n")

    def cclib_init(self, file, file_name):
        """
        Determine termination and error types (initial determination), create json files
        with cclib and load the data in the cclib json files
        """

        # cclib generation of json files with ccwrite
        termination, errortype, cclib_data, file = self.json_gen(file, file_name)
        outlines = []

        if errortype == "no_data":
            return termination, errortype, None, None, file 

        # calculations with 1 atom
        if cclib_data["properties"]["number of atoms"] == 1:
            cclib_data["vibrations"] = {"frequencies": [], "displacement": []}
            if not "energy" in cclib_data["properties"]:
                termination = "other"
                errortype = "not_specified"
            elif not "free energy" in cclib_data["properties"]["energy"]:
                errortype = "sp_calc"

        # general errors
        elif "vibrations" not in cclib_data:
            termination = "other"
            errortype = "not_specified"
            if "optimization" in cclib_data:
                # if the optimization finished, only a freq job is required
                if (
                    "done" in cclib_data["optimization"]
                    and cclib_data["optimization"]["done"]
                ):
                    errortype = "no_freq"

            # use very short reversed loop to find basis set incompatibilities and SCF errors
            outlines = read_file(os.getcwd(), self.args.w_dir_main, file)
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

        # normal terminations
        if "vibrations" in cclib_data or errortype == "sp_calc":
            # spin contamination analysis using user-defined thresholds
            if "S2 after annihilation" in cclib_data["properties"]:
                unpaired_e = cclib_data["properties"]["multiplicity"] - 1
                # this first part accounts for singlet diradicals (threshold is 10% of the spin before annihilation)
                if unpaired_e == 0:
                    if (
                        float(cclib_data["properties"]["S2 after annihilation"])
                        > abs(float(self.args.s2_threshold) / 100)
                        * cclib_data["properties"]["S2 before annihilation"]
                    ):
                        errortype = "spin_contaminated"
                else:
                    spin = unpaired_e * 0.5
                    s2_expected_value = spin * (spin + 1)
                    spin_diff = abs(
                        float(cclib_data["properties"]["S2 after annihilation"])
                        - s2_expected_value
                    )
                    if (
                        spin_diff
                        > abs(float(self.args.s2_threshold) / 100) * s2_expected_value
                    ):
                        errortype = "spin_contaminated"

        return termination, errortype, cclib_data, outlines, file 

    def analyze_normal(self, duplicate_data, errortype, cclib_data, file_name):
        """
        Analyze errors from normally terminated calculations
        """
        atom_types, cartesians = cclib_atoms_coords(cclib_data)
        dup_off = None
        if errortype == "none":
            # in eV, converted to hartree using the conversion factor from cclib
            E_dup = cclib_data["properties"]["energy"]["total"]
            E_dup = cclib.parser.utils.convertor(E_dup, "eV", "hartree")
            # in hartree
            try:
                H_dup = cclib_data["properties"]["enthalpy"]
                G_dup = cclib_data["properties"]["energy"]["free energy"]
            except (AttributeError, KeyError):
                if cclib_data["properties"]["number of atoms"] == 1:
                    if cclib_data["metadata"]["keywords line"].find("freq") == -1:
                        errortype = "sp_calc"
                        cclib_data["metadata"][
                            "ground or transition state"
                        ] = "SP calculation"
                    H_dup = E_dup
                    G_dup = E_dup
            try:
                ro_dup = cclib_data["properties"]["rotational"]["rotational constants"]
                if len(ro_dup) != 3:
                    ro_dup = None
            except:
                ro_dup = None

            # detects if this calculation is a duplicate
            for i, _ in enumerate(duplicate_data["Energies"]):
                E_diff = abs(E_dup - duplicate_data["Energies"][i])
                H_diff = abs(H_dup - duplicate_data["Enthalpies"][i])
                G_diff = abs(G_dup - duplicate_data["Gibbs"][i])
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

        if errortype == "none":
            duplicate_data["File"].append(file_name)
            duplicate_data["Energies"].append(E_dup)
            duplicate_data["Enthalpies"].append(H_dup)
            duplicate_data["Gibbs"].append(G_dup)
            duplicate_data["RO_constant"].append(ro_dup)

            initial_ifreqs = 0
            for freq in cclib_data["vibrations"]["frequencies"]:
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

            if len(atom_types) in [3, 4]:
                errortype = detect_linear(errortype, atom_types, cclib_data)

            # detects no convergence issues during freq calcs
            if self.args.freq_conv is not None:
                if (
                    errortype == "none"
                    and cclib_data["optimization"]["times converged"] == 1
                ):
                    errortype = "freq_no_conv"

        if errortype in ["extra_imag_freq", "freq_no_conv", "linear_mol_wrong"]:
            if errortype == "extra_imag_freq":
                cartesians = self.fix_imag_freqs(cclib_data, cartesians)

            # in case no previous OPT was done (only works if it's not a TS)
            opt_found = False
            for keyword in cclib_data["metadata"]["keywords line"].split():
                if keyword.lower().startswith("opt"):
                    opt_found = True

            if not opt_found:
                cclib_data["metadata"]["keywords line"] += " opt"

            if errortype == "freq_no_conv":
                # adjust the keywords so only FREQ is calculated
                new_keywords_line = ""
                for keyword in cclib_data["metadata"]["keywords line"].split():
                    if keyword.lower().startswith("opt"):
                        keyword = self.args.freq_conv
                        if (
                            cclib_data["metadata"]["ground or transition state"]
                            == "transition_state"
                        ):
                            keyword = keyword.replace("=(", "=(ts,noeigen,")
                    new_keywords_line += keyword
                    new_keywords_line += " "
                cclib_data["metadata"]["keywords line"] = new_keywords_line

            elif errortype == "linear_mol_wrong":
                cclib_data["metadata"]["keywords line"] += " symmetry=(PG=Cinfv)"

        return atom_types, cartesians, duplicate_data, errortype, cclib_data, dup_off

    def analyze_abnormal(self, errortype, cclib_data, outlines):
        """
        Analyze errors from calculations that did not finish normally
        """
        # for calcs with finished OPT but no freqs, adjust the keywords so only FREQ is calculated
        if errortype == "no_freq":
            new_keywords_line = ""
            for keyword in cclib_data["metadata"]["keywords line"].split():
                if keyword.lower().startswith("opt"):
                    keyword = ""
                else:
                    new_keywords_line += keyword
                    new_keywords_line += " "
            cclib_data["metadata"]["keywords line"] = new_keywords_line
            atom_types, cartesians = cclib_atoms_coords(cclib_data)
        else:
            # help to fix SCF convergence errors
            if errortype == "SCFerror":
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

            if errortype in ["not_specified", "SCFerror"]:
                if "geometric values" in cclib_data["optimization"]:
                    RMS_forces = [
                        row[1] for row in cclib_data["optimization"]["geometric values"]
                    ]
                    # cclib uses None when the values are corrupted in the output files, replace None for a large number
                    RMS_forces = [10000 if val is None else val for val in RMS_forces]
                    min_RMS = RMS_forces.index(min(RMS_forces))
                else:
                    # for optimizations that fail in the first step
                    min_RMS = 0

                atom_types, cartesians = QM_coords(
                    outlines,
                    min_RMS,
                    cclib_data["properties"]["number of atoms"],
                    "gaussian",
                    cclib_data["metadata"]["keywords line"],
                )

        return atom_types, cartesians, cclib_data

    def analyze_isom(self, file, cartesians, atom_types, errortype):
        """
        Check if the initial structure isomerized during QM geometry optimization
        """

        isomerized = False
        isom_valid = True
        init_csv = pd.DataFrame()
        try:
            os.chdir(self.args.isom_inputs)
        except FileNotFoundError:
            self.args.log.write("x  The PATH specified in isom_inputs doesn't exist!")
            isom_valid = False

        if not isom_valid:
            os.chdir(self.args.initial_dir)
            self.args.log.finalize()
            sys.exit()

        try:
            atoms_com, coords_com, atoms_and_coords = [], [], []
            if len(self.args.isom_type.split(".")) == 1:
                atoms_and_coords, _, _ = get_info_input(
                    f'{file.split(".")[0]}.{self.args.isom_type}'
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

            isom_data = {
                "Coords input": coords_com,
                "Coords output": cartesians,
                "Atoms input": atoms_com,
                "Atoms output": atom_types,
                "VdW radii fraction": self.args.vdwfrac,
                "Covalent radii fraction": self.args.covfrac,
                "Initial csv": init_csv,
            }

            isomerized = check_isomerization(isom_data, file)

        except FileNotFoundError:
            self.args.log.write(f"x  No com file were found for {os.path.basename(file)}, the check_geom test will be disabled for this calculation")

        if isomerized:
            errortype = "isomerization"

        os.chdir(self.args.w_dir_main)

        return errortype

    def qcorr_fixing(self, cclib_data, file, atom_types, cartesians):
        """
        Create com files for resubmission with the suggested protocols to correct the errors
        """

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

        if self.args.resume_qcorr:
            destination_fix = Path(
                f"{self.args.w_dir_main}/../../run_{self.args.round_num}/fixed_QM_inputs"
            )
        else:
            destination_fix = Path(
                f"{self.args.w_dir_main}/failed/run_{self.args.round_num}/fixed_QM_inputs"
            )

        if cclib_data["metadata"]["QM program"].lower().find("gaussian") > -1:
            program = "gaussian"
        elif cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
            program = "orca"

        if program in ["gaussian", "orca"]:
            qprep(
                destination=destination_fix,
                w_dir_main=self.args.w_dir_main,
                files=os.path.basename(file),
                charge=cclib_data["properties"]["charge"],
                mult=cclib_data["properties"]["multiplicity"],
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
        else:
            self.args.log.write(f"x  Couldn't create an input file to fix {os.path.basename(file)} (compatible programs: Gaussian and ORCA)\n")

    def json_gen(self, file, file_name):
        """
        Create a json file with cclib and load a dictionary
        """

        termination, errortype = "normal", "none"

        command_run_1 = ["ccwrite", "json", file]
        subprocess.run(command_run_1, capture_output=True)

        cclib_data = {}
        try:
            with open(file_name + ".json") as json_file:
                cclib_data = json.load(json_file)
        except FileNotFoundError:
            try:
                # this part avoids problems when using cclib from command lines (not complete file PATH)
                file = f'{self.args.initial_dir}/{file}'
                command_run_2 = ["ccwrite", "json", file]
                subprocess.run(command_run_2, capture_output=True)
                with open(file_name + ".json") as json_file:
                    cclib_data = json.load(json_file)
            except FileNotFoundError:
                termination = "other"
                errortype = "no_data"

        # add parameters that might be missing from cclib (depends on the version)
        if not hasattr(cclib_data, "metadata") and errortype != "no_data":
            cclib_data = get_json_data(self, file, cclib_data)

        # this is just a "dirty hack" until cclib is updated to be compatible for print mini in ORCA
        if hasattr(cclib_data, "metadata"):
            if cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
                if "final single point energy" in cclib_data["properties"]["energy"]:
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
        for mode, _ in enumerate(cclib_data["vibrations"]["frequencies"]):
            # moves along all imaginary freqs (ignoring the TS imag freq, assumed to be the most negative freq)
            if (
                mode == 0
                and cclib_data["metadata"]["ground or transition state"]
                == "transition_state"
            ):
                shift.append(0.0)
            else:
                if cclib_data["vibrations"]["frequencies"][mode] < 0.0:
                    shift.append(float(self.args.amplitude_ifreq))
                else:
                    shift.append(0.0)

            # The starting geometry is displaced along each normal mode according to the random shift
            for atom in range(0, cclib_data["properties"]["number of atoms"]):
                for coord in range(0, 3):
                    cartesians[atom][coord] = (
                        cartesians[atom][coord]
                        + cclib_data["vibrations"]["displacement"][mode][atom][coord]
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

        if errortype == "none" and termination == "normal":
            destination = destination_normal
            file_terms["finished"] += 1

        elif errortype == "sp_calc" and termination == "normal":
            destination = destination_normal.joinpath("SP_calcs/")
            file_terms["sp_calcs"] += 1

        elif errortype == "extra_imag_freq":
            destination = destination_error.joinpath("extra_imag_freq/")
            file_terms["extra_imag_freq"] += 1

        elif errortype == "ts_no_imag_freq":
            destination = destination_error.joinpath("ts_no_imag_freq/")
            file_terms["ts_no_imag_freq"] += 1

        elif errortype == "spin_contaminated":
            destination = destination_error.joinpath("spin_contaminated/")
            file_terms["spin_contaminated"] += 1

        elif errortype == "duplicate_calc":
            destination = destination_error.joinpath("duplicates/")
            file_terms["duplicate_calc"] += 1

        elif errortype == "atomicbasiserror":
            destination = destination_error.joinpath("error/basis_set_error/")
            file_terms["atom_error"] += 1

        elif errortype == "SCFerror":
            destination = destination_error.joinpath("error/scf_error/")
            file_terms["scf_error"] += 1

        elif errortype == "no_data":
            destination = destination_error.joinpath("error/no_data/")
            file_terms["no_data"] += 1

        elif errortype == "fail_geom_rules":
            destination = destination_error.joinpath("geom_rules_filter/")
            file_terms["geom_rules_qcorr"] += 1

        elif errortype == "isomerization":
            destination = destination_error.joinpath("isomerization/")
            file_terms["isomerized"] += 1

        elif errortype == "freq_no_conv":
            destination = destination_error.joinpath("freq_no_conv/")
            file_terms["freq_no_conv"] += 1

        elif errortype == "linear_mol_wrong":
            destination = destination_error.joinpath("linear_mol_wrong/")
            file_terms["linear_mol_wrong"] += 1

        else:
            destination = destination_error.joinpath("error/not_specified_error/")
            file_terms["not_specified"] += 1

        move_file(destination, self.args.w_dir_main, os.path.basename(file))

        return file_terms, destination

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
        if len(self.args.geom_rules) >= 1:
            ana_data.at[0, "geom_rules filter"] = file_terms["geom_rules_qcorr"]
        if self.args.isom_type is not None:
            ana_data.at[0, "Isomerization"] = file_terms["isomerized"]
        path_as_str = self.args.initial_dir.as_posix()
        csv_qcorr = path_as_str + f"/QCORR-run_{self.args.round_num}-stats.csv"
        ana_data.to_csv(csv_qcorr, index=False)

        return csv_qcorr
