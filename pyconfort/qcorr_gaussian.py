######################################################.
#        This file stores all the functions          #
#          used in the LOG file analyzer             #
######################################################.
import os
import glob
import sys
import subprocess
import numpy as np
import pandas as pd
import json
# the import for check and write genecp is probably inside the genecp function of qprep now, FIX!
from pyconfort.utils import periodic_table
from pyconfort.filter import exp_rules_output
from pyconfort.utils import move_file, bondi, rcov, get_info_com
from pyconfort.qprep_gaussian import write_qm_input_files,GaussianTemplate


def read_log_file(w_dir,file):
    """
    Reads through the QM output files and retrieves a list with all the lines
    """

    break_loop = False
    os.chdir(w_dir)
    try:
        outfile = open(file,"r")
        outlines = outfile.readlines()
    except FileNotFoundError:
        break_loop = True
        outfile, outlines = None, None

    return outlines, outfile, break_loop


def get_qcorr_params(outlines,preset_keywords_line,ts_opt,sp_calcs,charge_sp,mult_sp,isom_filt,no_check):
    """
    Retrieves termination types, error types, original number of imaginary frequencies, keywords
    line (input), number of atoms, charge, multiplicity, culation type (ground_state or 
    transition_state) and line number where NImag= (nimag_line) is located from QM files. 
    
    A combination of for forward and reverse for loops with different starting and ending points
    is used to speed up data reading and avoid reading the same lines multiple lines.

    Parameters
    ----------
    outlines : list of str
        Lines of the QM output files
    preset_keywords_line : str
        If a keyword line is specified in qm_input
    ts_opt : str
        OPT keywords section specified in the ts_input option
    sp_calcs : str
        Program (if any) selected for single-point calculations after processing output files
    charge_sp : int
        Charge for the following single point calculation
    mult_sp :
        Multiplicity for the following single point calculation
    isom_filt : str
        Type of isomerization filter selected in the --isom option (disable with isom = None)
    no_check : bool
        If True, skips analysis for QM output files (useful for single-point energy calculations)

    Returns
    -------
    TERMINATION : string
        Type of termination of the QM output file (i.e. normal, error, unfinished)
    ERRORTYPE : string
        Error type of the QM output file
    initial_IM_FREQS : int
        Initial amount of imaginary frequencies in the QM file adapted for transition state 
        calculations (if the file is a TS, the original number is reduced by 1)
    NATOMS : int
        Number of atoms in the calculation
    CHARGE : int
        Charge of the calculations used in the following input files
    MULT : int
        Multiplicity of the calculations used in the following input files
    keywords_line : string
        Keyword line (input) to use in subsequent input files
    calc_type : str
        Type of the QM calculation (ground_state or transition_state)
    nimag_line : int
        Line number of the line containing the amount of imaginary frequencies in the 
        QM calculation. This is tracked to avoid reading the same lines multiple times
    """    

    TERMINATION,ERRORTYPE  = 'unfinished','unknown'
    stop_term, initial_IM_FREQS = 0, 0
    
    # use reversed loops to find type of termination (faster than forward loops)
    for i in reversed(range(len(outlines)-15,len(outlines))):
        if stop_term == 1:
            break
        # Determine the kind of job termination
        if outlines[i].find("Normal termination") > -1:
            TERMINATION = "normal"
            ERRORTYPE = None
            stop_term += 1
        elif outlines[i].find("Error termination") > -1:
            TERMINATION = "error"
            if outlines[i-1].find("Atomic number out of range") > -1 or outlines[i-1].find("basis sets are only available") > -1 :
                ERRORTYPE = "atomicbasiserror"
            if outlines[i-3].find("SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error") > -1:
                ERRORTYPE = "SCFerror"
            stop_term += 1
    
    # if the calculation ends normally with no imaginary freqs, no more read_lines are necessary (saves time)
    # the nimag_line is defined so that part of the file is not read again
    nimag_line = len(outlines)
    if TERMINATION == 'normal':
        CHARGE, MULT, NATOMS, keywords_line = None, None, 0, ''
        calc_type = None
        for i in reversed(range(60,len(outlines)-5)):
            if outlines[i].find('NImag=') > -1 or outlines[i].find('PG=') > -1:
                nimag_line = i
                # merge a few lines since 'NIMag=' sometimes appears in two different lines
                line_with_ifreq = outlines[i-1].rstrip("\n")+outlines[i].rstrip("\n")+outlines[i+1].rstrip("\n")
                for part_line in line_with_ifreq.split('\\'):
                    if part_line.find('NImag=') > -1:
                        initial_IM_FREQS = int(part_line.split('=')[1])
                        break
                break

    # get keywords line (input), number of atoms, charge and multiplicity
    if TERMINATION != "normal" or initial_IM_FREQS > 0  or sp_calcs != None or isom_filt != None or no_check:
        NATOMS,CHARGE,MULT,keywords_line,initial_IM_FREQS,calc_type = get_init_info(outlines,preset_keywords_line,charge_sp,mult_sp,sp_calcs,initial_IM_FREQS,ts_opt)
        
        # detect calcs that finish before the functional or basis set was included
        if keywords_line == '':
            ERRORTYPE = 'before_E_calculation'

        # helps to fix SCF convergence errors
        if ERRORTYPE == 'SCFerror':
            if keywords_line.find(' scf=qc') > -1:
                pass
            else:
                keywords_line += ' scf=qc'
    
    return TERMINATION,ERRORTYPE,initial_IM_FREQS,NATOMS,CHARGE,MULT,keywords_line,calc_type,nimag_line


def get_init_info(outlines,preset_keywords_line,charge_sp,mult_sp,sp_calcs,initial_IM_FREQS,ts_opt):
    """
    Retrieves keywords line (input), number of atoms, charge and multiplicity from QM files. 
    
    This functions looks for the keywords line used originally in the out QM files, but if a 
    keyword line is specified in qm_input, this preset line is used for the generated input files. The
    preset line is modified if the QM file is a TS calculation. The new line uses the option 
    from ts_input for the OPT keyword section (i.e. ts_input = 'opt=(calcfc,noeigen,ts,maxstep=5)').
    If ts_input = None, the OPT section is not modified.

    Parameters
    ----------
    outlines : list of str
        Lines of the QM output files
    preset_keywords_line : str
        If a keyword line is specified in qm_input
    charge_sp : int
        Charge for the following single point calculation
    mult_sp : int
        Multiplicity for the following single point calculation
    sp_calcs : str
        Program (if any) selected for single-point calculations after processing output files
    initial_IM_FREQS : int
        Initial amount of imaginary frequencies in the QM file
    ts_opt : str
        OPT keywords section specified in the ts_input option

    Returns
    -------
    NATOMS : int
        Number of atoms in the calculation
    CHARGE : int
        Charge of the calculations used in the following input files
    MULT : int
        Multiplicity of the calculations used in the following input files
    keywords_line : string
        Original keyword line (input) from the output QM file
    initial_IM_FREQS : int
        Initial amount of imaginary frequencies in the QM file adapted for transition state 
        calculations (if the file is a TS, the original number is reduced by 1)
    calc_type : str
        Type of the QM calculation (ground_state or transition_state)
    """    

    end_keywords,keywords_line = False,''
    symbol_z_line, gradgrad_line = None, None
    CHARGE, MULT, NATOMS = None, None, 0
    
    # range optimized to skip Gaussian information (short for loop just to get charge, multiplicity and n of atoms)
    for i in range(60,len(outlines)):
        # retrieve the input line
        if outlines[i].startswith(' #'):
            keywords_line += outlines[i].rstrip("\n")
            for j in range(i,i+5):
                if outlines[j].find('----------'):
                    end_keywords = True
                if not end_keywords:
                    keywords_line += outlines[j].rstrip("\n")
            keywords_line = keywords_line[3:]

        # get number of atoms
        if outlines[i].find('Symbolic Z-matrix:') > -1:
            symbol_z_line = i
        elif outlines[i].find('GradGrad') > -1:
            gradgrad_line = i
        if symbol_z_line != None and gradgrad_line != None:
            NATOMS = gradgrad_line-symbol_z_line-4

        # Determine charge and multiplicity
        if outlines[i].find("Charge = ") > -1:
            if charge_sp != None and sp_calcs != None:
                CHARGE = charge_sp
            else:
                CHARGE = int(outlines[i].split()[2])
            if mult_sp != None and sp_calcs != None:
                MULT = mult_sp
            else:
                MULT = int(outlines[i].split()[5].rstrip("\n"))

        if NATOMS != 0:
            break

    # find if the calculation was for a ground or transition state
    calc_type = 'ground_state'
    calcfc_found, ts_found = False, False
    for keyword in keywords_line.split():
        if keyword.lower().find('calcfc') > -1:
            calcfc_found = True
        if keyword.lower().find('ts') > -1:
            ts_found = True

    if calcfc_found and ts_found:
        calc_type = 'transition_state'

    # substract the imaginary frequency from the TS in normal terminations
    if calc_type == 'transition_state':
        initial_IM_FREQS -= 1

    # if a preset input line was defined in args.qm_input, this overwrites the automatic protocol
    if preset_keywords_line != None:
        keywords_line = preset_keywords_line
        # add the ts options if an external input line was defined
        if calc_type == 'transition_state' and ts_opt != None:
            new_preset_keywords_line = '# '
            for keyword in preset_keywords_line.split():
                if keyword.lower().startswith('opt') > -1:
                    if keyword.lower().find('ts') > -1:
                        pass
                    else:
                        keyword = ts_opt
                new_preset_keywords_line += keyword
                new_preset_keywords_line += ' '
            keywords_line = new_preset_keywords_line

    return NATOMS,CHARGE,MULT,keywords_line,initial_IM_FREQS,calc_type


def get_input_geom(file, outlines, keywords_line, NATOMS, MULT, TERMINATION, ERRORTYPE, initial_IM_FREQS, nimag_line, calc_type, duplicate_filter, ifreq_threshold, amplitude_ifreq, s2_threshold, file_list, E_dup, H_dup, G_dup,no_check):
    """
    Analyzes QM output files with normal and not normal terminations and retrieves information
    to create subsequent input files.

    Parameters
    ----------
    Explained individually in its two sub-functions analysis_normal and analysis_not_normal

    Returns
    -------
    file_list : list of str
        Collects file names for the duplicate analysis
    E_dup : list of float
        Collects (electronic energy + ZPE) values for the duplicate analysis
    H_dup : list of float
        Collects enthalpy values for the duplicate analysis
    G_dup : list of float
        Collects Gibbs free energy values for the duplicate analysis 
    ERRORTYPE : string
        Error type of the QM output file
    keywords_line : string
        Line for the input file that will be generated
    ATOMTYPES : list of strings
        List containing the atoms of the system
    CARTESIANS : list of lists
        Cartesian coordinates used for further processing
    """    

    # this parameter sets the line to retrieve the geometry for normally and not normally terminated calcs
    stand_or = 0
        
    if ERRORTYPE == 'before_E_calculation':
        pass

    else:
        if TERMINATION == "normal":
            stand_or,ERRORTYPE,FREQS,NORMALMODE,file_list,E_dup,H_dup,G_dup = analysis_normal(file,outlines,ERRORTYPE,stand_or,nimag_line,NATOMS,no_check,duplicate_filter,file_list,E_dup,H_dup,G_dup,initial_IM_FREQS,ifreq_threshold,s2_threshold,MULT,calc_type)

        else:
            stand_or,keywords_line = analysis_not_normal(outlines,stand_or,no_check,keywords_line)

        # get atom types from the calculation and the cartesian coordinates to use in each case
        ATOMTYPES, CARTESIANS = [],[]
        per_tab = periodic_table()
        for i in range(stand_or+5,stand_or+5+NATOMS):
            massno = int(outlines[i].split()[1])
            if massno < len(per_tab):
                atom_symbol = per_tab[massno]
            else:
                atom_symbol = "XX"
            ATOMTYPES.append(atom_symbol)
            CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

        if ERRORTYPE == 'extra_imag_freq':
            CARTESIANS = fix_imag_freqs(NATOMS, CARTESIANS, FREQS, NORMALMODE, amplitude_ifreq)
        
        return file_list, E_dup, H_dup, G_dup, ERRORTYPE, keywords_line, ATOMTYPES, CARTESIANS


def analysis_normal(file,outlines,stand_or,ERRORTYPE,nimag_line,NATOMS,no_check,duplicate_filter,file_list,E_dup,H_dup,G_dup,initial_IM_FREQS,ifreq_threshold,s2_threshold,MULT,calc_type):
    """
    Analyzes QM output files with normal terminations. The coordinates of the last geometry
    of the file is targeted. This function also detects duplicates, undesired imaginary frequencies
    and systems with spin contamination. 
    
    For duplicates (--dup option), this function adds different thermochemistry values that 
    will be analyzed further ahead.
    
    For imaginary frequencies, this function will use the user defined threshold to collect 
    imaginary frequencies (--ifreq_cutoff). If the absolute value of a frequency is bigger than
    the absolute value of the threshold, the frequency is considered as negative (i.e. if 
    ifreq_cutoff = 50, all the frequencies < -50 cm-1 would be considered as negative frequencies). 
    
    For spin contamination, if the <S**2> operator deviates more than the user defined threshold
    (--s2_threshold) compared to the expected value, the calculation contains spin contamination. 
    For example, in a triplet with a expected value of <S**2> = 2, the acceptable range will be
    1.80-2.20 when s2_threshold = 10%.
    
    If the --nocheck option is set to True, all the calculations will be considered to have
    normal terminations with no error types and no imaginary frequencies.

    Parameters
    ----------
    file : string
        Filename
    outlines : list of str
        Lines of the QM output files
    stand_or : int
        Starts at 0 before the analysis
    ERRORTYPE : string
        Error type of the QM output file
    nimag_line : int
        Line number of the line containing the amount of imaginary frequencies in the 
        QM calculation. This is tracked to avoid reading the same lines multiple times
    NATOMS : int
        Number of atoms in the calculation
    no_check : bool
        If this is set to True, the analysis is skipped and the final geometry of
        the file is used as the target geometry
    duplicate_filter : bool
        If set to True, thermochemistry data is collected for a subsequent duplicate analysis
    file_list : list of str
        Collects file names for the duplicate analysis
    E_dup : list of float
        Collects (electronic energy + ZPE) values for the duplicate analysis
    H_dup : list of float
        Collects enthalpy values for the duplicate analysis
    G_dup : list of float
        Collects Gibbs free energy values for the duplicate analysis 
    initial_IM_FREQS : int
        Initial amount of imaginary frequencies in the QM file
    ifreq_threshold : float
        Threshold to consider imaginary frequencies
    s2_threshold : float
        Threshold to consider spin contamination
    MULT : int
        Multiplicity of the QM calculation
    calc_type : str
        Type of the QM calculation (ground_state or transition_state)

    Returns
    -------
    stand_or : int
        Line number containing the standard orientation of the target geometry
    ERRORTYPE : string
        Error type of the QM output file
    FREQS : list of float
        List containing the frequencies as floats
    NORMALMODE : list of matrixes
        Contains the normal modes for each frequencies (including all atoms)
    file_list : list of str
        Collects file names for the duplicate analysis
    E_dup : list of float
        Collects (electronic energy + ZPE) values for the duplicate analysis
    H_dup : list of float
        Collects enthalpy values for the duplicate analysis
    G_dup : list of float
        Collects Gibbs free energy values for the duplicate analysis 
    """

    IM_FREQS, start_dup = 0, False
    spin_contamination = False
    FREQS, NORMALMODE  = [],[]

    # reverse loop to speed up the reading of the output files
    # (also skips final part for normal termination which was already read plus dipole/multipole information)
    for i in reversed(range(0,nimag_line-(10*NATOMS))):
        # Sets where the final coordinates are inside the file
        if (outlines[i].find("Standard orientation:") > -1 or outlines[i].find("Input orientation:") > -1):
            stand_or = i
            break

        if not no_check:
            # analyze duplicates
            if duplicate_filter:
                if not start_dup:
                    if outlines[i].find('E (Thermal)'):
                        start_dup = True    
                else:
                    if outlines[i].find('Sum of electronic and zero-point Energies=') > -1:
                        E_dup.append(float(outlines[i].split()[-1]))
                        file_list.append(file)
                        
                    elif outlines[i].find('Sum of electronic and thermal Enthalpies=') > -1:
                        H_dup.append(float(outlines[i].split()[-1]))
                    elif outlines[i].find('Sum of electronic and thermal Free Energies=') > -1:
                        G_dup.append(float(outlines[i].split()[-1]))

            if initial_IM_FREQS != 0:
                # Get the negative frequencies and their modes (up to three negative frequencies)
                if outlines[i].find(" Frequencies -- ") > -1:
                    nfreqs = len(outlines[i].split())
                    for j in range(2, nfreqs):
                        if float(outlines[i].split()[j]) < 0 and abs(float(outlines[i].split()[j])) > abs(ifreq_threshold):
                            IM_FREQS += 1
                            FREQS.append(float(outlines[i].split()[j]))
                            NORMALMODE.append([])
                            for k in range(0,NATOMS):
                                NORMALMODE[(j-2)].append([float(outlines[i+5+k].split()[3*(j-2)+2]), float(outlines[i+5+k].split()[3*(j-2)+3]), float(outlines[i+5+k].split()[3*(j-2)+4])])

            # analyze spin contamination
            elif s2_threshold > 0.0:
                if outlines[i].find('S**2 before annihilation') > -1:
                    s2_value = outlines[i].split()[-1].rstrip("\n")
                    unpaired_e = MULT-1
                    spin = unpaired_e*0.5
                    s2_expected_value = spin*(spin+1)
                    spin_diff = abs(s2_value-s2_expected_value)
                    if spin_diff > abs(s2_threshold/100)*s2_expected_value:
                        spin_contamination = True

            # exclude TS imag frequency if needed
            if calc_type == 'transition_state':
                IM_FREQS -= 1
                FREQS.pop(0)
                NORMALMODE.pop(0)

            if IM_FREQS > 0:
                ERRORTYPE = 'extra_imag_freq'

            elif IM_FREQS < 0:
                ERRORTYPE = 'ts_no_imag_freq'

            if spin_contamination:
                ERRORTYPE = 'spin_contaminated'
            
    return stand_or,ERRORTYPE,FREQS,NORMALMODE,file_list,E_dup,H_dup,G_dup


def analysis_not_normal(outlines,stand_or,no_check,keywords_line):
    """
    Analyzes QM output files with error terminations or unfinished. The geometry
    with the lowest gradient during optimization is targeted since this geometry
    is the closest point to convergence during the optimization process. If the 
    --nocheck option is set to True, the last geometry of the file will be used instead.

    Parameters
    ----------
    outlines : list of str
        Lines of the QM output files
    stand_or : int
        Starts at 0 before the analysis
    no_check : bool
        If this is set to True, the analysis is skipped and the final geometry of
        the file is used as the target geometry
    keywords_line : string
        Line for the input file that will be generated. It might change 
        during the analysis if the optimization finished normally but the following
        frequency calculation failed (the opt keyword is removed in this case)

    Returns
    -------
    stand_or : int
        Line number containing the standard orientation of the target geometry
    keywords_line : string
        Line for the input file that will be generated
    """

    rms,stop_rms = 10000,0
    freq_only =  False
    normal_term_found, opt_found = False, False
    if not no_check:
        for i in reversed(range(60,len(outlines))):
            if outlines[i].find('Cartesian Forces:  Max') > -1:
                try:
                    if float(outlines[i].split()[5]) < rms:
                        rms = float(outlines[i].split()[5])
                        stop_rms = i
                    if normal_term_found and opt_found:
                        freq_only = True
                        break
                # sometimes Gaussian does not print valid numbers in this section
                except ValueError:
                    pass

            if not normal_term_found and not opt_found:
                if outlines[i].find("Normal termination") > -1:
                    normal_term_found = True
                elif outlines[i].strip().find('\\FOpt\\') > -1 or outlines[i].strip().find('\\FTS\\') > -1:
                    opt_found = True

        # only include Freq calculations to molecules that are already optimized (no OPT)
        if freq_only:
            new_keywords_line = '# '
            for keyword in keywords_line.split():
                if keyword.lower().startswith('opt') > -1:
                    keyword = ''
                new_keywords_line += keyword
                new_keywords_line += ' '
            keywords_line = new_keywords_line

    if stop_rms == 0:
        last_line = len(outlines)
    else:
        last_line = stop_rms
    
    for i in reversed(range(0,last_line)):
        # Sets where the final coordinates are inside the file
        if outlines[i].find("Standard orientation") > -1:
            stand_or = i
            break
    
    return stand_or,keywords_line


def fix_imag_freqs(NATOMS, CARTESIANS, FREQS, NORMALMODE, amplitude):
    """
    Fixes undersired (extra) imaginary frequencies from QM calculations.
    This function multiplies the imaginary normal mode vectors by the selected amplitude 
    (0.2 is the default amplitude in the pyQRC script from GitHub, user: bobbypaton).
    By default, all the extra imaginary modes are used (i.e. in calculations with three
    extra imaginary frequencies, all the three modes will be used to displace the atoms).
    This can be tuned with the --ifreq_cutoff option (i.e. only use freqs lower than -50 cm-1).

    Parameters
    ----------
    NATOMS : int
        Number of atoms in the calculation
    CARTESIANS : list of lists
        List of lists containing the molecular coordinates as floats
    FREQS : list of float
        List containing the frequencies as floats
    NORMALMODE : list of matrixes
        Contains the normal modes for each frequencies (including all atoms)
    amplitude : float
        Magnitude to displace the normal vector (commonly from -1 to 1)

    Returns
    -------
    CARTESIANS : list of lists
        New set of cartesian coordinates generated after displacing the original
        coordinates along the normal modes of the corresponding imaginary frequencies
    """

    # Multiplies the imaginary normal mode vector by this amount (from -1 to 1).
    # 0.2 is the default amplitude in the pyQRC script (GitHub, user: bobbypaton)
    shift = []

    # Save the original Cartesian coordinates before they are altered
    orig_carts = []
    for atom in range(0,NATOMS):
        orig_carts.append([CARTESIANS[atom][0], CARTESIANS[atom][1], CARTESIANS[atom][2]])

    # could get rid of atomic units here, if zpe_rat definition is changed
    for mode,_ in enumerate(FREQS):
        # Either moves along all imaginary freqs, or a specific mode requested by the user
        if FREQS[mode] < 0.0:
            shift.append(amplitude)
        else:
            shift.append(0.0)

        # The starting geometry is displaced along each normal mode according to the random shift
        for atom in range(0,NATOMS):
            for coord in range(0,3):
                CARTESIANS[atom][coord] = CARTESIANS[atom][coord] + NORMALMODE[mode][atom][coord] * shift[mode]

    return CARTESIANS

def check_for_final_folder(w_dir):
    """
    # Determines the folder where input files are gonna be generated in QCORR.
    """
    input_folder = w_dir+'\input_files'
    folder_count = 0
    
    if os.path.exists(input_folder):
        dir_list = os.listdir(input_folder)
        for folder in dir_list:
            if folder.find('run_') > -1:
                folder_count += 1

    if folder_count == 0:
        return 0
    else:
        num_com_folder = sum([len(d) for r, d, folder in os.walk(w_dir+'/input_files')])
        w_dir = w_dir+'/input_files/run_'+str(num_com_folder)
        return folder_count

def create_folder_and_com(com_type,w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,file,keywords_line,CHARGE,MULT):

    # creating new folder with new input gaussian files
    if com_type == 'analysis':
        new_gaussian_input_files = w_dir_main+'/input_files/run_'+str(round_num+1)

    try:
        os.makedirs(new_gaussian_input_files)
    except OSError:
        if  os.path.isdir(new_gaussian_input_files):
            os.chdir(new_gaussian_input_files)
        else:
            raise
        
    log.write('-> Creating new gaussian input file for {0} in {new_gaussian_input_files}')

    #error if both genecp and gen are specified
    if args.genecp_atoms != [] and args.gen_atoms != []:
        sys.exit("x  ERROR: Can't use Gen and GenECP at the same time")


    if com_type == 'sp':
        if args.suffix_sp == 'None':
            file_name = file.split(".")[0]
        else:
            file_name = file.split(".")[0]+'_'+args.suffix_sp
        type_com = 'sp'

    elif com_type == 'analysis':
        file_name = file.split(".")[0]
        type_com = 'analysis'
    
    write_qm_input_files(args, new_gaussian_input_files, file_name,
                        None, CHARGE, MULT, type_com, NATOMS, ATOMTYPES, CARTESIANS, keywords_line)


def organize_outputs(w_dir_main,round_num,file,TERMINATION,ERRORTYPE,file_terms):
    """
    1. Moves the QM output files to their corresponding folders after the analysis. 
    2. Keeps track of the number of calculations with the different types 
    of terminations and error types

    Parameters
    ----------
    w_dir_main : string
        Working directory
    round_num : string
        Round of analysis and input file generation in QCORR
    file : string
        Filename
    TERMINATION : string
        Type of termination of the QM output file (i.e. normal, error, unfinished)
    ERRORTYPE : string
        Type of error type of the QM output file (i.e. None, unknown, extra_imag_freq, etc)
    format : string
        File format
    file_terms : dictionary
        Keeps track of the number of calculations for each termination and error type
    """

    if ERRORTYPE == None and TERMINATION == "normal":
        destination = w_dir_main+'/success/output_files'
        file_terms['finished'] += 1

    elif ERRORTYPE == 'extra_imag_freq':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/imag_freq/'
        file_terms['imag_freq'] += 1

    elif ERRORTYPE == 'ts_no_imag_freq':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/ts_no_imag_freq/'
        file_terms['ts_no_imag_freq'] += 1

    elif ERRORTYPE == 'spin_contaminated':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/spin_contaminated/'
        file_terms['spin_contaminated'] += 1

    elif ERRORTYPE == 'duplicate_calc':
        destination = w_dir_main+'/duplicates/run_'+str(round_num)
        file_terms['duplicate_calc'] += 1

    elif TERMINATION == "error":
        if ERRORTYPE == "atomicbasiserror":
            destination = w_dir_main +'/failed/run_'+str(round_num)+'/error/basis_set_error'
            file_terms['atom_error'] += 1
        elif ERRORTYPE == "SCFerror":
            destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/scf_error'
            file_terms['scf_error'] += 1
        elif ERRORTYPE == "before_E_calculation":
            destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/before_E_calculation'
            file_terms['before_E_error'] += 1
        else:
            destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/unknown_error'
            file_terms['other_error'] += 1

    elif TERMINATION == "unfinished":
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/unfinished/'
        file_terms['unfinished'] += 1

    elif ERRORTYPE == 'fail_exp_rules':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/exp_rules_filter/'
        file_terms['exp_rules_qcorr'] += 1

    elif ERRORTYPE == 'isomerization':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/isomerization/'
        file_terms['check_geom_qcorr'] += 1
    
    if not os.path.isdir(destination):
        os.makedirs(destination)
    move_file(file, w_dir_main, destination)


def output_to_mol(file,format,log):
    #!!!DONT WE HAVE A VERY SIMILAR FUNCTION IN QPREP?
    """
    Input an XYZ, LOG or OUT file from QM calculations and converts it into
    a mol object.

    Parameters
    ----------
    file : string
        Filename
    format : string
        File format
    log : Logger object
        Writes data regarding the results from QCORR

    Returns
    -------
    mol
        Mol object
    ob_compat
        True if openbabel is installed correctly
    rdkit_compat
        True if RDKit is installed correctly
    """

    ob_compat = True
    rdkit_compat = True
    try:
        import openbabel as ob
    except (ModuleNotFoundError,AttributeError):
        log.write('\nx  Open Babel is not installed correctly, the exp_rules filter will be disabled')
        ob_compat = False
    try:
        from rdkit.Chem import AllChem as Chem
    except (ModuleNotFoundError,AttributeError):
        log.write('\nx  RDKit is not installed correctly, the exp_rules and check_geom filters will be disabled')
        rdkit_compat = False

    # transforms output file into mol object
    # for input (from com to xyz to mol)
    if format == 'xyz':
        cmd_obabel = ['obabel', '-ixyz', os.path.splitext(file)[0]+'.xyz', '-omol', '-O', os.path.splitext(file)[0]+'.mol']
    # for output (from log to mol)
    if format in ['log','out']:
        cmd_obabel = ['obabel', '-ilog', os.path.splitext(file)[0]+'.'+format, '-omol', '-O', os.path.splitext(file)[0]+'.mol']
    subprocess.run(cmd_obabel)
    mol = Chem.MolFromMolFile(file.split('.')[0]+'.mol')

    return mol,ob_compat,rdkit_compat


def output_processing(log_files, w_dir_main, args, log, round_num):
    """
    General function of the QCORR module that:
    1. Analyzes the QM output files and moves output files with normal termination 
    and no extra imaginary frequencies to the same folder
    2. Generates input files to fix errors and extra imaginary frequencies
    3. Generates input files with new keywords lines from the normally terminated 
    files from point 1 (i.e. single-point energy corrections). Optionally, 
    the analysis from points 1 and 2  might be disabled with the nocheck=True or --nocheck option.

    Parameters
    ----------
    log_files : list 
        Contains the filenames of QM output files to analyze
    w_dir_main : string
        Working directory
    args : class
        Class containing the different arguments
    log : Logger object
        Writes data regarding the results from QCORR
    round_num : string
        Round of analysis and input file generation in QCORR
    """

    # these lists will track duplicates if args.dup is activated
    file_list, E_dup, H_dup, G_dup = [],[],[],[]
    file_terms = {'finished': 0, 'imag_freq': 0, 'ts_no_imag_freq': 0, 
                'spin_contaminated': 0, 'duplicate_calc': 0, 'atom_error': 0,
                'scf_error': 0, 'before_E_error': 0, 'other_error': 0,
                'unfinished': 0, 'exp_rules_qcorr': 0, 'check_geom_qcorr': 0}

    for file in log_files:

        # read the file
        log.write(file)
        outlines, outfile, break_loop = read_log_file(w_dir_main,file)

        if break_loop:
            break

        # get parameters from the calculations (i.e. termination and error types, n of atoms, charge, mult, etc)
        TERMINATION, ERRORTYPE, initial_IM_FREQS, NATOMS, CHARGE, MULT, keywords_line, calc_type, nimag_line = get_qcorr_params(outlines,args.qm_input,args.ts_input,args.sp,args.charge_sp,args.mult_sp,args.isom,args.nocheck)
        
        # get new coordinates as input files to fix error (unless the errors cannot be fixed automatically
        # such as in atomic basis errors or errors happening too early on the calculation)
        if ERRORTYPE not in ['before_E_calculation',"atomicbasiserror"] or args.nocheck:

            # get geometry parameters and frequency information
            if TERMINATION != "normal" or initial_IM_FREQS != 0 or args.s2_threshold > 0.0 or args.sp != None or args.dup or args.isom != None:
                if args.nocheck:
                    TERMINATION = "normal"
                file_list, E_dup, H_dup, G_dup, ERRORTYPE, keywords_line, ATOMTYPES, CARTESIANS = get_input_geom(file, outlines, keywords_line, NATOMS, MULT, TERMINATION, ERRORTYPE, initial_IM_FREQS, nimag_line, calc_type, args.dup, args.ifreq_cutoff, args.amplitude_ifreq, args.s2_threshold, file_list, E_dup, H_dup, G_dup,args.nocheck)
            
            if args.isom != None:
                isomerized = False
                init_csv = pd.DataFrame()
                log.write("  ----- Geometrical check will be applied to the output file -----\n")
                try:
                    atoms_com, coords_com, atoms_and_coords = [],[],[]
                    if round_num != 0:
                        os.chdir(w_dir_main +'/input_files/run_'+str(round_num))
                    else:
                        pass
                    if args.isom == 'com':
                        atoms_and_coords,_ = get_info_com(file.split('.')[0]+'.com')
                    elif args.isom == 'gjf':
                        atoms_and_coords,_ = get_info_com(file.split('.')[0]+'.gjf')
                    elif args.isom.split('.')[1] == 'csv':
                        init_csv = pd.read_csv(args.isom)

                    for line in atoms_and_coords:
                        atoms_com.append(line.split()[0])
                        coords_com.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
                    
                    isomerized = check_isomerization(coords_com, CARTESIANS, atoms_com, ATOMTYPES, args.vdwfrac, args.covfrac, init_csv, file)
                
                except FileNotFoundError:
                    log.write("x  No com file were found for "+file+", the check_geom test will be disabled for this calculation")
                
                if isomerized:
                    ERRORTYPE = 'isomerization'

                os.chdir(w_dir_main)
        
            if len(args.exp_rules) >= 1:
                passing_rules = True
                valid_mol_gen = True
                log.write("  ----- Exp_rules filter(s) will be applied to the output file -----\n")
                try:
                    mol,ob_compat,rdkit_compat = output_to_mol(file,format_file,log)
                    print_error_exp_rules=False
                    if ob_compat and rdkit_compat:
                        passing_rules = exp_rules_output(mol,args,log,file,print_error_exp_rules)
                        if not passing_rules:
                            ERRORTYPE = 'fail_exp_rules'
                    os.remove(file.split('.')[0]+'.mol')
                except AttributeError:
                    valid_mol_gen = False
                    os.remove(file.split('.')[0]+'.mol')
                    log.write("The file could not be converted into a mol object, exp_rules filter(s) will be disabled\n")

        # is it necessary to close the file?!
        # close the file
        outfile.close()

        format_file = file.split('.')[1]
        
        if TERMINATION == "normal" and ERRORTYPE == None:
            # detects if this calculation is a duplicate
            if args.dup:
                for i,energy_value in enumerate(E_dup):
                    if i != len(E_dup)-1:
                        if abs(energy_value - E_dup[file_list.index(file)]) < abs(args.dup_threshold):
                            if abs(H_dup[i] - H_dup[file_list.index(file)]) < abs(args.dup_threshold):
                                if abs(G_dup[i] - G_dup[file_list.index(file)]) < abs(args.dup_threshold):                        
                                    ERRORTYPE = 'duplicate_calc'
        
        # This part places the calculations in different folders depending on the type of termination
        organize_outputs(w_dir_main,round_num,file,TERMINATION,ERRORTYPE,file_terms)

        # check if gen or genecp are active
        # right now, QCORR is only working with Gaussian output files
        ecp_list,orca_aux_section = [],[]
        if args.genecp_atoms != [] or args.gen_atoms != []:
            ecp_list = check_for_gen_or_genecp(ATOMTYPES,args,'analysis','gaussian')
        
        # create folders and set level of theory in COM files to fix imaginary freqs or not normal terminations
        if TERMINATION != "normal" or ERRORTYPE not in [None,'isomerization']:
            create_folder_and_com('analysis',w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,file,keywords_line,CHARGE, MULT)

        # this part creates input files for single-point energy calcs after reading from normally finished log files
        if TERMINATION == "normal" and ERRORTYPE == None:
            if args.sp != None or args.nics:

                if args.sp == 'gaussian':
                    # creating new folder with new input Gaussian files
                    single_point_input_files = w_dir_main+'/success/G16-SP_input_files'

                elif args.sp == 'orca':
                    # creating new folder with new input ORCA files
                    single_point_input_files = w_dir_main+'/success/ORCA-SP_input_files'

                elif args.sp == 'turbomole':
                  # creating new folder with new input Turbomole files
                  single_point_input_files = w_dir_main+'/success/output_files/TURBOMOLE-SP_input_files'

                if args.nics:
                    nics_input_files = w_dir_main+'/success/output_files/G16-NICS_input_files'

                # Options for genecp
                if args.sp == 'gaussian':
                    ecp_list = check_for_gen_or_genecp(ATOMTYPES,args,'sp','gaussian')

                elif args.sp == 'orca':
                    orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','orca')

                elif args.sp == 'turbomole':
                    # Check for gen or genecp
                    pass

                # Sets the folder and find the log files to analyze
                for lot_sp,bs_sp,bs_gcp_sp in zip(args.level_of_theory_sp,args.basis_set_sp,basis_set_for_genecp):
                    if args.sp == 'gaussian' or args.sp == 'orca' or args.sp == 'turbomole':
                        if str(bs_sp).find('/') > -1:
                            log.write('-> Creating new single point files for {0} in {1}/{2}-{3}'.format(file,single_point_input_files,lot_sp,bs_sp.split('/')[0]))
                        else:
                            log.write('-> Creating new single point files for {0} in {1}/{2}-{3}'.format(file,single_point_input_files,lot_sp,bs_sp))
                    if args.nics:
                        log.write('-> Creating NICS input files for {0} in {1}/{2}-{3}'.format(file,nics_input_files,lot_sp,bs_sp))

                    # eliminates * from the name (since folders cannot contain * in their names)
                    if bs_sp.find('**') > -1:
                        bs_sp = bs_sp.replace('**','(d,p)')
                    elif bs_sp.find('*') > -1:
                        bs_sp = bs_sp.replace('*','(d)')

                    if str(bs_sp).find('/') > -1:
                        dir_name = str(lot_sp) + '-' + str(bs_sp.split('/')[0])
                    else:
                        dir_name = str(lot_sp) + '-' + str(bs_sp)

                    keywords_opt = ''
                    if args.sp == 'gaussian':
                        if genecp == 'genecp' or  genecp == 'gen':
                            keywords_opt = lot_sp + '/' + genecp
                        else:
                            keywords_opt = lot_sp + '/' + bs_sp
                        if args.qm_input_sp != 'None':
                            keywords_opt += ' {0}'.format(args.qm_input_sp)
                        if args.empirical_dispersion_sp != 'None':
                            keywords_opt += ' empiricaldispersion={0}'.format(args.empirical_dispersion_sp)
                        if args.solvent_model_sp != 'gas_phase':
                            keywords_opt += ' scrf=({0},solvent={1})'.format(args.solvent_model_sp,args.solvent_name_sp)

                    if args.charge_sp != 'None':
                        CHARGE = args.charge_sp

                    if args.mult_sp != 'None':
                        MULT = args.mult_sp

                    if args.sp != None:
                        dir_for_sp = single_point_input_files+'/'+dir_name
                        if not os.path.isdir(dir_for_sp):
                            os.makedirs(dir_for_sp)
                        ('analysis',w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,file,keywords_line,CHARGE, MULT)
                        create_folder_and_com('sp',dir_for_sp,file,args,keywords_opt,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS)

                    if args.nics:
                        dir_for_nics = nics_input_files+'/'+dir_name
                        if not os.path.isdir(dir_for_nics):
                            os.makedirs(dir_for_nics)
                        create_folder_and_com('nics',dir_for_nics,file,args,keywords_opt,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS)

    os.chdir(w_dir_main)

    # moves the input files to respective folder
    destination = w_dir_main +'/input_files/run_'+str(round_num)
    for comfile in glob.glob('*.com'):
        move_file(comfile, w_dir_main, destination)
    for gjffile in glob.glob('*.gjf'):
        move_file(gjffile, w_dir_main, destination)

    #write to csv ana_data
    ana_data = pd.DataFrame()
    ana_data.at[0,'Total files'] = len(log_files)
    ana_data.at[0,'Normal termination'] = file_terms['finished']
    ana_data.at[0,'Imaginary frequencies'] = file_terms['imag_freq']
    ana_data.at[0,'TS with no imag. freq.'] = file_terms['ts_no_imag_freq']
    ana_data.at[0,'SCF error'] = file_terms['scf_error']
    ana_data.at[0,'Error before SCF'] = file_terms['before_E_error']
    ana_data.at[0,'Basis set error'] =  file_terms['atom_error']
    ana_data.at[0,'Other errors'] = file_terms['other_error']
    ana_data.at[0,'Unfinished'] = file_terms['unfinished']
    if args.s2_threshold > 0.0:
        ana_data.at[0,'Spin contamination'] = file_terms['spin_contaminated']
    if args.dup:
        ana_data.at[0,'Duplicates'] = file_terms['duplicate_calc']
    if len(args.exp_rules) >= 1:
        ana_data.at[0,'Exp_rules filter'] = file_terms['exp_rules_qcorr']
    if args.isom != None:
        ana_data.at[0,'Isomerization'] = file_terms['check_geom_qcorr']

    if not os.path.isdir(w_dir_main+'/csv_files/'):
        os.makedirs(w_dir_main+'/csv_files/')
    ana_data.to_csv(w_dir_main+'/csv_files/Analysis-Data-QCORR-run_'+str(round_num)+'.csv',index=False)


def check_isomerization(COORDINATES_com, COORDINATES_log, ATOMTYPES_com, ATOMTYPES_log, vdwfrac, covfrac, init_csv, file):
    """
    Inputs two molecules with the atoms in the same order and checks if any bond
    is too different between them. 
    
    Bonds are considered when the distance between two atoms is smaller than
    either the sum of their adjusted VDW radii or covalent radii (dist = n*R1 + n*R2, where
    n is a user defined parameter). 

    Bonds forming part of TSs are removed.

    Parameters
    ----------
    COORDINATES_com : list of lists containing atomic coordinates
        Molecule 1 (Theoretically, non-optimized)
    COORDINATES_log : list of lists containing atomic coordinates
        Molecule 2 (Theoretically, optimized)
    ATOMTYPES_com : list of atoms
        Molecule 1 (Theoretically, non-optimized)
    ATOMTYPES_log : list of atoms
        Molecule 2 (Theoretically, optimized)
    vdwfrac : float
        Fraction of the summed VDW radii (default is 0.5)
    covfrac : float
        Fraction of the summed covalent radii (default is 1.10)
    init_csv : dataframe
        Contains connectivity from the original non-optimized molecules (i.e. saved from CSEARCH)
    file : string
        Filename

    Returns
    -------
    bool
        True there is a clearly distorted bond within the geometries.
    """
    
    isomerized, diff = None, None

    # load connectivity matrix from the starting points and convert string into matrix
    if not init_csv.empty:
        filename = file.replace('_'+file.split('_')[-1],'')
        init_connectivity_string = init_csv[init_csv['code_name'] == filename]['initial_connectiv'][0]
        init_connectivity = json.loads(init_connectivity_string.replace('.',',').replace(',]','],').replace('],]',']]')) 
        ATOMTYPES_com = init_connectivity[0]

    else:
        init_connectivity = gen_connectivity(ATOMTYPES_com, COORDINATES_com, vdwfrac, covfrac)

    # in case the systems are not the same
    if len(ATOMTYPES_log) != len(ATOMTYPES_com):
        isomerized = True
    
    else:
        final_connectivity = gen_connectivity(ATOMTYPES_log, COORDINATES_log, vdwfrac, covfrac)
                
        # check connectivity differences from initial structure
        diff = final_connectivity - init_connectivity

        # remove bonds involved in TSs from connectivity matrixes
        if not init_csv.empty:
            if 'TS_atom_idx' in init_csv.columns:
                ts_atoms = init_csv[init_csv['code_name'] == filename]['TS_atom_idx'][0].split(',')
                for i,ts_idx in enumerate(ts_atoms):
                    for j,ts_idx_2 in enumerate(ts_atoms):
                        if j>i:
                            diff[int(ts_idx)][int(ts_idx_2)] = 0
                            diff[int(ts_idx_2)][int(ts_idx)] = 0

        isomerized = np.any(diff)

    return isomerized

def gen_connectivity(ATOMTYPES_conn, COORDINATES_conn, vdwfrac, covfrac):
    """
    Use VDW radii to infer a connectivity matrix
    """

    conn_mat = np.zeros((len(ATOMTYPES_conn), len(ATOMTYPES_conn)))
    for i, elem_i in enumerate(ATOMTYPES_conn):
        for j, elem_j in enumerate(ATOMTYPES_conn):
            if j > i:
                vdw_ij = bondi()[elem_i] + bondi()[elem_j]
                rcov_ij = rcov()[elem_i] + rcov()[elem_j]
                dist_ij = np.linalg.norm(np.array(COORDINATES_conn[i])-np.array(COORDINATES_conn[j]))
                if dist_ij / vdw_ij < vdwfrac or dist_ij / rcov_ij < covfrac:
                    conn_mat[i][j] = 1
                else: pass

    return conn_mat
