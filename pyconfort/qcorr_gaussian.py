######################################################.
#        This file stores all the functions          #
#          used in the LOG file analyzer             #
######################################################.
import os
import sys
import subprocess
from pyconfort.qcorr_rework import input_route_line,write_genecp,orca_file_gen
# the import for check and write genecp is probably inside the genecp function of qprep now, FIX!
from pyconfort.utils import periodic_table
from pyconfort.filter import exp_rules_output,check_geom_filter
from pyconfort.nics_conf import update_coord
from pyconfort.utils import move_file, com_2_xyz_2_sdf
from pyconfort.qprep_gaussian import write_gaussian_input_file

def write_header_and_coords(fileout,args,keywords_opt,file_name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,w_dir_initial,log,com_type=None):
    if com_type == 'nics':
        NATOMS,ATOMTYPES,CARTESIANS = update_coord(NATOMS,ATOMTYPES,CARTESIANS,args,log,file_name,w_dir_initial,'write')
    fileout.write("%mem="+str(args.mem)+"\n")
    fileout.write("%nprocshared="+str(args.nprocs)+"\n")
    fileout.write("# "+keywords_opt+"\n")
    fileout.write("\n")
    fileout.write(file_name.split('.')[0]+"\n")
    fileout.write(str(CHARGE)+' '+str(MULT)+'\n')
    for atom in range(0,NATOMS):
        fileout.write('{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}'.format(ATOMTYPES[atom], CARTESIANS[atom][0],  CARTESIANS[atom][1],  CARTESIANS[atom][2]))
        fileout.write("\n")
    fileout.write("\n")

# CREATION OF COM FILES
def new_com_file(com_type,w_dir_initial,log,new_gaussian_input_files,file,args,keywords_opt,file_name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,bs_com,lot_com,bs_gcp_com,orca_aux_section):

    if com_type == 'sp':
        if args.suffix_sp == 'None':
            file_name = file.split(".")[0]+'.com'

        else:
            file_name = file.split(".")[0]+'_'+args.suffix_sp+'.com'

    elif com_type == 'analysis':
        file_name = file.split(".")[0]+'.com'

    elif com_type == 'nics':
        file_name = file.split(".")[0]+'_nics.com'

    fileout = open(file_name, "w")

    write_header_and_coords(fileout,args,keywords_opt,file_name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,w_dir_initial,log,com_type)

    # write genecp/gen part
    if bs_com.lower() == 'genecp' or  bs_com.lower() == 'gen':
        if com_type == 'sp':
            type_gen = 'sp'
        elif com_type == 'analysis':
            type_gen = 'qcorr'

        write_genecp(type_gen,fileout,ecp_list,ecp_genecp_atoms,bs_com,bs_gcp_com,args,w_dir_initial,new_gaussian_input_files)

    if args.sp == 'gaussian' and com_type == 'sp':
        # final line for SP
        if args.qm_input_end_sp != 'None':
            fileout.write(args.qm_input_end_sp)
            fileout.write('\n\n')

    if args.QCORR == 'gaussian' and com_type == 'analysis':
        # final line for analysis
        if args.qm_input_end != 'None':
            fileout.write(args.qm_input_end)
            fileout.write('\n\n')

    fileout.close()

    if args.sp == 'orca' and com_type == 'sp':

        read_lines = open(file_name,"r").readlines()

        #create input file
        orca_file_gen(read_lines,file_name.split('.')[0]+'.inp',bs_com,lot_com,genecp,args.aux_atoms_orca_sp,args.aux_basis_set_genecp_atoms_sp,args.aux_fit_genecp_atoms_sp,CHARGE,MULT,orca_aux_section,args,args.qm_input_sp,args.solvent_model_sp,args.solvent_name_sp,args.cpcm_input_sp,args.orca_scf_iters_sp,args.mdci_orca_sp,args.print_mini_orca_sp)

        # removes the initial com file
        os.remove(file_name)

    if args.sp == 'turbomole' and com_type == 'sp':
	    # Do stuff, MISSING PART
        pass


def read_log_file(w_dir,file):
    break_loop = False
    os.chdir(w_dir)
    try:
        outfile = open(file,"r")
        outlines = outfile.readlines()
    except FileNotFoundError:
        break_loop = True
        outfile, outlines = None, None

    return outlines, outfile, break_loop

def get_name_charge_multiplicity(outlines):
    name, CHARGE, MULT = '', None, None
    stop_name = 0
    # only for name an and charge
    for i,outline in enumerate(outlines):
        if stop_name == 2:
            break
        # Get the name of the compound (specified in the title)
        if outline.find('Symbolic Z-matrix:') > -1:
            name = outlines[i-2]
            stop_name += 1
        # Determine charge and multiplicity
        if outline.find("Charge = ") > -1:
            CHARGE = int(outline.split()[2])
            MULT = int(outline.split()[5].rstrip("\n"))
            stop_name += 1

    return name, CHARGE, MULT

def get_qcorr_params(outlines,default_lot,default_bs,preset_keywords_line,ts_opt,sp_calcs,charge_sp,mult_sp):
    TERMINATION,ERRORTYPE  = 'unfinished','unknown'
    stop_term, initial_IM_FREQS, NATOMS = 0, 0, 0
    symbol_z_line, gradgrad_line = None, None

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
        for i in reversed(range(60,len(outlines)-5)):
            if outlines[i].finds('NImag=') or outlines[i].finds('PG='):
                nimag_line = i
                # merge a few lines since 'NIMag=' sometimes appears in two different lines
                line_with_ifreq = outlines[i-1].rstrip("\n")+outlines[i].rstrip("\n")+outlines[i+1].rstrip("\n")
                for part_line in line_with_ifreq.split('\\'):
                    if part_line.find('NImag=') > -1:
                        initial_IM_FREQS = part_line.split('=')[1]
                        break

    # get functional, basis set and NATOMS
    if TERMINATION != "normal" or initial_IM_FREQS > 0  or sp_calcs != None:
        end_keywords = False
        keywords_line = ''

        # range optimized to skip Gaussian information (short for loop just to get charge, multiplicity and n of atoms)
        for i in range(60,len(outlines)):
            # retrieve the input line
            if outlines[i].startswith('#'):
                keywords_line += outlines[i].rstrip("\n")
                for j in range(i,i+5):
                    if outlines[j].finds('----------'):
                        end_keywords = True
                    if not end_keywords:
                        keywords_line += outlines[j].rstrip("\n")

            # get number of atoms
            if outlines[i].find('Symbolic Z-matrix:') > -1:
                symbol_z_line = i
            elif outlines[i].find('GradGrad') > -1:
                gradgrad_line = i
            if symbol_z_line != None and gradgrad_line != None:
                NATOMS = gradgrad_line-symbol_z_line-4

            # Determine charge and multiplicity
            if outlines[i].find("Charge = ") > -1:
                if charge_sp != None:
                    CHARGE = charge_sp
                else:
                    CHARGE = int(outlines[i].split()[2])
                if mult_sp != None:
                    MULT = mult_sp
                else:
                    MULT = int(outlines[i].split()[5].rstrip("\n"))

            if NATOMS != 0:
                break

        # find if the calculation was for a ground or transition state
        calcfc_found, ts_found = False, False
        for keyword in keywords_line.split():
            if keyword.lower().find('calcfc') > -1:
                calcfc_found = True
            if keyword.lower().find('ts') > -1:
                ts_found = True

        calc_type = 'ground_state'
        if calcfc_found and ts_found:
            calc_type = 'transition_state'

        # substract the imaginary frequency from the TS in normal terminations
        if calc_type == 'transition_state':
            initial_IM_FREQS -= 1

        # if a preset input line was defined in args.qm_input, this overwrites the automatic protocol
        if preset_keywords_line != None:
            keywords_line = preset_keywords_line
            # add the ts options if an external input line was defined
            if calc_type == 'transition_state':
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

        # detect calcs that finish before the functional or basis set was included
        if default_lot == None or default_bs == None:
            if ERRORTYPE == 'unknown':
                ERRORTYPE = 'before_E_calculation'

        # helps to fix SCF convergence errors
        if ERRORTYPE == 'SCFerror':
            if keywords_line.find(' scf=qc') > -1:
                pass
            else:
                keywords_line += ' scf=qc'

    return TERMINATION,ERRORTYPE,initial_IM_FREQS,NATOMS,CHARGE,MULT,keywords_line,calc_type,nimag_line


def get_input_geom(file, outlines, keywords_line, NATOMS, MULT, TERMINATION, ERRORTYPE, initial_IM_FREQS, nimag_line, calc_type, duplicate_filter, single_point, ifreq_threshold, s2_threshold, amplitude_ifreq, file_list, E_dup, H_dup, G_dup):

    stand_or,NATOMS,IM_FREQS = 0,0,0
    rms,stop_rms = 10000,0
    spin_contamination,freq_only = False, False

    # reverse loop to speed up the reading of the output files
    # (also skips final part for normal termination which was already read plus dipole/multipole information)
    if duplicate_filter:
        file_list.append(file)

    if TERMINATION == "normal":
        for i in reversed(range(0,nimag_line-(10*NATOMS))):
            # Sets where the final coordinates are inside the file
            if (outlines[i].find("Standard orientation:") > -1 or outlines[i].find("Input orientation:") > -1) and stop_get_details_stand_or !=1 :
                stand_or = i
                break

            if duplicate_filter:
                start_dup = False
                if not start_dup:
                    if outlines[i].find('E (Thermal)'):
                        start_dup = True
                else:
                    if outlines[i].find('Sum of electronic and zero-point Energies=') > -1:
                        E_dup.append()
                    elif outlines[i].find('Sum of electronic and thermal Enthalpies=') > -1:
                        H_dup.append()
                    elif outlines[i].find('Sum of electronic and thermal Free Energies=') > -1:
                        G_dup.append()

            if initial_IM_FREQS != 0:
                # Get the negative frequencies and their modes (up to three negative frequencies)
                if outlines[i].find(" Frequencies -- ") > -1:
                    FREQS, NORMALMODE  = [],[]
                    nfreqs = len(outlines[i].split())
                    for j in range(2, nfreqs):
                        if float(outlines[i].split()[j]) < 0 and abs(float(outlines[i].split()[j])) > abs(ifreq_threshold):
                            IM_FREQS += 1
                            FREQS.append(float(outlines[i].split()[j]))
                            NORMALMODE.append([])
                            for k in range(0,NATOMS):
                                NORMALMODE[(j-2)].append([float(outlines[i+5+k].split()[3*(j-2)+2]), float(outlines[i+5+k].split()[3*(j-2)+3]), float(outlines[i+5+k].split()[3*(j-2)+4])])

            # analyze spin contamination. If the <S**2> operator deviates more than (s2_threshold) %
            # compared to the expected value, the calculation contains spin contamination
            elif s2_threshold > 0.0:
                if outlines[i].find('S**2 before annihilation') > -1:
                    s2_value = outlines[i].split()[-1].rstrip("\n")
                    unpaired_e = MULT-1
                    spin = unpaired_e*0.5
                    s2_expected_value = spin*(spin+1)
                    spin_diff = abs(s2_value-s2_expected_value)
                    if spin_diff > abs(s2_threshold/100)*s2_expected_value:
                        spin_contamination = True

        # exclude TS imag frequency
        if calc_type == 'transition_state':
            IM_FREQS -= 1
            FREQS.pop(0)
            NORMALMODE.pop(0)

        if IM_FREQS > 0:
            ERRORTYPE = 'extra_imag_freq'

        if IM_FREQS < 0:
            ERRORTYPE = 'ts_no_imag_freq'

        if spin_contamination:
            ERRORTYPE = 'spin_contaminated'

    else:
        normal_term_found, opt_found = False, False
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

    # get atom types fro the calculation and the cartesian coordinates to use in each case
    ATOMTYPES, CARTESIANS = [],[]
    for i in range(stand_or+5,stand_or+5+NATOMS):
        massno = int(outlines[i].split()[1])
        if massno < len(periodic_table):
            atom_symbol = periodic_table[massno]
        else:
            atom_symbol = "XX"
        ATOMTYPES.append(atom_symbol)
        CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

    if ERRORTYPE == 'extra_imag_freq':
        CARTESIANS = fix_imag_freqs(NATOMS, CARTESIANS, FREQS, NORMALMODE, amplitude_ifreq)

    return file_list, E_dup, H_dup, G_dup, ERRORTYPE, keywords_line, ATOMTYPES, CARTESIANS


def fix_imag_freqs(NATOMS, CARTESIANS, FREQS, NORMALMODE, amplitude):
    # Multiplies the imaginary normal mode vector by this amount (from -1 to 1).
    amplitude = args.amplitude_ifreq # 0.2 is the default in the pyQRC script (GitHub, user: bobbypaton)
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

def create_folder_and_com(w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE,MULT,orca_aux_section):
    # creating new folder with new input gaussian files
    new_gaussian_input_files = w_dir_main+'/input_files/run_'+str(round_num+1)

    try:
        os.makedirs(new_gaussian_input_files)
    except OSError:
        if  os.path.isdir(new_gaussian_input_files):
            os.chdir(new_gaussian_input_files)
        else:
            raise
    os.chdir(new_gaussian_input_files)
    log.write('-> Creating new gaussian input file for {0} in {1}/{2}'.format(file,lot,bs))

    #error if both genecp and gen are
    if args.genecp_atoms != [] and args.gen_atoms != []:
        sys.exit("x  ERROR: Can't use Gen and GenECP at the same time")

    com_type = 'analysis'
    # adapt this for write_gaussian_input_file, also change the function in qprep_gaussian
    new_com_file(com_type,w_dir_initial,log,new_gaussian_input_files,file,args,keywords_line,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,bs,lot,bs_gcp,orca_aux_section)

def organize_outputs(w_dir,w_dir_main,round_num,file,TERMINATION,ERRORTYPE,w_dir_fin,check_geom_qcorr):
    source = w_dir+'/'+file

    finished,imag_freq,ts_no_imag_freq,spin_contaminated,duplicate_calc = 0,0,0,0,0
    atom_error, scf_error, before_E_error, other_error = 0,0,0,0
    unfinished, exp_rules_qcorr, check_geom_qcorr = 0,0,0

    if ERRORTYPE == None and TERMINATION == "normal":
        destination = w_dir_fin
        move_file(source, destination)
        finished += 1

    elif ERRORTYPE == 'extra_imag_freq':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/imag_freq/'
        move_file(source, destination)
        imag_freq += 1

    elif ERRORTYPE == 'ts_no_imag_freq':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/ts_no_imag_freq/'
        move_file(source, destination)
        ts_no_imag_freq += 1

    elif ERRORTYPE == 'spin_contaminated':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/spin_contaminated/'
        move_file(source, destination)
        spin_contaminated += 1

    elif ERRORTYPE == 'duplicate_calc':
        destination = w_dir_main+'/duplicates/run_'+str(round_num)
        move_file(source, destination)
        duplicate_calc += 1

    elif TERMINATION == "error":
        if ERRORTYPE == "atomicbasiserror":
            destination = w_dir_main +'/failed/run_'+str(round_num)+'/error/basis_set_error'
            atom_error += 1
        elif ERRORTYPE == "SCFerror":
            destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/scf_error'
            scf_error += 1
        elif ERRORTYPE == "before_E_calculation":
            destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/before_E_calculation'
            before_E_error += 1
        else:
            destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/unknown_error'
            other_error += 1
        move_file(source, destination)

    elif TERMINATION == "unfinished":
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/unfinished/'
        move_file(source, destination)
        unfinished += 1

    elif ERRORTYPE == 'fail_exp_rules':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/exp_rules_filter/'
        move_file(source, destination)
        exp_rules_qcorr += 1

    elif ERRORTYPE == 'isomerization':
        destination = w_dir_main+'/failed/run_'+str(round_num)+'/geometry_changed/'
        move_file(source, destination)
        check_geom_qcorr += 1

    return finished,unfinished,atom_error,scf_error,before_E_error,imag_freq,ts_no_imag_freq,spin_contaminated,duplicate_calc,other_error,exp_rules_qcorr,check_geom_qcorr

# Output file to mol converter
def output_to_mol(file,format,log):
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

# DEFINTION OF OUTPUT ANALYSER and NMR FILES CREATOR
def output_analyzer(log_files,com_files, w_dir, w_dir_main, args, w_dir_fin, w_dir_initial, log, ana_data, round_num):

    input_route = input_route_line(args)

    if round_num == 1:
        # moves the com files to respective folder
        for file in com_files:
            source = w_dir+'/'+file
            destination = w_dir_main +'/input_files/run_'+str(round_num)
            move_file(source, destination)

    # these lists will track duplicates if args.dup is activated
    file_list, E_dup, H_dup, G_dup = [],[],[],[]

    for file in log_files:

        # read the file
        log.write(file)
        outlines, outfile, break_loop = read_log_file(w_dir,file)

        if break_loop:
            break

        # get parameters from the calculations (i.e. termination and error types, n of atoms, charge, mult, etc)
        TERMINATION, ERRORTYPE, initial_IM_FREQS, NATOMS, CHARGE, MULT, keywords_line, calc_type, nimag_line = get_qcorr_params(outlines,args.level_of_theory,args.basis_set,args.qm_input,args.ts_input,args.sp,args.charge_sp,args.mult_sp)

        # get new coordinates as input files to fix error (unless the errors cannot be fixed automatically
        # such as in atomic basis errors or errors happening too early on the calculation)
        if ERRORTYPE not in ['before_E_calculation',"atomicbasiserror"]:

            # whats trhe difference between this function and the 2 functions below?!
            # get geometry parameters and frequency information
            if TERMINATION != "normal" or initial_IM_FREQS != 0 or args.s2_threshold > 0.0 or args.sp != None or args.dup:
                file_list, E_dup, H_dup, G_dup, ERRORTYPE, keywords_line, ATOMTYPES, CARTESIANS = get_input_geom(file, outlines, keywords_line, NATOMS, MULT, TERMINATION, ERRORTYPE, initial_IM_FREQS, nimag_line, calc_type, args.dup, args.sp, args.ifreq_cutoff, args.amplitude_ifreq, args.s2_threshold, file_list, E_dup, H_dup, G_dup)

        # is it necessary to close the file?!
        # close the file
        outfile.close()

        # this part filters off conformers based on user-defined exp_rules
        passing_rules = True
        valid_mol_gen = True
        passing_geom = True
        format_file = file.split('.')[1]
        if TERMINATION == "normal" and ERRORTYPE == None:
            if len(args.exp_rules) >= 1:
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

            if args.check_geom and ERRORTYPE == None:
                log.write("  ----- Geometrical check will be applied to the output file -----\n")
                # this creates a mol object from the optimized log file
                mol,ob_compat,rdkit_compat = output_to_mol(file,format_file,log)
                # this creates a mol object from the input file
                try:
                    os.chdir(w_dir_main +'/input_files/run_'+str(round_num))
                    com_2_xyz_2_sdf(args.input,args.default_charge,os.path.splitext(file)[0]+'.com')
                    mol2,ob_compat,rdkit_compat = output_to_mol(file,'xyz',log)
                    passing_geom = check_geom_filter(mol,mol2,args.length_criteria)
                    if not passing_geom:
                        ERRORTYPE = 'isomerization'
                    # remove created files
                    os.remove(file.split('.')[0]+'.xyz')
                    os.remove(file.split('.')[0]+'.sdf')
                    os.remove(file.split('.')[0]+'.mol')

                except FileNotFoundError:
                    log.write("x  No com file were found for "+file+", the check_geom test will be disabled for this calculation")

                os.chdir(w_dir)
                os.remove(file.split('.')[0]+'.mol')

        elif args.check_geom and not valid_mol_gen:
            log.write("The file could not be converted into a mol object, check_geom test will be disabled\n")

        # detects if this calculation is a duplicate
        if args.dup:
            for i,energy_value in enumerate(E_dup):
                if abs(energy_value - E_dup[file_list.index(file)]) < abs(args.dup_threshold):
                    if abs(H_dup[i] - H_dup[file_list.index(file)]) < abs(args.dup_threshold):
                        if abs(G_dup[i] - G_dup[file_list.index(file)]) < abs(args.dup_threshold):
                            ERRORTYPE = 'duplicate_calc'

        # This part places the calculations in different folders depending on the type of termination
        finished,unfinished,atom_error,scf_error,before_E_error,imag_freq,ts_no_imag_freq,spin_contaminated,duplicate_calc,other_error,exp_rules_qcorr,check_geom_qcorr = organize_outputs(w_dir,w_dir_main,round_num,file,TERMINATION,ERRORTYPE,w_dir_fin,check_geom_qcorr)

        # check if gen or genecp are active
        # right now, QCORR is only working with Gaussian output files
        if args.genecp_atoms != [] or args.gen_atoms != []:
            ecp_list = check_for_gen_or_genecp(ATOMTYPES,args,'analysis','gaussian')

        # create folders and set level of theory in COM files to fix imaginary freqs or not normal terminations
        if TERMINATION != "normal" or ERRORTYPE != None:
            create_folder_and_com(w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE, MULT, orca_aux_section)

        # this part creates input files for single-point energy calcs after reading from normally finished log files
        if TERMINATION == "normal" and ERRORTYPE == None:
            if args.sp != None or args.nics:

                if args.sp == 'gaussian':
                    # creating new folder with new input Gaussian files
                    single_point_input_files = w_dir_fin+'/../G16-SP_input_files'

                elif args.sp == 'orca':
                    # creating new folder with new input ORCA files
                    single_point_input_files = w_dir_fin+'/../ORCA-SP_input_files'

                elif args.sp == 'turbomole':
                  # creating new folder with new input Turbomole files
                  single_point_input_files = w_dir_fin+'/../TURBOMOLE-SP_input_files'

                if args.nics:
                    nics_input_files = w_dir_fin+'/../G16-NICS_input_files'

                # Options for genecp
                if args.sp == 'gaussian':
                    ecp_list = check_for_gen_or_genecp(ATOMTYPES,args,'sp','gaussian')

                elif args.sp == 'orca':
                    orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','orca')

                elif args.sp == 'turbomole':
                    # Check for gen or genecp
                    pass

                # NEED TO FIX THIS PART FOR GENECP WITH SP!
                if genecp == None:
                    basis_set_for_genecp = args.basis_set_sp
                elif genecp == 'genecp' or genecp == 'gen':
                    basis_set_for_genecp = args.basis_set_genecp_atoms_sp

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
                        if not os.path.isdir(single_point_input_files+'/'+dir_name):
                            os.makedirs(single_point_input_files+'/'+dir_name)
                        os.chdir(single_point_input_files+'/'+dir_name)
                        new_com_file('sp',w_dir_initial,log,single_point_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

                    if args.nics:
                        if not os.path.isdir(nics_input_files+'/'+dir_name):
                            os.makedirs(nics_input_files+'/'+dir_name)
                        os.chdir(nics_input_files+'/'+dir_name)
                        new_com_file('nics',w_dir_initial,log,nics_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

    #write to csv ana_data
    ana_data.at[0,'Total files'] = len(log_files)
    ana_data.at[0,'Normal termination'] = finished
    ana_data.at[0,'Imaginary frequencies'] = imag_freq
    ana_data.at[0,'Spin contamination'] = spin_contaminated
    ana_data.at[0,'TS with no imag. freq.'] = ts_no_imag_freq
    ana_data.at[0,'SCF error'] = scf_error
    ana_data.at[0,'Error before SCF'] = before_E_error
    ana_data.at[0,'Basis set error'] =  atom_error
    ana_data.at[0,'Other errors'] = other_error
    ana_data.at[0,'Unfinished'] = unfinished
    if args.dup:
        ana_data.at[0,'Duplicates'] = duplicate_calc
    if len(args.exp_rules) >= 1:
        ana_data.at[0,'Exp_rules filter'] = exp_rules_qcorr
    if args.check_geom:
        ana_data.at[0,'Geometry changed'] = check_geom_qcorr

    if not os.path.isdir(w_dir_main+'/csv_files/'):
        os.makedirs(w_dir_main+'/csv_files/')
    ana_data.to_csv(w_dir_main+'/csv_files/Analysis-Data-QCORR-run_'+str(round_num)+'.csv',index=False)


 # DETECTION AND LISTING OF GEN/GENECP FROM COM FILES
def check_for_gen_or_genecp(ATOMTYPES,args,type_of_check,program_gen):

    ecp_list,ecp_atoms_include,orca_aux_section = [],[],False
    if type_of_check == 'analysis':
        if program_gen == 'gaussian':
            if args.genecp_atoms != []:
                ecp_atoms_include = args.genecp_atoms
            elif args.gen_atoms != []:
                ecp_atoms_include = args.gen_atoms
        elif program_gen == 'orca':
            if args.aux_atoms_orca != []:
                aux_atoms_include = args.aux_atoms_orca

    elif type_of_check == 'sp':
        if program_gen == 'gaussian':
            if args.genecp_atoms_sp != []:
                ecp_atoms_include = args.genecp_atoms_sp
            elif args.gen_atoms_sp != []:
                ecp_atoms_include = args.gen_atoms_sp
        elif program_gen == 'orca':
            if args.aux_atoms_orca_sp != []:
                aux_atoms_include = args.aux_atoms_orca_sp

    for _,atomtype in enumerate(ATOMTYPES):
        if program_gen == 'gaussian':
            if atomtype not in ecp_list and atomtype in periodic_table and atomtype in ecp_atoms_include:
                ecp_list.append(atomtype)
        elif program_gen == 'orca':
            if atomtype in aux_atoms_include:
                orca_aux_section = True

    if program_gen == 'gaussian':
    	return ecp_list
    elif program_gen == 'orca':
        return orca_aux_section


# CHECKS THE FOLDER OF FINAL LOG FILES
def check_for_final_folder(w_dir):
    ini_com_folder = sum(dirs.count('input_files') for _, dirs, _ in os.walk(w_dir))
    if ini_com_folder == 0:
        return w_dir, 1
    else:
        num_com_folder = sum([len(d) for r, d, folder in os.walk(w_dir+'/input_files')])
        w_dir = w_dir+'/input_files/run_'+str(num_com_folder)
        return w_dir, num_com_folder

def check_for_final_folderv2(w_dir):
	base_folder = Path(w_dir)
	inputs_folder = base_folder/'input_files'
	if not inputs_folder.exists():
		return w_dir, 1
	folders = [item for item in inputs_folder.iterdir() if item.isdir() and 'run_' in item.stem]
	get_number = lambda x: int(x.stem.rsplit('_',1)[1])
	last_folder = sorted(folders,key=get_number)[-1]
	last_number = get_number(last_folder)
	return str(last_folder), last_number
