#!/usr/bin/env python

######################################################.
# 		  This file stores all the functions 	     #
# 	   	    used in the LOG file analyzer	    	 #
######################################################.

import os
import sys
import subprocess
import shutil
import numpy as np
from pyconfort.qprep_gaussian import input_route_line,check_for_gen_or_genecp,write_genecp
from pyconfort.argument_parser import possible_atoms

possible_atoms = possible_atoms()

def moving_files(source, destination):
	if not os.path.isdir(destination):
		os.makedirs(destination)
	try:
		shutil.move(source, destination)
	except (FileExistsError,shutil.Error):
		os.chdir(destination)
		os.remove(source.split('/')[-1])
		shutil.move(source, destination)

def write_header_and_coords(fileout,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS):
	fileout.write("%mem="+str(args.mem)+"\n")
	fileout.write("%nprocshared="+str(args.nprocs)+"\n")
	fileout.write("# "+keywords_opt+"\n")
	fileout.write("\n")
	fileout.write(name+"\n")
	fileout.write(str(CHARGE)+' '+str(MULT)+'\n')
	for atom in range(0,NATOMS):
		fileout.write('{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}'.format(ATOMTYPES[atom], CARTESIANS[atom][0],  CARTESIANS[atom][1],  CARTESIANS[atom][2]))
		fileout.write("\n")
	fileout.write("\n")

# CREATION OF COM FILES
def new_com_file(com_type,w_dir,w_dir_initial,new_gaussian_input_files,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_com,lot_com,bs_gcp_com):

	if com_type == 'sp':
		if args.suffix_sp == 'None':
			file_name = file.split(".")[0]+'.com'
		else:
			file_name = file.split(".")[0]+'_'+args.suffix_sp+'.com'
	elif com_type == 'analysis':
		file_name = file.split(".")[0]+'.com'

	fileout = open(file_name, "w")

	write_header_and_coords(fileout,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS)

	if args.sp and TERMINATION == "normal" and IM_FREQS == 0 :
		fileout.write(args.last_line_for_sp)
		fileout.write('\n\n')

	# write genecp/gen part
	if genecp == 'genecp' or  genecp == 'gen':
		if com_type == 'sp':
			type_gen = 'sp'
		elif com_type == 'analysis':
			type_gen = 'qcorr'

		write_genecp(type_gen,fileout,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,bs_com,lot_com,bs_gcp_com,args,w_dir_initial,new_gaussian_input_files)

	fileout.close()

	# #submitting the gaussian file on summit
	if args.qsub:
		cmd_qsub = [args.submission_command, file_name]
		subprocess.call(cmd_qsub)

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

def get_initial_variables():
	rms = 10000
	stop_rms,stand_or,NATOMS,IM_FREQS,freqs_so_far,stop_name,stop_term,nfreqs,dist_rot_or = 0,0,0,0,0,0,0,0,0
	ATOMTYPES, CARTESIANS, FREQS, READMASS, FORCECONST, NORMALMODE = [],[],[],[],[],[]
	TERMINATION,ERRORTYPE  = 'unfinished','unknown'

	return rms,stop_rms,stand_or,NATOMS,IM_FREQS,freqs_so_far,stop_name,stop_term,nfreqs,ATOMTYPES,CARTESIANS,FREQS,READMASS,FORCECONST,NORMALMODE,TERMINATION,ERRORTYPE,dist_rot_or

def get_name_charge_multiplicity(outlines,stop_name):
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

def get_termination_type(outlines,stop_term,TERMINATION,ERRORTYPE):
	# use reversed loops to find type of termination (faster than forward loops)
	for i in reversed(range(len(outlines)-15,len(outlines))):
		if stop_term == 1:
			break
		# Determine the kind of job termination
		if outlines[i].find("Normal termination") > -1:
			TERMINATION = "normal"
			stop_term += 1
		elif outlines[i].find("Error termination") > -1:
			TERMINATION = "error"
			if outlines[i-1].find("Atomic number out of range") > -1 or outlines[i-1].find("basis sets are only available") > -1 :
				ERRORTYPE = "atomicbasiserror"
			if outlines[i-3].find("SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error") > -1:
				ERRORTYPE = "SCFerror"
			stop_term += 1
	return TERMINATION,ERRORTYPE

def get_geom_and_freq_for_normal(outlines, args, TERMINATION, NATOMS, FREQS, NORMALMODE, IM_FREQS, READMASS, FORCECONST, nfreqs, freqs_so_far, rms, stop_rms, dist_rot_or, stand_or):
	stop_get_details_stand_or, stop_get_details_dis_rot,finding_freq_line,stop_finding_freq_line = 0,0,0,0
	# reverse loop to speed up the reading of the output files
	for i in reversed(range(0,len(outlines))):
		if TERMINATION == "normal":
			if not args.frequencies:
				if stop_get_details_stand_or == 1 and stop_get_details_dis_rot == 1:
					break
			else:
				if stop_get_details_stand_or == 1 and stop_get_details_dis_rot == 1 and stop_finding_freq_line ==1 :
					break
			# Sets where the final coordinates are inside the file
			if stop_get_details_dis_rot !=1 and (outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") >-1) :
				if outlines[i-1].find("-------") > -1:
					dist_rot_or = i
					stop_get_details_dis_rot += 1
			if outlines[i].find("Standard orientation") > -1 and stop_get_details_stand_or !=1 :
				stand_or = i
				NATOMS = dist_rot_or-i-6
				stop_get_details_stand_or += 1
			if args.frequencies:
				if outlines[i].find(" Harmonic frequencies") > -1 and stop_finding_freq_line !=1 :
					finding_freq_line = i

	if args.frequencies:
		for i in range(finding_freq_line,len(outlines)):
			# Get the frequencies and identifies negative frequencies
			if outlines[i].find(" Frequencies -- ") > -1:
				nfreqs = len(outlines[i].split())
				for j in range(2, nfreqs):
					FREQS.append(float(outlines[i].split()[j]))
					NORMALMODE.append([])
					if float(outlines[i].split()[j]) < 0 and np.absolute(float(outlines[i].split()[j])) > args.ifreq_cutoff:
						IM_FREQS += 1
				for j in range(3, nfreqs+1):
					READMASS.append(float(outlines[i+1].split()[j]))
				for j in range(3, nfreqs+1):
					FORCECONST.append(float(outlines[i+2].split()[j]))
				for j in range(0,NATOMS):
					for k in range(0, nfreqs-2):
						NORMALMODE[(freqs_so_far + k)].append([float(outlines[i+5+j].split()[3*k+2]), float(outlines[i+5+j].split()[3*k+3]), float(outlines[i+5+j].split()[3*k+4])])
				freqs_so_far = freqs_so_far + nfreqs - 2
			if TERMINATION != "normal":
				if outlines[i].find('Cartesian Forces:  Max') > -1:
					try:
						if float(outlines[i].split()[5]) < rms:
							rms = float(outlines[i].split()[5])
							stop_rms = i
					except ValueError:
						rms = 10000
	return TERMINATION, NATOMS, FREQS, NORMALMODE, IM_FREQS, READMASS, FORCECONST, nfreqs, freqs_so_far, rms, stop_rms, dist_rot_or, stand_or

def get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS):
	for i in range(stand_or+5,stand_or+5+NATOMS):
		massno = int(outlines[i].split()[1])
		if massno < len(possible_atoms):
			atom_symbol = possible_atoms[massno]
		else:
			atom_symbol = "XX"
		ATOMTYPES.append(atom_symbol)
		CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

	return ATOMTYPES, CARTESIANS

def get_coords_not_normal(outlines, stop_rms, stand_or, dist_rot_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS):
	if stop_rms == 0:
		last_line = len(outlines)
	else:
		last_line = stop_rms
	stop_get_details_stand_or = 0
	stop_get_details_dis_rot = 0
	for i in reversed(range(0,last_line)):
		if stop_get_details_stand_or == 1 and stop_get_details_dis_rot == 1:
			break
		# Sets where the final coordinates are inside the file
		if outlines[i].find("Standard orientation") > -1 and stop_get_details_stand_or != 1:
			stand_or = i
			NATOMS = dist_rot_or-i-6
			stop_get_details_stand_or += 1
		if stop_get_details_stand_or != 1 and (outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") >-1):
			if outlines[i-1].find("-------") > -1:
				dist_rot_or = i
				stop_get_details_dis_rot += 1

	ATOMTYPES, CARTESIANS = get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)

	return ATOMTYPES, CARTESIANS, NATOMS

def fix_imag_freqs(NATOMS, CARTESIANS, args, FREQS, NORMALMODE):
	# Multiplies the imaginary normal mode vector by this amount (from -1 to 1).
	amplitude = args.amplitude_ifreq # 0.2 is the default in the pyQRC script (GitHub, user: bobbypaton)
	shift = []

	# Save the original Cartesian coordinates before they are altered
	orig_carts = []
	for atom in range(0,NATOMS):
		orig_carts.append([CARTESIANS[atom][0], CARTESIANS[atom][1], CARTESIANS[atom][2]])

	# could get rid of atomic units here, if zpe_rat definition is changed
	for mode,_ in enumerate(FREQS):
		# Either moves along any and all imaginary freqs, or a specific mode requested by the user
		if FREQS[mode] < 0.0:
			shift.append(amplitude)
		else:
			shift.append(0.0)

		# The starting geometry is displaced along each normal mode according to the random shift
		for atom in range(0,NATOMS):
			for coord in range(0,3):
				CARTESIANS[atom][coord] = CARTESIANS[atom][coord] + NORMALMODE[mode][atom][coord] * shift[mode]

	return CARTESIANS

def create_folder_and_com(w_dir,w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE,MULT):
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

	ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp = check_for_gen_or_genecp(ATOMTYPES,args)

	#error if both genecp and gen are
	if ecp_genecp_atoms and ecp_gen_atoms:
		sys.exit("x  ERROR: Can't use Gen and GenECP at the same time")

	if ERRORTYPE == 'SCFerror':
		input_route += ' scf=qc'
	if genecp == 'genecp' or  genecp == 'gen':
		keywords_opt = lot +'/'+ genecp+' '+ input_route
	else:
		keywords_opt = lot +'/'+ bs +' '+ input_route

	com_type = 'analysis'
	new_com_file(com_type, w_dir,w_dir_initial,new_gaussian_input_files,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs,lot,bs_gcp)

def create_folder_move_log_files(w_dir,w_dir_main,round_num,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin,finished,unfinished,atom_error,scf_error,imag_freq,other_error):
	source = w_dir+'/'+file

	if IM_FREQS == 0 and TERMINATION == "normal":
		destination = w_dir_fin
		moving_files(source, destination)
		finished += 1

	if IM_FREQS > 0:
		destination = w_dir_main+'/failed/run_'+str(round_num)+'/imag_freq/'
		moving_files(source, destination)
		imag_freq += 1

	if IM_FREQS == 0 and TERMINATION == "error":
		if ERRORTYPE == "atomicbasiserror":
			destination = w_dir_main +'/failed/run_'+str(round_num)+'/error/basis_set_error'
			atom_error += 1
		elif ERRORTYPE == "SCFerror":
			destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/scf_error'
			scf_error += 1
		else:
			destination = w_dir_main+'/failed/run_'+str(round_num)+'/error/unknown_error'
			other_error += 1
		moving_files(source, destination)

	elif IM_FREQS == 0 and TERMINATION == "unfinished":
		destination = w_dir_main+'/failed/run_'+str(round_num)+'/unfinished/'
		moving_files(source, destination)
		unfinished += 1

	return finished,unfinished,atom_error,scf_error,imag_freq,other_error

# DEFINTION OF OUTPUT ANALYSER and NMR FILES CREATOR
def output_analyzer(duplicates,log_files,com_files, w_dir, w_dir_main,lot, bs, bs_gcp, args, w_dir_fin, w_dir_initial, log, ana_data, round_num):

	input_route = input_route_line(args)
	finished,unfinished,atom_error,scf_error,imag_freq,other_error = 0,0,0,0,0,0

	if round_num == 1:
		#moves the comfiles to respective folder
		for file in com_files:
			source = w_dir+'/'+file
			destination = w_dir_main +'/input_files/run_'+str(round_num)
			moving_files(source, destination)

	for file in log_files:
		# read the file
		print(file)
		outlines, outfile, break_loop = read_log_file(w_dir,file)

		if break_loop:
			break

		# get initial parameters
		rms,stop_rms,stand_or,NATOMS,IM_FREQS,freqs_so_far,stop_name,stop_term,nfreqs,ATOMTYPES,CARTESIANS,FREQS,READMASS,FORCECONST,NORMALMODE,TERMINATION,ERRORTYPE,dist_rot_or = get_initial_variables()

		# get name, charge and multiplicity
		name, CHARGE, MULT = get_name_charge_multiplicity(outlines,stop_name)

		# get termination type
		TERMINATION,ERRORTYPE = get_termination_type(outlines,stop_term,TERMINATION,ERRORTYPE)

		# get geometry parameters and frequency information
		TERMINATION, NATOMS, FREQS, NORMALMODE, IM_FREQS, READMASS, FORCECONST, nfreqs, freqs_so_far, rms, stop_rms, dist_rot_or, stand_or = get_geom_and_freq_for_normal(outlines, args, TERMINATION, NATOMS, FREQS, NORMALMODE, IM_FREQS, READMASS, FORCECONST, nfreqs, freqs_so_far, rms, stop_rms, dist_rot_or, stand_or)

		# Get the coordinates for jobs that finished well with and without imag. freqs
		if TERMINATION == "normal" and IM_FREQS>0:
			ATOMTYPES, CARTESIANS = get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)

		# Get he coordinates for jobs that did not finished or finished with an error
		if TERMINATION != "normal":
			ATOMTYPES, CARTESIANS,NATOMS = get_coords_not_normal(outlines, stop_rms, stand_or, dist_rot_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)
		# This part fixes imaginary freqs (if any)
		if IM_FREQS > 0:
			CARTESIANS = fix_imag_freqs(NATOMS, CARTESIANS, args, FREQS, NORMALMODE)

		#close the file
		outfile.close()

		# This part places the calculations in different folders depending on the type of termination
		finished,unfinished,atom_error,scf_error,imag_freq,other_error = create_folder_move_log_files(w_dir,w_dir_main,round_num,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin,finished,unfinished,atom_error,scf_error,imag_freq,other_error)

		# check if gen or genecp are active
		ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp =  check_for_gen_or_genecp(ATOMTYPES,args)

		# create folders and set level of theory in COM files to fix imaginary freqs or not normal terminations
		if IM_FREQS > 0 or TERMINATION != "normal" and not os.path.exists(w_dir_main+'/failed/run_'+str(round_num)+'/error/basis_set_error/'+file):
			create_folder_and_com(w_dir,w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE, MULT)

		# adding in the NMR componenet only to the finished files after reading from normally finished log files
		if args.sp and TERMINATION == "normal" and IM_FREQS == 0:
			#get coordinates
			ATOMTYPES, CARTESIANS = get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)
			# creating new folder with new input gaussian files
			single_point_input_files = w_dir_fin+'/../G16-SP_input_files'
			# Options for genecp
			ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp = check_for_gen_or_genecp(ATOMTYPES,args)

			# Sets the folder and find the log files to analyze
			for lot_sp in args.level_of_theory_sp:
				for bs_sp in args.basis_set_sp:
					for bs_gcp_sp in args.basis_set_genecp_atoms_sp:
						log.write('-> Creating new single point files files for {0} in {1}/{2}-{3}\n'.format(file,single_point_input_files,lot_sp,bs_sp))
						dir_name = str(lot_sp) + '-' + str(bs_sp)
						if not os.path.isdir(single_point_input_files+'/'+dir_name):
							os.makedirs(single_point_input_files+'/'+dir_name)
						os.chdir(single_point_input_files+'/'+dir_name)

						if genecp == 'genecp' or  genecp == 'gen':
							keywords_opt = lot_sp+'/'+ genecp+' '+ args.input_for_sp
						else:
							keywords_opt = lot_sp+'/'+ bs_sp+' '+ args.input_for_sp
						if args.empirical_dispersion_sp != 'None':
							keywords_opt += ' empiricaldispersion={0}'.format(args.empirical_dispersion_sp)
						if args.solvent_model_sp != 'gas_phase':
							keywords_opt += ' scrf=({0},solvent={1})'.format(args.solvent_model_sp,args.solvent_name_sp)

						if args.charge_sp != 'None':
							CHARGE = args.charge_sp
						if args.mult_sp != 'None':
							MULT = args.mult_sp
						com_type = 'sp'
						new_com_file(com_type, w_dir,w_dir_initial,single_point_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp)

	#write to csv ana_data
	ana_data.at[0,'Total files'] = len(log_files)
	ana_data.at[0,'Normal termination'] = finished
	ana_data.at[0,'Imaginary frequencies'] = imag_freq
	ana_data.at[0,'SCF error'] = scf_error
	ana_data.at[0,'Basis set error'] =  atom_error
	ana_data.at[0,'Other errors'] = other_error
	ana_data.at[0,'Unfinished'] = unfinished
	if duplicates != False:
		ana_data.at[0,'Duplicates'] = duplicates

	if not os.path.isdir(w_dir_main+'/csv_files/'):
		os.makedirs(w_dir_main+'/csv_files/')
	ana_data.to_csv(w_dir_main+'/csv_files/Analysis-Data-QCORR-run_'+str(round_num)+'.csv',index=False)

# CHECKS THE FOLDER OF FINAL LOG FILES
def check_for_final_folder(w_dir):
	ini_com_folder = sum(dirs.count('input_files') for _, dirs, _ in os.walk(w_dir))
	if ini_com_folder == 0:
		return w_dir, 1
	else:
		num_com_folder = sum([len(d) for r, d, folder in os.walk(w_dir+'/input_files')])
		w_dir = w_dir+'/input_files/run_'+str(num_com_folder)
		return w_dir, num_com_folder

# CHECKING FOR DUPLICATES
def dup_calculation(val, w_dir,w_dir_main, args, log,round_num):

	# GoodVibes must be installed as a module (through pip or conda)
	cmd_dup = ['python', '-m', 'goodvibes', '--dup']
	for file in val:
		cmd_dup.append(file)
	subprocess.call(cmd_dup)

	#reading the txt files to get the DUPLICATES
	dup_file_list = []
	duplines = open('Goodvibes_output.dat',"r").readlines()

	for i,_ in enumerate(duplines):
		if duplines[i].find('duplicate') > -1:
			dup_file_list.append(duplines[i].split(' ')[2])

	duplicates = len(dup_file_list)

	#move the files to specific directory
	destination = w_dir_main+'/duplicates/run_'+str(round_num)
	source = w_dir_main+'/Goodvibes_output.dat'
	moving_files(source, destination)

	for source in dup_file_list:
		#finding the extension
		for file in val:
			if file.split('.')[0] == source:
				ext=file.split('.')[1]
		source=source+'.'+ext
		moving_files(source, destination)

	return duplicates
