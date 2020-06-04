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
import pandas as pd
from pyconfort.writer_functions import input_route_line
from pyconfort.argument_parser import possible_atoms

possible_atoms = possible_atoms()

def moving_log_files(source, destination, file):
	try:
		os.makedirs(destination)
		shutil.move(source, destination)
	except OSError:
		if  os.path.isdir(destination) and not os.path.exists(destination+file):
			shutil.move(source, destination)
		else:
			raise

# DETECTION OF GEN/GENECP
def check_for_gen_or_genecp(ATOMTYPES,args):
	# Options for genecp
	ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp = [],False,False,None

	for _,atomtype in enumerate(ATOMTYPES):
		if atomtype not in ecp_list and atomtype in possible_atoms:
			ecp_list.append(atomtype)
		if atomtype in args.genecp_atoms:
			ecp_genecp_atoms = True
		if atomtype in args.gen_atoms:
			ecp_gen_atoms = True
	if ecp_gen_atoms:
		genecp = 'gen'
	if ecp_genecp_atoms:
		genecp = 'genecp'

	return ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp

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

def write_genecp(fileout,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,bs_com,lot_com,bs_gcp_com,args,w_dir_initial,w_dir):
	for _,element_ecp in enumerate(ecp_list):
		if element_ecp not in (args.genecp_atoms or args.gen_atoms):
			fileout.write(element_ecp+' ')
	fileout.write('0\n')
	fileout.write(bs_com+'\n')
	fileout.write('****\n')

	if len(bs_gcp_com.split('.')) > 1:
		if bs_gcp_com.split('.')[1] == 'txt' or bs_gcp_com.split('.')[1] == 'yaml':
			os.chdir(w_dir_initial)
			read_lines = open(bs_gcp_com,"r").readlines()
			os.chdir(w_dir)
			#getting the title line
			for line in read_lines:
				fileout.write(line)
			fileout.write('\n\n')
	else:
		for _,element_ecp in enumerate(ecp_list):
			if element_ecp in args.genecp_atoms :
				fileout.write(element_ecp+' ')
			elif element_ecp in args.gen_atoms :
				fileout.write(element_ecp+' ')
		fileout.write('0\n')
		fileout.write(bs_gcp_com+'\n')
		fileout.write('****\n\n')
		if ecp_genecp_atoms:
			for _,element_ecp in enumerate(ecp_list):
				if element_ecp in args.genecp_atoms:
					fileout.write(element_ecp+' ')
			fileout.write('0\n')
			fileout.write(bs_gcp_com+'\n\n')

# CREATION OF COM FILES
def new_com_file(w_dir,w_dir_initial,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_com,lot_com,bs_gcp_com):
	fileout = open(file.split(".")[0]+'.com', "w")

	write_header_and_coords(fileout,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS)

	if genecp == 'genecp' or  genecp == 'gen':
		write_genecp(fileout,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,bs_com,lot_com,bs_gcp_com,args,w_dir_initial,w_dir)

	if args.sp and TERMINATION == "normal" and IM_FREQS == 0 :
		fileout.write(args.last_line_for_sp)
		fileout.write('\n\n')

	fileout.close()

def read_log_file(w_dir,file):
	break_loop = False
	os.chdir(w_dir)
	try:
		outfile = open(file,"r")
	except FileNotFoundError:
		break_loop = True
	outlines = outfile.readlines()

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
			ERRORTYPE = ERRORTYPE
			stop_term += 1
		elif outlines[i].find("Error termination") > -1:
			TERMINATION = "error"
			if outlines[i-1].find("Atomic number out of range") > -1 or outlines[i-1].find("basis sets are only available") > -1 :
				ERRORTYPE = "atomicbasiserror"
			if outlines[i-3].find("SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error") > -1:
				ERRORTYPE = "SCFerror"
			stop_term += 1
		else:
			TERMINATION = TERMINATION
			ERRORTYPE = ERRORTYPE
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
					if float(outlines[i].split()[j]) < 0.0:
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
					if float(outlines[i].split()[5]) < rms:
						rms = float(outlines[i].split()[5])
						stop_rms = i
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

def create_folder_and_com(w_dir,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE,MULT):
	# creating new folder with new input gaussian files
	new_gaussian_input_files = w_dir+'/new_gaussian_input_files/'

	try:
		os.makedirs(new_gaussian_input_files)
	except OSError:
		if  os.path.isdir(new_gaussian_input_files):
			os.chdir(new_gaussian_input_files)
		else:
			raise
	os.chdir(new_gaussian_input_files)
	log.write('-> Creating new gaussian input file for {0} in {1}/{2}'.format(file,lot,bs))

	ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp =  check_for_gen_or_genecp(ATOMTYPES,args)

	#error if both genecp and gen are
	if ecp_genecp_atoms and ecp_gen_atoms:
		sys.exit("ERROR: Can't use Gen and GenECP at the same time")

	if ERRORTYPE == 'SCFerror':
		input_route += ' scf=qc'
	if genecp == 'genecp' or  genecp == 'gen':
		keywords_opt = lot +'/'+ genecp+' '+ input_route
	else:
		keywords_opt = lot +'/'+ bs +' '+ input_route

	new_com_file(w_dir,w_dir_initial,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs,lot,bs_gcp)

def create_folder_move_log_files(w_dir,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin):
	source = w_dir+'/'+file

	if IM_FREQS == 0 and TERMINATION == "normal":
		destination = w_dir_fin
		moving_log_files(source, destination, file)

	if IM_FREQS > 0:
		destination = w_dir+'/imaginary_frequencies/'
		moving_log_files(source, destination, file)

	if IM_FREQS == 0 and TERMINATION == "error":
		if ERRORTYPE == "atomicbasiserror":
			destination = w_dir+'/failed_error/atomic_basis_error'
		elif ERRORTYPE == "SCFerror":
			destination = w_dir+'/failed_error/SCF_error'
		else:
			destination = w_dir+'/failed_error/unknown_error'
		moving_log_files(source, destination, file)

	elif IM_FREQS == 0 and TERMINATION == "unfinished":
		destination = w_dir+'/failed_unfinished/'
		moving_log_files(source, destination, file)

# DEFINTION OF OUTPUT ANALYSER and NMR FILES CREATOR
def output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin, w_dir_initial, log):

	input_route = input_route_line(args)

	for file in log_files:
		# read the file
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
		if TERMINATION == "normal":
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
		create_folder_move_log_files(w_dir,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin)

		# check if gen or genecp are active
		ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp =  check_for_gen_or_genecp(ATOMTYPES,args)

		# create folders and set level of theory in COM files to fix imaginary freqs or not normal terminations
		if IM_FREQS > 0 or TERMINATION != "normal" and not os.path.exists(w_dir+'/failed_error/atomic_basis_error/'+file):
			create_folder_and_com(w_dir,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE, MULT)

		# adding in the NMR componenet only to the finished files after reading from normally finished log files
		if args.sp and TERMINATION == "normal" and IM_FREQS == 0:
			# creating new folder with new input gaussian files
			single_point_input_files = w_dir_fin+'/single_point_input_files'
			# Options for genecp
			ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp =  check_for_gen_or_genecp(ATOMTYPES,args)

			# Sets the folder and find the log files to analyze
			for lot_sp in args.level_of_theory_sp:
				for bs_sp in args.basis_set_sp:
					for bs_gcp_sp in args.basis_set_genecp_atoms_sp:
						log.write('-> Creating new single point files files for {0} in {1}/{2}\n'.format(file,lot_sp,bs_sp))
						folder = single_point_input_files + '/' + str(lot_sp) + '-' + str(bs_sp)
						try:
							os.makedirs(folder)
							os.chdir(folder)
						except OSError:
							if os.path.isdir(folder):
								os.chdir(folder)
							else:
								raise
						if genecp == 'genecp' or  genecp == 'gen':
							if args.dispersion_correction_sp:
								if args.solvent_model_sp == 'gas_phase':
									keywords_opt = lot_sp+'/'+ genecp+' '+ args.input_for_sp + ' empiricaldispersion={0}'.format(args.empirical_dispersion_sp)
								else:
									keywords_opt = lot_sp+'/'+ genecp+' '+ args.input_for_sp + ' scrf=({0},solvent={1}) empiricaldispersion={2}  '.format(args.solvent_model_sp,args.solvent_name_sp,args.empirical_dispersion_sp)
							else:
								if args.solvent_model_sp == 'gas_phase':
									keywords_opt = lot_sp+'/'+ genecp+' '+ args.input_for_sp
								else:
									keywords_opt = lot_sp+'/'+ genecp+' '+ args.input_for_sp + ' scrf=({0},solvent={1}) '.format(args.solvent_model_sp,args.solvent_name_sp)
						else:
							if args.dispersion_correction_sp:
								if args.solvent_model_sp == 'gas_phase':
									keywords_opt = lot_sp+'/'+ bs_sp+' '+ args.input_for_sp + ' empiricaldispersion={0}'.format(args.empirical_dispersion_sp)
								else:
									keywords_opt = lot_sp+'/'+ bs_sp+' '+ args.input_for_sp + ' scrf=({0},solvent={1}) empiricaldispersion={2}  '.format(args.solvent_model_sp,args.solvent_name_sp,args.empirical_dispersion_sp)
							else:
								if args.solvent_model_sp == 'gas_phase':
									keywords_opt = lot_sp+'/'+ bs_sp+' '+ args.input_for_sp
								else:
									keywords_opt = lot_sp+'/'+ bs_sp+' '+ args.input_for_sp + ' scrf=({0},solvent={1}) '.format(args.solvent_model_sp,args.solvent_name_sp)
						new_com_file(w_dir,w_dir_initial,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp)

# CHECKS THE FOLDER OF FINAL LOG FILES
def check_for_final_folder(w_dir,log):
	dir_found = False
	while not dir_found:
		temp_dir = w_dir+'/new_gaussian_input_files'
		if os.path.isdir(temp_dir):
			w_dir = temp_dir
		else:
			dir_found =True
	return w_dir

# CALCULATION OF BOLTZMANN FACTORS
def boltz_calculation(val,i,log):
	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--csv', '--boltz', '--output', str(i), val]
	subprocess.call(cmd_boltz)

# CHECKING FOR DUPLICATES
def dup_calculation(val,w_dir, agrs,log):
	# GoodVibes must be installed as a module (through pip or conda)
	cmd_dup = ['python', '-m', 'goodvibes', '--dup', val, '>', 'duplicate_files_checked.txt']
	subprocess.call(cmd_dup)

	#reading the txt files to get the DUPLICATES
	dup_file_list = []
	dupfile = open('duplicate_files_checked.txt',"r")
	duplines = dupfile.readlines()

	for i,_ in enumerate(duplines):
		if duplines[i].find('duplicate') > -1:
			dup_file_list.append(duplines[i].split(' ')[1])

	#move the files to specific directory
	destination = w_dir+'Duplicates/'
	for source in dup_file_list:
		try:
			os.makedirs(destination)
			shutil.move(source, destination)
		except OSError:
			if  os.path.isdir(destination) and not os.path.exists(destination):
				shutil.move(source, destination)
			else:
				raise

# COMBINING FILES FOR DIFFERENT MOLECULES
def combine_files(csv_files, lot, bs, args,log):
	columns = ['Structure', 'E', 'ZPE', 'H', 'T.S', 'T.qh-S', 'G(T)', 'qh-G(T)']
	#final dataframe with only the boltzmann averaged values
	final_file_avg_thermo_data = pd.DataFrame(columns=columns)
	compare_G = pd.DataFrame(columns=['Structure_of_min_conf','min_qh-G(T)','boltz_avg_qh-G(T)'])

	files = []
	#combine all the csv_files

	for f in csv_files:

		log.write(f)

		df = pd.read_csv(f, skiprows = 16)
		# df['Structure']= df['Structure'].astype(str)
		df = df.rename(columns={"   Structure": "Structure"})

		#dropping the ************* line
		df = df.drop(df.index[0])
		df.iloc[-1] = np.nan

		for col in columns:
			if col == 'Structure':
				#identifyin the minmum energy if the conformers
				min_G = df['qh-G(T)'].min()
				#getting the name of the structure of the min G
				idx_name_of_min_conf = df['qh-G(T)'].idxmin() - 1
				name_of_min_conf = df.iloc[idx_name_of_min_conf]['Structure']

			elif col != 'Structure':
				boltz_avg = np.sum(df[col] * df['Boltz'])
				df.at[df.index[-1], col] = boltz_avg
				if col == 'qh-G(T)':
					compare_G = compare_G.append({'Structure_of_min_conf': name_of_min_conf,'min_qh-G(T)': min_G,'boltz_avg_qh-G(T)': boltz_avg}, ignore_index=True)

		final_file_avg_thermo_data = final_file_avg_thermo_data.append({'Structure':name_of_min_conf , 'E': df.iloc[-1]['E'] , 'ZPE': df.iloc[-1]['ZPE'], 'H':df.iloc[-1]['H'] , 'T.S':df.iloc[-1]['T.S'] , 'T.qh-S':df.iloc[-1]['T.qh-S'] , 'G(T)': df.iloc[-1]['G(T)'], 'qh-G(T)':df.iloc[-1]['qh-G(T)'] },ignore_index=True)

		files.append(df)

	final_file_all_data = pd.concat(files, axis=0, ignore_index=True)

	#combined_csv = pd.concat([pd.read_csv(f, skiprows = 14, skipfooter = 1) for f in csv_files ])
	#change directory to write all files in one place
	destination = args.path+'All csv files/'+ str(lot)+ '-'+ str(bs)
	try:
		os.makedirs(destination)
	except OSError:
		if  os.path.isdir(destination):
			pass
		else:
			raise
	os.chdir(destination)

	#export to csv
	final_file_all_data.to_csv( str(lot) + '-' + str(bs) + '_all_molecules_all data.csv', index=False, encoding='utf-8-sig')
	final_file_avg_thermo_data.to_csv( str(lot) + '-' + str(bs) + '_all_molecules_avg_thermo_data.csv', index=False, encoding='utf-8-sig')
	compare_G.to_csv( str(lot) + '-' + str(bs) + '_all_molecules_compare_G(T).csv', index=False, encoding='utf-8-sig')
