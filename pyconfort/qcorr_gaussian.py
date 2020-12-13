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
from pyconfort.qprep_gaussian import input_route_line,check_for_gen_or_genecp,write_genecp,orca_file_gen
from pyconfort.argument_parser import possible_atoms
from pyconfort.filter import exp_rules_output,check_geom_filter
from pyconfort.csearch import com_2_xyz_2_sdf
from pyconfort.nics_conf import update_coord,get_coords_normal

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

def write_header_and_coords(fileout,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,w_dir_initial,log,com_type=None):
	if com_type == 'nics':
		NATOMS,ATOMTYPES,CARTESIANS = update_coord(NATOMS,ATOMTYPES,CARTESIANS,args,log,name,w_dir_initial,'write')
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
def new_com_file(com_type,w_dir_initial,log,new_gaussian_input_files,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_com,lot_com,bs_gcp_com,orca_aux_section):

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

	write_header_and_coords(fileout,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,w_dir_initial,log,com_type)

	# write genecp/gen part
	if genecp == 'genecp' or  genecp == 'gen':
		if com_type == 'sp':
			type_gen = 'sp'
		elif com_type == 'analysis':
			type_gen = 'qcorr'

		write_genecp(type_gen,fileout,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,bs_com,lot_com,bs_gcp_com,args,w_dir_initial,new_gaussian_input_files)

	if args.sp == 'gaussian' and com_type == 'sp':
		# final line for SP
		if args.last_line_for_sp != 'None':
			fileout.write(args.last_line_for_sp)
			fileout.write('\n\n')

	if args.QCORR == 'gaussian' and com_type == 'analysis':
		# final line for analysis
		if args.last_line_for_input != 'None':
			fileout.write(args.last_line_for_input)
			fileout.write('\n\n')

	fileout.close()

	if args.sp == 'orca' and com_type == 'sp':

		read_lines = open(file_name,"r").readlines()

		#create input file
		orca_file_gen(read_lines,file_name.split('.')[0]+'.inp',bs_com,lot_com,genecp,args.aux_atoms_orca_sp,args.aux_basis_set_genecp_atoms_sp,args.aux_fit_genecp_atoms_sp,CHARGE,MULT,orca_aux_section,args,args.set_input_line_sp,args.solvent_model_sp,args.solvent_name_sp,args.cpcm_input_sp,args.orca_scf_iters_sp,args.mdci_orca_sp,args.print_mini_orca_sp)

		# removes the initial com file
		os.remove(file_name)

	#submitting the gaussian file on summit
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
			if (outlines[i].find("Standard orientation") > -1 or outlines[i].find("Input orientation") > -1) and stop_get_details_stand_or !=1 :
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
	ATOMTYPES, CARTESIANS,stand_or = get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)
	return ATOMTYPES, CARTESIANS, NATOMS,stand_or

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

def create_folder_and_com(w_dir,w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE,MULT,orca_aux_section):
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
	if ecp_genecp_atoms and ecp_gen_atoms:
		sys.exit("x  ERROR: Can't use Gen and GenECP at the same time")

	if ERRORTYPE == 'SCFerror':
		input_route += ' scf=qc'
	if genecp == 'genecp' or  genecp == 'gen':
		keywords_opt = lot +'/'+ genecp+' '+ input_route
	else:
		keywords_opt = lot +'/'+ bs +' '+ input_route

	com_type = 'analysis'
	new_com_file(com_type,w_dir_initial,log,new_gaussian_input_files,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs,lot,bs_gcp,orca_aux_section)

def create_folder_move_log_files(w_dir,w_dir_main,round_num,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin,finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,passing_rules,passing_geom,check_geom_qcorr):
	source = w_dir+'/'+file
	if IM_FREQS == 0 and TERMINATION == "normal" and passing_rules and passing_geom:
		destination = w_dir_fin
		moving_files(source, destination)
		finished += 1

	elif IM_FREQS > 0:
		destination = w_dir_main+'/failed/run_'+str(round_num)+'/imag_freq/'
		moving_files(source, destination)
		imag_freq += 1

	elif IM_FREQS == 0 and TERMINATION == "error":
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

	elif not passing_rules:
		destination = w_dir_main+'/failed/run_'+str(round_num)+'/exp_rules_filter/'
		moving_files(source, destination)
		exp_rules_qcorr += 1

	elif not passing_geom:
		destination = w_dir_main+'/failed/run_'+str(round_num)+'/geometry_changed/'
		moving_files(source, destination)
		check_geom_qcorr += 1

	return finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,check_geom_qcorr

# Output file to mol converter
def output_to_mol(file,format,mol_name):
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
	if format == 'log':
		cmd_obabel = ['obabel', '-ilog', os.path.splitext(file)[0]+'.log', '-omol', '-O', os.path.splitext(file)[0]+'.mol']
	subprocess.run(cmd_obabel)
	mol = Chem.MolFromMolFile(mol_name.split('.')[0]+'.mol')

	return mol,ob_compat,rdkit_compat

# DEFINTION OF OUTPUT ANALYSER and NMR FILES CREATOR
def output_analyzer(duplicates,log_files,com_files, w_dir, w_dir_main,lot, bs, bs_gcp, args, w_dir_fin, w_dir_initial, log, ana_data, round_num):

	input_route = input_route_line(args)
	finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,check_geom_qcorr = 0,0,0,0,0,0,0,0

	if round_num == 1:
		#moves the comfiles to respective folder
		for file in com_files:
			source = w_dir+'/'+file
			destination = w_dir_main +'/input_files/run_'+str(round_num)
			moving_files(source, destination)

	for file in log_files:
		# read the file
		log.write(file)
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
			ATOMTYPES, CARTESIANS, stand_or = get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)

		# Get he coordinates for jobs that did not finished or finished with an error
		if TERMINATION != "normal":
			ATOMTYPES, CARTESIANS,NATOMS, stand_or = get_coords_not_normal(outlines, stop_rms, stand_or, dist_rot_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)
		# This part fixes imaginary freqs (if any)
		if IM_FREQS > 0:
			CARTESIANS = fix_imag_freqs(NATOMS, CARTESIANS, args, FREQS, NORMALMODE)

		#close the file
		outfile.close()

		# this part filters off conformers based on user-defined exp_rules
		passing_rules = True
		valid_mol_gen = True
		if len(args.exp_rules) >= 1:
			if TERMINATION == "normal" and IM_FREQS == 0:
				log.write("  ----- Exp_rules filter(s) will be applied to the output file -----\n")
				try:
					mol,ob_compat,rdkit_compat = output_to_mol(file,'log',file)
					print_error_exp_rules=False
					if ob_compat and rdkit_compat:
						passing_rules = exp_rules_output(mol,args,log,file,print_error_exp_rules,ob_compat,rdkit_compat)
					os.remove(file.split('.')[0]+'.mol')
				except AttributeError:
					valid_mol_gen = False
					os.remove(file.split('.')[0]+'.mol')
					log.write("The file could not be converted into a mol object, exp_rules filter(s) will be disabled\n")

		passing_geom = True
		if args.check_geom and passing_rules and valid_mol_gen:
			if TERMINATION == "normal" and IM_FREQS == 0:
				log.write("  ----- Geometrical check will be applied to the output file -----\n")
				# this creates a mol object from the optimized log file
				mol,ob_compat,rdkit_compat = output_to_mol(file,'log',file)
				# this creates a mol object from the input file
				try:
					os.chdir(w_dir_main +'/input_files/run_'+str(round_num))
					com_2_xyz_2_sdf(args,os.path.splitext(file)[0]+'.com')
					mol2,ob_compat,rdkit_compat = output_to_mol(file,'xyz',file)
					passing_geom = check_geom_filter(mol,mol2,args)
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

		# This part places the calculations in different folders depending on the type of termination
		finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,check_geom_qcorr = create_folder_move_log_files(w_dir,w_dir_main,round_num,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin,finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,passing_rules,passing_geom,check_geom_qcorr)

		# check if gen or genecp are active
		# right now, QCORR is only working with Gaussian output files

		ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'analysis','gaussian')

		# create folders and set level of theory in COM files to fix imaginary freqs or not normal terminations
		if IM_FREQS > 0 or TERMINATION != "normal" and not os.path.exists(w_dir_main+'/failed/run_'+str(round_num)+'/error/basis_set_error/'+file):
			create_folder_and_com(w_dir,w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE, MULT, orca_aux_section)

		# adding in the NMR componenet only to the finished files after reading from normally finished log files
		if TERMINATION == "normal" and IM_FREQS == 0 and passing_rules and passing_geom:
			if args.sp == 'gaussian' or args.sp == 'orca' or args.nics:

				#get coordinates
				ATOMTYPES, CARTESIANS,stand_or = get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)

				if args.sp == 'gaussian':
					# creating new folder with new input gaussian files
					single_point_input_files = w_dir_fin+'/../G16-SP_input_files'

				elif args.sp == 'orca':
					# creating new folder with new input gaussian files
					single_point_input_files = w_dir_fin+'/../ORCA-SP_input_files'

				if args.nics:
					nics_input_files = w_dir_fin+'/../G16-NICS_input_files'

				# Options for genecp
				if args.sp == 'gaussian':
					ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','gaussian')

				elif args.sp == 'orca':
					ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','orca')

				# this avoids problems related to genecp
				if genecp == None:
					basis_set_for_genecp = args.basis_set_sp
				elif genecp == 'genecp' or genecp == 'gen':
					basis_set_for_genecp = args.basis_set_genecp_atoms_sp

				# Sets the folder and find the log files to analyze
				for lot_sp,bs_sp,bs_gcp_sp in zip(args.level_of_theory_sp,args.basis_set_sp,basis_set_for_genecp):
					if args.sp == 'gaussian' or args.sp == 'orca':
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
						if args.set_input_line_sp != 'None':
							keywords_opt += ' {0}'.format(args.set_input_line_sp)
						if args.empirical_dispersion_sp != 'None':
							keywords_opt += ' empiricaldispersion={0}'.format(args.empirical_dispersion_sp)
						if args.solvent_model_sp != 'gas_phase':
							keywords_opt += ' scrf=({0},solvent={1})'.format(args.solvent_model_sp,args.solvent_name_sp)

					if args.charge_sp != 'None':
						CHARGE = args.charge_sp

					if args.mult_sp != 'None':
						MULT = args.mult_sp

					if args.sp == 'gaussian' or args.sp == 'orca':
						if not os.path.isdir(single_point_input_files+'/'+dir_name):
							os.makedirs(single_point_input_files+'/'+dir_name)
						os.chdir(single_point_input_files+'/'+dir_name)
						new_com_file('sp',w_dir_initial,log,single_point_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

					if args.nics:
						if not os.path.isdir(nics_input_files+'/'+dir_name):
							os.makedirs(nics_input_files+'/'+dir_name)
						os.chdir(nics_input_files+'/'+dir_name)
						new_com_file('nics',w_dir_initial,log,nics_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

	#write to csv ana_data
	if duplicates=='None':
		ana_data.at[0,'Total files'] = len(log_files)
	else:
		ana_data.at[0,'Total files'] = len(log_files)+int(duplicates) # since duplicates are moved before anything else
	ana_data.at[0,'Normal termination'] = finished
	ana_data.at[0,'Imaginary frequencies'] = imag_freq
	ana_data.at[0,'SCF error'] = scf_error
	ana_data.at[0,'Basis set error'] =  atom_error
	ana_data.at[0,'Other errors'] = other_error
	ana_data.at[0,'Unfinished'] = unfinished
	if args.dup:
		ana_data.at[0,'Duplicates'] = duplicates
	if len(args.exp_rules) >= 1:
		ana_data.at[0,'Exp_rules filter'] = exp_rules_qcorr
	if args.check_geom:
		ana_data.at[0,'Geometry changed'] = check_geom_qcorr

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
	source = 'Goodvibes_output.dat'
	moving_files(source, destination)

	for source in dup_file_list:
		#finding the extension
		for file in val:
			if file.split('.')[0] == source:
				ext=file.split('.')[1]
		source=source+'.'+ext
		moving_files(source, destination)

	return duplicates
