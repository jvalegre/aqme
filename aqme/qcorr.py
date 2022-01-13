######################################################.
#        This file stores all the functions          #
#          used in the LOG file analyzer             #
######################################################.
import os
import sys
import glob
from numpy import disp
import pandas as pd
import json
import cclib
import subprocess

from aqme.utils import periodic_table
from aqme.filter import geom_rules_output

from aqme.utils import (
    move_file,
    get_info_com,
    Logger,
    load_from_yaml,
    check_isomerization,
    read_file,
    output_to_mol,
    cclib_atoms_coords,
)
from aqme.qprep import qprep
from aqme.argument_parser import set_options


def check_for_final_folder(w_dir):
	"""
	# Determines the folder where input files are gonna be generated in QCORR.
	"""
	input_folder = w_dir+'/fixed_QM_input_files'
	folder_count = 0
	
	if os.path.exists(input_folder):
		dir_list = os.listdir(input_folder)
		for folder in dir_list:
			if folder.find('run_') > -1:
				folder_count += 1

	if folder_count == 0:
		return 0
	else:
		num_com_folder = sum([len(d) for r, d, folder in os.walk(w_dir+'/fixed_QM_input_files')])
		w_dir = w_dir+'/fixed_QM_input_files/run_'+str(num_com_folder)
		return folder_count


class qcorr():
	"""
	Class containing all the functions from the QCORR module.

	Parameters
	----------
	qm_files : list 
		Contains the filenames of QM output files to analyze
	round_num : int
		Round of analysis from QCORR
	w_dir_main : str
		Working directory
	mem : str
		Memory used in the calculations
	nprocs : int
		Number of processors used in the calculations
	chk : bool
		Include the chk input line in new input files
	yaml_file : str
		Option to parse the variables using a yaml file (specify the filename)
	qm_input : str
		Keywords line for new input files
	s2_threshold :  float
		Cut off for spin contamination during analysis in \%\ of the expected value 
		(i.e. multiplicity 3 has an the expected <S**2> of 2.0, if s2_threshold = 10,
		the <S**2> value is allowed to be 2.0 +- 0.2). Set s2_threshold = 0 to deactivate this option.
	bs_gen : str
		Basis set used for gen(ECP) atoms	
	bs : str
		Basis set used for non gen(ECP) atoms in gen(ECP) calculations
	gen_atoms : list of str
		Atoms included in the gen(ECP) basis set
	qm_end : str
		Final line in the new input files
	amplitude_ifreq : float
		Amplitude used to displace the imaginary frequencies to fix
	ifreq_cutoff : float
		Cut off for to consider whether a frequency is imaginary (absolute of the specified value is used)
	fullcheck : bool
		Perform an analysis to detect whether the calculations were done homogeneously
		(i.e. same level of theory, solvent, grid size, etc)
	author : str
		Author of the calculations
	program : str
		Program required to create the new input files
	kwargs : argument class
		Specify any arguments from the QCORR module
	"""
	
	def __init__(self, qm_files=[], round_num=0, w_dir_main=os.getcwd(), 
				mem='', nprocs=0, chk=False, yaml_file=None, qm_input='', s2_threshold=10.0,
				bs_gen='', bs='', gen_atoms=[], qm_end='', amplitude_ifreq=0.2,
				ifreq_cutoff=0.0, fullcheck=False, author='', program='gaussian', **kwargs):
		
		if qm_files in ['*.log','*.out']:
			self.qm_files = glob.glob(qm_files.lower())
			for new_file in glob.glob(qm_files.upper()):
				if new_file not in self.qm_files:
					self.qm_files.append(new_file)
		else:
			self.qm_files = qm_files
		self.w_dir_main = w_dir_main
		self.mem = mem
		self.nprocs = nprocs
		self.chk = chk
		self.amplitude_ifreq = amplitude_ifreq
		self.ifreq_cutoff = ifreq_cutoff
		self.qm_input = qm_input
		self.s2_threshold = s2_threshold
		self.fullcheck = fullcheck
		self.bs_gen = bs_gen
		self.bs = bs
		if gen_atoms != []:
			gen_atoms_list = gen_atoms.split(',')
			self.gen_atoms = gen_atoms_list
		else:
			self.gen_atoms = gen_atoms
		self.qm_end = qm_end
		self.author = author
		self.program = program
		
		if 'options' in kwargs:
			self.args = kwargs['options']
		else:
			self.args = set_options(kwargs)
		
		self.args.varfile = yaml_file
		
		if yaml_file is not None:
			self.args, self.log = load_from_yaml(self.args, self.log)
		
		# detects cycle of analysis (0 represents the starting point)
		self.round_num = round_num 

		# start a log file to track the QCORR module
		if not os.path.isdir(self.w_dir_main+'/dat_files/'):
			os.makedirs(self.w_dir_main+'/dat_files/')
		self.log = Logger(self.w_dir_main+'/dat_files/pyCONFORT',f'QCORR-run_{str(self.round_num)}')
		self.log.write("\no  Analyzing output files in {}\n".format(self.w_dir_main))
		
		if len(qm_files) == 0:
			self.log.write('x  There are no output files in this folder.')
			sys.exit('x  There are no output files in this folder.')
		self.qcorr_processing()


	def qcorr_processing(self):
		"""
		General function of the QCORR module that:
		1. Analyzes the QM output files and moves output files with normal termination 
		and no extra imaginary frequencies to the same folder
		2. Generates input files to fix errors and extra imaginary frequencies
		3. Generates input files with new keywords lines from the normally terminated 
		files from point 1 (i.e. single-point energy corrections). Optionally, 
		the analysis from points 1 and 2  might be disabled with the nocheck=True or --nocheck option.
		"""

		file_terms = {'finished': 0, 'extra_imag_freq': 0, 'ts_no_imag_freq': 0, 
			'spin_contaminated': 0, 'duplicate_calc': 0, 'atom_error': 0,
			'scf_error': 0, 'before_E_error': 0, 'not_specified': 0,
			'geom_rules_qcorr': 0, 'check_geom_qcorr': 0}   

		for file in self.qm_files:

			file_name = file.split('.')[0]

			# create a json file with cclib and load a dictionary. This protocol is
			# favored over the traditional ccread since more data will be add to
			# the json file at the end
			command_run_1 = ['ccwrite', 'json', file]
			subprocess.run(command_run_1)
			with open(file_name+'.json') as json_file:
				cclib_data = json.load(json_file)

			# get number of atoms, multiplicity and number of imaginary freqs
			n_atoms = cclib_data['properties']['number of atoms']
			charge = cclib_data['properties']['charge']
			mult = cclib_data['properties']['multiplicity']

			# determine the calculation type (ground or transition state), keyword line,
			# grid size, and more parameters reading part of the output QM file
			outlines, outfile = read_file(self.w_dir_main,file)
			
			keywords_line,calc_type,mem,nprocs,qm_program,author,qm_solv,qm_emp = self.get_init_info(outlines)

			grid_size,s2_operator,level_of_theory = self.grid_s2_info(outlines,qm_program)

			# determine error/unfinished vs normal terminations
			try:
				# freq for normal terminations
				if n_atoms == 1:
					freqs = []
					freq_displacements = []
				else:
					freqs = cclib_data['vibrations']['frequencies']
					freq_displacements = cclib_data['vibrations']['displacement']
				termination = 'normal'
				errortype = 'none'
				
				# spin contamination analysis using user-defined thresholds
				unpaired_e = mult-1
				spin = unpaired_e*0.5
				s2_expected_value = spin*(spin+1)
				spin_diff = abs(float(s2_operator)-s2_expected_value)
				if spin_diff > abs(self.s2_threshold/100)*s2_expected_value:
					errortype = 'spin_contaminated'

			except (AttributeError,KeyError):
				try:
					# if the optimization finished, only a freq job is required
					if cclib_data['optimization']['done']:
						termination = 'other'
						errortype = 'no_freq'
				except (AttributeError,KeyError):
					try:
						# general errors
						atoms_dummy = cclib_data['atoms']['elements']['number']
						termination = 'other'
						errortype = 'not_specified'
					except (AttributeError,KeyError):
						# the program lost the molecular information (i.e. crashed at the start)
						errortype = 'before_E_error'

			# use very short reversed loop to find basis set incompatibilities
			if termination == 'other' and errortype == 'not_specified':
				for i in reversed(range(len(outlines)-15,len(outlines))):
					if outlines[i].find("Error termination") > -1:
						if outlines[i-1].find("Atomic number out of range") > -1 or outlines[i-1].find("basis sets are only available") > -1:
							errortype = "atomicbasiserror"	
							break

			# check for undesired imaginary freqs and data used by GoodVibes
			if termination == 'normal':
				symmno,point_group,roconst,rotemp = self.symm_rot_data(outlines,qm_program)
				
				if errortype != 'spin_contaminated':
					initial_ifreqs = 0
					for freq in freqs:
						if float(freq) < 0 and abs(float(freq)) > abs(self.ifreq_cutoff):
							initial_ifreqs += 1

					# exclude TS imag frequency
					if calc_type == 'transition_state':
						initial_ifreqs -= 1

					# gives new coordinates by displacing the normal mode(s) of the negative freq(s)
					if initial_ifreqs > 0:
						errortype = 'extra_imag_freq'
						atom_types,cartesians = cclib_atoms_coords(cclib_data)
						cartesians = self.fix_imag_freqs(n_atoms, cartesians, freqs, freq_displacements, calc_type)

						qcorr_calcs = qprep(destination=f'{os.getcwd()}/qm_input/run_{self.round_num}',
											molecule=file_name, charge=charge, mult=mult,
											program=self.program, atom_types=atom_types,
											cartesians=cartesians, qm_input=keywords_line,
											mem=mem, nprocs=nprocs, chk=self.chk, qm_end=self.qm_end,
											bs_gen=self.bs_gen, bs=self.bs, gen_atoms=self.gen_atoms)

					elif initial_ifreqs < 0:
						errortype = 'ts_no_imag_freq'

			# for calcs with finished OPT but no freqs
			elif termination != 'normal' and errortype not in ['ts_no_imag_freq','atomicbasiserror','before_E_error']:
				if errortype == 'no_freq':
					# adjust the keywords so only FREQ is calculated
					new_keywords_line = ''
					for keyword in keywords_line.split():
						if keyword.lower().startswith('opt'):
							keyword = ''
						else:
							new_keywords_line += keyword
							new_keywords_line += ' '
					keywords_line = new_keywords_line

					atom_types,cartesians = cclib_atoms_coords(cclib_data)

				else:
					# helps to fix SCF convergence errors
					if errortype == 'SCFerror':
						if keywords_line.find(' scf=qc') > -1:
							pass
						else:
							keywords_line += ' scf=qc'
					
					if errortype in ['not_specified','SCFerror']:
						RMS_forces = [row[1] for row in cclib_data['optimization']['geometric values']]
						min_RMS = RMS_forces.index(min(RMS_forces))
						
						atom_types,cartesians = [],[]
						per_tab = periodic_table()
						count_RMS = -1
						for i in range(60,len(outlines)):
							if outlines[i].find('Standard orientation:') > -1:
								count_RMS += 1
							if count_RMS == min_RMS:
								for j in range(i+5,i+5+n_atoms):
									massno = int(outlines[j].split()[1])
									if massno < len(per_tab):
										atom_symbol = per_tab[massno]
									else:
										atom_symbol = "XX"
									atom_types.append(atom_symbol)
									cartesians.append([float(outlines[j].split()[3]), float(outlines[j].split()[4]), float(outlines[j].split()[5])])
								break
					
				qcorr_calcs = qprep(destination=f'{os.getcwd()}/qm_input/run_{self.round_num}',
							molecule=file_name, charge=charge, mult=mult,
							program=self.program, atom_types=atom_types,
							cartesians=cartesians, qm_input=keywords_line,
							mem=mem, nprocs=nprocs, chk=self.chk, qm_end=self.qm_end,
							bs_gen=self.bs_gen, bs=self.bs, gen_atoms=self.gen_atoms)

			print(f'{file}: Termination = {termination}, Error type = {errortype}')
			self.log.write(f'{file}: Termination = {termination}, Error type = {errortype}')

			# close the file
			outfile.close()

			# This part places the calculations in different folders depending on the type of termination
			file_terms,destination = self.organize_outputs(file,termination,errortype,file_terms)

			# write information about the QCORR analysis in a csv
			self.write_qcorr_csv(file_terms)

			# append new data from AQME into the json file from cclib and move json files
			# into the corresponding folders
			aqme_data = {'AQME data': {}}
			aqme_data['AQME data']['keywords line'] = keywords_line
			aqme_data['AQME data']['calculation type'] = calc_type
			aqme_data['AQME data']['memory'] = mem
			aqme_data['AQME data']['number of procs'] = nprocs
			aqme_data['AQME data']['QM program'] = qm_program
			aqme_data['AQME data']['author'] = author
			aqme_data['AQME data']['solvation'] = qm_solv
			aqme_data['AQME data']['dispersion'] = qm_emp
			aqme_data['AQME data']['grid'] = grid_size
			aqme_data['AQME data']['<S**2>'] = s2_operator
			aqme_data['AQME data']['level of theory'] = level_of_theory
			if termination == 'normal':
				aqme_data['AQME data']['symmetry number'] = symmno
				aqme_data['AQME data']['point group'] = point_group
				aqme_data['AQME data']['rotational constant'] = roconst
				aqme_data['AQME data']['rotational temperature'] = rotemp
			
			# update and move the json files to their corresponding destination
			with open(file_name+'.json', "r+") as json_file:
				data = json.load(json_file)
				data.update(aqme_data)
				json_file.seek(0)
				json.dump(data, json_file, indent=4, separators=(", ", ": "))

			
			destination_json = destination+'json_files/'
			if not os.path.isdir(destination_json):
				os.makedirs(destination_json)
			move_file(file_name+'.json', self.w_dir_main, destination_json)

			
	# make a function for self.fullcheck que te chequee todo:
	# *funcion check dentro de qcorr aparte*
	# 	mismo grid size
	# 	no spin contamination
	# 	no isomeriz
	# 	no dups
	# 	same solvent
	# 	same disp
	# 	same leveloftheory
	# 	in the additional check, take the input file in G16 and look for SCRF y emp/empirical, luego split con (, ), ,, = y ordena la lista, tiene que ser igual en todos los casos

	# include geom filters and isomeriz filters too
	# 			if self.args.isom != None:
	# 				isomerized = False
	# 				init_csv = pd.DataFrame()
	# 				folder_dir_isom = ''
	# 				self.log.write("  ----- Geometrical check will be applied to the output file -----\n")
	# 				try:
	# 					atoms_com, coords_com, atoms_and_coords = [],[],[]
	# 					if self.round_num != 0:
	# 						folder_dir_isom = self.w_dir_main +'/fixed_QM_input_files/run_'+str(self.round_num)
	# 					else:
	# 						pass
	# 					if self.args.isom == 'com':
	# 						atoms_and_coords,_ = get_info_com(folder_dir_isom+file.split('.')[0]+'.com')
	# 					elif self.args.isom == 'gjf':
	# 						atoms_and_coords,_ = get_info_com(folder_dir_isom+file.split('.')[0]+'.gjf')
	# 					elif self.args.isom.split('.')[1] == 'csv':
	# 						init_csv = pd.read_csv(self.args.isom)

	# 					for line in atoms_and_coords:
	# 						atoms_com.append(line.split()[0])
	# 						coords_com.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
						
	# 					isomerized = check_isomerization(coords_com, cartesians, atoms_com, atom_types, self.args.vdwfrac, self.args.covfrac, init_csv, file)
					
	# 				except FileNotFoundError:
	# 					self.log.write("x  No com file were found for "+file+", the check_geom test will be disabled for this calculation")
					
	# 				if isomerized:
	# 					errortype = 'isomerization'

	# 			if len(self.args.geom_rules) >= 1:
	# 				passing_rules = True
	# 				valid_mol_gen = True
	# 				self.log.write("  ----- geom_rules filter(s) will be applied to the output file -----\n")
	# 				try:
	# 					format_file = file.split('.')[1]
	# 					mol,ob_compat,rdkit_compat = output_to_mol(file,format_file,self.log)
	# 					print_error_geom_rules=False
	# 					if ob_compat and rdkit_compat:
	# 						passing_rules = geom_rules_output(mol,self.args,self.log,file,print_error_geom_rules)
	# 						if not passing_rules:
	# 							errortype = 'fail_geom_rules'
	# 					os.remove(file.split('.')[0]+'.mol')
	# 				except AttributeError:
	# 					valid_mol_gen = False
	# 					os.remove(file.split('.')[0]+'.mol')
	# 					self.log.write("The file could not be converted into a mol object, geom_rules filter(s) will be disabled\n")



	def get_init_info(self,outlines):
		"""
		Retrieves information from QM files that is not included in the cclib analysis. 


		Parameters
		----------
		outlines : list of str
			Lines of the QM output files

		Returns
		-------
		keywords_line : string
			Original keyword line (input) from the output QM file
		calc_type : str
			Type of the QM calculation (ground_state or transition_state)
		mem : str
			Memory used in the calculations
		nprocs : int
			Number of processors used in the calculations
		qm_program : str
			Program and version used in the calculation
		author : str
			Author of the calculations
		qm_solv : str
			Solvation model used in the calculations
		qm_emp : str
			Empirical model (if any) used in the calculations
		"""    

		end_keywords,keywords_line,author = False,'',''
		mem,nprocs,qm_program,skip_lines = '',0,'',0
		qm_solv,qm_emp = 'gas_phase','no_dispersion'
		orca_solv_1, orca_solv_2, orca_solv_3 = 'gas phase','',''
		
		for i in range(0,len(outlines)):
			if outlines[i].find(' Gaussian ') > -1 and outlines[i].find('Revision') > -1:
				qm_program = outlines[i][1:-1]

			elif outlines[i].find('* O   R   C   A *') > -1:
				qm_program = 'ORCA'

			if qm_program.lower().find('gaussian') > -1:
				# retrieve the multiple parameters of the calculation
				if outlines[i].startswith(' #'):
					keywords_line += outlines[i].rstrip("\n")
					for j in range(i+1,i+100):
						if outlines[j].find('----------') > -1:
							end_keywords = True
							break
						if not end_keywords:
							keywords_line += outlines[j].rstrip("\n")[1:]
					keywords_line = keywords_line[3:]

				elif outlines[i].find('%mem') > -1:
					mem = outlines[i].strip().split('=')[-1]

				elif outlines[i].find('%nprocs') > -1:
					nprocs = outlines[i].strip().split('=')[-1]
				
				if keywords_line != '':
					break
			
			elif qm_program.lower().find('orca') > -1:
				# get program version
				if outlines[i].find('Program Version') > -1:
					qm_program = "ORCA, version " + outlines[i].split()[2]
				
				# retrieve the input lines, mem and nprocs
				if outlines[i].startswith('NAME ='):
					for j in range(i+1,i+100):
						if outlines[j].find('xyz') > -1 and outlines[j].find('*') > -1:
							end_keywords = True
							break

						if not end_keywords:
							if outlines[j].lower().find('%maxcore') > -1:
								mem = outlines[j].split()[-1]+'MB'

							# multiple ways to describe %pal
							elif outlines[j].lower().find('%pal') > -1 and outlines[j].lower().find('nprocs') and outlines[j].lower().find('end'):
								nprocs = int(outlines[j].split()[-2])
							elif outlines[j].lower().find('%pal') > -1 and outlines[j].lower().find('nprocs') and outlines[j+1].lower().find('end'):
								nprocs = int(outlines[j].split()[-1])
								skip_lines = 1
							elif outlines[j].lower().find('%pal') > -1 and outlines[j+1].lower().find('nprocs') and outlines[j+1].lower().find('end'):
								nprocs = int(outlines[j+1].split()[-2])
								skip_lines = 1
							elif outlines[j].lower().find('%pal') > -1 and outlines[j+1].lower().find('nprocs') > -1 and outlines[j+2].lower().find('end'):
								nprocs = int(outlines[j+1].split()[-1])
								skip_lines = 2
							elif skip_lines != 0:
								skip_lines -= 1

							else:
								keywords_line += outlines[j][6:]

				# ORCA parsing for solvation model
				elif outlines[i].find('CPCM SOLVATION MODEL') > -1:
					orca_solv_1 = 'CPCM, '
				elif outlines[i].find('SMD-CDS solvent descriptors:') > -1:
					orca_solv_2 = 'SMD, '
				elif outlines[i].find('Solvent:       ') > -1:
					orca_solv_3 = outlines[i].strip().split()[-1].lower()

				if outlines[i].find('SCF ITERATIONS') > -1:
					break
		
		if self.author != '':
			author = self.author

		# user-defined keywords line, mem and nprics overwrites previously used parameters
		if self.qm_input != '':
			keywords_line = self.qm_input

		if self.mem != '':
			mem = self.mem
		elif mem == '':
			mem = '8GB'

		if self.nprocs != 0:
			nprocs = self.nprocs
		elif nprocs == 0:
			nprocs = 4

		# find if the calculation was for a ground or transition state and get individual keywords
		calc_type = 'ground_state'
		calcfc_found, ts_found = False, False

		if qm_program.lower().find('gaussian') > -1:
			for keyword in keywords_line.split():
				if keyword.lower().find('opt') > -1:
					if keyword.lower().find('calcfc') > -1:
						calcfc_found = True
					if keyword.lower().find('ts') > -1:
						ts_found = True
				elif keyword.lower().startswith('scrf'):
					qm_solv = keyword
				elif keyword.lower().startswith('emp'):
					qm_emp = keyword
			if calcfc_found and ts_found:
				calc_type = 'transition_state'

		
		elif qm_program.lower().find('orca') > -1:
			for keyword in keywords_line.split():
				if keyword.lower() in ['optts','scants']:
					ts_found = True
					if keyword.lower() == 'scants':
						keyword = 'OptTS'
			if ts_found:
				calc_type = 'transition_state'
			if orca_solv_1 == 'CPCM, ' and orca_solv_2 == 'SMD, ':
				qm_solv = orca_solv_2 + orca_solv_3
			else:		
				qm_solv = orca_solv_1 + orca_solv_2 + orca_solv_3
	
		return keywords_line,calc_type,mem,nprocs,qm_program,author,qm_solv,qm_emp


	def symm_rot_data(self,outlines,qm_program):
		"""
		Retrieves information from QM files regarding symmetry and rotational parameters used by GoodVibes.


		Parameters
		----------
		outlines : list of str
			Lines of the QM output files
		qm_program : str
			Program and version used in the calculation

		Returns
		-------
		symmno : int
			Rotational symmetry number
		point_group : int
			Symmetry point group
		roconst : list of float
			Rotational constants in GHz
		rotemp : list of float
			Rotational temperatures
		"""  
		
		symmno, point_group, roconst, rotemp = '','',[],[]
		
		for i in reversed(range(0,len(outlines))):
			if qm_program.lower().find('gaussian') > -1:
				if outlines[i].strip().startswith('Rotational symmetry number'):
					symmno = int((outlines[i].strip().split()[3]).split(".")[0])

				elif outlines[i].strip().startswith('Full point group'):
					point_group = outlines[i].strip().split()[3]
					break

				elif outlines[i].strip().startswith('Rotational constants (GHZ):'):
					try:
						roconst = [float(outlines[i].strip().replace(':', ' ').split()[3]),
										float(outlines[i].strip().replace(':', ' ').split()[4]),
										float(outlines[i].strip().replace(':', ' ').split()[5])]
					except ValueError:
						if outlines[i].strip().find('********'):
							roconst = [float(outlines[i].strip().replace(':', ' ').split()[4]),
											float(outlines[i].strip().replace(':', ' ').split()[5])]
				
				elif outlines[i].strip().startswith('Rotational temperature '):
					rotemp = [float(outlines[i].strip().split()[3])]
				elif outlines[i].strip().startswith('Rotational temperatures'):
					try:
						rotemp = [float(outlines[i].strip().split()[3]), float(outlines[i].strip().split()[4]),
								  float(outlines[i].strip().split()[5])]
					except ValueError:
						if outlines[i].strip().find('********'):
							rotemp = [float(outlines[i].strip().split()[4]), float(outlines[i].strip().split()[5])]
	
		return symmno,point_group,roconst,rotemp


	def grid_s2_info(self,outlines,qm_program):
		"""
		Checks multiple parameters from QM files that are not included in the cclib analysis. 


		Parameters
		----------
		outlines : list of str
			Lines of the QM output files
		qm_program : str
			Program and version used in the calculation

		Returns
		-------
		grid_size : str
			Grid size used in the calculations
		s2_operator : float
			<S**2> value from open-shell QM calculations
		level_of_theory : str
			Functional and basis set used in the QM calculations
		"""  

		grid_size,s2_operator,level_of_theory = '',0,''
		grid_lookup = {1: 'sg1', 2: 'coarse', 4: 'fine', 5: 'ultrafine', 7: 'superfine'}
		found_theory,level,bs = False,'',''
		line_options = ['\\Freq\\','\\SP\\']

		for i in reversed(range(0,len(outlines)-15)):
			# get grid size
			if qm_program.lower().find('gaussian') > -1:
				if outlines[i].strip().startswith('ExpMin='):
					IRadAn = int(outlines[i].strip().split()[-3])
					grid_size = grid_lookup[IRadAn]
				elif outlines[i].find('S**2 before annihilation') > -1:
					s2_operator = float(outlines[i].strip().split()[-1])

				# get functional and basis set
				if not found_theory:
					if outlines[i].strip().find('External calculation') > -1:
						level, bs = 'ext', 'ext'
						level_of_theory = '/'.join([level, bs])
						found_theory = True
					for option in line_options:
						if option in (outlines[i]+outlines[i+1]).strip():
							level, bs = ((outlines[i]+outlines[i+1]).strip().split("\\")[4:6])
							found_theory = True
					# Remove the restricted R or unrestricted U label
					if level != '':
						if level[0] in ('R', 'U'):
							level = level[1:]
						level_of_theory = '/'.join([level, bs])

			if level_of_theory != '':
				found_theory = True
			
			elif qm_program.lower().find('orca') > -1:
				if outlines[i].lower().find('grid') > -1 and outlines[i].startswith('|  '):
					for keyword in outlines[i].lower().split():
						if keyword.find('grid') > -1:
							if grid_size != '':
								grid_size += ' '
							grid_size += keyword
				elif outlines[i].find('xyz') > -1 and outlines[i].find('*') > -1:
					if grid_size == '':
						grid_size = 'default multigrid'

			if grid_size != '':
				break

		return grid_size,s2_operator,level_of_theory


	def fix_imag_freqs(self, n_atoms, cartesians, freqs, freq_displacements, calc_type):
		"""
		Fixes undersired (extra) imaginary frequencies from QM calculations.
		This function multiplies the imaginary normal mode vectors by the selected amplitude 
		(0.2 is the default amplitude in the pyQRC script from GitHub, user: bobbypaton).
		By default, all the extra imaginary modes are used (i.e. in calculations with three
		extra imaginary frequencies, all the three modes will be used to displace the atoms).
		This can be tuned with the --ifreq_cutoff option (i.e. only use freqs lower than -50 cm-1).

		Parameters
		----------
		n_atoms : int
			Number of atoms in the calculation
		cartesians : list of lists
			List of lists containing the molecular coordinates as floats
		freqs : list of float
			List containing the frequencies as floats
		freq_displacements : list of matrixes
			Contains the normal modes for each frequencies (including all atoms)
		calc_type : str
			Type of the QM calculation (ground_state or transition_state)

		Returns
		-------
		cartesians : list of lists
			New set of cartesian coordinates generated after displacing the original
			coordinates along the normal modes of the corresponding imaginary frequencies
		"""

		shift = []

		# could get rid of atomic units here, if zpe_rat definition is changed
		for mode,_ in enumerate(freqs):
			# moves along all imaginary freqs (ignoring the TS imag freq, assumed to be the most negative freq)
			if mode == 0 and calc_type == 'transition_state':
				shift.append(0.0)
			else:
				if freqs[mode] < 0.0:
					shift.append(self.amplitude_ifreq)
				else:
					shift.append(0.0)

			# The starting geometry is displaced along each normal mode according to the random shift
			for atom in range(0,n_atoms):
				for coord in range(0,3):
					cartesians[atom][coord] = cartesians[atom][coord] + freq_displacements[mode][atom][coord] * shift[mode]

		return cartesians


	def organize_outputs(self,file,termination,errortype,file_terms):
		"""
		1. Moves the QM output files to their corresponding folders after the analysis. 
		2. Keeps track of the number of calculations with the different types 
		of terminations and error types

		Parameters
		----------
		file : str
			Output file
		termination : string
			Type of termination of the QM output file (i.e. normal, error, unfinished)
		errortype : string
			Type of error type of the QM output file (i.e. None, not_specified, extra_imag_freq, etc)
		file_terms : dictionary
			Keeps track of the number of calculations for each termination and error type

		Returns
		-------
		destination : str
			Path to store the file
		"""
		
		if errortype == 'none' and termination == "normal":
			destination = self.w_dir_main+'/successful_QM_outputs/'
			file_terms['finished'] += 1

		elif errortype == 'extra_imag_freq':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/extra_imag_freq/'
			file_terms['extra_imag_freq'] += 1

		elif errortype == 'ts_no_imag_freq':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/ts_no_imag_freq/'
			file_terms['ts_no_imag_freq'] += 1

		elif errortype == 'spin_contaminated':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/spin_contaminated/'
			file_terms['spin_contaminated'] += 1

		elif errortype == 'duplicate_calc':
			destination = self.w_dir_main+'/duplicates/run_'+str(self.round_num)+'/'
			file_terms['duplicate_calc'] += 1

		elif errortype == "atomicbasiserror":
			destination = self.w_dir_main +'/failed/run_'+str(self.round_num)+'/error/basis_set_error/'
			file_terms['atom_error'] += 1
		elif errortype == "SCFerror":
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/error/scf_error/'
			file_terms['scf_error'] += 1
		elif errortype == "before_E_calculation":
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/error/before_E_calculation/'
			file_terms['before_E_error'] += 1

		elif errortype == 'fail_geom_rules':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/geom_rules_filter/'
			file_terms['geom_rules_qcorr'] += 1

		elif errortype == 'isomerization':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/isomerization/'
			file_terms['check_geom_qcorr'] += 1

		else:
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/error/not_specified_error/'
			file_terms['not_specified'] += 1
		
		if not os.path.isdir(destination):
			os.makedirs(destination)
		move_file(file, self.w_dir_main, destination)

		return file_terms,destination

					
	def write_qcorr_csv(self,file_terms):
		"""
		Write information about the QCORR analysis in a csv
		"""
		ana_data = pd.DataFrame()
		ana_data.at[0,'Total files'] = len(self.qm_files)
		ana_data.at[0,'Normal termination'] = file_terms['finished']
		ana_data.at[0,'Extra imag. freq.'] = file_terms['extra_imag_freq']
		ana_data.at[0,'TS with no imag. freq.'] = file_terms['ts_no_imag_freq']
		ana_data.at[0,'SCF error'] = file_terms['scf_error']
		ana_data.at[0,'Error before SCF'] = file_terms['before_E_error']
		ana_data.at[0,'Basis set error'] =  file_terms['atom_error']
		ana_data.at[0,'Other errors'] = file_terms['not_specified']
		if self.args.s2_threshold > 0.0:
			ana_data.at[0,'Spin contamination'] = file_terms['spin_contaminated']
		if self.args.dup:
			ana_data.at[0,'Duplicates'] = file_terms['duplicate_calc']
		if len(self.args.geom_rules) >= 1:
			ana_data.at[0,'geom_rules filter'] = file_terms['geom_rules_qcorr']
		if self.args.isom != None:
			ana_data.at[0,'Isomerization'] = file_terms['check_geom_qcorr']

		if not os.path.isdir(self.w_dir_main+'/csv_files/'):
			os.makedirs(self.w_dir_main+'/csv_files/')
		ana_data.to_csv(self.w_dir_main+'/csv_files/Analysis-Data-QCORR-run_'+str(self.round_num)+'.csv',index=False)


def json2input(json_files=[], source=os.getcwd(), destination=os.getcwd(), suffix='', charge=None, mult=None,
				mem='8GB', nprocs=4, chk=False, yaml_file=None, qm_input='', bs_gen='', 
				bs='', gen_atoms=[], qm_end='', program='gaussian'):
	'''
	Reads a json file and use QPREP to generate input files.

	Parameters
	----------
	qm_files : list 
		Contains the filenames of QM output files to analyze
	source : str
		Folder with the json files to process
	destination : str
		Destination to create the new input files
	suffix : str
		Suffix for the new input files
	charge : int
		Charge of the calculations used in the following input files
	mult : int
		Multiplicity of the calculations used in the following input files
	mem : str
		Memory used in the calculations
	nprocs : int
		Number of processors used in the calculations
	chk : bool
		Include the chk input line in new input files
	yaml_file : str
		Option to parse the variables using a yaml file (specify the filename)
	qm_input : str
		Keywords line for new input files
	bs_gen : str
		Basis set used for gen(ECP) atoms	
	bs : str
		Basis set used for non gen(ECP) atoms in gen(ECP) calculations
	gen_atoms : list of str
		Atoms included in the gen(ECP) basis set
	qm_end : str
		Final line in the new input files
	program : str
		Program required to create the new input files
	kwargs : argument class
		Specify any arguments from the QCORR module
	'''
	
	w_dir_main=os.getcwd()

	os.chdir(source)

	if json_files == '*.json':
		json_files = glob.glob('*.json')
	
	for file in json_files:
		file_name = file.split('.')[0]
		with open(file) as json_file:
			cclib_data = json.load(json_file)
		try:
			atom_types,cartesians = cclib_atoms_coords(cclib_data)
		except (AttributeError,KeyError):
			print('x  The json files do not contain coordinates and/or atom type information')
		
		# if no charge and multiplicity are specified, they are read from the json file
		if charge == None:
			charge = cclib_data['properties']['charge']
		if mult == None:
			mult = cclib_data['properties']['multiplicity']
		
		if charge == None:
			print('x  No charge was specified in the json file or function input (i.e. json2input(charge=0) )')
		elif mult == None:
			print('x  No multiplicity was specified in the json file or function input (i.e. json2input(mult=1) )')
		
		json_calcs = qprep(destination=destination,
					molecule=file_name, charge=charge, mult=mult,
					program=program, atom_types=atom_types,
					cartesians=cartesians, qm_input=qm_input,
					mem=mem, nprocs=nprocs, chk=chk, qm_end=qm_end,
					bs_gen=bs_gen, bs=bs, gen_atoms=gen_atoms, suffix=suffix)

	os.chdir(w_dir_main)
	
