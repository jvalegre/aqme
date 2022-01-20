######################################################.
#        This file stores all the functions          #
#          used in the LOG file analyzer             #
######################################################.
import os
import sys
import glob
import pandas as pd
import json
import cclib
import subprocess
from pathlib import Path

from aqme.utils import periodic_table,get_info_input
from aqme.filter import geom_rules_output

from aqme.utils import (
	move_file,
    Logger,
    load_from_yaml,
    check_isomerization,
    read_file,
    cclib_atoms_coords,
)

from aqme.qprep import qprep
from aqme.argument_parser import set_options


class qcorr():
	"""
	Class containing all the functions from the QCORR module.

	Parameters
	----------
	qm_files : list 
		Filenames of QM output files to analyze
	w_dir_main : str
		Working directory
	dup_threshold : float
		Energy (in hartree) used as the energy difference in E, H and G to detect duplicates
	mem : str
		Memory for the QM calculations (i) Gaussian: total memory; (ii) ORCA: memory per processor
	nprocs : int
		Number of processors used in the QM calculations
	chk : bool
		Include the chk input line in new input files for Gaussian calculations
	qm_input : str
		Keywords line for new input files
	s2_threshold :  float
		Cut off for spin contamination during analysis in \%\ of the expected value 
		(i.e. multiplicity 3 has an the expected <S**2> of 2.0, if s2_threshold = 10,
		the <S**2> value is allowed to be 2.0 +- 0.2). Set s2_threshold = 0 to deactivate this option.
	isom : str
		Check for isomerization from the initial input file to the resulting output files. 
		It requires the extension of the initial input files (i.e. isom='com') and the folder of
		the input files must be added in the isom_inputs option
	isom_inputs : str
		Folder containing the initial input files to check for isomerization
	vdwfrac : float
		Fraction of the summed VDW radii that constitutes a bond between two atoms in the isomerization filter
	covfrac : float
		Fraction of the summed covalent radii that constitutes a bond between two atoms in the isomerization filter
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
	freq_conv : str
		If a string is defined, it will remove calculations that converged during optimization but did not
		convergence in the subsequent frequency calculation. Options: opt sections as strings i.e. 
		(opt=(calcfc,maxstep=5)). If readfc is specified in the string, the chk option must be included as well.
		Turn this option off by using freq_conv=False.
	ifreq_cutoff : float
		Cut off for to consider whether a frequency is imaginary (absolute of the specified value is used)
	fullcheck : bool
		Perform an analysis to detect whether the calculations were done homogeneously
		(i.e. same level of theory, solvent, grid size, etc)
	author : str
		Author of the calculations
	program : str
		Program required to create the new input files
	varfile : str
		Option to parse the variables using a yaml file (specify the filename)
	kwargs : argument class
		Specify any arguments from the QCORR module
	"""
	
	def __init__(self, qm_files=[], w_dir_main=os.getcwd(), dup_threshold=0.0001,
				mem='4GB', nprocs=2, chk=False, qm_input='', s2_threshold=10.0, 
				isom=False, isom_inputs=os.getcwd(), vdwfrac=0.50, covfrac=1.10, 
				bs_gen='', bs='', gen_atoms=[], qm_end='', amplitude_ifreq=0.2, freq_conv='opt=(calcfc,maxstep=5)',
				ifreq_cutoff=0.0, fullcheck=True, author='', program='gaussian', varfile=None, **kwargs):
		
		self.initial_dir = Path(os.getcwd())
		self.w_dir_main = Path(w_dir_main)
		self.dup_threshold = dup_threshold
		self.mem = mem
		self.nprocs = nprocs
		self.chk = chk
		self.amplitude_ifreq = amplitude_ifreq
		self.freq_conv = freq_conv
		self.ifreq_cutoff = ifreq_cutoff
		self.qm_input = qm_input
		self.s2_threshold = s2_threshold
		self.fullcheck = fullcheck
		self.bs_gen = bs_gen
		self.bs = bs
		self.gen_atoms = gen_atoms
		self.qm_end = qm_end
		self.author = author
		self.program = program
		self.isom = isom
		self.isom_inputs = Path(isom_inputs)
		self.vdwfrac = vdwfrac
		self.covfrac = covfrac
		
		if "options" in kwargs:
			self.args = kwargs["options"]
		else:
			self.args = set_options(kwargs)

		self.args.varfile = varfile

		if varfile is not None:
			self.args, self.log = load_from_yaml(self.args, self.log)
			self.w_dir_main = Path(self.args.w_dir_main)
			self.qm_files = self.args.qm_files
			self.dup_threshold = self.args.dup_threshold
			self.mem = self.args.mem
			self.nprocs = self.args.nprocs
			self.chk = self.args.chk
			self.amplitude_ifreq = self.args.amplitude_ifreq
			self.freq_conv = self.args.freq_conv
			self.ifreq_cutoff = self.args.ifreq_cutoff
			self.qm_input = self.args.qm_input
			self.s2_threshold = self.args.s2_threshold
			self.fullcheck = self.args.fullcheck
			self.gen_atoms = self.args.gen_atoms
			self.bs_gen = self.args.bs_gen
			self.bs = self.args.bs
			self.qm_end = self.args.qm_end
			self.program = self.args.program
			self.author = self.args.author
			self.isom = self.args.isom
			self.isom_inputs = Path(self.args.isom_inputs)
			self.vdwfrac = self.args.vdwfrac
			self.covfrac = self.args.covfrac

		# go to working folder and detect QM output files
		os.chdir(self.w_dir_main)
		if isinstance(qm_files, list): 
			self.qm_files = qm_files
		else:
			self.qm_files = glob.glob(qm_files)

		# detects cycle of analysis (0 represents the starting point)
		self.round_num = check_run(w_dir_main)

		# start a log file to track the QCORR module
		error_setup = False
		try:
			self.log = Logger(self.w_dir_main / 'QCORR-run',f'{str(self.round_num)}')
		except FileNotFoundError:
			print('x  The PATH specified as input in the w_dir_main option might be invalid!')
			error_setup = True

		if len(self.qm_files) == 0 and not error_setup:
			print(f'x  There are no output files in {self.w_dir_main}.')
			self.log.write(f'x  There are no output files in {self.w_dir_main}.')
			self.log.finalize()
			move_file(self.w_dir_main.joinpath('dat_files/'), self.w_dir_main,f'QCORR-run_{str(self.round_num)}.dat')
			error_setup = True

		if error_setup:
			# this is added to avoid path problems in jupyter notebooks
			os.chdir(self.initial_dir)
			sys.exit()

		self.log.write(f"o  Analyzing output files in {self.w_dir_main}\n")
		self.qcorr_processing()

		self.log.finalize()
		move_file(self.w_dir_main.joinpath('dat_files/'), self.w_dir_main,f'QCORR-run_{str(self.round_num)}.dat')
		
		# this is added to avoid path problems in jupyter notebooks
		os.chdir(self.initial_dir)

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

		file_terms = {'finished': 0, 'sp_calcs' : 0, 'extra_imag_freq': 0, 'ts_no_imag_freq': 0, 
			'freq_no_conv': 0, 'spin_contaminated': 0, 'duplicate_calc': 0, 'atom_error': 0,
			'scf_error': 0, 'before_E_error': 0, 'not_specified': 0,
			'geom_rules_qcorr': 0, 'isomerized': 0}   
		
		E_dup_list, H_dup_list, G_dup_list = [],[],[]

		for file in self.qm_files:

			file_name = file.split('.')[0]

			# create a json file with cclib and load a dictionary. This protocol is
			# favored over the traditional ccread since more data will be added to
			# the json file at the end
			command_run_1 = ['ccwrite', 'json', file]
			subprocess.run(command_run_1)
			cclib_incomp = False
			try:
				with open(file_name+'.json') as json_file:
					cclib_data = json.load(json_file)
			except FileNotFoundError:
				cclib_incomp = True
			if cclib_incomp:
				print(f'x  Potential cclib compatibility problem with file {file}')
				self.log.write(f'x  Potential cclib compatibility problem with file {file}')
				self.log.finalize()
				move_file(self.w_dir_main.joinpath('dat_files/'), self.w_dir_main,f'QCORR-run_{str(self.round_num)}.dat')

				# this is added to avoid path problems in jupyter notebooks
				os.chdir(self.initial_dir)
				sys.exit()

			# get number of atoms, multiplicity and number of imaginary freqs
			n_atoms = cclib_data['properties']['number of atoms']
			charge = cclib_data['properties']['charge']
			mult = cclib_data['properties']['multiplicity']

			# determine the calculation type (ground or transition state), keyword line,
			# grid size, and more parameters reading part of the output QM file
			outlines, outfile = read_file(self.w_dir_main,file)
			
			keywords_line,calc_type,mem,nprocs,qm_program,author,qm_solv,qm_emp = self.get_init_info(outlines)

			grid_size,s2_operator,s2_preannhi,level_of_theory,creation_date = self.grid_s2_info(outlines,qm_program)

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
				# this first part accounts for singlet diradicals (threshold is 10% of the spin before annihilation)
				if unpaired_e == 0 and s2_operator != 0:
					if s2_operator > abs(self.s2_threshold/100)*s2_preannhi:
						errortype = 'spin_contaminated'
				else:
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
						# single-point calcs (normal terminations with no opt)
						cclib_data['optimization']['geometric values']
						termination = 'other'
						errortype = 'not_specified'
					except (AttributeError,KeyError):
						try:
							# general errors
							cclib_data['atoms']['elements']['number']
							termination = 'normal'
							errortype = 'sp_calc'
							calc_type = 'SP calculation'
						except (AttributeError,KeyError):
							# the program lost the molecular information (i.e. crashed at the start)
							termination = 'other'
							errortype = 'before_E_error'
							
			# use very short reversed loop to find basis set incompatibilities
			if errortype in ['not_specified','sp_calc']:
				for i in reversed(range(len(outlines)-15,len(outlines))):
					if outlines[i-1].find('Atomic number out of range') > -1 or outlines[i-1].find('basis sets are only available') > -1:
						errortype = 'atomicbasiserror'
						break	
					if outlines[i].find('SCF Error') > -1:
						errortype = 'SCFerror'
						break
					if outlines[i].find('Normal termination') > -1:
						break

			# check for undesired imaginary freqs and data used by GoodVibes
			if termination == 'normal':
				atom_types,cartesians = cclib_atoms_coords(cclib_data)
				if errortype == 'none':
					# in eV, converted to hartree using the conversion factor from cclib
					E_dup = cclib_data['properties']['energy']['total']
					E_dup = E_dup/27.21138505
					# in hartree
					H_dup = cclib_data['properties']['enthalpy']
					G_dup = cclib_data['properties']['energy']['free energy']

					# detects if this calculation is a duplicate
					for i,_ in enumerate(E_dup_list):
						if abs(E_dup - E_dup_list[i]) < abs(self.dup_threshold):
							if abs(H_dup - H_dup_list[i]) < abs(self.dup_threshold):
								if abs(G_dup - G_dup_list[i]) < abs(self.dup_threshold):                        
									errortype = 'duplicate_calc'

				if errortype == 'none': 
					E_dup_list.append(E_dup)
					H_dup_list.append(H_dup)
					G_dup_list.append(G_dup)

					symmno,point_group,roconst,rotemp,errortype = self.symm_rot_data(outlines,qm_program,errortype)
					
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
					
					elif initial_ifreqs < 0:
						errortype = 'ts_no_imag_freq'

					if errortype in ['extra_imag_freq','freq_no_conv']:
						if errortype == 'extra_imag_freq':
							cartesians = self.fix_imag_freqs(n_atoms, cartesians, freqs, freq_displacements, calc_type)
							# in case no previous OPT was done (only works if it's not a TS)
							opt_found = False
							for keyword in keywords_line.split():
								if keyword.lower().startswith('opt'):
									opt_found = True
							if not opt_found:
								keywords_line += ' opt'
								
						elif errortype == 'freq_no_conv':
							# adjust the keywords so only FREQ is calculated
							new_keywords_line = ''
							for keyword in keywords_line.split():
								if keyword.lower().startswith('opt'):
									keyword = self.freq_conv
								new_keywords_line += keyword
								new_keywords_line += ' '
							keywords_line = new_keywords_line

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
						if keywords_line.find(' scf=xqc') > -1 or keywords_line.find(' scf=qc') > -1:
							new_keywords_line = ''
							for keyword in keywords_line.split():
								if keyword == 'scf=xqc':
									keyword = 'scf=qc'
								new_keywords_line += keyword
								new_keywords_line += ' '
							keywords_line = new_keywords_line
								
						else:
							keywords_line += ' scf=xqc'
					
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

			# check for isomerization
			if self.isom != False and errortype not in ['ts_no_imag_freq','atomicbasiserror','before_E_error']:
				isomerized = False
				init_csv = pd.DataFrame()
				os.chdir(self.isom_inputs)
				try:
					atoms_com, coords_com, atoms_and_coords = [],[],[]
					if len(self.isom.split('.')) == 1:
						atoms_and_coords,_ = get_info_input(f'{file_name}.{self.isom}')

					elif self.isom.split('.')[1] != 'csv':
						init_csv = pd.read_csv(self.args.isom)

					for line in atoms_and_coords:
						atoms_com.append(line.split()[0])
						coords_com.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])

					isomerized = check_isomerization(coords_com, cartesians, atoms_com, atom_types, self.vdwfrac, self.covfrac, init_csv, file)
				
				except FileNotFoundError:
					print(f'x  No com file were found for {file}, the check_geom test will be disabled for this calculation')
					self.log.write(f'x  No com file were found for {file}, the check_geom test will be disabled for this calculation')
				
				if isomerized:
					errortype = 'isomerization'
				os.chdir(self.w_dir_main)

			if errortype not in ['ts_no_imag_freq','atomicbasiserror','before_E_error','isomerization','duplicate_calc','spin_contaminated','none','sp_calc']:
				qcorr_calcs = qprep(destination=Path(f'{os.getcwd()}/unsuccessful_QM_outputs/run_{self.round_num}/fixed_QM_inputs'), w_dir_main=self.w_dir_main,
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
			aqme_data['AQME data']['creation date'] =creation_date
			if termination == 'normal' and errortype == 'none':
				aqme_data['AQME data']['symmetry number'] = symmno
				aqme_data['AQME data']['point group'] = point_group
				aqme_data['AQME data']['rotational constants'] = roconst
				aqme_data['AQME data']['rotational temperatures'] = rotemp
			
			# update and move the json files to their corresponding destination
			with open(file_name+'.json', "r+") as json_file:
				data = json.load(json_file)
				data.update(aqme_data)
				json_file.seek(0)
				json.dump(data, json_file, indent=4, separators=(", ", ": "))

			destination_json = destination.joinpath('json_files/')
			move_file(destination_json, self.w_dir_main, file_name+'.json')

		if self.fullcheck:
			try: 
				destination_fullcheck = self.w_dir_main.joinpath('successful_QM_outputs/json_files/')
				json_files = glob.glob(f'{destination_fullcheck}/*.json')
			
				full_check(w_dir_main=destination_fullcheck,destination_fullcheck=destination_fullcheck,json_files=json_files)
			
			except FileNotFoundError:
				print('x  No normal terminations with no errors to run the full check analysis')
				self.log.write('x  No normal terminations with no errors to run the full check analysis')	
		
			
	# include geom filters (ongoing work)

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
		qm_solv,qm_emp = 'gas_phase','none'
		orca_solv_1, orca_solv_2, orca_solv_3 = 'gas phase','',''
		
		for i in range(0,len(outlines)):
			if outlines[i].find(' Gaussian ') > -1 and outlines[i].find('Revision') > -1:
				qm_program = outlines[i][1:-2]

			elif outlines[i].find('* O   R   C   A *') > -1:
				qm_program = 'ORCA'

			if qm_program.lower().find('gaussian') > -1:
				# retrieve multiple parameters of the calculation
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
							elif outlines[j].lower().find('%pal') > -1 and outlines[j].lower().find('nprocs') > -1 and outlines[j].lower().find('end') > -1:
								nprocs = int(outlines[j].split()[-2])
							elif outlines[j].lower().find('%pal') > -1 and outlines[j].lower().find('nprocs') > -1 and outlines[j+1].lower().find('end') > -1:
								nprocs = int(outlines[j].split()[-1])
								skip_lines = 1
							elif outlines[j].lower().find('%pal') > -1 and outlines[j+1].lower().find('nprocs') > -1 and outlines[j+1].lower().find('end') > -1:
								nprocs = int(outlines[j+1].split()[-2])
								skip_lines = 1
							elif outlines[j].lower().find('%pal') > -1 and outlines[j+1].lower().find('nprocs') > -1 and outlines[j+2].lower().find('end') > -1:
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


	def symm_rot_data(self,outlines,qm_program,errortype):
		"""
		Retrieves information from QM files regarding symmetry and rotational parameters used by GoodVibes.


		Parameters
		----------
		outlines : list of str
			Lines of the QM output files
		qm_program : str
			Program and version used in the calculation
		errortype : string
			Type of error type of the QM output file (i.e. None, not_specified, extra_imag_freq, etc)

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
		for i in reversed(range(0,len(outlines)-15)):
			if qm_program.lower().find('gaussian') > -1:
				if outlines[i].find('Rotational symmetry number') > -1:
					symmno = int((outlines[i].strip().split()[3]).split(".")[0])

				elif outlines[i].find('Full point group') > -1:
					point_group = outlines[i].strip().split()[3]

				elif outlines[i].find('Rotational constants (GHZ):') > -1:
					try:
						roconst = [float(outlines[i].strip().replace(':', ' ').split()[3]),
										float(outlines[i].strip().replace(':', ' ').split()[4]),
										float(outlines[i].strip().replace(':', ' ').split()[5])]
					except ValueError:
						if outlines[i].find('********') > -1:
							roconst = [float(outlines[i].strip().replace(':', ' ').split()[4]),
											float(outlines[i].strip().replace(':', ' ').split()[5])]
				
				elif outlines[i].find('Rotational temperature ') > -1:
					rotemp = [float(outlines[i].strip().split()[3])]
				elif outlines[i].find('Rotational temperatures') > -1:
					try:
						rotemp = [float(outlines[i].strip().split()[3]), float(outlines[i].strip().split()[4]),
								  float(outlines[i].strip().split()[5])]
					except ValueError:
						if outlines[i].find('********') > -1:
							rotemp = [float(outlines[i].strip().split()[4]), float(outlines[i].strip().split()[5])]
							
				if point_group != '' and symmno != '' and roconst != [] and rotemp != []:
					break

				elif self.freq_conv and outlines[i].find('Converged?') > -1:
					for j in range(i+1,i+5):
						if outlines[j].strip().split()[-1] == 'NO':
							errortype = 'freq_no_conv'
					if errortype == 'freq_no_conv':
						break
				
				elif outlines[i].find('Normal termination') > -1:
					break

		return symmno,point_group,roconst,rotemp,errortype


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
		creation_date : str
			Date of creation of the QM calculations (when they finish)
		"""  

		grid_size,s2_operator,s2_preannhi,s2_found,level_of_theory,creation_date = '',0,0,False,'',''
		grid_lookup = {1: 'sg1', 2: 'coarse', 4: 'fine', 5: 'ultrafine', 7: 'superfine'}
		found_theory,level,bs = False,'',''
		line_options = ['\\freq\\','|freq|','\\sp\\','|sp|']

		for i in reversed(range(0,len(outlines))):
			# get grid size
			if qm_program.lower().find('gaussian') > -1:
				if outlines[i].find('ExpMin=') > -1:
					IRadAn = int(outlines[i].strip().split()[-3])
					grid_size = grid_lookup[IRadAn]
				elif outlines[i].find('S**2 before annihilation') > -1 and not s2_found:
					s2_operator = float(outlines[i].strip().split()[-1])
					s2_preannhi = float(outlines[i].strip().split()[-3][:-1])
					s2_found = True

				elif outlines[i].find('Normal termination') > -1:
					creation_date = ' '.join(item for item in outlines[i].strip().split()[-4:])
					creation_date = creation_date[:-1]

				# get functional and basis set
				if not found_theory:
					if outlines[i].find('External calculation') > -1:
						level, bs = 'ext', 'ext'
						level_of_theory = '/'.join([level, bs])
						found_theory = True
					if i < len(outlines)-1:
						for option in line_options:
							if option in (outlines[i].lower()+outlines[i+1].lower()):
								level, bs = ((outlines[i]+outlines[i+1]).strip().replace('|','\\').split("\\")[4:6])
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

		# designed to detect G4 calcs
		if level_of_theory == 'HF/GFHFB2':
			level_of_theory = 'G4'

		return grid_size,s2_operator,s2_preannhi,level_of_theory,creation_date


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
			destination = self.w_dir_main.joinpath('successful_QM_outputs/')
			file_terms['finished'] += 1

		elif errortype == 'sp_calc' and termination == "normal":
			destination = self.w_dir_main.joinpath('successful_QM_outputs/SP_calcs/')
			file_terms['sp_calcs'] += 1

		elif errortype == 'extra_imag_freq':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/extra_imag_freq/')
			file_terms['extra_imag_freq'] += 1

		elif errortype == 'ts_no_imag_freq':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/ts_no_imag_freq/')
			file_terms['ts_no_imag_freq'] += 1

		elif errortype == 'spin_contaminated':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/spin_contaminated/')
			file_terms['spin_contaminated'] += 1

		elif errortype == 'duplicate_calc':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/duplicates/')
			file_terms['duplicate_calc'] += 1

		elif errortype == "atomicbasiserror":
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/error/basis_set_error/')
			file_terms['atom_error'] += 1

		elif errortype == "SCFerror":
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/error/scf_error/')
			file_terms['scf_error'] += 1

		elif errortype == "before_E_calculation":
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/error/before_E_calculation/')
			file_terms['before_E_error'] += 1

		elif errortype == 'fail_geom_rules':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/geom_rules_filter/')
			file_terms['geom_rules_qcorr'] += 1

		elif errortype == 'isomerization':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/isomerization/')
			file_terms['isomerized'] += 1

		elif errortype == 'freq_no_conv':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/freq_no_conv/')
			file_terms['freq_no_conv'] += 1

		else:
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/error/not_specified_error/')
			file_terms['not_specified'] += 1
		
		move_file(destination, self.w_dir_main, file)

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
		ana_data.at[0,'Freq not converged'] = file_terms['freq_no_conv']
		ana_data.at[0,'SCF error'] = file_terms['scf_error']
		ana_data.at[0,'Error before SCF'] = file_terms['before_E_error']
		ana_data.at[0,'Basis set error'] =  file_terms['atom_error']
		ana_data.at[0,'Other errors'] = file_terms['not_specified']
		if self.args.s2_threshold > 0.0:
			ana_data.at[0,'Spin contamination'] = file_terms['spin_contaminated']
		ana_data.at[0,'Duplicates'] = file_terms['duplicate_calc']
		if len(self.args.geom_rules) >= 1:
			ana_data.at[0,'geom_rules filter'] = file_terms['geom_rules_qcorr']
		if self.isom != False:
			ana_data.at[0,'Isomerization'] = file_terms['isomerized']
		path_as_str = self.w_dir_main.as_posix()
		ana_data.to_csv(path_as_str+f'/QCORR-run_{self.round_num}-stats.csv',index=False)

		move_file(self.w_dir_main.joinpath('dat_files/'), self.w_dir_main,f'QCORR-run_{self.round_num}-stats.csv')


def full_check(w_dir_main=os.getcwd(),destination_fullcheck='',json_files=glob.glob('*.json')):
	"""
	Checks that multiple calculations were done following the same protocols, including
	program and version, grid size, level of theory, dispersion and solvation model.

	Parameters
	----------
	w_dir_main : str
		Working directory
	destination_fullcheck : str
		Destination to create the file with the full check
	json_files : list of str
		json files to compare (glob.glob('*.json') and '*.json are both valid inputs to
		include all the json files from a folder)	
	"""

	initial_dir = os.getcwd()
	w_dir_main = Path(w_dir_main)
	os.chdir(w_dir_main)

	if json_files == '*.json' or json_files == '*json':
		json_files=glob.glob('*.json')

	df_fullcheck = pd.DataFrame(columns=['file', 'program', 'grid', 'lot', 'dispersion', 'solvation'])

	for file in json_files:
		file_name = file.split('.')[0]
		with open(file) as json_file:
			cclib_data = json.load(json_file)
		program = cclib_data['AQME data']['QM program']
		grid = cclib_data['AQME data']['grid']
		lot = cclib_data['AQME data']['level of theory']
		dispersion = cclib_data['AQME data']['dispersion']
		solvation = cclib_data['AQME data']['solvation']
		df_fullcheck.loc[len(df_fullcheck.index)] = [file_name, program, grid, lot, dispersion, solvation] 
	
	fullcheck_file = '--QCORR_Fullcheck_Analysis--.dat'
	fullcheck_txt = '-- Full check analysis --'

	for prop in df_fullcheck.columns:
		if prop != 'file':
			unique_props = df_fullcheck[prop].unique()
			if len(unique_props) > 1:
				fullcheck_txt += f'\nx  Different {prop} used in the calculations:'
				for unique_prop in unique_props:
					file_names = df_fullcheck["file"].loc[df_fullcheck[prop] == unique_prop]
					fullcheck_txt += f'\n     * {unique_prop} in:'
					for file_name in file_names:
						adapted_name = file_name.replace('/','\\').split("\\")[-1]
						fullcheck_txt += f'\n       - {adapted_name}'
			else:
				fullcheck_txt += f'\no  Same {prop} ({unique_props[0]}) used in all the calculations'
			
	fullcheck_analysis = open(fullcheck_file, 'w')
	print(fullcheck_txt)
	fullcheck_analysis.write(fullcheck_txt)	
	fullcheck_analysis.close()

	if destination_fullcheck == '':
		destination_fullcheck = w_dir_main.joinpath('successful_QM_outputs/json_files/')
	else:
		destination_fullcheck = Path(destination_fullcheck)
	move_file(destination_fullcheck, w_dir_main, fullcheck_file)
	
	os.chdir(initial_dir)


def json2input(json_files=[], w_dir_main=os.getcwd(), destination=None, suffix='', 
				charge=None, mult=None,	mem='8GB', nprocs=4, chk=False, qm_input='', bs_gen='', 
				bs='', gen_atoms=[], qm_end='', program='gaussian'):
	'''
	Reads a json file and use QPREP to generate input files.

	Parameters
	----------
	json_files : list 
		Filenames of json files to analyze
	w_dir_main : str
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
	'''
	
	w_dir_initial=os.getcwd()
	os.chdir(w_dir_main)

	charge_initial = charge
	mult_initial = mult

	if destination is None:
		destination = Path(w_dir_main).joinpath("QCALC")
	else:
		destination = Path(destination)

	if not isinstance(json_files, list): 
		json_files = glob.glob(json_files)

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
		if charge_initial == None:
			charge = cclib_data['properties']['charge']
		if mult_initial == None:
			mult = cclib_data['properties']['multiplicity']
		if charge == None:
			print('x  No charge was specified in the json file or function input (i.e. json2input(charge=0) )')
		elif mult == None:
			print('x  No multiplicity was specified in the json file or function input (i.e. json2input(mult=1) )')

		json_calcs = qprep(destination=destination, w_dir_main=w_dir_main,
					molecule=file_name, charge=charge, mult=mult,
					program=program, atom_types=atom_types,
					cartesians=cartesians, qm_input=qm_input,
					mem=mem, nprocs=nprocs, chk=chk, qm_end=qm_end,
					bs_gen=bs_gen, bs=bs, gen_atoms=gen_atoms, suffix=suffix)
	
	print(f'o  Final input files were generated in {destination}')

	os.chdir(w_dir_initial)
	

def check_run(w_dir):
	"""
	# Determines the folder where input files are gonna be generated in QCORR.
	"""

	input_folder = w_dir+'/unsuccessful_QM_outputs/'
	folder_count = 0

	if os.path.exists(input_folder):
		dir_list = os.listdir(input_folder)
		for folder in dir_list:
			if folder.find('run_') > -1:
				folder_count += 1

	return folder_count