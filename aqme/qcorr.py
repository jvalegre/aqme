######################################################.
#        This file stores all the functions          #
#          used in the LOG file analyzer             #
######################################################.
import os
import sys
import glob
import pandas as pd
import json
import subprocess
from pathlib import Path

from aqme.utils import (
	move_file,
    Logger,
    load_from_yaml,
    check_isomerization,
    read_file,
    cclib_atoms_coords,
	periodic_table,
	get_info_input,
	get_metadata,
	get_s2,
	detect_linear
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
		keywords_line for new input files
	s2_threshold :  float
		Cut off for spin contamination during analysis in \%\ of the expected value (i.e. multiplicity 3 has an the expected <S**2> of 2.0, if s2_threshold = 10, the <S**2> value is allowed to be 2.0 +- 0.2). Set s2_threshold = 0 to deactivate this option.
	isom : str
		Check for isomerization from the initial input file to the resulting output files. It requires the extension of the initial input files (i.e. isom='com') and the folder of the input files must be added in the isom_inputs option
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
		If a string is defined, it will remove calculations that converged during optimization but did not convergence in the subsequent frequency calculation. Options: opt sections as strings i.e. (opt=(calcfc,maxstep=5)). If readfc is specified in the string, the chk option must be included as well.
	ifreq_cutoff : float
		Cut off for to consider whether a frequency is imaginary (absolute of the specified value is used)
	fullcheck : bool
		Perform an analysis to detect whether the calculations were done homogeneously (i.e. same level of theory, solvent, grid size, etc)
	program : str
		Program required to create the new input files
	varfile : str
		Option to parse the variables using a yaml file (specify the filename)
	kwargs : argument class
		Specify any arguments from the QCORR module
	"""

	def __init__(self, qm_files='', w_dir_main=os.getcwd(), dup_threshold=0.0001,
				mem='4GB', nprocs=2, chk=False, qm_input='', s2_threshold=10.0,
				isom=None, isom_inputs=os.getcwd(), vdwfrac=0.50, covfrac=1.10,
				bs_gen='', bs='', gen_atoms=[], qm_end='', amplitude_ifreq=0.2, freq_conv=None,
				ifreq_cutoff=0.0, fullcheck=True, program='gaussian', varfile=None, **kwargs):

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

		1. Analyzes the QM output files and moves output files with normal termination and no extra imaginary frequencies to the same folder
		2. Generates input files to fix errors and extra imaginary frequencies
		3. Generates input files with new keywords_lines from the normally terminated files from point 1 (i.e. single-point energy corrections). Optionally, the analysis from points 1 and 2  might be disabled with the nocheck=True or --nocheck option.
		"""

		file_terms = {'finished': 0, 'sp_calcs' : 0, 'extra_imag_freq': 0, 'ts_no_imag_freq': 0,
			'freq_no_conv': 0, 'spin_contaminated': 0, 'duplicate_calc': 0, 'atom_error': 0,
			'scf_error': 0, 'no_data': 0, 'linear_mol_wrong': 0, 'not_specified': 0,
			'geom_rules_qcorr': 0, 'isomerized': 0}
		
		E_dup_list, H_dup_list, G_dup_list = [],[],[]

		for file in self.qm_files:
			# get initial cclib data and termination/error types
			file_name = file.split('.')[0]
			termination,errortype,cclib_data,n_atoms,freqs,freq_displacements,charge,mult,keywords_line,calc_type,mem,nprocs = self.cclib_init(file,file_name)

			if errortype == 'no_data':
				file_terms,_ = self.organize_outputs(file,termination,errortype,file_terms)
				continue

			# use very short reversed loop to find basis set incompatibilities and SCF errors
			if errortype in ['not_specified','sp_calc']:
				outlines = read_file(self.w_dir_main,file)
				for i in reversed(range(len(outlines)-15,len(outlines))):
					if outlines[i-1].find('Atomic number out of range') > -1 or outlines[i-1].find('basis sets are only available') > -1:
						errortype = 'atomicbasiserror'
						break
					if outlines[i].find('SCF Error') > -1:
						errortype = 'SCFerror'
						break
					if outlines[i].find('Normal termination') > -1:
						break

			# check for duplicates and wrong number of freqs
			if termination == 'normal':
				atom_types,cartesians = cclib_atoms_coords(cclib_data)

			if errortype == 'none':
				# in eV, converted to hartree using the conversion factor from cclib
				E_dup = cclib_data['properties']['energy']['total']
				E_dup = E_dup/27.21138505
				# in hartree
				try:
					H_dup = cclib_data['properties']['enthalpy']
					G_dup = cclib_data['properties']['energy']['free energy']
				except (AttributeError,KeyError):
					if n_atoms == 1:
						if keywords_line.find('freq') == -1:
							errortype = 'sp_calc'
							calc_type = 'SP calculation'
						H_dup = E_dup
						G_dup = E_dup
				# detects if this calculation is a duplicate
				for i,_ in enumerate(E_dup_list):
					E_diff = abs(E_dup - E_dup_list[i])
					H_diff = abs(H_dup - H_dup_list[i])
					G_diff = abs(G_dup - G_dup_list[i])
					if max([E_diff,H_diff,G_diff]) < abs(self.dup_threshold):
						errortype = 'duplicate_calc'

			if errortype == 'none':
				E_dup_list.append(E_dup)
				H_dup_list.append(H_dup)
				G_dup_list.append(G_dup)

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
				
				if len(atom_types) in [3,4]:
					errortype = detect_linear(errortype,atom_types,freqs)

			if errortype in ['extra_imag_freq','freq_no_conv','linear_mol_wrong']:
				if errortype == 'extra_imag_freq':
					cartesians = self.fix_imag_freqs(n_atoms, cartesians, freqs, freq_displacements, calc_type)

				# in case no previous OPT was done (only works if it's not a TS)
				opt_found = False
				for keyword in keywords_line.split():
					if keyword.lower().startswith('opt'):
						opt_found = True
				if not opt_found:
					keywords_line += ' opt'

				if errortype == 'linear_mol_wrong':
					keywords_line += ' symmetry=(PG=Cinfv)'

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
			elif termination != 'normal' and errortype not in ['ts_no_imag_freq','atomicbasiserror','no_data']:
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
						for _,line in enumerate(outlines):
							if line.find('Standard orientation:') > -1:
								count_RMS += 1
							if count_RMS == min_RMS:
								range_lines = i+5,i+5+n_atoms
								break
						for i in range(range_lines):
							massno = int(outlines[i].split()[1])
							if massno < len(per_tab):
								atom_symbol = per_tab[massno]
							else:
								atom_symbol = "XX"
							atom_types.append(atom_symbol)
							cartesians.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

			# check for isomerization
			if self.isom is not None and errortype not in ['ts_no_imag_freq','atomicbasiserror','no_data']:
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

			# move initial QM input files (if the files are placed in the same folder as the output files)
			if os.path.exists(f'{self.w_dir_main}/{file_name}.com'):
				move_file(self.w_dir_main.joinpath('initial_QM_inputs/'), self.w_dir_main,f'{file_name}.com')

			if errortype not in ['ts_no_imag_freq','atomicbasiserror','no_data','isomerization','duplicate_calc','spin_contaminated','none','sp_calc']:
				# user-defined keywords_line, mem and nprics overwrites previously used parameters
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

				qprep(destination=Path(f'{self.w_dir_main}/unsuccessful_QM_outputs/run_{self.round_num}/fixed_QM_inputs'), w_dir_main=self.w_dir_main,
					molecule=file_name, charge=charge, mult=mult,
					program=self.program, atom_types=atom_types,
					cartesians=cartesians, qm_input=keywords_line,
					mem=mem, nprocs=nprocs, chk=self.chk, qm_end=self.qm_end,
					bs_gen=self.bs_gen, bs=self.bs, gen_atoms=self.gen_atoms)

			print(f'{file}: Termination = {termination}, Error type = {errortype}')
			self.log.write(f'{file}: Termination = {termination}, Error type = {errortype}')

			# This part places the calculations in different folders depending on the type of termination
			file_terms,destination = self.organize_outputs(file,termination,errortype,file_terms)
			
			destination_json = destination.joinpath('json_files/')
			move_file(destination_json, self.w_dir_main, file_name+'.json')

			# write information about the QCORR analysis in a csv
			self.write_qcorr_csv(file_terms)

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
	# 					mol = output_to_mol(file,format_file)
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


	def cclib_init(self,file,file_name):
		'''
		Determine termination and error types (initial determination), create json files
		with cclib and load the data in the cclib json files
		'''
		# create a json file with cclib and load a dictionary. This protocol is
		# favored over the traditional ccread since more data will be added to
		# the json file at the end
		command_run_1 = ['ccwrite', 'json', file]
		subprocess.run(command_run_1)
		
		try:
			with open(file_name+'.json') as json_file:
				cclib_data = json.load(json_file)
		except FileNotFoundError:
			print(f'x  Potential cclib compatibility problem with file {file}')
			self.log.write(f'x  Potential cclib compatibility problem with file {file}')
			termination = 'other'
			errortype = 'no_data'
			
		# retrieves information from the cclib parser
		n_atoms = cclib_data['properties']['number of atoms']
		charge = cclib_data['properties']['charge']
		mult = cclib_data['properties']['multiplicity']
		keywords_line,calc_type,mem,nprocs = get_metadata(cclib_data)
		s2_after_anni,s2_before_anni = get_s2(cclib_data)

		# determine error/unfinished vs normal terminations and freqs for normal terminations
		try:
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
			if unpaired_e == 0 and s2_after_anni != 0:
				if s2_after_anni > abs(self.s2_threshold/100)*s2_before_anni:
					errortype = 'spin_contaminated'
			else:
				spin = unpaired_e*0.5
				s2_expected_value = spin*(spin+1)
				spin_diff = abs(float(s2_after_anni)-s2_expected_value)
				if spin_diff > abs(self.s2_threshold/100)*s2_expected_value:
					errortype = 'spin_contaminated'
			
			if errortype == 'none' and self.freq_conv is not None:
				freq_conv = cclib_data['optimization']['geometric values'][-1]
				freq_conv_targets = cclib_data['optimization']['geometric targets']
				for i,conv in enumerate(freq_conv):
					if conv > freq_conv_targets[i]:
						errortype = 'freq_no_conv'

		except (AttributeError,KeyError):
			# if the optimization finished, only a freq job is required
			try:
				if cclib_data['optimization']['done']:
					termination = 'other'
					errortype = 'no_freq'
			except (AttributeError,KeyError):
				# single-point calcs (normal terminations with no opt)
				try:
					cclib_data['optimization']['geometric values']
					termination = 'other'
					errortype = 'not_specified'
				except (AttributeError,KeyError):
					# general errors
					cclib_data['atoms']['elements']['number']
					termination = 'normal'
					errortype = 'sp_calc'
					calc_type = 'SP calculation'

		return termination,errortype,cclib_data,n_atoms,freqs,freq_displacements,charge,mult,keywords_line,calc_type,mem,nprocs


	def fix_imag_freqs(self, n_atoms, cartesians, freqs, freq_displacements, calc_type):
		"""
		Fixes undersired (extra) imaginary frequencies from QM calculations. This function multiplies the imaginary normal mode vectors by the selected amplitude (0.2 is the default amplitude in the pyQRC script from GitHub, user: bobbypaton).	By default, all the extra imaginary modes are used (i.e. in calculations with three	extra imaginary frequencies, all the three modes will be used to displace the atoms). This can be tuned with the --ifreq_cutoff option (i.e. only use freqs lower than -50 cm-1).

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

		elif errortype == "no_data":
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/error/no_data/')
			file_terms['no_data'] += 1

		elif errortype == 'fail_geom_rules':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/geom_rules_filter/')
			file_terms['geom_rules_qcorr'] += 1

		elif errortype == 'isomerization':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/isomerization/')
			file_terms['isomerized'] += 1

		elif errortype == 'freq_no_conv':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/freq_no_conv/')
			file_terms['freq_no_conv'] += 1

		elif errortype == 'linear_mol_wrong':
			destination = self.w_dir_main.joinpath(f'unsuccessful_QM_outputs/run_{self.round_num}/linear_mol_wrong/')
			file_terms['linear_mol_wrong'] += 1

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
		ana_data.at[0,'Single-point calcs'] = file_terms['sp_calcs']
		ana_data.at[0,'Extra imag. freq.'] = file_terms['extra_imag_freq']
		ana_data.at[0,'TS with no imag. freq.'] = file_terms['ts_no_imag_freq']
		ana_data.at[0,'Freq not converged'] = file_terms['freq_no_conv']
		ana_data.at[0,'Linear mol with wrong n of freqs'] = file_terms['linear_mol_wrong']
		ana_data.at[0,'SCF error'] = file_terms['scf_error']
		ana_data.at[0,'Error before SCF'] = file_terms['no_data']
		ana_data.at[0,'Basis set error'] =  file_terms['atom_error']
		ana_data.at[0,'Other errors'] = file_terms['not_specified']
		if self.args.s2_threshold > 0.0:
			ana_data.at[0,'Spin contamination'] = file_terms['spin_contaminated']
		ana_data.at[0,'Duplicates'] = file_terms['duplicate_calc']
		if len(self.args.geom_rules) >= 1:
			ana_data.at[0,'geom_rules filter'] = file_terms['geom_rules_qcorr']
		if self.isom is not None:
			ana_data.at[0,'Isomerization'] = file_terms['isomerized']
		path_as_str = self.w_dir_main.as_posix()
		ana_data.to_csv(path_as_str+f'/QCORR-run_{self.round_num}-stats.csv',index=False)

		move_file(self.w_dir_main.joinpath('dat_files/'), self.w_dir_main,f'QCORR-run_{self.round_num}-stats.csv')


def full_check(w_dir_main=os.getcwd(),destination_fullcheck='',json_files='*.json'):
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

	df_fullcheck = pd.DataFrame(columns=['file', 'program', 'grid_type', 'level_of_theory', 'dispersion', 'solvation'])

	for file in json_files:
		file_name = file.split('.')[0]
		with open(file) as json_file:
			cclib_data = json.load(json_file)

		program = cclib_data['metadata']['QM program']
		solvation = cclib_data['metadata']['solvation']
		dispersion = cclib_data['metadata']['dispersion model']			
		grid_type = cclib_data['metadata']['grid type']
		functional = cclib_data['metadata']['functional']
		bs = cclib_data['metadata']['basis set']
		if functional != '' or bs != '':
			level_of_theory = '/'.join([functional, bs])
		else:
			level_of_theory = ''
		# designed to detect G4 calcs
		if level_of_theory == 'HF/GFHFB2':
			level_of_theory = 'G4'
		df_fullcheck.loc[len(df_fullcheck.index)] = [file_name, program, grid_type, level_of_theory, dispersion, solvation]
	
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


def json2input(json_files='', w_dir_main=os.getcwd(), destination=None, suffix='',
				charge=None, mult=None,	mem='8GB', nprocs=4, chk=False, qm_input='', bs_gen='',
				bs='', gen_atoms=[], qm_end='', program='gaussian'):
	'''
	Reads a json file and use QPREP to generate input files.

	Parameters
	----------
	json_files : list of str
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
		keywords_line for new input files
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
		if charge_initial is None:
			charge = cclib_data['properties']['charge']
		if mult_initial is None:
			mult = cclib_data['properties']['multiplicity']
		if charge is None:
			print('x  No charge was specified in the json file or function input (i.e. json2input(charge=0) )')
		elif mult is None:
			print('x  No multiplicity was specified in the json file or function input (i.e. json2input(mult=1) )')

		qprep(destination=destination, w_dir_main=w_dir_main,
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
	folder_count = 1

	if os.path.exists(input_folder):
		dir_list = os.listdir(input_folder)
		for folder in dir_list:
			if folder.find('run_') > -1:
				folder_count += 1

	return folder_count
	