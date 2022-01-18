######################################################.
#        This file stores all the functions          #
#          used in the LOG file analyzer             #
######################################################.
import os
import sys
import glob
import pandas as pd
import json

# the import for check and write genecp is probably inside the genecp function of qprep now, FIX!
from aqme.utils import periodic_table
from aqme.filter import geom_rules_output
from aqme.utils import (move_file, get_info_input, Logger, load_from_yaml, 
							check_isomerization, read_file, output_to_mol)
from aqme.argument_parser import set_options


def check_run(w_dir):
	"""
	# Determines the folder where input files are gonna be generated in QCORR.
	"""
	input_folder = w_dir+'/fixed_QM_inputs_files'
	folder_count = 0
	
	if os.path.exists(input_folder):
		dir_list = os.listdir(input_folder)
		for folder in dir_list:
			if folder.find('run_') > -1:
				folder_count += 1

	if folder_count == 0:
		return 0
	else:
		num_com_folder = sum([len(d) for r, d, folder in os.walk(w_dir+'/fixed_QM_inputs_files')])
		w_dir = w_dir+'/fixed_QM_inputs_files/run_'+str(num_com_folder)
		return folder_count


class qcorr():
	"""
	Class containing all the functions from the QCORR module.

	Parameters
	----------
	qm_files : list 
		Contains the filenames of QM output files to analyze
	w_dir_main : str
		Working directory
	yaml_file : str
		Option to parse the variables using a yaml file (specify the filename)
	kwargs : argument class
		Specify any arguments from the QCORR module
	"""
	
	def __init__(self, qm_files=[], round_num=0, w_dir_main=os.getcwd(), 
				yaml_file=None, **kwargs):
		
		self.qm_files = qm_files
		self.w_dir_main = w_dir_main
		
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
		if not os.path.isdir(self.w_dir_main+'/dat_json_files/'):
			os.makedirs(self.w_dir_main+'/dat_json_files/')
		self.log = Logger(self.w_dir_main+'/dat_json_files/aqme',f'QCORR-run_{str(self.round_num)}')
		self.log.write("\no  Analyzing output files in {}\n".format(self.w_dir_main))
		
		if len(qm_files) == 0:
			self.log.write('x  There are no output files in this folder.')
			sys.exit('x  There are no output files in this folder.')


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

		# these lists will contain information returned by QCORR
		file_name_list, charge_list, mult_list, molecular_mass_list = [],[],[],[]
		atom_types_list, cartesians_list = [],[]
		freqs_list, freq_displacements_list = [],[]
		keywords_line_list, grid_size_list = [],[]
		termination_list, errortype_list = [],[]

		# these lists will track duplicates if args.dup is activated
		file_list, E_dup, H_dup, G_dup = [],[],[],[]
		file_terms = {'finished': 0, 'imag_freq': 0, 'ts_no_imag_freq': 0, 
					'spin_contaminated': 0, 'duplicate_calc': 0, 'atom_error': 0,
					'scf_error': 0, 'before_E_error': 0, 'other_error': 0,
					'unfinished': 0, 'geom_rules_qcorr': 0, 'check_geom_qcorr': 0}   

		for file in self.qm_files:

			file_name = file.split('.')[0]

			# read the file
			self.log.write(file)
			outlines, outfile, break_loop = read_file(self.w_dir_main,file)

			if break_loop:
				break

			# get parameters from the calculations (i.e. termination and error types, n of atoms, charge, mult, etc)
			termination, errortype, initial_ifreqs, n_atoms, charge, mult, keywords_line, calc_type, nimag_line, grid_size = self.get_qcorr_params(outlines)
			
			# get new coordinates as input files to fix error (unless the errors cannot be fixed automatically
			# such as in atomic basis errors or errors happening too early on the calculation)
			if errortype not in ['before_E_calculation',"atomicbasiserror"] or self.args.nocheck:

				# get geometry parameters and frequency information
				file_list, E_dup, H_dup, G_dup, errortype, keywords_line, atom_types, cartesians, freqs, freq_displacements, molecular_mass = self.get_input_geom(file, outlines, keywords_line, n_atoms, mult, termination, errortype, initial_ifreqs, nimag_line, calc_type, file_list, E_dup, H_dup, G_dup, self.args.nocheck)
					
				if self.args.isom != None:
					isomerized = False
					init_csv = pd.DataFrame()
					folder_dir_isom = ''
					self.log.write("  ----- Geometrical check will be applied to the output file -----\n")
					try:
						atoms_com, coords_com, atoms_and_coords = [],[],[]
						if self.round_num != 0:
							folder_dir_isom = self.w_dir_main +'/fixed_QM_inputs_files/run_'+str(self.round_num)
						else:
							pass
						if self.args.isom == 'com':
							atoms_and_coords,_ = get_info_input(folder_dir_isom+file.split('.')[0]+'.com')
						elif self.args.isom == 'gjf':
							atoms_and_coords,_ = get_info_input(folder_dir_isom+file.split('.')[0]+'.gjf')
						elif self.args.isom.split('.')[1] == 'csv':
							init_csv = pd.read_csv(self.args.isom)

						for line in atoms_and_coords:
							atoms_com.append(line.split()[0])
							coords_com.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
						
						isomerized = check_isomerization(coords_com, cartesians, atoms_com, atom_types, self.args.vdwfrac, self.args.covfrac, init_csv, file)
					
					except FileNotFoundError:
						self.log.write("x  No com file were found for "+file+", the check_geom test will be disabled for this calculation")
					
					if isomerized:
						errortype = 'isomerization'

				if len(self.args.geom_rules) >= 1:
					passing_rules = True
					valid_mol_gen = True
					self.log.write("  ----- geom_rules filter(s) will be applied to the output file -----\n")
					try:
						format_file = file.split('.')[1]
						mol,ob_compat,rdkit_compat = output_to_mol(file,format_file,self.log)
						print_error_geom_rules=False
						if ob_compat and rdkit_compat:
							passing_rules = geom_rules_output(mol,self.args,self.log,file,print_error_geom_rules)
							if not passing_rules:
								errortype = 'fail_geom_rules'
						os.remove(file.split('.')[0]+'.mol')
					except AttributeError:
						valid_mol_gen = False
						os.remove(file.split('.')[0]+'.mol')
						self.log.write("The file could not be converted into a mol object, geom_rules filter(s) will be disabled\n")

			# is it necessary to close the file?!
			# close the file
			outfile.close()
			
			if termination == "normal" and errortype == None:
				# detects if this calculation is a duplicate
				if self.args.dup:
					for i,energy_value in enumerate(E_dup):
						if i != len(E_dup)-1:
							if abs(energy_value - E_dup[file_list.index(file)]) < abs(self.args.dup_threshold):
								if abs(H_dup[i] - H_dup[file_list.index(file)]) < abs(self.args.dup_threshold):
									if abs(G_dup[i] - G_dup[file_list.index(file)]) < abs(self.args.dup_threshold):                        
										errortype = 'duplicate_calc'
			
			# This part places the calculations in different folders depending on the type of termination
			file_terms = self.organize_outputs(file,termination,errortype,file_terms)

			file_name_list.append(file_name)
			charge_list.append(charge)
			mult_list.append(mult)
			molecular_mass_list.append(molecular_mass)
			atom_types_list.append(atom_types)
			cartesians_list.append(cartesians)
			keywords_line_list.append(keywords_line)
			freqs_list.append(freqs)
			freq_displacements_list.append(freq_displacements)
			termination_list.append(termination)
			errortype_list.append(errortype)
			grid_size_list.append(grid_size)
			
		# write information about the QCORR analysis in a csv
		self.write_qcorr_csv(file_terms)
		
		# exports all the data as a DataFrame and splits it into successful and unsuccessful terminations
		qcorr_data_dict = {'File_name': file_name_list, 'Grid_size': grid_size_list, 'Charge': charge_list, 'Mult': mult_list,
					'MW': molecular_mass_list, 'Atom_types': atom_types_list, 'Cartesians': cartesians_list,
					'keywords_line': keywords_line_list, 'Freqs': freqs_list, 'Normal_mode': freq_displacements_list,
					'Termination': termination_list, 'Error_type': errortype_list}
		
		df_qcorr = pd.DataFrame(qcorr_data_dict)
		df_qcorr_success = df_qcorr[df_qcorr['Error_type'].isna()]
		df_qcorr_success = df_qcorr_success.drop(['Error_type','Termination'], axis=1)

		df_qcorr_error = df_qcorr[df_qcorr['Error_type'].notna()]
		df_qcorr_error = df_qcorr_error.drop(['MW'], axis=1)

		# loads previous json objects
		os.chdir(self.w_dir_main+'/dat_json_files/')
		json_files = glob.glob('*.json')
		df_success_json, df_error_json = pd.DataFrame(),pd.DataFrame()
		found_success_json, found_error_json = False, False
		round_success_global = -1
		round_error_global = -1
		for json_file in json_files:
			if json_file.find('qcorr_success_run_') > -1:
				round_success = json_file.split('.json')[0].split('_')[-1]
				if int(round_success) > round_success_global:
					round_success_global = int(round_success)
			elif json_file.find('qcorr_error_run_') > -1:
				round_error = json_file.split('.json')[0].split('_')[-1]
				if int(round_error) > round_error_global:
					round_error_global = int(round_error)

		if f'qcorr_success_run_{round_success_global}.json' in json_files:
			df_success_json = pd.read_json(f'qcorr_success_run_{round_success_global}.json')
			found_success_json = True
		if f'qcorr_error_run_{round_error_global}.json' in json_files:
			df_error_json = pd.read_json(f'qcorr_error_run_{round_error_global}.json')
			found_error_json= True

		os.chdir(self.w_dir_main)

		# merges new Normal termination results with previous results
		if found_success_json:
			for new_qcorr_name in df_qcorr_success['File_name']:
				if not new_qcorr_name in df_success_json.values:
					df_success_json = df_success_json.append(df_qcorr_success[df_qcorr_success['File_name'] == new_qcorr_name])
					# removes this error from past runs (since it was fixed in the new run)
					if new_qcorr_name in df_error_json.values:
						df_error_json = df_error_json[df_error_json['File_name'] != new_qcorr_name]

		else:
			df_success_json = df_qcorr_success

		if found_error_json:
			for new_qcorr_name in df_qcorr_error['File_name']:
				if not new_qcorr_name in df_error_json.values:
					df_error_json = df_error_json.append(df_qcorr_error[df_qcorr_error['File_name'] == new_qcorr_name])
		else:
			df_error_json = df_qcorr_error

		# sort values and save json files
		df_success_json.sort_values(by=['File_name'])
		df_success_json = df_success_json.reset_index(drop=True)
		df_success_json = df_success_json.to_dict()

		df_error_json.sort_values(by=['File_name'])
		df_error_json = df_error_json.reset_index(drop=True)
		df_error_json = df_error_json.to_dict()

		with open(f'{self.w_dir_main}/dat_json_files/qcorr_success_run_{int(round_error_global)+1}.json', 'w') as f:
			json.dump(df_success_json, f, indent=4, separators=(", ", ": ")) 
		with open(f'{self.w_dir_main}/dat_json_files/qcorr_error_run_{int(round_error_global)+1}.json', 'w') as f:
			json.dump(df_error_json, f, indent=4, separators=(", ", ": ")) 

		return df_qcorr_success,df_qcorr_error

	
	def get_qcorr_params(self,outlines):
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

		Returns
		-------
		termination : string
			Type of termination of the QM output file (i.e. normal, error, unfinished)
		errortype : string
			Error type of the QM output file
		initial_ifreqs : int
			Initial amount of imaginary frequencies in the QM file adapted for transition state 
			calculations (if the file is a TS, the original number is reduced by 1)
		n_atoms : int
			Number of atoms in the calculation
		charge : int
			Charge of the calculations used in the following input files
		mult : int
			Multiplicity of the calculations used in the following input files
		keywords_line : string
			Keyword line (input) to use in subsequent input files
		calc_type : str
			Type of the QM calculation (ground_state or transition_state)
		nimag_line : int
			Line number of the line containing the amount of imaginary frequencies in the 
			QM calculation. This is tracked to avoid reading the same lines multiple times
		"""    

		termination,errortype  = 'unfinished','unknown'
		stop_term, initial_ifreqs = 0, 0
		
		# use reversed loops to find type of termination (faster than forward loops)
		for i in reversed(range(len(outlines)-15,len(outlines))):
			if stop_term == 1:
				break
			# Determine the kind of job termination
			if outlines[i].find("Normal termination") > -1:
				termination = "normal"
				errortype = None
				stop_term += 1
			elif outlines[i].find("Error termination") > -1:
				termination = "error"
				if outlines[i-1].find("Atomic number out of range") > -1 or outlines[i-1].find("basis sets are only available") > -1 :
					errortype = "atomicbasiserror"
				if outlines[i-3].find("SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error SCF Error") > -1:
					errortype = "SCFerror"
				stop_term += 1
		
		# if the calculation ends normally with no imaginary freqs, no more read_lines are necessary (saves time)
		# the nimag_line is defined so that part of the file is not read again
		nimag_line = len(outlines)
		if termination == 'normal':
			charge, mult, n_atoms, keywords_line = None, None, 0, ''
			calc_type = None
			for i in reversed(range(60,len(outlines)-5)):
				if outlines[i].find('NImag=') > -1 or outlines[i].find('PG=') > -1:
					nimag_line = i
					# merge a few lines since 'NIMag=' sometimes appears in two different lines
					line_with_ifreq = outlines[i-1].rstrip("\n")+outlines[i].rstrip("\n")+outlines[i+1].rstrip("\n")
					for part_line in line_with_ifreq.split('\\'):
						if part_line.find('NImag=') > -1:
							initial_ifreqs = int(part_line.split('=')[1])
							break
					break
			
		# get keywords line (input), number of atoms, charge and multiplicity
		n_atoms,charge,mult,keywords_line,initial_ifreqs,calc_type,grid_size = self.get_init_info(outlines,initial_ifreqs)
		
		# detect calcs that finish before the functional or basis set was included
		if keywords_line == '':
			errortype = 'before_E_calculation'
		
		# helps to fix SCF convergence errors
		if errortype == 'SCFerror':
			if keywords_line.find(' scf=qc') > -1:
				pass
			else:
				keywords_line += ' scf=qc'
		
		return termination,errortype,initial_ifreqs,n_atoms,charge,mult,keywords_line,calc_type,nimag_line,grid_size


	def get_init_info(self,outlines,initial_ifreqs):
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
		initial_ifreqs : int
			Initial amount of imaginary frequencies in the QM file

		Returns
		-------
		n_atoms : int
			Number of atoms in the calculation
		charge : int
			Charge of the calculations used in the following input files
		mult : int
			Multiplicity of the calculations used in the following input files
		keywords_line : string
			Original keyword line (input) from the output QM file
		initial_ifreqs : int
			Initial amount of imaginary frequencies in the QM file adapted for transition state 
			calculations (if the file is a TS, the original number is reduced by 1)
		calc_type : str
			Type of the QM calculation (ground_state or transition_state)
		"""    

		end_keywords,keywords_line,grid_size = False,'',''
		symbol_z_line, gradgrad_line = None, None
		charge, mult, n_atoms = None, None, 0
		
		# range optimized to skip Gaussian information (short for loop just to get charge, multiplicity and n of atoms)
		for i in range(60,len(outlines)):
			# retrieve the input line
			if outlines[i].startswith(' #'):
				keywords_line += outlines[i].rstrip("\n")
				for j in range(i+1,i+10):
					if outlines[j].find('----------') > -1:
						end_keywords = True
						break
					if not end_keywords:
						keywords_line += outlines[j].rstrip("\n")[1:]
				keywords_line = keywords_line[3:]
			
			# get number of atoms
			if n_atoms == 0:
				if outlines[i].find('Symbolic Z-matrix:') > -1:
					symbol_z_line = i
				elif outlines[i].find('GradGrad') > -1:
					gradgrad_line = i

			# Determine charge and multiplicity
			if outlines[i].find("Charge = ") > -1:
				if self.args.charge_sp != None and self.args.sp != None:
					charge = self.args.charge_sp
				else:
					charge = int(outlines[i].split()[2])
				if self.args.mult_sp != None and self.args.sp != None:
					mult = self.args.mult_sp
				else:
					mult = int(outlines[i].split()[5].rstrip("\n"))
				
			# Get grid size
			elif outlines[i].strip().startswith('ExpMin='):
				IRadAn = int(outlines[i].strip().split()[-3])
				if IRadAn == 1:
					grid_size = 'sg1'
				elif IRadAn == 2:
					grid_size = 'coarse'
				elif IRadAn == 4:
					grid_size = 'fine'		
				elif IRadAn == 5:
					grid_size = 'ultrafine'	
				elif IRadAn == 7:
					grid_size = 'superfine'	

			if symbol_z_line != None and gradgrad_line != None:
				n_atoms = gradgrad_line-symbol_z_line-4
				
			if n_atoms != 0 and grid_size != '':
				break

		# find if the calculation was for a ground or transition state
		calc_type = 'ground_state'
		calcfc_found, ts_found = False, False
		for keyword in keywords_line.split():
			if keyword.lower().find('opt') > -1:
				if keyword.lower().find('calcfc') > -1:
					calcfc_found = True
				if keyword.lower().find('ts') > -1:
					ts_found = True
		
		if calcfc_found and ts_found:
			calc_type = 'transition_state'
		
		# substract the imaginary frequency from the TS in normal terminations
		if calc_type == 'transition_state':
			initial_ifreqs -= 1
		
		return n_atoms,charge,mult,keywords_line,initial_ifreqs,calc_type,grid_size


	def get_input_geom(self, file, outlines, keywords_line, n_atoms, mult, termination, errortype, initial_ifreqs, nimag_line, calc_type, file_list, E_dup, H_dup, G_dup, nocheck):
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
		errortype : string
			Error type of the QM output file
		keywords_line : string
			Line for the input file that will be generated
		atom_types : list of strings
			List containing the atoms of the system
		cartesians : list of lists
			Cartesian coordinates used for further processing
		molecular_mass: str
			Molecular weight of the molecule in amu read from the output QM file
		"""    

		# this parameter sets the line to retrieve the geometry for normally and not normally terminated calcs
		stand_or = 0
			
		if errortype == 'before_E_calculation':
			pass

		else:
			if termination == "normal" or nocheck:
				stand_or,errortype,freqs,freq_displacements,file_list,E_dup,H_dup,G_dup,molecular_mass = self.analysis_normal(file,outlines,stand_or,errortype,nimag_line,n_atoms,file_list,E_dup,H_dup,G_dup,initial_ifreqs,mult,calc_type)

			else:
				freqs, freq_displacements, molecular_mass = [],[],None
				stand_or,keywords_line = self.analysis_not_normal(outlines,stand_or,keywords_line)

			# get atom types from the calculation and the cartesian coordinates to use in each case
			atom_types, cartesians = [],[]
			
			per_tab = periodic_table()
			for i in range(stand_or+5,stand_or+5+n_atoms):
				massno = int(outlines[i].split()[1])
				if massno < len(per_tab):
					atom_symbol = per_tab[massno]
				else:
					atom_symbol = "XX"
				atom_types.append(atom_symbol)
				cartesians.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

			if errortype == 'extra_imag_freq':
				cartesians = self.fix_imag_freqs(n_atoms, cartesians, freqs, freq_displacements)
			
			# if a preset input line was defined in args.qm_input, this overwrites the automatic protocol
			if self.args.qm_input != '':
				if termination != "normal" or errortype == 'extra_imag_freq':
					keywords_line = self.args.qm_input
					# add the ts options if an external input line was defined
					if calc_type == 'transition_state' and self.args.ts_input != None:
						new_preset_keywords_line = '# '
						for keyword in self.args.qm_input.split():
							if keyword.lower().startswith('opt') > -1:
								if keyword.lower().find('ts') > -1:
									pass
								else:
									keyword = self.args.ts_input
							new_preset_keywords_line += keyword
							new_preset_keywords_line += ' '
						keywords_line = new_preset_keywords_line

			return file_list, E_dup, H_dup, G_dup, errortype, keywords_line, atom_types, cartesians, freqs, freq_displacements, molecular_mass


	def analysis_normal(self,file,outlines,stand_or,errortype,nimag_line,n_atoms,file_list,E_dup,H_dup,G_dup,initial_ifreqs,mult,calc_type):
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
		errortype : string
			Error type of the QM output file
		nimag_line : int
			Line number of the line containing the amount of imaginary frequencies in the 
			QM calculation. This is tracked to avoid reading the same lines multiple times
		n_atoms : int
			Number of atoms in the calculation
		file_list : list of str
			Collects file names for the duplicate analysis
		E_dup : list of float
			Collects (electronic energy + ZPE) values for the duplicate analysis
		H_dup : list of float
			Collects enthalpy values for the duplicate analysis
		G_dup : list of float
			Collects Gibbs free energy values for the duplicate analysis 
		initial_ifreqs : int
			Initial amount of imaginary frequencies in the QM file
		mult : int
			Multiplicity of the QM calculation
		calc_type : str
			Type of the QM calculation (ground_state or transition_state)

		Returns
		-------
		stand_or : int
			Line number containing the standard orientation of the target geometry
		errortype : string
			Error type of the QM output file
		freqs : list of float
			List containing the frequencies as floats
		freq_displacements : list of matrixes
			Contains the normal modes for each frequencies (including all atoms)
		file_list : list of str
			Collects file names for the duplicate analysis
		E_dup : list of float
			Collects (electronic energy + ZPE) values for the duplicate analysis
		H_dup : list of float
			Collects enthalpy values for the duplicate analysis
		G_dup : list of float
			Collects Gibbs free energy values for the duplicate analysis 
		molecular_mass: str
			Molecular weight of the molecule in amu read from the output QM file
		"""

		ifreqs, molecular_mass, freqs_so_far, start_dup = 0, 0, 0, False
		spin_contamination = False
		freqs, freq_displacements  = [],[]

		# reverse loop to speed up the reading of the output files
		# (also skips final part for normal termination which was already read plus dipole/multipole information)
		for i in reversed(range(0,nimag_line-(10*n_atoms))):
			# Sets where the final coordinates are inside the file
			if (outlines[i].find("Standard orientation:") > -1 or outlines[i].find("Input orientation:") > -1):
				stand_or = i
				break

			if not self.args.nocheck:
				# analyze duplicates using thermodata
				if self.args.dup:
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

				# Get frequencies and their modes
				if outlines[i].find(" Frequencies -- ") > -1:
					nfreqs = len(outlines[i].split())
					for j in range(2, nfreqs):
						if float(outlines[i].split()[j]) < 0 and abs(float(outlines[i].split()[j])) > abs(self.args.ifreq_cutoff):
							ifreqs += 1
						freqs.append(float(outlines[i].split()[j]))
						freq_displacements.append([])
					for j in range(0,n_atoms):
						for k in range(0, nfreqs-2):
							freq_displacements[(freqs_so_far + k)].append([float(outlines[i+5+j].split()[3*k+2]), float(outlines[i+5+j].split()[3*k+3]), float(outlines[i+5+j].split()[3*k+4])])
					freqs_so_far = freqs_so_far + nfreqs - 2
				
				# Get molecular mass
				elif outlines[i].strip().startswith('Molecular mass:'):
					molecular_mass = float(outlines[i].strip().split()[2])																			

				# analyze spin contamination
				elif self.args.s2_threshold > 0.0:
					if outlines[i].find('S**2 before annihilation') > -1:
						s2_value = outlines[i].split()[-1].rstrip("\n")
						unpaired_e = mult-1
						spin = unpaired_e*0.5
						s2_expected_value = spin*(spin+1)
						spin_diff = abs(float(s2_value)-s2_expected_value)
						if spin_diff > abs(self.args.s2_threshold/100)*s2_expected_value:
							spin_contamination = True

		# exclude TS imag frequency if needed
		if calc_type == 'transition_state':
			ifreqs -= 1

		if ifreqs > 0:
			errortype = 'extra_imag_freq'

		elif ifreqs < 0:
			errortype = 'ts_no_imag_freq'

		if spin_contamination:
			errortype = 'spin_contaminated'
		
		# sorts the frequencies as they are shown in Gaussian (correction from using a reversed read lines loop)
		zipped_lists = zip(freqs, freq_displacements)
		sorted_pairs = sorted(zipped_lists, key=lambda x: x[0])

		tuples = zip(*sorted_pairs)
		freqs, freq_displacements = [ list(tuple) for tuple in  tuples]
	
		return stand_or,errortype,freqs,freq_displacements,file_list,E_dup,H_dup,G_dup,molecular_mass


	def analysis_not_normal(self,outlines,stand_or,keywords_line):
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
		if not self.args.nocheck:
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


	def fix_imag_freqs(self, n_atoms, cartesians, freqs, freq_displacements):
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

		Returns
		-------
		cartesians : list of lists
			New set of cartesian coordinates generated after displacing the original
			coordinates along the normal modes of the corresponding imaginary frequencies
		"""

		shift = []

		# could get rid of atomic units here, if zpe_rat definition is changed
		for mode,_ in enumerate(freqs):
			# Either moves along all imaginary freqs, or a specific mode requested by the user
			if freqs[mode] < 0.0:
				shift.append(self.args.amplitude_ifreq)
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
		file : string
			Filename
		termination : string
			Type of termination of the QM output file (i.e. normal, error, unfinished)
		errortype : string
			Type of error type of the QM output file (i.e. None, unknown, extra_imag_freq, etc)
		file_terms : dictionary
			Keeps track of the number of calculations for each termination and error type
		"""
		
		if errortype == None and termination == "normal":
			destination = self.w_dir_main+'/successful_QM_outputs'
			file_terms['finished'] += 1

		elif errortype == 'extra_imag_freq':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/imag_freq/'
			file_terms['imag_freq'] += 1

		elif errortype == 'ts_no_imag_freq':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/ts_no_imag_freq/'
			file_terms['ts_no_imag_freq'] += 1

		elif errortype == 'spin_contaminated':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/spin_contaminated/'
			file_terms['spin_contaminated'] += 1

		elif errortype == 'duplicate_calc':
			destination = self.w_dir_main+'/duplicates/run_'+str(self.round_num)
			file_terms['duplicate_calc'] += 1

		elif termination == "error":
			if errortype == "atomicbasiserror":
				destination = self.w_dir_main +'/failed/run_'+str(self.round_num)+'/error/basis_set_error'
				file_terms['atom_error'] += 1
			elif errortype == "SCFerror":
				destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/error/scf_error'
				file_terms['scf_error'] += 1
			elif errortype == "before_E_calculation":
				destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/error/before_E_calculation'
				file_terms['before_E_error'] += 1
			else:
				destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/error/unknown_error'
				file_terms['other_error'] += 1

		elif termination == "unfinished":
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/unfinished/'
			file_terms['unfinished'] += 1

		elif errortype == 'fail_geom_rules':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/geom_rules_filter/'
			file_terms['geom_rules_qcorr'] += 1

		elif errortype == 'isomerization':
			destination = self.w_dir_main+'/failed/run_'+str(self.round_num)+'/isomerization/'
			file_terms['check_geom_qcorr'] += 1
		
		if not os.path.isdir(destination):
			os.makedirs(destination)
		# move_file(file, self.w_dir_main, destination)

		return file_terms

					
	def write_qcorr_csv(self,file_terms):
		"""
		Write information about the QCORR analysis in a csv
		"""
		ana_data = pd.DataFrame()
		ana_data.at[0,'Total files'] = len(self.qm_files)
		ana_data.at[0,'Normal termination'] = file_terms['finished']
		ana_data.at[0,'Imaginary frequencies'] = file_terms['imag_freq']
		ana_data.at[0,'TS with no imag. freq.'] = file_terms['ts_no_imag_freq']
		ana_data.at[0,'SCF error'] = file_terms['scf_error']
		ana_data.at[0,'Error before SCF'] = file_terms['before_E_error']
		ana_data.at[0,'Basis set error'] =  file_terms['atom_error']
		ana_data.at[0,'Other errors'] = file_terms['other_error']
		ana_data.at[0,'Unfinished'] = file_terms['unfinished']
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
