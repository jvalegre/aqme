#!/usr/bin/env python

#####################################################.
# 		   This file stores all the main functions 	#
#####################################################.

import glob
import sys
import os
import time
import pandas as pd
import subprocess
import yaml
import concurrent.futures as futures
from progress.bar import IncrementalBar
from rdkit.Chem import AllChem as Chem
from pyconfort.csearch import check_for_pieces, check_charge_smi, clean_args, compute_confs, com_2_xyz_2_sdf, mol_from_sdf_or_mol_or_mol2,creation_of_dup_csv,Logger
from pyconfort.filter import exp_rules_output

from pyconfort.qprep_gaussian import read_energies, write_gaussian_input_file, moving_files, convert_xyz_to_sdf, get_name_and_charge
from pyconfort.qcorr_gaussian import output_analyzer, check_for_final_folder, dup_calculation

from pyconfort.grapher import graph
from pyconfort.descp import calculate_parameters
from pyconfort.nmr import calculate_boltz_and_nmr
from pyconfort.energy import calculate_boltz_and_energy,calculate_avg_and_energy
from pyconfort.dbstep_conf import calculate_db_parameters,calculate_boltz_and_dbstep
from pyconfort.nics_conf import calculate_nics_parameters,calculate_boltz_for_nics,calculate_avg_nics
from pyconfort.cclib_conf import calculate_cclib,calcualte_average_cclib_parameter,calculate_boltz_for_cclib
from pyconfort.cmin import mult_min
#need to and in energy

#load paramters from yaml file
def load_from_yaml(args,log):
	# Variables will be updated from YAML file
	try:
		if args.varfile is not None:
			if os.path.exists(args.varfile):
				if os.path.splitext(args.varfile)[1] == '.yaml':
					log.write("\no  Importing pyCONFORT parameters from " + args.varfile)
					with open(args.varfile, 'r') as file:
						try:
							param_list = yaml.load(file, Loader=yaml.SafeLoader)
						except yaml.scanner.ScannerError:
							log.write("\nx  Error while reading " + args.varfile+ ". Edit the yaml file and try again (i.e. use ':' instead of '=' to specify variables)")
							sys.exit()
			for param in param_list:
				if hasattr(args, param):
					if getattr(args, param) != param_list[param]:
						log.write("o  RESET " + param + " from " + str(getattr(args, param)) + " to " + str(param_list[param]))
						setattr(args, param, param_list[param])
					else:
						log.write("o  DEFAULT " + param + " : " + str(getattr(args, param)))
	except UnboundLocalError:
		log.write("\no  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code.\n")


#creation of csv for QCORR
def creation_of_ana_csv(args,duplicates):

	columns_list = ['Total files', 'Normal termination', 'Imaginary frequencies', 'SCF error', 'Basis set error', 'Other errors', 'Unfinished']
	if args.dup:
		columns_list.append('Duplicates')
	if len(args.exp_rules) >= 1:
		columns_list.append('Exp_rules filter')
	if args.check_geom:
		columns_list.append('Geometry changed')
	ana_data = pd.DataFrame(columns = columns_list)

	return ana_data

# main function to generate conformers
def csearch_main(w_dir_initial,args,log_overall):
	# input file format specified
	file_format = os.path.splitext(args.input)[1]

	if file_format not in ['.smi', '.sdf', '.cdx', '.csv','.com','.gjf','.mol','.mol2','.xyz','.txt','.yaml','.yml','.rtf']:
		log_overall.write("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!")
		sys.exit()

	if not os.path.exists(args.input):
		log_overall.write("\nx  INPUT FILE NOT FOUND!")
		sys.exit()

	# sets up the chosen force field (this fixes some problems in case MMFF is replaced by UFF)
	ori_ff = args.ff
	ori_charge = args.charge_default

	# if large system increase stack size
	# if args.STACKSIZE != '1G':
	# 	os.environ['OMP_STACKSIZE'] = args.STACKSIZE

	with futures.ProcessPoolExecutor(max_workers=args.cpus) as executor:
		# Submit a set of asynchronous jobs
		jobs = []
		count_mol = 0
		# SMILES input specified
		smi_derivatives = ['.smi', '.txt', '.yaml', '.yml', '.rtf']
		if file_format in smi_derivatives:
			smifile = open(args.input)
			for i, line in enumerate(smifile):
				if line == '\n':
					pass
				else:
					toks = line.split()
					#editing part
					smi = toks[0]
					smi = check_for_pieces(smi)
					mol = Chem.MolFromSmiles(smi)
					try:
						clean_args(args,ori_ff,mol,ori_charge)
						if args.charge_default == 'auto':
							if not args.metal_complex:
								args.charge_default = check_charge_smi(smi)
						if args.prefix == 'None':
							name = ''.join(toks[1:])
						else:
							name = str(args.prefix)+str(i)+'_'+''.join(toks[1:])

						job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
						jobs.append(job)
						count_mol +=1

						# compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)
					except AttributeError:
						log_overall.write("\nx  Wrong SMILES string ("+smi+") found (not compatible with RDKit or ANI/xTB if selected)! This compound will be omitted\n")


		# CSV file with one columns SMILES and code_name
		elif os.path.splitext(args.input)[1] == '.csv':
			csv_smiles = pd.read_csv(args.input)
			for i in range(len(csv_smiles)):
				#assigning names and smi i  each loop
				smi = csv_smiles.loc[i, 'SMILES']
				smi = check_for_pieces(smi)
				mol = Chem.MolFromSmiles(smi)
				clean_args(args,ori_ff,mol,ori_charge)
				if args.charge_default == 'auto':
					if not args.metal_complex:
						args.charge_default = check_charge_smi(smi)
				if args.prefix == 'None':
					name = csv_smiles.loc[i, 'code_name']
				else:
					name = 'comp_'+str(i)+'_'+csv_smiles.loc[i, 'code_name']
				job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
				jobs.append(job)
				count_mol +=1
				# compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)

		# CDX file
		elif os.path.splitext(args.input)[1] == '.cdx':
			#converting to smiles from chemdraw
			cmd_cdx = ['obabel', '-icdx', args.input, '-osmi', '-Ocdx.smi']
			subprocess.call(cmd_cdx)
			smifile = open('cdx.smi',"r")

			for i, smi in enumerate(smifile):
				smi = check_for_pieces(smi)
				mol = Chem.MolFromSmiles(smi)
				clean_args(args,ori_ff,mol,ori_charge)
				if args.charge_default == 'auto':
					if not args.metal_complex:
						args.charge_default = check_charge_smi(smi)
				name = args.input.split('.')[0] +'_'+ str(i)
				job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
				jobs.append(job)
				count_mol +=1

			smifile.close()
			os.remove('cdx.smi')

		# COM file
		elif os.path.splitext(args.input)[1] == '.gjf' or os.path.splitext(args.input)[1] == '.com' or  os.path.splitext(args.input)[1] == '.xyz' :
			#converting to sdf from comfile to preserve geometry
			charge_com = com_2_xyz_2_sdf(args)
			sdffile = os.path.splitext(args.input)[0]+'.sdf'
			if os.path.splitext(args.input)[1] == '.gjf' or os.path.splitext(args.input)[1] == '.com':
				suppl, _, _ = mol_from_sdf_or_mol_or_mol2(sdffile)
			if  os.path.splitext(args.input)[1] == '.xyz':
				suppl, _, charge_com_list = mol_from_sdf_or_mol_or_mol2(sdffile)
				charge_com = charge_com_list[0]
			name = os.path.splitext(args.input)[0]
			counter_for_template = 0
			i=0
			for mol in suppl:
				clean_args(args,ori_ff,mol,ori_charge)
				if args.charge_default == 'auto':
					args.charge_default = charge_com
				job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
				jobs.append(job)
				count_mol +=1
				i += 1

		# SDF file
		elif os.path.splitext(args.input)[1] == '.sdf' or os.path.splitext(args.input)[1] == '.mol' or os.path.splitext(args.input)[1] == '.mol2':
			suppl, IDs, charges = mol_from_sdf_or_mol_or_mol2(args.input)
			i=0
			if os.path.splitext(args.input)[1] == '.sdf':
				for mol,name,charge_sdf in zip(suppl,IDs,charges):
					clean_args(args,ori_ff,mol,ori_charge)
					if args.charge_default == 'auto':
						args.charge_default = charge_sdf
					job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
					jobs.append(job)
					count_mol +=1
					i += 1
			elif os.path.splitext(args.input)[1] == '.mol' or os.path.splitext(args.input)[1] == '.mol2':
				if args.charge_default == 'auto':
					args.charge_default = charges[0]
				name = IDs[0]
				mol = suppl
				job = executor.submit(compute_confs,w_dir_initial,mol,name,args,i)
				jobs.append(job)
				count_mol +=1
				i += 1

		final_dup_data = creation_of_dup_csv(args)
		bar = IncrementalBar('o  Number of finished jobs from CSEARCH', max = count_mol)
		# Process the job results (in submission order) and save the conformers.
		for i,job in enumerate(jobs):
			total_data = job.result()
			frames = [final_dup_data, total_data]
			final_dup_data = pd.concat(frames,ignore_index=True,sort=True)
			bar.next()
		bar.finish()

		# removing temporary files
		temp_files = ['gfn2.out', 'xTB_opt.traj', 'ANI1_opt.traj', 'wbo', 'xtbrestart','ase.opt','xtb.opt','gfnff_topo']
		for file in temp_files:
			if os.path.exists(file):
				os.remove(file)

	return final_dup_data


def cmin_main(w_dir_initial,args,log_overall,dup_data):
	try:
		os.makedirs(w_dir_initial+'/CMIN/dat_files')
	except OSError:
		if os.path.isdir(w_dir_initial+'/CMIN/dat_files'):
			pass
	bar = IncrementalBar('o  Number of finished jobs from CMIN', max = len(dup_data))
	for dup_data_idx in range(len(dup_data)):
		start_time = time.time()
		update_to_rdkit = dup_data.at[dup_data_idx,'update_to_rdkit']
		name = dup_data.at[dup_data_idx,'Molecule']
		log = Logger(w_dir_initial+'/CMIN/dat_files/'+name, args.output_name)
		if dup_data.at[dup_data_idx,'status'] != -1:
			try:
				if args.CMIN=='ani':
					min_suffix = 'ani'
				elif args.CMIN=='xtb':
					min_suffix = 'xtb'
				if args.CSEARCH=='rdkit':
					mult_min(name+'_'+'rdkit', args, min_suffix, log, dup_data, dup_data_idx)
				elif args.CSEARCH=='summ' and not update_to_rdkit :
					mult_min(name+'_'+'summ', args, min_suffix, log, dup_data, dup_data_idx)
				elif args.CSEARCH=='summ' and update_to_rdkit :
					mult_min(name+'_'+'summ', args, min_suffix, log, dup_data, dup_data_idx)
				elif args.CSEARCH=='fullmonte' and not update_to_rdkit:
					mult_min(name+'_'+'fullmonte', args, min_suffix, log, dup_data, dup_data_idx)
				elif args.CSEARCH=='fullmonte' and update_to_rdkit:
					mult_min(name+'_'+'fullmonte', args, min_suffix, log, dup_data, dup_data_idx)
				dup_data.at[dup_data_idx, 'CMIN time (seconds)'] = round(time.time() - start_time,2)
			except:
				pass
		bar.next()
	bar.finish()

	# removing temporary files
	temp_files = ['gfn2.out', 'xTB_opt.traj', 'ANI1_opt.traj', 'wbo', 'xtbrestart','ase.opt','xtb.opt','gfnff_topo']
	for file in temp_files:
		if os.path.exists(file):
			os.remove(file)

	return dup_data

#writing gauss main
def qprep_main(w_dir_initial,args,log):

	if len(args.exp_rules) >= 1:
		conf_files =  glob.glob('*_rules.sdf')
	# define the SDF files to convert to COM Gaussian files
	elif not args.CMIN=='xtb' and not args.CMIN=='ani' and args.CSEARCH=='rdkit':
		conf_files =  glob.glob('*_rdkit.sdf')
	elif not args.CMIN=='xtb' and not args.CMIN=='ani' and args.CSEARCH=='summ':
		conf_files =  glob.glob('*_summ.sdf')
	elif not args.CMIN=='xtb' and not args.CMIN=='ani' and args.CSEARCH=='fullmonte':
		conf_files =  glob.glob('*_fullmonte.sdf')
	elif args.CMIN=='xtb' and not args.CMIN=='ani':
		conf_files =  glob.glob('*_xtb.sdf')
	elif args.CMIN=='ani' and not args.CMIN=='xtb':
		conf_files =  glob.glob('*_ani.sdf')
	elif args.CMIN=='ani' and args.CMIN=='xtb':
		conf_files =  glob.glob('*_ani.sdf') + glob.glob('*_xtb.sdf')
	else:
		conf_files =  glob.glob('*.sdf')

	if args.com_from_xyz:
		xyz_files =  glob.glob('*.xyz')
		convert_xyz_to_sdf(xyz_files,args,log)
		conf_files =  glob.glob('*.sdf')

	# names for directories created
	if args.QPREP == 'gaussian':
		g_dir = 'QMCALC/G16'
	elif args.QPREP == 'orca':
		g_dir = 'QMCALC/ORCA'

	if len(conf_files) != 0:
		#read in dup_data to get the overall charge of MOLECULES
		invalid_files = []
		try:
			charge_data = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv', usecols=['Molecule','Overall charge'])
		except:
			charge_data = pd.DataFrame()
			for i,conf_file in enumerate(conf_files):
				try:
					suppl = Chem.SDMolSupplier(conf_file, removeHs=False)
				except OSError:
					suppl = False
				if suppl:
					mol = suppl[0]
					charge = mol.GetProp('Real charge')
					charge_data.at[i,'Overall charge'] = charge
				else:
					invalid_files.append(conf_file)
					charge_data.at[i,'Overall charge'] = 'Invalid'

				name = get_name_and_charge(conf_file.split('.')[0],None)
				charge_data.at[i,'Molecule'] = name

		for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
			# only create this directory if single point calculation is requested
			if bs.find('**') > -1:
				bs = bs.replace('**','(d,p)')
			elif bs.find('*') > -1:
				bs = bs.replace('*','(d)')

			if str(bs).find('/') > -1:
				folder = g_dir + '/' + str(lot) + '-' + str(bs).split('/')[0]
			else:
				folder = g_dir + '/' + str(lot) + '-' + str(bs)
			log.write("\no  Preparing QM input files in {}".format(folder))

			# this variable keeps trakc of folder creation
			folder_error = False
			try:
				os.makedirs(folder)
			except OSError:
				if os.path.isdir(folder):
					pass
				else:
					log.write('\nx  The QMCALC folder from QPREP could not be created, probably due to incompatible characters')
					folder_error = True

			# writing the com files
			# check conf_file exists, parse energies and then write DFT input
			for file in conf_files:
				if file not in invalid_files:
					if os.path.exists(file):
						if args.verbose:
							if not folder_error:
								log.write("   -> Converting from {}".format(file))
						energies = read_energies(file,log)
						name = file.split('.')[0]

						write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args, log, charge_data, w_dir_initial)
	else:
		log.write('\nx  No SDF files detected to convert to gaussian COM files')

# moving files after compute and/or write_gauss
def move_sdf_main(args):
	src = os.getcwd()
	if len(args.exp_rules) >= 1:
		exp_rules_files = glob.glob('*_filter_exp_rules.sdf')
	if args.CMIN=='xtb':
		all_xtb_conf_files = glob.glob('*_xtb.sdf')
		destination_xtb = src +'/CMIN/xtb'
		for file in all_xtb_conf_files:
			moving_files(destination_xtb,src,file)
		all_xtb_conf_files_all = glob.glob('*_xtb_all_confs.sdf')
		destination_xtb_all = src +'/CMIN/xtb_all_confs'
		for file in all_xtb_conf_files_all:
			moving_files(destination_xtb_all,src,file)
		if len(args.exp_rules) >= 1:
			destination_exp_rules = src +'/CMIN/xtb/filter_exp_rules/'
			for file in exp_rules_files:
				moving_files(destination_exp_rules,src,file)
	if args.CMIN=='ani':
		all_ani_conf_files = glob.glob('*_ani.sdf')
		destination_ani = src +'/CMIN/ani'
		for file in all_ani_conf_files:
			moving_files(destination_ani,src,file)
		all_ani_conf_files_all = glob.glob('*_ani_all_confs.sdf')
		destination_ani_all = src +'/CMIN/ani_all_confs'
		for file in all_ani_conf_files_all:
			moving_files(destination_ani_all,src,file)
		if len(args.exp_rules) >= 1:
			destination_exp_rules = src +'/CMIN/ani/filter_exp_rules/'
			for file in exp_rules_files:
				moving_files(destination_exp_rules,src,file)
	if args.CSEARCH=='rdkit':
		all_name_conf_files = glob.glob('*_rdkit.sdf')
		destination_rdkit = src+ '/CSEARCH/rdkit'
		for file in all_name_conf_files:
			moving_files(destination_rdkit,src,file)
		if len(args.exp_rules) >= 1:
			destination_exp_rules = src +'/CSEARCH/rdkit/filter_exp_rules/'
			for file in exp_rules_files:
				moving_files(destination_exp_rules,src,file)

	if args.CSEARCH=='summ':
		all_name_conf_files = glob.glob('*_summ.sdf')
		destination_rdkit = src+ '/CSEARCH/summ'
		for file in all_name_conf_files:
			moving_files(destination_rdkit,src,file)
		if len(args.exp_rules) >= 1 and args.CMIN is None:
			destination_exp_rules = src +'/CSEARCH/summ/filter_exp_rules/'
			for file in exp_rules_files:
				moving_files(destination_exp_rules,src,file)

	if args.CSEARCH=='fullmonte':
		all_name_conf_files = glob.glob('*_fullmonte.sdf')
		destination_rdkit = src+ '/CSEARCH/fullmonte'
		for file in all_name_conf_files:
			moving_files(destination_rdkit,src,file)
		if len(args.exp_rules) >= 1 and args.CMIN is None:
			destination_exp_rules = src +'/CSEARCH/fullmonte/filter_exp_rules/'
			for file in exp_rules_files:
				moving_files(destination_exp_rules,src,file)

	if args.CSEARCH is None:
		if len(args.exp_rules) >= 1:
			destination_exp_rules = src +'/QMCALC/SDF_input/filter_exp_rules/'
			for file in exp_rules_files:
				moving_files(destination_exp_rules,src,file)

		all_conf_files = glob.glob('*.sdf')
		destination = src +'/QMCALC/SDF_input/'
		for file in all_conf_files:
			moving_files(destination,src,file)


	if args.com_from_xyz:
		all_xyz_conf_files = glob.glob('*.xyz')+glob.glob('*.sdf')
		destination_xyz = 'QMCALC/xyz_and_sdf'
		for file in all_xyz_conf_files:
			moving_files(destination_xyz,src,file)

#finding the file type to move for analysis
def get_com_or_log_out_files(type,name):
	files = []
	if type =='output':
		formats = ['*.log','*.LOG','*.out','*.OUT','*json']
	elif type =='input':
		formats =['*.com','*.gjf']
	for _,format in enumerate(formats):
		if name is None:
			all_files = enumerate(glob.glob(format))
		else:
			all_files = enumerate(glob.glob(name+format))
		for _,file in all_files:
			if file not in files:
				files.append(file)
	return files

# main part of the analysis functions
def qcorr_gaussian_main(duplicates,w_dir_initial,args,log):
	# when you run analysis in a folder full of output files
	if not os.path.exists(w_dir_initial+'/QMCALC'):
		w_dir_main = os.getcwd()
		w_dir_fin = w_dir_main+'/success/output_files'
		for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
			if not os.path.isdir(w_dir_main+'/dat_files/'):
				os.makedirs(w_dir_main+'/dat_files/')
			w_dir,round_num = check_for_final_folder(w_dir_main)
			os.chdir(w_dir)
			log = Logger(w_dir_main+'/dat_files/pyCONFORT-QCORR-run_'+str(round_num), args.output_name)
			ana_data = creation_of_ana_csv(args,duplicates)
			log.write("\no  Analyzing output files in {}\n".format(w_dir))
			log_files = get_com_or_log_out_files('output',None)
			if len(log_files) == 0:
				log.write('x  There are no output files in this folder.')
			com_files = get_com_or_log_out_files('input',None)
			output_analyzer(duplicates,log_files, com_files, w_dir, w_dir_main, lot, bs, bs_gcp, args, w_dir_fin, w_dir_initial, log, ana_data, round_num)
		os.chdir(w_dir_main)
	# when you specify multiple levels of theory
	else:
		args.path = w_dir_initial+'/QMCALC/G16/'
		# Sets the folder and find the log files to analyze
		for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
			#for main folder
			if str(bs).find('/') > -1:
				w_dir_main = args.path + str(lot) + '-' + str(bs).split('/')[0]
			else:
				w_dir_main = args.path + str(lot) + '-' + str(bs)
			if not os.path.isdir(w_dir_main+'/dat_files/'):
				os.makedirs(w_dir_main+'/dat_files/')
			#check if New_Gaussian_Input_Files folder exists
			w_dir,round_num = check_for_final_folder(w_dir_main)
			log = Logger(w_dir_main+'/dat_files/pyCONFORT-QCORR-run_'+str(round_num), args.output_name)
			#assign the path to the finished directory.
			if str(bs).find('/') > -1:
				w_dir_fin = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files'
			else:
				w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
			os.chdir(w_dir)
			ana_data = creation_of_ana_csv(args,duplicates)
			log.write("\no  Analyzing output files in {}\n".format(w_dir))
			log_files = get_com_or_log_out_files('output',None)
			com_files = get_com_or_log_out_files('input',None)
			output_analyzer(duplicates,log_files, com_files, w_dir, w_dir_main , lot, bs, bs_gcp, args, w_dir_fin, w_dir_initial, log, ana_data, round_num)
		os.chdir(args.path)
	os.chdir(w_dir_initial)

#removing the duplicates
def dup_main(args,log,w_dir_initial):
	if not os.path.exists(w_dir_initial+'/QMCALC'):
		w_dir_main = os.getcwd()
		w_dir,round_num = check_for_final_folder(w_dir_main)
		os.chdir(w_dir)
		log_files = get_com_or_log_out_files('output',None)
		if len(log_files) != 0:
			duplicates = dup_calculation(log_files,w_dir,w_dir_main,args,log,round_num)
		else:
			log.write('x  There are no output files in this folder.')
			duplicates = 'None'
		os.chdir(w_dir_main)
	else:
		if args.QCORR=='gaussian':
			args.path = w_dir_initial+'/QMCALC/G16/'
		# Sets the folder and find the log files to analyze
		for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
			if str(bs).find('/') > -1:
				w_dir_main = args.path + str(lot) + '-' + str(bs).split('/')[0]
			else:
				w_dir_main = args.path + str(lot) + '-' + str(bs)
			os.chdir(w_dir_main)
			w_dir,round_num = check_for_final_folder(w_dir_main)
			os.chdir(w_dir)
			log_files = get_com_or_log_out_files('output',None)
			# change molecules to a range as files will have codes in a continous manner
			if len(log_files) != 0:
				duplicates = dup_calculation(log_files,w_dir,w_dir_main,args,log,round_num)
			else:
				log.write('x  There are no output files in this folder.')
				duplicates = 'None'

	return duplicates

#getting descriptors
def geom_par_main(args,log,w_dir_initial):
	#get sdf FILES from csv
	pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no  Calculating paramters for molecule : {0} ".format(name))

		sdf_ani,sdf_xtb = None,None
		if os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name+'_rdkit.sdf'):
				sdf_rdkit =  w_dir_initial+'/CSEARCH/rdkit/'+name+'_rdkit.sdf'
		elif os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name+'_summ.sdf'):
			sdf_rdkit = w_dir_initial+'/CSEARCH/summ/'+name+'_summ.sdf'
		elif os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name+'_fullmonte.sdf'):
			sdf_rdkit = w_dir_initial+'/CSEARCH/fullmonte/'+name+'_fullmonte.sdf'
		if os.path.exists(w_dir_initial+'/CMIN/xtb/'+name+'_xtb.sdf'):
			sdf_xtb =  w_dir_initial+'/CMIN/xtb/'+name+'_xtb.sdf'
		if os.path.exists(w_dir_initial+'/CMIN/ani/'+name+'_ani.sdf'):
			sdf_ani = w_dir_initial+'/CMIN/ani/'+name+'_ani.sdf'
		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
			# Sets the folder and find the log files to analyze
			for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
				#assign the path to the finished directory.
				if str(bs).find('/') > -1:
					w_dir = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files'
				else:
					w_dir = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
				os.chdir(w_dir)
				log_files = get_com_or_log_out_files('output',name)
				calculate_parameters(sdf_rdkit,sdf_ani,sdf_xtb,log_files,args,log,w_dir_initial,name,lot,bs)
		else:
			calculate_parameters(sdf_rdkit,sdf_ani,sdf_xtb,None,args,log,w_dir_initial,name,None,None)

		os.chdir(w_dir_initial)

#function to plot graphs
def graph_main(args,log,w_dir_initial):
	#get sdf FILES from csv
	pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no  Plotting graphs for molecule : {0} ".format(name))

		sdf_ani,sdf_xtb = None,None
		if os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name+'_rdkit.sdf'):
			sdf_rdkit =  w_dir_initial+'/CSEARCH/rdkit/'+name+'_rdkit.sdf'
		elif os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name+'_summ.sdf'):
			sdf_rdkit = w_dir_initial+'/CSEARCH/summ/'+name+'_summ.sdf'
		elif os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name+'_fullmonte.sdf'):
			sdf_rdkit = w_dir_initial+'/CSEARCH/fullmonte/'+name+'_fullmonte.sdf'
		if os.path.exists(w_dir_initial+'/CMIN/xtb_all_confs/'+name+'_xtb_all_confs.sdf'):
			sdf_xtb =  w_dir_initial+'/CMIN/xtb_all_confs/'+name+'_xtb_all_confs.sdf'
		if os.path.exists(w_dir_initial+'/CMIN/ani_all_confs/'+name+'_ani_all_confs.sdf'):
			sdf_ani = w_dir_initial+'/CMIN/ani_all_confs/'+name+'_ani_all_confs.sdf'
		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
			# Sets the folder and find the log files to analyze
			for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
				#assign the path to the finished directory.
				w_dir = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
				os.chdir(w_dir)
				log_files = get_com_or_log_out_files('output',name)
				if os.path.exists(args.path + str(lot) + '-' + str(bs) +'/success/G16-SP_input_files'):
					for lot_sp,bs_sp,bs_gcp_sp in zip(args.level_of_theory_sp,args.basis_set_sp,args.basis_set_genecp_atoms_sp):
						w_dir_sp = args.path + str(lot) + '-' + str(bs) +'/success/G16-SP_input_files'+'/'+str(lot_sp)+'-'+str(bs_sp)
						sp_files = get_com_or_log_out_files('output',name)
						graph(sdf_rdkit,sdf_xtb,sdf_ani,log_files,sp_files,args,log,lot,bs,lot_sp,bs_sp,name,w_dir_initial,w_dir_sp,w_dir,'g16')
				if os.path.exists(args.path + str(lot) + '-' + str(bs) +'/success/ORCA-SP_input_files'):
					for lot_sp,bs_sp,bs_gcp_sp in zip(args.level_of_theory_sp,args.basis_set_sp,args.basis_set_genecp_atoms_sp):
						w_dir_sp = args.path + str(lot) + '-' + str(bs) +'/success/ORCA-SP_input_files'+'/'+str(lot_sp)+'-'+str(bs_sp.split('/')[0])
						os.chdir(w_dir_sp)
						sp_files = get_com_or_log_out_files('output',name)
						os.chdir(w_dir)
						graph(sdf_rdkit,sdf_xtb,sdf_ani,log_files,sp_files,args,log,lot,bs,lot_sp,bs_sp.split('/')[0],name,w_dir_initial,w_dir_sp,w_dir,'orca')
				else:
					graph(sdf_rdkit,sdf_xtb,sdf_ani,log_files,None,args,log,lot,bs,None,None,name,w_dir_initial,None,w_dir,None)

		else:
			graph(sdf_rdkit,sdf_xtb,sdf_ani,None,None,args,log,None,None,None,None,name,w_dir_initial,None,None,None)

	os.chdir(w_dir_initial)


#function for comparison of nmr
def nmr_main(args,log,w_dir_initial):

	if os.path.exists(w_dir_initial+'/QMCALC/G16'):
		args.path = w_dir_initial+'/QMCALC/G16/'
	else:
		log.write('\nx  The path for NMR analysis was not set up properly! (check the tutorials for more information)')

	#get sdf FILES from csv
	try:
		pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	except FileNotFoundError:
		# detects all the unique molecules from the success folder
		if str(args.basis_set[0]).find('/') > -1:
			w_dir_fin = args.path + str(args.level_of_theory[0]) + '-' + str(args.basis_set[0]).split('/')[0] +'/success/output_files'
		else:
			w_dir_fin = args.path + str(args.level_of_theory[0]) + '-' + str(args.basis_set[0]) +'/success/output_files'
		os.chdir(w_dir_fin)

		nmr_list = []
		standard_suffixes = ['xtb','ani']
		for name_nmr in glob.glob('*.*'):
			# discard_charact keeps track of the extra characters after the name of the molecule
			discard_charact = len(name_nmr.split('.')[1])+1
			potential_unique = name_nmr.split('.')[0].split('_')
			for i in reversed(range(len(potential_unique))):
				try:
					if potential_unique[i] not in standard_suffixes:
						 int(potential_unique[i])
					discard_charact += len(potential_unique[i])
					discard_charact += 1
				except ValueError:
						original_name = name_nmr[:-discard_charact]
						if original_name not in nmr_list:
							nmr_list.append(original_name)
		pd_name = pd.DataFrame(data=nmr_list, columns=['Molecule'])

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no NMR analysis for molecule : {0} ".format(name))

		# Sets the folder and find the log files to analyze
		for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
			#assign the path to the finished directory.
			if str(bs).find('/') > -1:
				w_dir_fin = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files'
			else:
				w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
			os.chdir(w_dir_fin)

			log_files = get_com_or_log_out_files('output',name)
			if len(log_files) != 0:
				calculate_boltz_and_nmr(log_files,args,log,name,w_dir_fin,w_dir_initial,lot,bs)
	os.chdir(w_dir_initial)

def energy_main(args,log,w_dir_initial):
	#get sdf FILES from csv
	pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no Boltzmann average energy analysis for molecule : {0} ".format(name))
		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
		# Sets the folder and find the log files to analyze
		for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
			#assign the path to the finished directory.
			if str(bs).find('/') > -1:
				w_dir_fin = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files/'
			else:
				w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/success/output_files/'
			os.chdir(w_dir_fin)
			log_files = get_com_or_log_out_files('output',name)
			if len(log_files) != 0:
				calculate_boltz_and_energy(log_files,args,log,name,w_dir_fin,w_dir_initial,lot,bs)

	#combining the combining all files in different folders
	w_dir_boltz = w_dir_initial+'/QPRED/energy/boltz/'

	for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
		#assign the path to the finished directory.
		if str(bs).find('/') > -1:
			w_dir_fin = w_dir_boltz + str(lot) + '-' + str(bs).split('/')[0]
		else:
			w_dir_fin = w_dir_boltz + str(lot) + '-' + str(bs)
		os.chdir(w_dir_fin)
		dat_files = glob.glob('*.dat')
		if len(dat_files) != 0:
			calculate_avg_and_energy(dat_files,args,log,name,w_dir_fin,w_dir_initial,w_dir_boltz,lot,bs)

	os.chdir(w_dir_initial)


def dbstep_par_main(args,log,w_dir_initial):
	#get sdf FILES from csv
	pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no  Calculating paramters for molecule : {0} ".format(name))

		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
			# Sets the folder and find the log files to analyze
			for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
				#assign the path to the finished directory.
				if str(bs).find('/') > -1:
					w_dir = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files'
				else:
					w_dir = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
				os.chdir(w_dir)
				log_files = get_com_or_log_out_files('output',name)
				calculate_db_parameters(log_files,args,log,w_dir_initial,name,lot,bs)
				calculate_boltz_and_dbstep(log_files,args,log,name,w_dir,w_dir_initial,lot,bs)
		os.chdir(w_dir_initial)

def nics_par_main(args,log,w_dir_initial):
	#get sdf FILES from csv
	pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no  Calculating nics for molecule : {0} ".format(name))

		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
			# Sets the folder and find the log files to analyze
			for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
				#assign the path to the finished directory.
				if str(bs).find('/') > -1:
					w_dir = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files'
				else:
					w_dir = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
				os.chdir(w_dir)
				log_files = get_com_or_log_out_files('output',name)
				#do boltz firsst
				calculate_boltz_for_nics(log_files,args,log,name,w_dir,w_dir_initial,lot,bs)
				for lot_sp,bs_sp,bs_gcp_sp in zip(args.level_of_theory_sp, args.basis_set_sp,args.basis_set_genecp_atoms_sp):
					if str(bs).find('/') > -1:
						w_dir_sp = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/G16-NICS_input_files/'+str(lot_sp)+'-'+str(bs_sp)
					else:
						w_dir_sp = args.path + str(lot) + '-' + str(bs) +'/success/G16-NICS_input_files/'+str(lot_sp)+'-'+str(bs_sp)
					os.chdir(w_dir_sp)
					log_files_sp = get_com_or_log_out_files('output',name)
					calculate_nics_parameters(log_files_sp,args,log,w_dir_initial,name,lot_sp,bs_sp)
					calculate_avg_nics(log_files_sp,args,log,name,w_dir_sp,w_dir_initial,lot_sp,bs_sp)
		os.chdir(w_dir_initial)

def cclib_main(args,log,w_dir_initial):
	#get sdf FILES from csv
	pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no  Calculating cclib paramters for molecule : {0} ".format(name))
		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
			# Sets the folder and find the log files to analyze
			for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
				#assign the path to the finished directory.
				if str(bs).find('/') > -1:
					w_dir = args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files'
				else:
					w_dir = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
				os.chdir(w_dir)
				log_files = get_com_or_log_out_files('output',name)
				#do boltz firsst
				calculate_cclib(log_files,args,log,name,w_dir,w_dir_initial,lot,bs)
				calculate_boltz_for_cclib(log_files,args,log,name,w_dir,w_dir_initial,lot,bs)
				if str(bs).find('/') > -1:
					os.chdir(w_dir_initial + '/QPRED/cclib-json/all_confs_cclib/'+str(lot)+'-'+str(bs).split('/')[0])
				else:
					os.chdir(w_dir_initial + '/QPRED/cclib-json/all_confs_cclib/'+str(lot)+'-'+str(bs))
				json_files = get_com_or_log_out_files('output',name)
				calcualte_average_cclib_parameter(json_files,args,log,name,w_dir,w_dir_initial,lot,bs)


# MAIN OPTION FOR DISCARDING MOLECULES BASED ON USER INPUT DATA (REFERRED AS EXPERIMENTAL RULES)
def exp_rules_main(args,log,exp_rules_active):
	if exp_rules_active:
		if args.verbose:
			log.write("\n   ----- Applying experimental rules to write the new confs file -----")
		# do 2 cases, for RDKit only and RDKIt+xTB
		if args.CMIN == 'xtb':
			conf_files =  glob.glob('*_xtb.sdf')
		if args.CMIN == 'ani':
			conf_files =  glob.glob('*_ani.sdf')
		if args.CMIN is None:
			if args.CSEARCH=='rdkit':
				conf_files =  glob.glob('*_rdkit.sdf')
			elif args.CSEARCH=='summ':
				conf_files =  glob.glob('*_summ.sdf')
			elif args.CSEARCH=='fullmonte':
				conf_files =  glob.glob('*_fullmonte.sdf')
			else:
				conf_files =  glob.glob('*.sdf')

		for file in conf_files:
			try:
				allmols = Chem.SDMolSupplier(file, removeHs=False)
			except OSError:
				pass
			if allmols:
				sdwriter = Chem.SDWriter(file.split('.')[0]+'_filter_exp_rules.sdf')
				print_error_exp_rules = 0
				for mol in allmols:
					check_mol = True
					ob_compat = True
					rdkit_compat = True
					check_mol = exp_rules_output(mol,args,log,file,print_error_exp_rules,ob_compat,rdkit_compat)
					print_error_exp_rules += 1
					if check_mol:
						sdwriter.write(mol)

				sdwriter.close()
