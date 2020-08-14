#!/usr/bin/env python

#####################################################.
# 		   This file stores all the main functions 	#
#####################################################.

import glob
import sys
import os
import pandas as pd
import subprocess
import yaml
from rdkit.Chem import AllChem as Chem
from pyconfort.csearch import check_for_pieces, check_charge_smi, clean_args, compute_confs, com_2_xyz_2_sdf, mol_from_sdf_or_mol_or_mol2
from pyconfort.filter import exp_rules_output

from pyconfort.qprep_gaussian import read_energies, write_gaussian_input_file, moving_files, convert_xyz_to_sdf
from pyconfort.qcorr_gaussian import output_analyzer, check_for_final_folder, dup_calculation

from pyconfort.grapher import graph
from pyconfort.descp import calculate_parameters
from pyconfort.nmr import calculate_boltz_and_nmr
from pyconfort.energy import calculate_boltz_and_energy,calculate_avg_and_energy
#need to and in energy

#class for logging
class Logger:
	# Class Logger to writargs.input.split('.')[0] output to a file
	def __init__(self, filein, append):
		# Logger to write the output to a file
		suffix = 'dat'
		self.log = open('{0}_{1}.{2}'.format(filein, append, suffix), 'w')

	def write(self, message):
		print(message, end='\n')
		self.log.write(message+ "\n")

	def fatal(self, message):
		print(message, end='\n')
		self.log.write(message + "\n")
		self.finalize()
		sys.exit(1)

	def finalize(self):
		self.log.close()

#load paramters from yaml file
def load_from_yaml(args,log):
	# Variables will be updated from YAML file
	try:
		if args.varfile is not None:
			if os.path.exists(args.varfile):
				if os.path.splitext(args.varfile)[1] == '.yaml':
					log.write("\no  IMPORTING VARIABLES FROM " + args.varfile)
					with open(args.varfile, 'r') as file:
						param_list = yaml.load(file, Loader=yaml.SafeLoader)
			for param in param_list:
				if hasattr(args, param):
					if getattr(args, param) != param_list[param]:
						log.write("o  RESET " + param + " from " + str(getattr(args, param)) + " to " + str(param_list[param]))
						setattr(args, param, param_list[param])
					else:
						log.write("o  DEFAULT " + param + " : " + str(getattr(args, param)))
	except UnboundLocalError:
		log.write("\no  The specified yaml file containing parameters was not found! Make sure that the valid params file is in the folder where you are running the code.\n")

#creation of csv for csearch
def creation_of_dup_csv(args):
	# writing the list of DUPLICATES
	if args.CSEARCH=='rdkit':
		if not args.CMIN=='xtb' and not args.CMIN=='ani1ccx':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','time (seconds)','Overall charge'])
		elif args.CMIN=='xtb' and not args.CMIN=='ani1ccx':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge'])
		elif args.CMIN=='ani1ccx' and not args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','time (seconds)','Overall charge'])
		elif args.CMIN=='ani1ccx' and args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time(seconds)','Overall charge'])
	elif args.CSEARCH=='fullmonte':
		if not args.CMIN=='xtb' and not args.CMIN=='ani1ccx':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-Unique-conformers','time (seconds)','Overall charge'])# ,'FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates',
		elif args.CMIN=='xtb' and not args.CMIN=='ani1ccx':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge']) #'FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates',
		elif args.CMIN=='ani1ccx' and not args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','time (seconds)','Overall charge'])#'FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates',
		elif args.CMIN=='ani1ccx' and args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','FullMonte-conformers','FullMonte-energy-window', 'FullMonte-initial_energy_threshold','FullMonte-RMSD-and-energy-duplicates','FullMonte-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time(seconds)','Overall charge'])
	elif args.CSEARCH=='summ':
		if not args.CMIN=='xtb' and not args.CMIN=='ani1ccx':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','time (seconds)','Overall charge'])
		elif args.CMIN=='xtb' and not args.CMIN=='ani1ccx':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window','RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge'])
		elif args.CMIN=='ani1ccx' and not args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','time (seconds)','Overall charge'])
		elif args.CMIN=='ani1ccx' and args.CMIN=='xtb':
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-initial_energy_threshold','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','summ-conformers','summ-energy-window', 'summ-initial_energy_threshold','summ-RMSD-and-energy-duplicates','summ-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge'])
	return dup_data

#creation of csv for qcorr
def creation_of_ana_csv(args):
	if args.QCORR=='gaussian':
		ana_data =  pd.DataFrame(columns = ['Total Files','Normal Termination', 'Imaginary frequencies', 'SCF Error','Basis Set Error','Other Errors','Unfinished'])
	return ana_data

# main function to generate conformers
def csearch_main(w_dir_initial,dup_data,args,log,start_time):
	# input file format specified
	file_format = os.path.splitext(args.input)[1]

	if file_format not in ['.smi', '.sdf', '.cdx', '.csv','.com','.gjf','.mol','.mol2','.xyz','.txt']:
		log.write("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!")
		sys.exit()

	if not os.path.exists(args.input):
		log.write("\nx  INPUT FILE NOT FOUND!")
		sys.exit()

	# sets up the chosen force field (this fixes some problems in case MMFF is replaced by UFF)
	ori_ff = args.ff
	ori_charge = args.charge_default

	# SMILES input specified
	if file_format == '.smi' or file_format =='.txt':
		smifile = open(args.input)
		#used only for template
		counter_for_template = 0

		for i, line in enumerate(smifile):
			toks = line.split()
			#editing part
			smi = toks[0]
			smi = check_for_pieces(smi)
			mol = Chem.MolFromSmiles(smi)
			clean_args(args,ori_ff,mol,ori_charge)
			if args.charge_default == 'auto':
				if not args.metal_complex:
					args.charge_default = check_charge_smi(smi)
			if args.prefix == 'None':
				name = ''.join(toks[1:])
			else:
				name = str(args.prefix)+str(i)+'_'+''.join(toks[1:])
			compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)

	# CSV file with one columns SMILES and code_name
	elif os.path.splitext(args.input)[1] == '.csv':
		csv_smiles = pd.read_csv(args.input)
		counter_for_template =0
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
			compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)

	# CDX file
	elif os.path.splitext(args.input)[1] == '.cdx':
		#converting to smiles from chemdraw
		cmd_cdx = ['obabel', '-icdx', args.input, '-osmi', '-O', 'cdx.smi']
		subprocess.call(cmd_cdx)
		smifile = open('cdx.smi',"r")

		counter_for_template = 0
		for i, smi in enumerate(smifile):
			smi = check_for_pieces(smi)
			mol = Chem.MolFromSmiles(smi)
			clean_args(args,ori_ff,mol,ori_charge)
			if args.charge_default == 'auto':
				if not args.metal_complex:
					args.charge_default = check_charge_smi(smi)
			name = 'comp' + str(i)+'_'
			compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)

		os.remove('cdx.smi')

	# COM file
	elif os.path.splitext(args.input)[1] == '.gjf' or os.path.splitext(args.input)[1] == '.com' or  os.path.splitext(args.input)[1] == '.xyz' :
		#converting to sdf from comfile to preserve geometry
		charge_com = com_2_xyz_2_sdf(args)
		sdffile = os.path.splitext(args.input)[0]+'.sdf'
		suppl = Chem.SDMolSupplier(sdffile, removeHs=False)
		name = os.path.splitext(args.input)[0]
		counter_for_template = 0
		i=0
		for mol in suppl:
			clean_args(args,ori_ff,mol,ori_charge)
			args.charge_default = charge_com
			compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)
			i += 1

	# SDF file
	elif os.path.splitext(args.input)[1] == '.sdf' or os.path.splitext(args.input)[1] == '.mol' or os.path.splitext(args.input)[1] == '.mol2':
		suppl, IDs, charges = mol_from_sdf_or_mol_or_mol2(args)
		counter_for_template = 0
		i=0
		if os.path.splitext(args.input)[1] == '.sdf':
			for mol,name,charge_sdf in zip(suppl,IDs,charges):
				clean_args(args,ori_ff,mol,ori_charge)
				args.charge_default = charge_sdf
				compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)
				i += 1
		elif os.path.splitext(args.input)[1] == '.mol' or os.path.splitext(args.input)[1] == '.mol2':
			args.charge_default = charges[0]
			name = IDs[0]
			mol = suppl
			compute_confs(w_dir_initial,mol,name,args,log,dup_data,counter_for_template,i,start_time)
			i += 1

	if not os.path.isdir(w_dir_initial+'/CSEARCH/csv_files'):
		os.makedirs(w_dir_initial+'/CSEARCH/csv_files/')
	dup_data.to_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv',index=False)

#writing gauss main
def qprep_gaussian_main(args,log):

	if args.exp_rules:
		conf_files =  glob.glob('*_rules.sdf')
	# define the SDF files to convert to COM Gaussian files
	elif not args.CMIN=='xtb' and not args.CMIN=='ani1ccx' and args.CSEARCH=='rdkit':
		conf_files =  glob.glob('*_rdkit.sdf')
	elif not args.CMIN=='xtb' and not args.CMIN=='ani1ccx' and args.CSEARCH=='summ':
		conf_files =  glob.glob('*_summ.sdf')
	elif not args.CMIN=='xtb' and not args.CMIN=='ani1ccx' and args.CSEARCH=='fullmonte':
		conf_files =  glob.glob('*_fullmonte.sdf')
	elif args.CMIN=='xtb' and not args.CMIN=='ani1ccx':
		conf_files =  glob.glob('*_xtb.sdf')
	elif args.CMIN=='ani1ccx' and not args.CMIN=='xtb':
		conf_files =  glob.glob('*_ani.sdf')
	elif args.CMIN=='ani1ccx' and args.CMIN=='xtb':
		conf_files =  glob.glob('*_ani.sdf') + glob.glob('*_xtb.sdf')
	else:
		conf_files =  glob.glob('*.sdf')

	if args.com_from_xyz:
		xyz_files =  glob.glob('*.xyz')
		convert_xyz_to_sdf(xyz_files,args,log)
		conf_files =  glob.glob('*.sdf')

	# names for directories created
	sp_dir = 'QMCALC/G16-SP'
	g_dir = 'QMCALC/G16'
	#fixing genecp to LAL2DZ if empty

	if len(conf_files) != 0:
		#read in dup_data to get the overall charge of MOLECULES
		try:
			charge_data = pd.read_csv(w_dir_initial+'/pyCONFORT_csv_files/confgen/'+args.input.split('.')[0]+'-Confgen-Data.csv', usecols=['Molecule','Overall charge'])
		except:
			charge_data = None
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					# only create this directory if single point calculation is requested
					if bs.find('**') > -1:
						bs = bs.replace('**','(d,p)')
					elif bs.find('*') > -1:
						bs = bs.replace('*','(d)')
					if args.single_point:
						folder = sp_dir + '/' + str(lot) + '-' + str(bs)
						log.write("\no  Preparing single-point inputs in {}".format(folder))
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
							print('\nx  The QMCALC folder from QPREP could not be created, probably due to incompatible characters')
							folder_error = True
					# writing the com files
					# check conf_file exists, parse energies and then write DFT input
					for file in conf_files:
						if os.path.exists(file):
							if args.verbose:
								if not folder_error:
									log.write("   -> Converting from {}".format(file))
							energies = read_energies(file,log)
							name = file.split('.')[0]

							write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args, log, charge_data)
	else:
		log.write('\nx  No SDF files detected to convert to gaussian COM files')

# moving files after compute and/or write_gauss
def move_sdf_main(args):
	src = os.getcwd()
	if args.CMIN=='xtb':
		all_xtb_conf_files = glob.glob('*_xtb.sdf')
		destination_xtb = src +'/CSEARCH/xtb'
		for file in all_xtb_conf_files:
			moving_files(destination_xtb,src,file)
	if args.CMIN=='ani1ccx':
		all_ani_conf_files = glob.glob('*_ani.sdf')
		destination_ani = src +'/CSEARCH/ani1ccx'
		for file in all_ani_conf_files:
			moving_files(destination_ani,src,file)
	if args.CSEARCH=='rdkit':
		all_name_conf_files = glob.glob('*_rdkit.sdf')
		destination_rdkit = src+ '/CSEARCH/rdkit'
		for file in all_name_conf_files:
			moving_files(destination_rdkit,src,file)
	if args.CSEARCH=='summ':
		all_name_conf_files = glob.glob('*_summ.sdf')
		destination_rdkit = src+ '/CSEARCH/summ'
		for file in all_name_conf_files:
			moving_files(destination_rdkit,src,file)
	if args.CSEARCH=='fullmonte':
		all_name_conf_files = glob.glob('*_fullmonte.sdf')
		destination_rdkit = src+ '/CSEARCH/fullmonte'
		for file in all_name_conf_files:
			moving_files(destination_rdkit,src,file)
	if args.com_from_xyz:
		all_xyz_conf_files = glob.glob('*.xyz')+glob.glob('*.sdf')
		destination_xyz = 'QMCALC/xyz_and_sdf'
		for file in all_xyz_conf_files:
			moving_files(destination_xyz,src,file)

#finding the file type to move for analysis
def get_com_or_log_out_files(type,name=None):
	files = []
	if type =='output':
		formats = ['*.log','*.LOG','*.out','*.OUT']
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
def qcorr_gaussian_main(w_dir_initial,args,log):
	# when you run analysis in a folder full of output files
	if not os.path.exists(w_dir_initial+'/QMCALC'):
		w_dir = os.getcwd()
		w_dir_fin = w_dir+'/success/output_files'
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					folder = w_dir_initial
					ana_data = creation_of_ana_csv(args)
					log.write("\no  ANALYZING OUTPUT FILES IN {}\n".format(folder))
					log_files = get_com_or_log_out_files('output')
					com_files = get_com_or_log_out_files('input')
					output_analyzer(log_files,com_files, w_dir,w_dir, lot, bs, bs_gcp, args, w_dir_fin,w_dir_initial,log,ana_data,1)
		os.chdir(w_dir)
	# when you specify multiple levels of theory
	else:
		if args.QCORR=='gaussian':
			args.path = w_dir_initial+'/QMCALC/G16/'
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					#for main folder
					w_dir_main = args.path + str(lot) + '-' + str(bs)
					if not os.path.isdir(w_dir_main+'/dat_files/'):
						os.makedirs(w_dir_main+'/dat_files/')
					#for currect w_dir folder
					w_dir = args.path + str(lot) + '-' + str(bs)
					#check if New_Gaussian_Input_Files folder exists
					w_dir,round_num = check_for_final_folder(w_dir)
					log = Logger(w_dir_main+'/dat_files/pyCONFORT-analysis-run_'+str(round_num), args.output_name)
					#assign the path to the finished directory.
					w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
					os.chdir(w_dir)
					ana_data = creation_of_ana_csv(args)
					log.write("\no  ANALYZING OUTPUT FILES IN {}\n".format(w_dir))
					log_files = get_com_or_log_out_files('output')
					com_files = get_com_or_log_out_files('input')
					output_analyzer(log_files, com_files, w_dir, w_dir_main , lot, bs, bs_gcp, args, w_dir_fin,w_dir_initial,log,ana_data,round_num)
		os.chdir(args.path)
	os.chdir(w_dir_initial)

#removing the duplicates
def dup_main(args,log,w_dir_initial):
	if not os.path.exists(w_dir_initial+'/QMCALC'):
		w_dir = os.getcwd()
		log_files = get_com_or_log_out_files('output')
		if len(log_files) != 0:
			dup_calculation(log_files,w_dir,w_dir,args,log,1)
		else:
			log.write(' There are no log or out files in this folder.')
	else:
		if args.QCORR=='gaussian':
			args.path = w_dir_initial+'/QMCALC/G16/'
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir_main = args.path + str(lot) + '-' + str(bs)
				w_dir = args.path + str(lot) + '-' + str(bs)
				os.chdir(w_dir)
				w_dir,round_num = check_for_final_folder(w_dir)
				os.chdir(w_dir)
				# change molecules to a range as files will have codes in a continous manner
				log_files = get_com_or_log_out_files('output')
				if len(log_files) != 0:
					dup_calculation(log_files,w_dir,w_dir_main,args,log,round_num)
				else:
					log.write(' There are no any log or out files in this folder.')

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
		if os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name+'_xtb.sdf'):
			sdf_xtb =  w_dir_initial+'/CSEARCH/xtb/'+name+'_xtb.sdf'
		if os.path.exists(w_dir_initial+'/CSEARCH/ani1ccx/'+name+'_ani.sdf'):
			sdf_ani = w_dir_initial+'/CSEARCH/ani1ccx/'+name+'_ani.sdf'
		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
			# Sets the folder and find the log files to analyze
			for lot in args.level_of_theory:
				for bs in args.basis_set:
					for bs_gcp in args.basis_set_genecp_atoms:
						#assign the path to the finished directory.
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
		if os.path.exists(w_dir_initial+'/CSEARCH/all_confs/'+name+'_xtb.sdf'):
			sdf_xtb =  w_dir_initial+'/CSEARCH/all_confs/'+name+'_xtb.sdf'
		if os.path.exists(w_dir_initial+'/CSEARCH/all_confs/'+name+'_ani.sdf'):
			sdf_ani = w_dir_initial+'/CSEARCH/all_confs/'+name+'_ani.sdf'
		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
			# Sets the folder and find the log files to analyze
			for lot in args.level_of_theory:
				for bs in args.basis_set:
					for bs_gcp in args.basis_set_genecp_atoms:
						#assign the path to the finished directory.
						w_dir = args.path + str(lot) + '-' + str(bs) +'/success/output_files'
						os.chdir(w_dir)
						log_files = get_com_or_log_out_files('output',name)
						graph(sdf_rdkit,sdf_xtb,sdf_ani,log_files,args,log,lot,bs,name,w_dir_initial)
		else:
			graph(sdf_rdkit,sdf_xtb,sdf_ani,None,args,log,None,None,name,w_dir_initial)

	os.chdir(w_dir_initial)


#function for compariosn of nmr
def nmr_main(args,log,w_dir_initial):

	#get sdf FILES from csv
	pd_name = pd.read_csv(w_dir_initial+'/CSEARCH/csv_files/'+args.input.split('.')[0]+'-CSEARCH-Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no NMR analysis for molecule : {0} ".format(name))
		if os.path.exists(w_dir_initial+'/QMCALC/G16'):
			args.path = w_dir_initial+'/QMCALC/G16/'
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					#assign the path to the finished directory.
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
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					#assign the path to the finished directory.
					w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/success/output_files/'
					os.chdir(w_dir_fin)
					log_files = get_com_or_log_out_files('output',name)
					if len(log_files) != 0:
						calculate_boltz_and_energy(log_files,args,log,name,w_dir_fin,w_dir_initial,lot,bs)

	#combining the combining all files in different folders
	w_dir_boltz = w_dir_initial+'/QPRED/energy/boltz/'

	for lot in args.level_of_theory:
		for bs in args.basis_set:
			for bs_gcp in args.basis_set_genecp_atoms:
				#assign the path to the finished directory.
				w_dir_fin = w_dir_boltz + str(lot) + '-' + str(bs)
				os.chdir(w_dir_fin)
				dat_files = glob.glob('*.dat')
				if len(dat_files) != 0:
					calculate_avg_and_energy(dat_files,args,log,name,w_dir_fin,w_dir_initial,w_dir_boltz,lot,bs)

	os.chdir(w_dir_initial)


# MAIN OPTION FOR DISCARDING MOLECULES BASED ON USER INPUT DATA (REFERRED AS EXPERIMENTAL RULES)
def exp_rules_main(args,log):
	if args.verbose:
		log.write("\n   ----- Applying experimental rules to write the new confs file -----")
	# do 2 cases, for RDKit only and RDKIt+xTB
	if not args.CMIN=='xtb':
		if args.CSEARCH=='rdkit':
			conf_files =  glob.glob('*_rdkit.sdf')
		elif args.CSEARCH=='summ':
			conf_files =  glob.glob('*_summ.sdf')
	else:
		conf_files =  glob.glob('*_xtb.sdf')

	for file in conf_files:
		allmols = Chem.SDMolSupplier(file, removeHs=False)
		if allmols is None:
			log.write("Could not open "+ file)
			sys.exit(-1)

		sdwriter = Chem.SDWriter(file.split('.')[0]+'_filter_exp_rules.sdf')

		for mol in allmols:
			check_mol = True
			check_mol = exp_rules_output(mol,args,log)
			if check_mol:
				sdwriter.write(mol)
		sdwriter.close()
