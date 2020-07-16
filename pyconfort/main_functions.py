#!/usr/bin/env python

#####################################################.
# 		   This file stores all the main functions 	#
#####################################################.

import glob
import sys
import os
import pandas as pd
import subprocess
from rdkit.Chem import AllChem as Chem
from pyconfort.confgen_functions import check_for_pieces, check_charge_smi, clean_args, compute_confs, com_2_xyz_2_sdf, mol_from_sdf_or_mol_or_mol2
from pyconfort.writer_functions import read_energies, write_gaussian_input_file, moving_sdf_files, convert_xyz_to_sdf
from pyconfort.filter_functions import exp_rules_output
from pyconfort.analyzer_functions import output_analyzer, check_for_final_folder, dup_calculation
from pyconfort.grapher import graph

# main function to generate conformers
def compute_main(w_dir_initial,dup_data,args,log,start_time):
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
		suppl = Chem.SDMolSupplier(sdffile)
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

	dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)

def write_gauss_main(args,log):

	if args.exp_rules:
		conf_files =  glob.glob('*_rules.sdf')
	# define the SDF files to convert to COM Gaussian files
	elif not args.xtb and not args.ANI1ccx and args.nodihedrals:
		conf_files =  glob.glob('*_rdkit.sdf')
	elif not args.xtb and not args.ANI1ccx and not args.nodihedrals:
		conf_files =  glob.glob('*_rdkit_rotated.sdf')
	elif args.xtb and not args.ANI1ccx:
		conf_files =  glob.glob('*_xtb.sdf')
	elif args.ANI1ccx and not args.xtb:
		conf_files =  glob.glob('*_ani.sdf')
	elif args.ANI1ccx and args.xtb:
		conf_files =  glob.glob('*_ani.sdf') + glob.glob('*_xtb.sdf')
	else:
		conf_files =  glob.glob('*.sdf')

	if args.com_from_xyz:
		xyz_files =  glob.glob('*.xyz')
		convert_xyz_to_sdf(xyz_files,args,log)
		conf_files =  glob.glob('*.sdf')
		print(conf_files)

	# names for directories created
	sp_dir = 'generated_sp_files'
	g_dir = 'generated_gaussian_files'
	#fixing genecp to LAL2DZ if empty

	if len(conf_files) != 0:
		#read in dup_data to get the overall charge of MOLECULES
		if not args.com_from_xyz:
			charge_data = pd.read_csv(args.input.split('.')[0]+'-Duplicates Data.csv', usecols=['Molecule','Overall charge'])
		else:
			charge_data = None
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					# only create this directory if single point calculation is requested
					if args.single_point:
						folder = sp_dir + '/' + str(lot) + '-' + str(bs)
						log.write("\no  PREPARING SINGLE POINT INPUTS in {}".format(folder))
					else:
						folder = g_dir + '/' + str(lot) + '-' + str(bs)
						log.write("\no  Preparing Gaussian COM files in {}".format(folder))
					try:
						os.makedirs(folder)
					except OSError:
						if os.path.isdir(folder):
							pass
						else:
							raise
					# writing the com files
					# check conf_file exists, parse energies and then write DFT input
					for file in conf_files:
						if os.path.exists(file):
							if args.verbose:
								log.write("   -> Converting from {}".format(file))
							energies = read_energies(file,log)
							name = file.split('.')[0]

							write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args, log, charge_data)
	else:
		log.write('\nx  No SDF files detected to convert to gaussian COM files')

# moving files after compute and/or write_gauss
def move_sdf_main(args):
	src = os.getcwd()
	if args.xtb:
		all_xtb_conf_files = glob.glob('*_xtb.sdf') + glob.glob('*_xtb_all_confs.sdf')
		destination_xtb = src +'/xtb_minimised_generated_sdf_files'
		for file in all_xtb_conf_files:
			moving_sdf_files(destination_xtb,src,file)
	if args.ANI1ccx:
		all_ani_conf_files = glob.glob('*_ani.sdf') + glob.glob('*_ani_all_confs.sdf')
		destination_ani = src +'/ani1ccx_minimised_generated_sdf_files'
		for file in all_ani_conf_files:
			moving_sdf_files(destination_ani,src,file)
	if args.compute or args.write_gauss:
		all_name_conf_files = glob.glob('*_rdkit*.sdf')
		destination_rdkit = 'rdkit_generated_sdf_files'
		for file in all_name_conf_files:
			moving_sdf_files(destination_rdkit,src,file)
	if args.com_from_xyz:
		all_xyz_conf_files = glob.glob('*.xyz')+glob.glob('*.sdf')
		destination_xyz = 'xyz_and_sdf_files_from_xyz'
		for file in all_xyz_conf_files:
			moving_sdf_files(destination_xyz,src,file)

def get_log_out_files():
	log_files = []
	output_formats = ['*.log','*.LOG','*.out','*.OUT']
	for _,output_format in enumerate(output_formats):
		for _,file_out in enumerate(glob.glob(output_format)):
			if file_out not in log_files:
				log_files.append(file_out)
	return log_files

# main part of the analysis functions
def analysis_main(w_dir_initial,args,log,ana_data):
	# when you run analysis in a folder full of output files
	if args.path == '':
		log_files = get_log_out_files()
		w_dir = os.getcwd()
		w_dir_fin = w_dir+'/finished'

		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					folder = w_dir_initial
					log.write("\no  ANALYZING OUTPUT FILES IN {}\n".format(folder))
					output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin,w_dir_initial,log,ana_data)
		os.chdir(w_dir)
	# when you specify multiple levels of theory
	else:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					w_dir = args.path + str(lot) + '-' + str(bs)
					#check if New_Gaussian_Input_Files folder exists
					w_dir = check_for_final_folder(w_dir,log)
					#assign the path to the finished directory.
					w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/finished'
					os.chdir(w_dir)
					log_files = get_log_out_files()
					folder = w_dir + '/' + str(lot) + '-' + str(bs)
					log.write("\no  ANALYZING OUTPUT FILES IN {}\n".format(folder))
					output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin,w_dir_initial,log,ana_data)
		os.chdir(args.path)

def dup_main(args,log):
	if args.path == '':
		w_dir = os.getcwd()
		log_files = get_log_out_files()
		if len(log_files) != 0:
			dup_calculation(log_files,w_dir,args,log)
		else:
			log.write(' There are not any log or out files in this folder.')
	else:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs)
				os.chdir(w_dir)
				# change molecules to a range as files will have codes in a continous manner
				log_files = get_log_out_files()
				if len(log_files) != 0:
					dup_calculation(log_files,w_dir,args,log)
				else:
					log.write(' There are not any log or out files in this folder.')

def qsub_main(args,log):
	#chceck if ech level of theory has a folder New gaussin FILES
	for lot in args.level_of_theory:
		for bs in args.basis_set:
			w_dir = args.path + str(lot) + '-' + str(bs) +'/'
			#check if New_Gaussian_Input_Files folder exists
			w_dir = check_for_final_folder(w_dir,log)
			os.chdir(w_dir)
			cmd_qsub = [args.submission_command, '*.com']
			subprocess.call(cmd_qsub)

def graph_main(args,log,w_dir_initial):

	#get sdf FILES from csv
	pd_name = pd.read_csv(args.input.split('.')[0]+'-Duplicates Data.csv')

	for i in range(len(pd_name)):
		name = pd_name.loc[i,'Molecule']

		log.write("\no  Plotting graphs for molecule : {0} ".format(name))

		sdf_ani,sdf_xtb = None,None
		if args.nodihedrals:
			sdf_rdkit =  w_dir_initial+'/rdkit_generated_sdf_files/'+name+'_rdkit.sdf'
		elif not args.nodihedrals:
			sdf_rdkit = w_dir_initial+'/rdkit_generated_sdf_files/'+name+'_rdkit_rotated.sdf'
		if args.xtb:
			sdf_xtb =  w_dir_initial+'/xtb_minimised_generated_sdf_files/'+name+'_xtb_all_confs.sdf'
		if args.ANI1ccx:
			sdf_ani = w_dir_initial+'/ani1ccx_minimised_generated_sdf_files/'+name+'_ani_all_confs.sdf'

		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					#assign the path to the finished directory.
					w_dir = args.path + str(lot) + '-' + str(bs) +'/finished'
					os.chdir(w_dir)
					log_files = glob.glob(name+'_*.log')
					graph(sdf_rdkit,sdf_xtb,sdf_ani,log_files,args,log,lot,bs,name,w_dir_initial)


# MAIN OPTION FOR DISCARDING MOLECULES BASED ON USER INPUT DATA (REFERRED AS EXPERIMENTAL RULES)
def exp_rules_main(args,log):
	if args.verbose:
		log.write("\n   ----- Applying experimental rules to write the new confs file -----")
	# do 2 cases, for RDKit only and RDKIt+xTB
	if not args.xtb:
		if args.nodihedrals:
			conf_files =  glob.glob('*_rdkit.sdf')
		else:
			conf_files =  glob.glob('*_rdkit_rotated.sdf')
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
