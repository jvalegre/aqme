#!/usr/bin/env python

""".####################################################################################.
#########################################################################################
###																					  ###
###  pyCONFORT is a tool that allows to carry out automated:						  ###
###  (1) Conformational searches and creation of COM files using RDKit, xTB and ANI1  ###
###  (2) LOG file processing (detects imaginary freqs and error terminations		  ###
###      and creates new COM files)													  ###
###  (3) Use LOG files to create new COM files with new keywords (i.e. single-point   ###
###      corrections after geometry optimization)									  ###
###  																				  ###
#########################################################################################
###  																				  ###
###  Version: v1.0.1, Release date: 22-May-2020								     	  ###
###  																				  ###
#########################################################################################
###  																				  ###
###  Authors: Shree Sowndarya S. V., Juan V. Alegre Requena, Robert S. Paton		  ###
###  																				  ###
###  Please, report any bugs or suggestions to:										  ###
###  svss@colostate.edu or juanvi89@hotmail.com  									  ###
###																					  ###
#########################################################################################
.####################################################################################."""

from __future__ import print_function
import os
import subprocess
import glob
import time
from pyconfort.argument_parser import parser_args
from pyconfort.confgen_functions import compute_main, clean_args, substituted_mol
from pyconfort.analyzer_functions import output_analyzer, check_for_final_folder, dup_calculation, combine_files, boltz_calculation
from pyconfort.writer_functions import creation_of_dup_csv, load_from_yaml, Logger, moving_sdf_files, write_gauss_main
from pyconfort.filter_functions import exp_rules_main

def main():

	args = parser_args()
	# Define the logging object
	log = Logger("pyCONFORT", args.output_name)
	#time
	start_time = time.time()
	#if needed to load from a yaml file
	load_from_yaml(args,log)

	#creation of csv to write dup data
	dup_data = creation_of_dup_csv(args)

	# this will perform conformational analysis and create inputs for Gaussian
	if args.compute:
		compute_main(dup_data,args,log,start_time)

	#applying rules to discard certain conformers based on rules that the user define
	if args.exp_rules:
		exp_rules_main(args,log)

	if args.write_gauss:
		write_gauss_main(args,log)

	#moving files after compute and write_gauss or only after compute
	#moving all the sdf files to a separate folder after writing gaussian files
	src = os.getcwd()
	if args.xtb:
		all_xtb_conf_files = glob.glob('*_xtb.sdf')
		destination_xtb = src +'/xtb_minimised_generated_sdf_files'
		for file in all_xtb_conf_files:
			moving_sdf_files(destination_xtb ,src,file)
	elif args.ANI1ccx:
		all_ani_conf_files = glob.glob('*_ani.sdf')
		destination_ani = src +'/ani1ccx_minimised_generated_sdf_files'
		for file in all_ani_conf_files:
			moving_sdf_files(destination_ani,src,file)
	else:
		all_name_conf_files = glob.glob('*rdkit*.sdf')
		destination_rdkit = 'rdkit_generated_sdf_files'
		for file in all_name_conf_files:
			moving_sdf_files(destination_rdkit,src,file)

	if args.analysis:
		# adding in for general analysis
		# need to specify the lot, bs as arguments for each analysis
		if args.path == '':
			log_files = glob.glob('*.LOG'.lower())+glob.glob('*.LOG')
			w_dir = os.getcwd()
			w_dir_fin = w_dir+'/finished'
			for lot in args.level_of_theory:
				for bs in args.basis_set:
					for bs_gcp in args.basis_set_genecp_atoms:
							output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin,log)

		#taking the path
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
						#log.write(w_dir)
						os.chdir(w_dir)
						#log.write(w_dir)
						log_files = glob.glob('*.log')+glob.glob('*.LOG')
						output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin,log)

	#adding the part to check for resubmission of the newly created gaussian files.
	if args.qsub:
		#chceck if ech level of theory has a folder New gaussin FILES
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'
				#check if New_Gaussian_Input_Files folder exists
				w_dir = check_for_final_folder(w_dir,log)
				os.chdir(w_dir)
				cmd_qsub = [args.submission_command, '*.com']
				subprocess.run(cmd_qsub)

	#once all files are finished are in the Finished folder
	if args.dup:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'+'finished'
				os.chdir(w_dir)
				#can change molecules to a range as files will have codes in a continous manner
				try:
					log_files = glob.glob('*.log')
					if len(log_files) != 0:
						val = ' '.join(log_files)
						dup_calculation(val,w_dir,args,log)
					else:
						log.write(' Files for are not there!')

				except:
					pass

	#once all files are finished are in the Finished folder
	if args.boltz:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'+'finished'
				os.chdir(w_dir)
				#can change molecules to a range as files will have codes in a continous manner
				for i in range(args.maxnumber):
					#grab all the corresponding files make sure to renamme prefix when working with differnet files
					try:
						log_files = glob.glob('RE' + '_' + str(i)+'_'+'confs_low.log')
						if len(log_files) != 0:
							val = ' '.join(log_files)
							boltz_calculation(val,i,log)
						else:
							log.write(' Files for {} are not there!'.format(i))
					except:
						pass

	if args.combine:
		#combines the files and gives the boltzmann weighted energies
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) + '/finished'
				os.chdir(w_dir)
				#read the csv log_files
				csv_files = glob.glob('Goodvibes*.csv')
				combine_files(csv_files, lot, bs, args,log)

if __name__ == "__main__":
	main()
