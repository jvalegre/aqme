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
###  Version: v1.0, Release date: 24-April-2020										  ###
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
import argparse
import os
import sys
import subprocess
import glob
import time
import pandas as pd
from rdkit.Chem import AllChem as Chem
from DBGEN.db_gen_functions import creation_of_dup_csv,load_from_yaml,compute_confs,clean_args, substituted_mol, Logger, read_energies, write_gaussian_input_file, exp_rules_output,moving_sdf_files
from DBGEN.db_gen_functions import output_analyzer, check_for_final_folder, dup_calculation, combine_files, boltz_calculation

def parser_args():
	parser = argparse.ArgumentParser(description="Generate conformers depending on type of optimization (change parameters in db_gen_PATHS.py file).")

	#Input details
	parser.add_argument("--varfile", dest="varfile", default=None, help="Parameters in YAML format")
	parser.add_argument("-i", "--input", help="File containing molecular structure(s)",dest="input", default=" ")
	parser.add_argument("--output_name", dest="output_name", default="output", metavar="output_name", help="Change output filename to DBGEN_\"output\".dat")
	parser.add_argument("--output", dest="output", default=".sdf", metavar="output", help="The extension of the SDF files written")

	#metal complex
	parser.add_argument("--metal_complex", action="store_true", default=False, help="Request metal complex with coord. no. 4, 5 or 6")
	parser.add_argument("--metal",  help="Specify metallic element", default=["Ir"], dest="metal", type=str)
	parser.add_argument("--complex_spin",  help="Multiplicity of metal complex", default="1", dest="complex_spin", type=int)
	parser.add_argument("--complex_coord", help="Coord. no. of metal complex (automatically updates)", default=[], dest="complex_coord", type=int)
	parser.add_argument("--complex_type",  help="Geometry about metal (e.g. octahedral)", default="octahedral", dest="complex_type", type=str)
	parser.add_argument("--m_oxi",  help="Metal oxidation state", default=["3"], dest="m_oxi", type=int)
	parser.add_argument("--metal_idx",  help="Metal index (automatically updates)", default=[], dest="metal_idx", type=int)
	parser.add_argument("--charge",  help="Charge of metal complex (automatically updates)", default=[], dest="charge", type=int)
	parser.add_argument("--charge_default",  help="Charge default to be considered", default="0", dest="charge_default", type=int)
	parser.add_argument("--metal_sym",  help="Symbols of metals to be considered from list (automatically updates)", default=[], dest="metal_sym", type=str)

	#NCI complex
	parser.add_argument("--nci_complex", action="store_true", default=False, help="Request NCI complexes")
	parser.add_argument("-m", "--maxnumber", help="Number of compounds", type=int, metavar="maxnumber")
	parser.add_argument("--prefix", help="Prefix for naming files", default=None, metavar="prefix")

	#work the script has to do
	parser.add_argument("-w", "--compute", action="store_true", default=False, help="Perform conformational analysis")
	parser.add_argument("--write_gauss", action="store_true", default=False, help="Create input files for Gaussian")
	parser.add_argument("-a", "--analysis", action="store_true", default=False, help="Fix and analyze Gaussian outputs")
	parser.add_argument("-r", "--resubmit", action="store_true", default=False, help="Resubmit Gaussian input files")
	parser.add_argument("--sp", action="store_true", default=False, help="Resubmit Gaussian single point input files")

	#Post analysis
	parser.add_argument("--dup",action="store_true",default=False, help="Remove Duplicates after DFT optimization")
	parser.add_argument("-n","--nmr",action="store_true",default=False, help="Create Files for Single Point which includes NMR calculation after DFT Optimization")
	parser.add_argument("-b","--boltz", action="store_true", default=False, help="Boltzmann factor for each conformers from Gaussian output files")
	parser.add_argument("-f","--combine", action="store_true", default=False, help="Combine files of differnt molecules including boltzmann weighted energies")

	#aaply exp rules
	parser.add_argument("--exp_rules", action="store_true", default=False, help="Experimental rules applied to make Gaussian input files")
	parser.add_argument("--angle_off", type=float,help="Any limit to set for check rules",default=30)

	#pass the argument for path for the gaussian folder.
	parser.add_argument("--path", help="Path for analysis/boltzmann factor/combining files where the gaussian folder created is present",dest="path", default="")
	parser.add_argument("-v","--verbose",action="store_true",default=False, help="verbose output")

	#argumets for conformer generation
	parser.add_argument("--ANI1ccx", "--ani", action="store_true",default=False, help="request ANI1ccx optimizations")
	parser.add_argument("--xtb", action="store_true",default=False, help="request xtb optimizations")
	parser.add_argument("--ewin_min", action="store",default=40.0, help="energy window to print conformers for minimization using xTB or ANI1ccx (kcal/mol)", type=float)
	parser.add_argument("--ewin_rdkit", action="store",default=40.0, help="energy window to print conformers for RDKit (kcal/mol)", type=float)
	parser.add_argument("--opt_fmax", action="store",default=0.05, help="fmax value used in xTB and AN1 optimizations", type=float)
	parser.add_argument("--opt_steps", action="store",default=1000, help="max cycles used in xTB and AN1 optimizations", type=int)
	parser.add_argument("--opt_steps_RDKit", action="store",default=1000, help="max cycles used in RDKit optimizations", type=int)
	parser.add_argument("--time","-t",action='store_true',default=False,help="request program runtime")
	parser.add_argument("--heavyonly", help="only consider torsion angles involving heavy (non H) elements (default=True)", default=True, metavar="heavyonly")
	parser.add_argument("--constraints", action="store_true",default=None, help="distance constraint")
	parser.add_argument("--nodihedrals", action="store_true", default=False, help="turn off dihedral scan")
	parser.add_argument("-d","--degree", type=float,help="Amount, in degrees, to enumerate torsions by (default 30.0)",default=30.0)
	parser.add_argument("--max_torsions",type=int,help="Skip any molecules with more than this many torsions (default 5)",default=5)
	parser.add_argument("--sample", help="number of conformers to sample to get non-torsional differences (default 100)", default=100, type=int, metavar="sample")
	parser.add_argument("--auto_sample", help="final factor to multiply in the auto mode for the sample option (default 20)", default=20, type=int, metavar="auto_sample")
	parser.add_argument("--ff", help="force field (MMFF or UFF)", default="MMFF", metavar="ff")
	parser.add_argument("--seed", help="random seed (default 062609)", default="062609", type=int, metavar="s")
	parser.add_argument("--rms_threshold", help="cutoff for considering sampled conformers the same (default 0.25)", default=0.25, type=float, metavar="R")
	parser.add_argument("--max_matches_RMSD", help="iteration cutoff for considering  matches in sampled conformers the same (default 1000000 )", default=1000000 , type=float, metavar="max_matches_RMSD")
	parser.add_argument("--energy_threshold", dest="energy_threshold",action="store",default=0.05, help="energy difference between unique conformers")
	parser.add_argument("--initial_energy_threshold", dest="initial_energy_threshold",action="store",default=0.01, help="energy difference between unique conformers for the first filter of only E")
	parser.add_argument("--max_MolWt", help="Max. molecular weight of molecule", default=1000, type=int, metavar="max_MolWt")
	parser.add_argument("--large_sys", action="store_true",default=False, help="Large systems for xtb optimizations")
	parser.add_argument("--STACKSIZE", help="STACKSIZE for optimization of large systems", default="500m")

	#arguments for gaussian files Creation
	parser.add_argument("-l", "--level_of_theory",help="Level of Theory", default=['wB97xd'], dest="level_of_theory", type=str, nargs='*')
	parser.add_argument("--basis_set",  help="Basis Set", default=['6-31g*'], dest="basis_set", type=str, nargs='*')
	parser.add_argument("--basis_set_genecp_atoms",default=['LANL2DZ'], help="Basis Set genecp/gen: Can specify only one as basis_set", dest="basis_set_genecp_atoms", type=str, nargs='?')
	parser.add_argument("--genecp_atoms",  help="genecp atoms",default=[], dest="genecp_atoms",type=str, nargs='*')
	parser.add_argument("--gen_atoms",  help="gen atoms",default=[], dest="gen_atoms",type=str, nargs='*')
	parser.add_argument("--max_cycle_opt", help="Number of cycles for DFT optimization", default="300", type=int, dest="max_cycle_opt")
	parser.add_argument("--frequencies",action="store_true", default=False, help="Request only optimization without any frequency calculation")
	parser.add_argument("--single_point",action="store_true", default=False, help="Request only single point calculation")
	parser.add_argument("--lowest_only", action="store_true", default=False, help="Lowest conformer to write for gaussian")
	parser.add_argument("--lowest_n", action="store_true", default=False, help="Lowest Number of conformers to write for gaussian")
	parser.add_argument("--energy_threshold_for_gaussian", help="cutoff for considering sampled conformers for gaussian input", default="4.0", type=float, dest="energy_threshold_for_gaussian")
	parser.add_argument("--dispersion_correction",action="store_true", default=False, help="Add Dispersion Correction")
	parser.add_argument("--empirical_dispersion",  help="Type of Dispersion ", default="D3BJ", dest="empirical_dispersion", type=str)
	parser.add_argument("--solvent_model",  help="Type of solvent model", default="gas_phase", dest="solvent_model", type=str)
	parser.add_argument("--solvent_name",  help="Name of Solvent", default="Acetonitrile", dest="solvent_name", type=str)
	parser.add_argument("--nprocs", help="Number of Processors", default="24", type=int, dest="nprocs")
	parser.add_argument("--mem", help="Memory", default="96GB", type=str, dest="mem")
	parser.add_argument("--chk", action="store_true", default=False, help="Create .chk files for Gaussian")


	#autoprep kind of single point inputs
	parser.add_argument("--level_of_theory_sp",help="Level of Theory for single point after optimization", default=['wB97xd'], dest="level_of_theory_sp", type=str, nargs='*')
	parser.add_argument("--basis_set_sp",  help="Basis Set for single point after optimization", default=['6-31g*'], dest="basis_set_sp", type=str, nargs='*')
	parser.add_argument("--basis_set_genecp_atoms_sp",default=['LANL2DZ'], help="Basis Set genecp/gen: Can specify only one for single point after optimization", dest="basis_set_genecp_atoms_sp", type=str, nargs='?')
	parser.add_argument("--dispersion_correction_sp",action="store_true", default=False, help="Add Dispersion Correction for single point after optimization")
	parser.add_argument("--empirical_dispersion_sp",  help="Type of Dispersion for single point after optimization", default="D3BJ", dest="empirical_dispersion_sp", type=str)
	parser.add_argument("--solvent_model_sp",  help="Type of solvent model for single point after optimization", default="gas_phase", dest="solvent_model_sp", type=str)
	parser.add_argument("--solvent_name_sp",  help="Name of Solvent for single point after optimization", default="Acetonitrile", dest="solvent_name_sp", type=str)
	parser.add_argument("--input_for_sp",  help="Input line for Single point after DFT optimization ", default="nmr=giao", dest="input_for_sp", type=str)
	parser.add_argument("--last_line_for_sp",  help="Last input line for Single point after DFT optimization ", default="", dest="last_line_for_sp", type=str)


	# submIssion of Gaussion files
	parser.add_argument("--qsub", action="store_true", default=False, help="Submit Gaussian files")
	parser.add_argument("--submission_command",  help="Queueing system that the submission is done on", default="qsub_summit", metavar="submission_command", type=str)

	args = parser.parse_args()

	return args

def main():

	args = parser_args()
	# Define the logging object
	log = Logger("DBGEN", args.output_name)
	#time
	start_time = time.time()
	#if needed to load from a yaml file
	load_from_yaml(args,log)

	#creation of csv to write dup data
	dup_data = creation_of_dup_csv(args)

	# this will perform conformational analysis and create inputs for Gaussian
	if args.compute:

		# input file format specified
		file_format = os.path.splitext(args.input)[1]

		if file_format not in ['.smi', '.sdf', '.cdx', '.csv','.com','.gjf']:
			log.write("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!")
			sys.exit()

		# sets up the chosen force field (this fixes some problems in case MMFF is replaced by UFF)
		ori_ff = args.ff

		# SMILES input specified
		if file_format == '.smi':
			smifile = open(args.input)
			#used only for template
			counter_for_template =0
			for i, line in enumerate(smifile):
				toks = line.split()
				#editing part
				smi = toks[0]
				clean_args(args,ori_ff,smi)
				if not args.prefix:
					name = ''.join(toks[1:])
				else:
					name = args.prefix+str(i)+'_'+''.join(toks[1:])

				compute_confs(smi,name,args,log,dup_data,counter_for_template,i,start_time)
			dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)

		# CSV file with one columns SMILES and code_name
		elif os.path.splitext(args.input)[1] == '.csv':
			csv_smiles = pd.read_csv(args.input)
			counter_for_template =0
			for i in range(len(csv_smiles)):
				#assigning names and smi i  each loop
				if not args.prefix:
					name = csv_smiles.loc[i, 'code_name']
				# else:
				# 	name = 'comp_'+str(m)+'_'+csv_smiles.loc[i, 'code_name']
				smi = csv_smiles.loc[i, 'SMILES']
				clean_args(args,ori_ff,smi)
				compute_confs(smi,name,args,log,dup_data,counter_for_template,i,start_time)
			dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)

		# CDX file
		elif os.path.splitext(args.input)[1] == '.cdx':
			#converting to smiles from chemdraw
			cmd_cdx = ['obabel', '-icdx', args.input, '-osmi', '-O', 'cdx.smi']
			subprocess.run(cmd_cdx)
			smifile = open('cdx.smi',"r")

			counter_for_template = 0
			for i, smi in enumerate(smifile):
				clean_args(args,ori_ff,smi)
				name = 'comp' + str(i)+'_'
				compute_confs(smi,name,args,log,dup_data,counter_for_template,i,start_time)
			dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)

	#applying rule to get the necessary conformers only
	if args.exp_rules:
		if args.verbose:
			log.write("   ----- Applying experimental rules to write the new confs file -----")
		### do 2 cases, for RDKit only and RDKIt+xTB
		#grab all the gaussian files
		if not args.xtb and not args.ANI1ccx:
			conf_files =  glob.glob('*_rdkit.sdf')
		elif args.xtb:
			conf_files =  glob.glob('*_xtb.sdf')
		elif args.ANI1ccx:
			conf_files =  glob.glob('*_ani.sdf')

		for file in conf_files:
			allmols = Chem.SDMolSupplier(file, removeHs=False)
			if allmols is None:
				log.write("Could not open "+ file)
				sys.exit(-1)

			sdwriter = Chem.SDWriter(file.split('.')[0]+args.exp_rules_output_ext)

			for mol in allmols:
				check_mol = True
				check_mol = exp_rules_output(mol,args,log)
				if check_mol:
					sdwriter.write(mol)
			sdwriter.close()

	if args.write_gauss:
		if args.exp_rules:
			conf_files =  glob.glob('*_rules.sdf')
		# define the SDF files to convert to COM Gaussian files
		elif not args.xtb and not args.ANI1ccx and args.nodihedrals:
				conf_files =  glob.glob('*_rdkit.sdf')
		elif not args.xtb and not args.ANI1ccx and not args.nodihedrals:
			conf_files =  glob.glob('*_rdkit_rotated.sdf')
		elif args.xtb:
			conf_files =  glob.glob('*_xtb.sdf')
		elif args.ANI1ccx:
			conf_files =  glob.glob('*_ani.sdf')
		else:
			conf_files =  glob.glob('*.sdf')

		# names for directories created
		sp_dir = 'generated_sp_files'
		g_dir = 'generated_gaussian_files'

		#read in dup_data to get the overall charge of MOLECULES
		charge_data = pd.read_csv(args.input.split('.')[0]+'-Duplicates Data.csv', usecols=['Molecule','Overall charge'])

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
					# check conf_file exists, parse energies and then write dft input
					for file in conf_files:
						if os.path.exists(file):
							if args.verbose:
								log.write("   -> Converting from {}".format(file))
							energies = read_energies(file,log)
							name = os.path.splitext(file)[0]

							write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args,log,charge_data)

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
	all_name_conf_files = glob.glob('*.sdf')
	destination_rdkit = 'rdkit_generated_sdf_files'
	for file in all_name_conf_files:
		moving_sdf_files(destination_rdkit,src,file)

	if args.analysis:
		#adding in for general analysis
		#need to specify the lot, bs as arguments for each analysis
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
