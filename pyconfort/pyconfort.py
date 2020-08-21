#!/usr/bin/env python

#########################################################################################.
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
#########################################################################################.

from __future__ import print_function
import os
import time
from pyconfort.argument_parser import parser_args
from pyconfort.mainf import csearch_main, exp_rules_main, qprep_gaussian_main, move_sdf_main, qcorr_gaussian_main,dup_main,graph_main,geom_par_main,nmr_main,energy_main,creation_of_dup_csv,load_from_yaml,Logger,creation_of_ana_csv

def main():
	# working directory and arguments
	w_dir_initial = os.getcwd()
	args = parser_args()

	log = Logger("pyCONFORT", args.output_name)
	#time
	start_time = time.time()
	#if needed to load from a yaml file
	load_from_yaml(args,log)

	#setting variable if needed
	for i,_ in enumerate(args.basis_set):
		if len(args.basis_set) != len(args.basis_set_genecp_atoms):
			args.basis_set_genecp_atoms.append('')
	for i,_ in enumerate(args.basis_set_sp):
		if len(args.basis_set_sp) != len(args.basis_set_genecp_atoms_sp):
			args.basis_set_genecp_atoms_sp.append('')

	#CSEARCH AND CMIN
	if args.CSEARCH=='rdkit' or args.CSEARCH=='summ' or args.CSEARCH=='fullmonte':
		#creation of csv to write dup data
		dup_data = creation_of_dup_csv(args)
		csearch_main(w_dir_initial,dup_data,args,log,start_time)
		os.chdir(w_dir_initial)

	#applying rules to discard certain conformers based on rules that the user define
	if args.exp_rules != False:
		exp_rules_main(args,log)
		os.chdir(w_dir_initial)

	#QPREP
	if args.QPREP=='gaussian':
		qprep_gaussian_main(w_dir_initial,args,log)
		os.chdir(w_dir_initial)

	if args.CSEARCH=='rdkit' or args.CSEARCH=='summ' or args.CSEARCH=='fullmonte':
		# moving files after compute and/or write_gauss
		move_sdf_main(args)
		os.chdir(w_dir_initial)

	#QCORR
	if args.QCORR=='gaussian':
		log.write("\no  Writing analysis of output files in respective folders in csv_files\n")
		# main part of the duplicate function
		if args.dup:
			duplicates = dup_main(args, log, w_dir_initial)
			os.chdir(w_dir_initial)
		else:
			duplicates = False

		# main part of the output file analyzer for errors/imag freqs
		qcorr_gaussian_main(duplicates,w_dir_initial,args,log)
		os.chdir(w_dir_initial)

	#QPRED
	if args.QPRED=='nmr':
		nmr_main(args,log,w_dir_initial)
	if args.QPRED=='energy':
		energy_main(args,log,w_dir_initial)
	os.chdir(w_dir_initial)

	#QSTAT
	if args.QSTAT=='descp':
		geom_par_main(args,log,w_dir_initial)
	if args.QSTAT=='graph':
		graph_main(args,log,w_dir_initial)
	os.chdir(w_dir_initial)

	log.finalize()

	try:
		os.rename('pyCONFORT_output.dat','pyCONFORT_{0}.dat'.format(args.output_name))
	except FileExistsError:
		os.remove('pyCONFORT_{0}.dat'.format(args.output_name))
		os.rename('pyCONFORT_output.dat','pyCONFORT_{0}.dat'.format(args.output_name))


if __name__ == "__main__":
	main()
