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
from pyconfort.main_functions import compute_main, exp_rules_main, write_gauss_main, move_sdf_main, analysis_main, dup_main, qsub_main,graph_main
from pyconfort.writer_functions import creation_of_dup_csv, load_from_yaml, Logger

def main():
	# working directory and arguments
	w_dir_initial = os.getcwd()
	args = parser_args()
	# Define the logging object
	log = Logger("pyCONFORT", args.output_name)
	#time
	start_time = time.time()
	#if needed to load from a yaml file
	load_from_yaml(args,log)

	#setting defaults back
	if len(args.basis_set_genecp_atoms) == 0:
		args.basis_set_genecp_atoms = ['LANL2DZ']

	#creation of csv to write dup data
	dup_data = creation_of_dup_csv(args)

	# this will perform conformational analysis and create inputs for Gaussian
	if args.compute:
		compute_main(w_dir_initial,dup_data,args,log,start_time)

	#applying rules to discard certain conformers based on rules that the user define
	if args.exp_rules:
		exp_rules_main(args,log)

	# main part for writing COM files from SDF files
	if args.write_gauss:
		write_gauss_main(args,log)

	# moving files after compute and/or write_gauss
	move_sdf_main(args)

	# main part of the duplicate function
	if args.dup:
		dup_main(args,log)

	# main part of the analysis functions
	if args.analysis:
		analysis_main(w_dir_initial,args,log)

	# main part of the automated workflow (submission of COM files and analyzer)
	if args.qsub:
		qsub_main(args,log)

	# main part of the automated workflow (submission of COM files and analyzer)
	if args.graph:
		graph_main(args,log,w_dir_initial)

if __name__ == "__main__":
	main()
