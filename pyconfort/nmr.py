#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	       used for genrating details for nmr       #
#####################################################.

from rdkit.Chem import AllChem as Chem
import numpy as np
import os,subprocess
import glob
import pandas as pd

from pyconfort.cheshire_lookup import cheshire
from pyconfort.qprep_gaussian import moving_files


def nmr_stats(y_exp,y_pred):
	sum=0
	for i,j in zip(y_exp,y_pred):
		sum += abs(float(i)-float(j))
	mae =sum/len(y_pred)
	y_sd = y_exp.astype(float).subtract(y_pred.astype(float))
	sd = y_sd.std()
	return mae,sd

def get_exp_data(args,name,w_dir_initial,final_shieldings,conf_idx,conf_sym):
	if args.nmr_exp == 'fromsdf':
		sdf_file = w_dir_initial+'/'+name+'.sdf'
		sdf_lines = open(sdf_file,'r').readlines()
		found_nmr = 0
		for i,line in enumerate(sdf_lines):
			if line.find('> <NMREDATA_ASSIGNMENT>') > -1:
				start = i+1
				found_nmr = 1
			if len(line.strip()) == 0 and found_nmr == 1 :
				stop = i
				break

		atom_num,atom_sheilding = [],[]
		for i in range(start,stop):
			split_line = sdf_lines[i].split(',')
			if len(split_line) == 3:
				if len(split_line[1].split('-'))>1:
					split_line_sheild = split_line[1].split('-')
					split_line_atom = split_line[2].split('-')
					for idx, shield in zip(conf_idx,final_shieldings):
						if idx in split_line_atom:
							split_line_sheild_min  = split_line_sheild -  shield
							idx_min = split_line_sheild_min.index(min(split_line_sheild_min))
							atom_num.append(split_line_atom[idx_min])
							atom_sheilding.append(split_line_sheild[idx_min])
				else:
					atom_num.append(split_line[2].strip().split('\\')[0])
					atom_sheilding.append(split_line[1].strip())
			if len(split_line) > 3:
				atom_shield_cal,atom_idx_cal,sym = [],[],None
				calcuate_idx="-"
				for j in range(2,len(split_line)):
					atom_idx_cal.append(j)
					atom_shield_cal.append(final_shieldings[conf_idx.index(j)])
					final_shieldings = final_shieldings.pop(conf_idx.index(j))
					conf_idx = conf_idx.pop(conf_idx.index(j))
					sym = conf_sym[(conf_idx.index(j))]
					conf_sym = conf_sym.pop(conf_idx.index(j))

				calculate_avg = np.average(atom_shield_cal)
				calcuate_idx = calcuate_idx.join(atom_idx_cal)
				final_shieldings.append(calculate_avg)
				conf_idx.append(calcuate_idx)
				conf_idx.append(sym)

				atom_num.append(calcuate_idx)
				atom_sheilding.append(split_line[1].strip())

		return atom_num,atom_sheilding

def calculate_nmr(nmr_log_files,args,log,name,w_dir_fin,w_dir_initial,lot_sp,bs_sp,lot,bs):

	#define the lot etc.,
	[opt_method, opt_basis, opt_solv] = [lot,bs,[args.solvent_model,args.solvent_name]]
	[nmr_method, nmr_basis, nmr_solv] = [lot_sp,bs_sp,[args.solvent_model_sp,args.solvent_name]]

	## Matching against the CHESHIRE database of scaling factors
	if args.nmr_online:
		slope, intercept,tms_ref = [],[], None
		for nuc in args.nmr_nucleus:
			scale = cheshire(args.nmr_online, nuc, opt_method, opt_basis, opt_solv, nmr_method, nmr_basis, nmr_solv, args.nmr_aos,log)
			try:
				slope.append(float(scale.split()[1]))
				intercept.append(float(scale.split()[3]))
			except AttributeError:
				log.write("   No scaling factors found for this level of theory! Input the values as arguments for pyCONFORT!"); exit()

			log.write("\no   The slope for nucleus {0} = {1}".format(nuc, slope))
			log.write("\no  The intercept for nucleus {0} = {1}".format(nuc, intercept))

	else:
		slope = args.nmr_slope
		intercept = args.nmr_intercept
		tms_ref = args.nmr_tms_ref

	final_shieldings = []

	for num_file,file in enumerate(nmr_log_files):

		# list of H and C shieldings for each individual conformer
		conf_shieldings, conf_idx,conf_sym = [],[],[]

		# read the file and detects start and stop point enclosing the NMR section of the LOG file
		start, stop = 0,0

		outfile = open(file,"r")
		outlines = outfile.readlines()

		for i in range(len(outlines)):
			if outlines[i].find('SCF GIAO Magnetic shielding tensor (ppm):') > -1:
				start = i
			if outlines[i].find('Population analysis using the SCF Density') > -1:
				stop = i
				break

		# stores H and C NMR shieldings within the NMR section of the LOG file
		for i in range(start,stop):
			try:
				if outlines[i].split()[1] in args.nmr_nucleus:
					atom_nuc = outlines[i].split()[1]
					#assigning values from arrays
					index = args.nmr_nucleus.index(atom_nuc)
					TMS_ref_nuc = tms_ref[index]
					slope_nuc = slope[index]
					intercept_nuc = intercept[index]

					conf_idx.append(outlines[i].split()[0])
					conf_sym.append(atom_nuc)

					if TMS_ref_nuc is not None:
						conf_shieldings.append(TMS_ref_nuc-float(outlines[i].split()[4]))
					else:
						conf_shieldings.append((intercept_nuc-float(outlines[i].split()[4])/(-slope_nuc)))
			except: pass

		outfile.close()

		# get the Boltzmann probabilities of each conformer from a GoodVibes file
		boltz_file = w_dir_initial+'/QPRED/nmr/boltz/'+str(lot)+'-'+str(bs) + '/Goodvibes_'+name+'.dat'

		boltz_outfile = open(boltz_file,"r")
		boltz_outlines = boltz_outfile.readlines()
		for i in range(len(boltz_outlines)):
			# I remove the NMR from the file names using [0:-4]
			if boltz_outlines[i].find(file.split('_NMR')[0]) > -1:
				boltz_factor = float(boltz_outlines[i].split()[-1])
				break

		# Multiply the shieldings by the corresponding Boltzmann factor
		conf_shieldings = [x * boltz_factor for x in conf_shieldings]
		final_shieldings.append(conf_shieldings)


	#final additons
	final_shieldings = np.array(final_shieldings)
	final_shieldings = np.sum(final_shieldings, axis=0)
	conf_idx = np.array(conf_idx)
	conf_sym = np.array(conf_sym)

	#getting experimental shiftd for get_atom
	atom_num_exp,atom_sheilding_exp = get_exp_data(args,name,w_dir_initial,final_shieldings,conf_idx)

	#make final array with atom_num, exp, dft
	df = pd.DataFrame({'Atom': atom_num_exp, 'Shielding-Exp': atom_sheilding_exp }, columns=['Atom', 'Shielding-Exp'])

	for i, atom_num in enumerate(atom_num_exp):
		for j,atom_cal in enumerate(conf_idx):
			if atom_num == atom_cal:
				df.at[i,'Shielding-Cal'] = final_shieldings[j]
				df.at[i,'Atom-Symbol'] = conf_sym[j]

	#calculation of MAE and SD
	for atom_nuc in args.nmr_nucleus:
		df_nuc = df.loc[df['Atom-Symbol'] == atom_nuc]
		mae_nuc,sd_nuc = nmr_stats(df_nuc['Shielding-Exp'],df_nuc['Shielding-Cal'])
		log.write("\no  The MAE and SD for molecule {0} is {1} and {2} respectively for {3}".format(name,mae_nuc,sd_nuc,atom_nuc))

	w_dir = os.getcwd()
	folder = w_dir_initial+'/QPRED/nmr/average-nmr/'+str(lot_sp)+'-'+str(bs_sp)
	log.write("\no  Preparing final NMR results in {}".format(folder))
	try:
		os.makedirs(folder)
		os.chdir(folder)
	except OSError:
		if os.path.isdir(folder):
			os.chdir(folder)
		else:
			pass
	export_param_excel = df.to_csv(name+'_'+lot_sp+'-'+bs_sp+'_predicted_shifts.csv', index = None)

def calculate_boltz_and_nmr(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):

	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
	for file in val:
		cmd_boltz.append(file)
	subprocess.call(cmd_boltz)

	#writing to coorect places
	destination = w_dir_initial+'/QPRED/nmr/boltz/'+str(lot)+'-'+str(bs)
	moving_files(destination, os.getcwd(),'Goodvibes_'+name+'.dat')


	for lot_sp in args.level_of_theory_sp:
		for bs_sp in args.basis_set_sp:
			for bs_gcp_sp in args.basis_set_genecp_atoms_sp:
				#performing the nmr calculations
				dir_sp_nmr =  w_dir_fin+'/../G16-SP_input_files/'+str(lot_sp)+'-'+str(bs_sp)
				os.chdir(dir_sp_nmr)
				#grabbing the respective NMR files for a given molecules
				nmr_log_files = glob.glob(name+'*_NMR.log') + glob.glob(name+'*_NMR_*.log')
				calculate_nmr(nmr_log_files,args,log,name,w_dir_fin,w_dir_initial,lot_sp,bs_sp,lot,bs)
