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
	sd = y_pred.std()
	return mae,sd

def get_exp_data(args,name,w_dir_initial):
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
			if len(sdf_lines[i].split(',')) == 3:
				atom_num.append(sdf_lines[i].split(',')[2].strip().split('\\')[0])
				atom_sheilding.append(sdf_lines[i].split(',')[1].strip())
			if len(sdf_lines[i].split(',')) > 3:
				for j in range(2,len(sdf_lines[i].split(','))):
					atom_num.append(sdf_lines[i].split(',')[j].strip().split('\\')[0])
					atom_sheilding.append(sdf_lines[i].split(',')[1].strip())

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

	final_H_shieldings,final_C_shieldings = [],[]

	for num_file,file in enumerate(nmr_log_files):

		# list of H and C shieldings for each individual conformer
		conf_H_shieldings, conf_H_idx, conf_C_shieldings, conf_C_idx, nmr_file_idx,conf_H_sym,conf_C_sym = [],[],[],[],[],[],[]

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
				if '1H' in args.nmr_nucleus:
					if outlines[i].split()[1] == 'H':
						#assigning values from arrays
						index = args.nmr_nucleus.index('1H')
						TMS_H_ref = tms_ref[index]
						slope_H = slope[index]
						intercept_H = intercept[index]

						#conf_H_idx.append('H'+str(outlines[i].split()[0]))
						conf_H_idx.append(outlines[i].split()[0])
						conf_H_sym.append('H')
						if TMS_H_ref is not None:
							conf_H_shieldings.append(TMS_H_ref-float(outlines[i].split()[4]))
						else:
							conf_H_shieldings.append((intercept_H-float(outlines[i].split()[4])/(-slope_H)))
				if '13C' in args.nmr_nucleus:
					if outlines[i].split()[1] == 'C':
						#assigning values from arrays
						index = args.nmr_nucleus.index('13C')
						TMS_C_ref = tms_ref[index]
						slope_C = slope[index]
						intercept_C = intercept[index]

						#conf_C_idx.append('C'+str(outlines[i].split()[0]))
						conf_C_idx.append(outlines[i].split()[0])
						conf_C_sym.append('C')
						if TMS_C_ref is not None:
							conf_C_shieldings.append(TMS_C_ref-float(outlines[i].split()[4]))
						else:
							conf_C_shieldings.append((intercept_C-float(outlines[i].split()[4])/(-slope_C)))
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
		if '1H' in args.nmr_nucleus:
			conf_H_shieldings = [x * boltz_factor for x in conf_H_shieldings]
			final_H_shieldings.append(conf_H_shieldings)

		if '13C' in args.nmr_nucleus:
			conf_C_shieldings = [x * boltz_factor for x in conf_C_shieldings]
			#conf_C_shieldings = np.asarray(conf_C_shieldings, dtype=np.float64)*boltz_factor
			final_C_shieldings.append(conf_C_shieldings)

	#final additons
	if '1H' in args.nmr_nucleus:
		final_H_shieldings = np.array(final_H_shieldings)
		final_H_shieldings = np.sum(final_H_shieldings, axis=0)
		conf_H_idx = np.array(conf_H_idx)
		conf_H_sym = np.array(conf_H_sym)


	if '13C' in args.nmr_nucleus:
		final_C_shieldings = np.array(final_C_shieldings)
		final_C_shieldings = np.sum(final_C_shieldings,axis=0)
		conf_C_idx = np.array(conf_C_idx)
		conf_C_sym = np.array(conf_C_sym)

	# concatenate H and C data one after the other (for printing in the CSV)
	if final_C_shieldings is not None:
		all_shieldings = np.concatenate((final_H_shieldings, final_C_shieldings))
		all_idx = np.concatenate((conf_H_idx, conf_C_idx))
		all_sym = np.concatenate((conf_H_sym, conf_C_sym))


	#getting experimental shiftd for get_atom
	atom_num_exp,atom_sheilding_exp = get_exp_data(args,name,w_dir_initial)

	#make final array with atom_num, exp, dft
	df = pd.DataFrame({'Atom': atom_num_exp, 'Shielding-Exp': atom_sheilding_exp }, columns=['Atom', 'Shielding-Exp'])

	for i, atom_num in enumerate(atom_num_exp):
		for j,atom_cal in enumerate(all_idx):
			if atom_num == atom_cal:
				df.at[i,'Shielding-Cal'] = all_shieldings[j]
				df.at[i,'Atom-Symbol'] = all_sym[j]

	#calculation of MAE and SD
	if '1H' in args.nmr_nucleus:
		df_H = df.loc[df['Atom-Symbol'] == 'H']
		mae_H,sd_H = nmr_stats(df_H['Shielding-Exp'],df_H['Shielding-Cal'])
		log.write("\no  The MAE and SD for molecule {0} is {1} and {2} respectively for 1H".format(name,mae_H,sd_H))
	if '13C' in args.nmr_nucleus:
		df_C = df.loc[df['Atom-Symbol'] == 'C']
		mae_C,sd_C = nmr_stats(df_C['Shielding-Exp'],df_C['Shielding-Cal'])
		log.write("\no  The MAE and SD for molecule {0} is {1} and {2} respectively for 13C".format(name,mae_C,sd_C))

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
